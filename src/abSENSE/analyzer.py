import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats


from abSENSE.parameters import AbsenseParameters
from abSENSE.exceptions import MissingGeneException, MissingSpeciesException
from abSENSE.results import FitResult, ErrorResult, NotEnoughDataResult, SampledResult


class AbsenseAnalyzer():
    def __init__(self, params: AbsenseParameters, seed=0):
        self.random = np.random.default_rng(seed)
        self.e_value = params.e_value
        self.distances = pd.read_csv(
            params.distances,
            delimiter='\t',
            comment='#',
            names=['specie', 'distance'],
            index_col=0,
        )

        self.bitscores = pd.read_csv(
            params.bitscores,
            delimiter='\t',
            comment='#',
            index_col=0,
        )

        self.species = self.distances.index
        self.genes = self.bitscores.index

        if params.include_only is not None:
            self.species = pd.Series(params.include_only.split(','))

        if params.gene_lengths is None:
            self.gene_lengths = pd.DataFrame(
                {'length': params.default_gene_length},
                index=self.genes,
            )
        else:
            self.gene_lengths = pd.read_csv(
                params.gene_lengths,
                delimiter='\t',
                comment='#',
                names=['gene', 'length'],
                index_col=0,
            )

        if params.db_lengths is None:
            self.db_lengths = pd.DataFrame(
                {'length': params.default_db_length},
                index=self.species,
            )
        else:
            self.db_lengths = pd.read_csv(
                params.db_lengths,
                delimiter='\t',
                comment='#',
                names=['specie', 'length'],
                index_col=0,
            )

        self.validate_species()
        self.validate_genes()

    def validate_species(self):
        """Checks that the requested species are present in the bitscores and lengths."""
        # all species in species are in bitscores
        self.bitscores = self.bitscores.loc[:, self.bitscores.columns.isin(self.species)]
        # species not found in column
        if not self.species.isin(self.bitscores.columns).all():
            missing = list(self.species[~self.species.isin(self.bitscores.columns)])
            raise MissingSpeciesException(
                'Unable to find all requested species in bitscores. '
                f'Missing: {missing}'
            )

        # reorder columns to match species
        self.bitscores = self.bitscores[self.species]

        if not self.species.isin(self.db_lengths.index).all():
            missing = list(self.species[~self.species.isin(self.db_lengths.index)])
            raise MissingSpeciesException(
                'Unable to find all requested species in database lengths. '
                f'Missing: {missing}'
            )

    def validate_genes(self):
        """Checks that the genes present in bitscores are found in gene lengths."""

        if not self.genes.isin(self.gene_lengths.index).all():
            missing = list(self.genes[~self.genes.isin(self.gene_lengths.index)])
            raise MissingGeneException(
                'Unable to find all requested genes in gene lengths. '
                f'Missing: {missing}'
            )

    def total_genes(self):
        return len(self.bitscores)

    def fit_genes(self):
        for row in self.bitscores.iterrows():
            yield self.fit_gene(*row)


    def fit_gene(self, gene, bitscore: pd.Series) -> FitResult:
        data = self.distances.merge(
            bitscore.rename('score'),
            left_index=True,
            right_index=True,
        )

        gene_length = self.gene_lengths.loc[gene]
        ambiguous = data['score'].isna()
        in_fit = ~data['score'].isna() & (data['score'] != 0)

        if in_fit.sum() <= 2:
            return NotEnoughDataResult(gene)

        try:
            (a_fit, b_fit), covariance = scipy.optimize.curve_fit(
                exponential,
                data.loc[in_fit, 'distance'],
                data.loc[in_fit, 'score'],
                bounds=((-np.inf, 0), (np.inf, np.inf))
            )
        except (RuntimeError, scipy.optimize.OptimizeWarning):
            return ErrorResult(gene)

        predictions = exponential(self.distances, a_fit, b_fit)

        if np.isinf(covariance).any():
            return ErrorResult(gene, predictions=predictions['distance'])

        intervals = self.find_confidence_interval(
            self.sample_parameters(a_fit, b_fit, covariance),
            self.bit_threshold(gene_length),
        )

        result = data.join(intervals).assign(
            ambiguous=ambiguous,
            in_fit=in_fit,
            prediction=predictions,
        ).round(2)

        return SampledResult(
            gene,
            result=result,
            a_fit=a_fit,
            b_fit=b_fit,
        )


    def sample_parameters(self, a_fit, b_fit, covariance):
        """
        function to, where possible, use maximum likelihood estimates of a and b
        parameter plus estimated covariance matrix to directly sample from the
        probability distribution of a and b (assume Gaussian with mean of max
        likelihood estimates and given covariance structure)
        """
        return self.random.multivariate_normal([a_fit, b_fit], covariance, size=200)

    def bit_threshold(self, gene_length):
        result = np.log2(gene_length) + np.log2(self.db_lengths) - np.log2(self.e_value)
        return result.rename(columns={'length': 'bit_threshold'})

    def find_confidence_interval(self, sampled_parameters, bit_threshold):
        """Gives an empirical estimate of the prediction interval.
        function to take each of the sampled a, b values and use them to sample
        directly from the distribution of scores taking into account the Gaussian
        noise (a function of distance, a, b)."""
        # species x samples
        point_estimates = exponential(
            self.distances.to_numpy(),
            sampled_parameters[:, 0],
            sampled_parameters[:, 1],
        )

        # species x samples
        exp = np.exp(-1*sampled_parameters[:, 1] * self.distances.to_numpy())
        variance = sampled_parameters[:, 0] * (1 - exp) * exp

        # if variance is negative, need to ignore those samples
        invalid = variance <= 0
        # set to 0 to handle issues with sqrt and random normal
        variance[invalid] = 0

        # with normal, can't broadcast a shape on matrix inputs, instead
        # perform a standard normal draw and do rescaling manually
        # species x samples x draws
        draws = (self.random.standard_normal(size=(*point_estimates.shape, 200)) *
                 np.sqrt(variance[..., None]) + point_estimates[..., None])

        # set invalid values to nan, keep first entry as point estimate
        # this works because the invalid entries have 0 variance.
        draws[invalid, 1:] = np.nan

        # now find mean and std, ignoring nans
        mean = np.nanmean(draws, axis=(1, 2))
        std = np.nanstd(draws, axis=(1, 2))

        # calculate p values analytically from std estimate
        p_values = scipy.stats.norm.cdf(bit_threshold['bit_threshold'], mean, std)
        low, high = scipy.stats.norm.interval(0.99, mean, std)

        return pd.DataFrame(
            {'p_values': p_values,
             'low_interval': low,
             'high_interval': high,
             },
            index=self.species)



def exponential(x, a, b):
    return a * np.exp(-b * x)
