"""Module for performing the abSENSE analysis."""
from __future__ import annotations

from typing import Callable, Generator

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats

from abSENSE.exceptions import MissingGeneException, MissingSpeciesException
from abSENSE.parameters import AbsenseParameters
from abSENSE.results import ErrorResult, FitResult, NotEnoughDataResult, SampledResult
from abSENSE.utilities import exponential, find_confidence_interval, sample_parameters


class AbsenseAnalyzer:
    """The analyzer object."""

    def __init__(self, params: AbsenseParameters, seed: int = 0):
        self.random = np.random.default_rng(seed)
        self.e_value = params.e_value
        self.distances = pd.read_csv(
            params.distances,
            delimiter="\t",
            comment="#",
            names=np.array(["specie", "distance"]),
            index_col=0,
        )

        self.bitscores = pd.read_csv(
            params.bitscores,
            delimiter="\t",
            comment="#",
            index_col=0,
        )

        self.species = self.distances.index
        self.genes = self.bitscores.index

        if params.include_only is not None:
            self.species = pd.Index(params.include_only.split(","))

        if params.gene_lengths is None:
            self.gene_lengths = pd.DataFrame(
                {"length": params.default_gene_length},
                index=self.genes,
            )
        else:
            self.gene_lengths = pd.read_csv(
                params.gene_lengths,
                delimiter="\t",
                comment="#",
                names=np.array(["gene", "length"]),
                index_col=0,
            )

        if params.db_lengths is None:
            self.db_lengths = pd.DataFrame(
                {"length": params.default_db_length},
                index=self.species,
            )
        else:
            self.db_lengths = pd.read_csv(
                params.db_lengths,
                delimiter="\t",
                comment="#",
                names=np.array(["specie", "length"]),
                index_col=0,
            )

        self.validate_species()
        self.validate_genes()

        self.gene_filter: Callable[[str], bool] = lambda x: True
        if params.need_plots():
            self.gene_filter = params.plot_test()

    def validate_species(self) -> None:
        """Checks that the requested species are present in the bitscores and lengths."""
        # all species in species are in bitscores
        self.bitscores = self.bitscores.loc[
            :, self.bitscores.columns.isin(self.species)
        ]
        # species not found in column
        if not self.species.isin(self.bitscores.columns).all():
            missing = list(self.species[~self.species.isin(self.bitscores.columns)])
            raise MissingSpeciesException(
                "Unable to find all requested species in bitscores. "
                f"Missing: {missing}"
            )

        # reorder columns to match species
        self.bitscores = self.bitscores[self.species]

        if not self.species.isin(self.db_lengths.index).all():
            missing = list(self.species[~self.species.isin(self.db_lengths.index)])
            raise MissingSpeciesException(
                "Unable to find all requested species in database lengths. "
                f"Missing: {missing}"
            )

    def validate_genes(self) -> None:
        """Checks that the genes present in bitscores are found in gene lengths."""

        if not self.genes.isin(self.gene_lengths.index).all():
            missing = list(self.genes[~self.genes.isin(self.gene_lengths.index)])
            raise MissingGeneException(
                "Unable to find all requested genes in gene lengths. "
                f"Missing: {missing}"
            )

    def total_genes(self) -> int:
        """The total number of genes in this analysis.

        Returns:
            Number of genes to analyze
        """
        return sum(self.gene_filter(gene) for gene in self.bitscores.index)

    def fit_genes(self) -> Generator[FitResult, None, None]:
        """Fit all genes in analysis.

        Returns:
            Generator of FitResults for each gene
        """
        for gene, bitscore in self.bitscores.iterrows():
            if self.gene_filter(gene):
                yield self.fit_gene(str(gene), bitscore)

    def fit_gene(self, gene: str, bitscore: pd.Series) -> FitResult:
        """Fit a gene.

        Args:
            gene: gene name
            bitscore: the bitscores for the gene
        """
        data = self.distances.merge(
            bitscore.rename("score"),
            left_index=True,
            right_index=True,
        )

        gene_length = self.gene_lengths.loc[gene]
        ambiguous = data["score"].isna()
        in_fit = ~data["score"].isna() & (data["score"] != 0)

        if in_fit.sum() <= 2:
            return NotEnoughDataResult(gene)

        try:
            fit_result = scipy.optimize.curve_fit(
                exponential,
                data.loc[in_fit, "distance"],
                data.loc[in_fit, "score"],
                bounds=((-np.inf, 0), (np.inf, np.inf)),
            )
            (a_fit, b_fit), covariance = fit_result[:2]
        except (RuntimeError, scipy.optimize.OptimizeWarning):
            return ErrorResult(gene)

        predictions = exponential(self.distances, a_fit, b_fit)

        if np.isinf(covariance).any():
            return ErrorResult(gene, predictions=predictions["distance"])

        bit_threshold = self.bit_threshold(gene_length)
        global_threshold = bit_threshold.mean()["bit_threshold"]
        correlation = self.correlation(data.score)

        intervals = find_confidence_interval(
            self.random,
            self.distances.to_numpy(),
            sample_parameters(self.random, a_fit, b_fit, covariance),
            bit_threshold,
            index=self.species,
        )

        result = data.join(intervals).assign(
            ambiguous=ambiguous,
            in_fit=in_fit,
            prediction=predictions,
        )

        return SampledResult(
            gene,
            result=result,
            a_fit=a_fit,
            b_fit=b_fit,
            bit_threshold=global_threshold,
            correlation=correlation,
            covariance=covariance,
        )

    def bit_threshold(self, gene_length: float) -> pd.DataFrame:
        """From the provided gene_length, estimate bit_thresholds for all species."""
        result = (
            self.db_lengths.apply(np.log2)
            + np.log2(gene_length)
            - np.log2(self.e_value)
        )
        return pd.DataFrame({"bit_threshold": result["length"]}, index=result.index)

    def correlation(self, scores: pd.Series) -> float:
        """For the given distance/scores, estimate the correlation."""
        data = self.distances.assign(score=scores)
        data = data[data.score > 0]
        _, _, correlation, _, _ = scipy.stats.linregress(
            data.distance, np.log(data.score)
        )
        return float(correlation)
