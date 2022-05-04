"""Recorder module for handling conversion and recording of absense analyses."""
from contextlib import contextmanager
import os
from datetime import datetime

from abSENSE.parameters import AbsenseParameters


# TODO once analysis is contained in it's own object, refactor to recorder subclasses
class FileRecorder:
    """Contains file handles to record absense analysis results as text files."""
    def __init__(self, params, species):
        os.makedirs(params.output_directory, exist_ok=True)
        self.output_dir = params.output_directory
        self.species = species
        self.predict_all = params.predict_all
        self.filenames = {
            'bitscores': 'predicted_bitscores.tsv',
            'low': '99PI_lower_prediction.tsv',
            'high': '99PI_high_prediction.tsv',
            'failure': 'failure_probabilities.tsv',
            'params': 'parameters.tsv',
        }
        self._files = {}
        self._info_file = None

    @contextmanager
    def open(self):
        """Context managed file opening."""
        for key, file in self.filenames.items():
            self._files[key] = open(f'{self.output_dir}/{file}', 'w')
        self._info_file = open(f'{self.output_dir}/run_info.txt', 'w')

        try:
            yield self

        finally:
            for handle in self._files.values():
                handle.close()
            self._info_file.close()

    def write_headers(self):
        """Write brief summary header for each file."""
        species_header = "Gene\t" + "\t".join(self.species) + "\n"
        self._files['bitscores'].write(
            "# maximum likelihood bitscore predictions "
            "for each tested gene in each species\n"
            + species_header
        )

        self._files['low'].write(
            "# the lower bound of the 99% bitscore prediction "
            "interval for each tested gene in each species\n"
            + species_header
        )

        self._files['high'].write(
            "# the upper bound of the 99% bitscore prediction "
            "interval for each tested gene in each species\n"
            + species_header
        )

        self._files['failure'].write(
            "# the probability of a homolog being undetected "
            "at the specified significance threshold (see run info file) in each "
            "tested gene in each species\n"
            + species_header
        )

        self._files['params'].write(
            "# the best-fit parameters "
            "(performed using only bitscores from species not omitted from the fit; "
            "see run info file) for a and b for each gene\n"
            "Gene\ta\tb\n"  # not species header
        )

    def write_info(
        self,
        params: AbsenseParameters,
    ):
        """Write the run information to the info file."""
        self._info_file.write(f"abSENSE analysis run on {params.start_time}\n")
        self._info_file.write(f"Input bitscore file: {params.bitscores.name}\n")
        self._info_file.write(f"Input distance file: {params.distances.name}\n")

        if params.gene_lengths is None:
            self._info_file.write(f"Gene length (for all genes): {params.default_gene_length} (default)\n")
        else:
            self._info_file.write(f"Gene length file: {params.gene_lengths}\n")

        if params.db_lengths is None:
            self._info_file.write(f"Database length (for all species): {params.default_db_length} (default)\n")
        else:
            self._info_file.write(f"Database length file: {params.db_lengths}\n")

        self._info_file.write("Species used in fit: ")
        if params.include_only is None:
            self._info_file.write("All (default)")
        else:
            self._info_file.write(' '.join(params.include_only.split(',')))
        self._info_file.write("\n")

        self._info_file.write(f"E-value threshold: {params.e_value}\n")

    def write_gene(self, gene):
        """Record the gene column for all files."""
        for file in self._files.values():
            file.write(gene)

    def analysis_error(self, predictions=None):
        """Record an analysis error for this gene."""
        self._write_str('analysis_error', predictions)

    def not_enough_data(self):
        """Record the gene does not have enough data."""
        self._write_str('not_enough_data')

    def _write_str(self, entry, bitscore_overrides=None):
        """Write the entry for each column in this row."""
        line = "\t" + "\t".join([entry] * len(self.species)) + "\n"
        for key, file in self._files.items():
            if key == 'params':
                file.write(f"\t{entry}\t{entry}\n")
            elif bitscore_overrides is not None and key == 'bitscores':
                file.write(
                    "\t" +
                    "\t".join(str(round(val, 2)) for val in bitscore_overrides) +
                    "\n")
            else:
                file.write(line)

    def write_result(
        self,
        prediction,
        high,
        low,
        pval,
        realscore,
        is_considered,
        is_ambiguous,
    ):
        """Record the fit values for this gene, species."""

        high = round(high, 2)
        low = round(low, 2)
        pval = round(pval, 2)

        site_type = ''
        additional_score = ''

        if is_considered:
            site_type = 'Ortholog_detected'
            additional_score = f':{realscore}'
        elif is_ambiguous:
            site_type = 'Homolog_detected(orthology_ambiguous)'

        if site_type == '':
            self._files['bitscores'].write(f"\t{prediction}")
            self._files['high'].write(f"\t{high}")
            self._files['low'].write(f"\t{low}")
            self._files['failure'].write(f"\t{pval}")

        else:
            if self.predict_all:
                self._files['bitscores'].write(f"\t{prediction}({site_type}{additional_score})")
                self._files['high'].write(f"\t{high}({site_type})")
                self._files['low'].write(f"\t{low}({site_type})")
                self._files['failure'].write(f"\t{pval}({site_type})")
            else:
                for key, file in self._files.items():
                    if key != 'params':
                        file.write(f'\t{site_type}')

    def write_params(
        self,
        a_prediction,
        b_prediction,
    ):
        """Record the param values for this gene."""
        self._files['params'].write(f"\t{a_prediction}\t{b_prediction}")

    def finalize_row(self):
        """Write newlines to all files."""
        for file in self._files.values():
            file.write('\n')
