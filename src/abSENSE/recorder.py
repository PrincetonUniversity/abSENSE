"""Recorder module for handling conversion and recording of absense analyses."""
from __future__ import annotations

import os
from abc import ABC, abstractmethod
from contextlib import contextmanager
from typing import TYPE_CHECKING, Generator, TextIO

import numpy as np
import pandas as pd

from abSENSE.parameters import AbsenseParameters
from abSENSE.plotting import FitPlot

if TYPE_CHECKING:
    from abSENSE.results import SampledResult


class Recorder(ABC):
    """Base class of an object for recording results."""

    def __init__(self, params: AbsenseParameters, species: pd.Index):
        os.makedirs(params.output_directory, exist_ok=True)
        self.output_dir = params.output_directory
        self.species = species
        self._info_file: TextIO | None = None

    @abstractmethod
    @contextmanager
    def open(self) -> Generator[Recorder, None, None]:
        """Context managed resources."""

    @abstractmethod
    def write_headers(self) -> None:
        """Called once prior to streaming data, setup headers."""

    def write_info(self, params: AbsenseParameters) -> None:
        """Write the run information to the info file."""
        if self._info_file is None:
            return

        self._info_file.write(f"abSENSE analysis run on {params.start_time}\n")
        self._info_file.write(f"Input bitscore file: {params.bitscores.name}\n")
        self._info_file.write(f"Input distance file: {params.distances.name}\n")

        if params.gene_lengths is None:
            self._info_file.write(
                f"Gene length (for all genes): {params.default_gene_length} (default)\n"
            )
        else:
            self._info_file.write(f"Gene length file: {params.gene_lengths.name}\n")

        self._info_file.write("Database length ")
        if params.db_lengths is None:
            self._info_file.write(
                f"(for all species): {params.default_db_length} (default)\n"
            )
        else:
            self._info_file.write(f"file: {params.db_lengths.name}\n")

        self._info_file.write("Species used in fit: ")
        if params.include_only is None:
            self._info_file.write("All (default)")
        else:
            self._info_file.write(" ".join(params.include_only.split(",")))
        self._info_file.write("\n")

        self._info_file.write(f"E-value threshold: {params.e_value}\n")

    @staticmethod
    def build_recorder(params: AbsenseParameters, species: pd.Index) -> Recorder:
        """Factory method to return Recorder subtype based on parameters."""
        if params.validate:
            return ValidationRecorder(params, species)
        return FileRecorder(params, species)


class FileRecorder(Recorder):
    """Contains file handles to record absense analysis results as text files."""

    def __init__(self, params: AbsenseParameters, species: pd.Index):
        super().__init__(params, species)
        self.predict_all = params.predict_all
        self.should_plot = params.plot_test()
        if params.need_plots():
            os.makedirs(f"{self.output_dir}/plots", exist_ok=True)
        self.filenames = {
            "bitscores": "predicted_bitscores.tsv",
            "low": "99PI_lower_prediction.tsv",
            "high": "99PI_high_prediction.tsv",
            "failure": "failure_probabilities.tsv",
            "params": "parameters.tsv",
        }
        self._files: dict[str, TextIO] = {}

    @contextmanager
    def open(self) -> Generator[Recorder, None, None]:
        """Context managed file opening."""
        for key, file in self.filenames.items():
            self._files[key] = open(f"{self.output_dir}/{file}", "w", encoding="utf8")
        self._info_file = open(f"{self.output_dir}/run_info.txt", "w", encoding="utf8")

        try:
            yield self

        finally:
            for handle in self._files.values():
                handle.close()
            self._info_file.close()

    def write_headers(self) -> None:
        """Write brief summary header for each file."""
        species_header = "Gene\t" + "\t".join(self.species) + "\n"
        self._files["bitscores"].write(
            "# maximum likelihood bitscore predictions "
            "for each tested gene in each species\n" + species_header
        )

        self._files["low"].write(
            "# the lower bound of the 99% bitscore prediction "
            "interval for each tested gene in each species\n" + species_header
        )

        self._files["high"].write(
            "# the upper bound of the 99% bitscore prediction "
            "interval for each tested gene in each species\n" + species_header
        )

        self._files["failure"].write(
            "# the probability of a homolog being undetected "
            "at the specified significance threshold (see run info file) in each "
            "tested gene in each species\n" + species_header
        )

        self._files["params"].write(
            "# the best-fit parameters "
            "(performed using only bitscores from species not omitted from the fit; "
            "see run info file) for a and b for each gene\n"
            "Gene\ta\tb\n"  # not species header
        )

    def write_gene(self, gene: str) -> None:
        """Record the gene column for all files."""
        for file in self._files.values():
            file.write(gene)

    def analysis_error(self, predictions: list[float] | None = None) -> None:
        """Record an analysis error for this gene."""
        self._write_str("analysis_error", predictions)

    def not_enough_data(self) -> None:
        """Record the gene does not have enough data."""
        self._write_str("not_enough_data")

    def _write_str(
        self,
        entry: str,
        bitscore_overrides: list[float] | None = None,
    ) -> None:
        """Write the entry for each column in this row."""
        line = "\t" + "\t".join([entry] * len(self.species)) + "\n"
        for key, file in self._files.items():
            if key == "params":
                file.write(f"\t{entry}\t{entry}\n")
            elif bitscore_overrides is not None and key == "bitscores":
                file.write(
                    "\t"
                    + "\t".join(str(round(val, 2)) for val in bitscore_overrides)
                    + "\n"
                )
            else:
                file.write(line)

    def write_result(
        self,
        prediction: float,
        high: float,
        low: float,
        pval: float,
        realscore: float,
        is_considered: bool,
        is_ambiguous: bool,
    ) -> None:
        """Record the fit values for this gene, species."""

        high = round(high, 2)
        low = round(low, 2)
        pval = round(pval, 2)

        site_type = ""
        additional_score = ""

        if is_considered:
            site_type = "Ortholog_detected"
            additional_score = f":{realscore}"
        elif is_ambiguous:
            site_type = "Homolog_detected(orthology_ambiguous)"

        if site_type == "":
            self._files["bitscores"].write(f"\t{prediction}")
            self._files["high"].write(f"\t{high}")
            self._files["low"].write(f"\t{low}")
            self._files["failure"].write(f"\t{pval}")

        else:
            if self.predict_all:
                self._files["bitscores"].write(
                    f"\t{prediction}({site_type}{additional_score})"
                )
                self._files["high"].write(f"\t{high}({site_type})")
                self._files["low"].write(f"\t{low}({site_type})")
                self._files["failure"].write(f"\t{pval}({site_type})")
            else:
                for key, file in self._files.items():
                    if key != "params":
                        file.write(f"\t{site_type}")

    def write_params(
        self,
        a_prediction: float,
        b_prediction: float,
    ) -> None:
        """Record the param values for this gene."""
        self._files["params"].write(f"\t{a_prediction}\t{b_prediction}")

    def finalize_row(self) -> None:
        """Write newlines to all files."""
        for file in self._files.values():
            file.write("\n")

    def plot(
        self,
        result: SampledResult,
        outfile: str | None = None,
    ) -> None:
        """Produce a plot from the provided SampledResult.

        Args:
            result: the sampled result
            outfile: if set, will save to the specified file.  Default to result.gene
        """
        if not self.should_plot(result.gene):
            return

        if outfile is None:
            outfile = f"{self.output_dir}/plots/{result.gene}.svg"

        plot = FitPlot()

        plot.generate_plot(result)

        plot.save(outfile)


class ValidationRecorder(Recorder):
    """Contains file handles to record absense analysis results as text files."""

    def __init__(self, params: AbsenseParameters, species: pd.Index):
        super().__init__(params, species)
        self.species = species
        self.results = pd.DataFrame(
            np.zeros((len(self.species), 3)),
            index=self.species,
            columns=["count", "mre", "p>0.05"],
        )

        self.p_values = pd.Series(
            "",
            index=self.species,
        )
        self._result_file: None | TextIO = None

    @contextmanager
    def open(self) -> Generator[Recorder, None, None]:
        """Context managed file opening."""
        self._info_file = open(f"{self.output_dir}/run_info.txt", "w", encoding="utf8")
        self._result_file = open(
            f"{self.output_dir}/validation.tsv", "w", encoding="utf8"
        )

        try:
            yield self
            # on return, write output
            self._result_file.write(
                "# count: number of valid estimates\n"
                "# mre: mean, percent relative error\n"
                "# p>0.05: number of estimates with p(not detectable) > 0.05\n"
                "# p_values: comma separated list of p(not detectable)\n"
            )

            self.results["mre"] = self.results.apply(
                lambda x: x["mre"] / x["count"] * 100 if x["count"] != 0 else np.nan,
                axis=1,
            )
            # remove any repeated commas (\3) or at front (\1) or back (\2)
            self.results["p_values"] = self.p_values.str.replace(
                r"(^,+)|(,+$)|(,)+", r"\3", regex=True
            )
            # loc maintains original order
            self.results.loc[self.species].to_csv(
                self._result_file,
                sep="\t",
                float_format="%.2f",
                columns=["count", "mre", "p>0.05", "p_values"],
            )

        finally:
            self._info_file.close()
            self._result_file.close()

    def write_info(self, params: AbsenseParameters) -> None:
        """Write the run information to the info file."""
        super().write_info(params)
        if self._info_file is None:
            return
        self._info_file.write("Run in validation mode\n")

    def write_headers(self) -> None:
        pass

    def aggregate(self, result: pd.DataFrame) -> None:
        """Record the result from a single gene to this validation recorder."""
        self.results = self.results.combine(
            result.loc[:, result.columns != "p_values"],  # type: ignore[index]
            np.add,
            fill_value=0,
        )

        # fails if no results were SampledResults
        if "p_values" in result:
            self.p_values = self.p_values.combine(
                result["p_values"],
                lambda a, b: str(a) + str(b) + ",",
                fill_value="",
            )
