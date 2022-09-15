"""Wrapper for parameters from click."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, TextIO


@dataclass
class AbsenseParameters:
    """Collection of parameters from the command line."""

    bitscores: TextIO
    db_lengths: TextIO | None
    distances: TextIO
    e_value: float
    gene_lengths: TextIO | None
    include_only: str | None
    out_dir: str | None
    plot: str
    plot_all: bool
    predict_all: bool
    start_time: str
    validate: bool
    default_gene_length: float = 400
    default_db_length: float = 8_000_000

    @property
    def output_directory(self) -> str:
        """Get the output directory supplied or the default."""
        if self.out_dir is None:
            return f"abSENSE_results_{self.start_time}"

        return self.out_dir

    def need_plots(self) -> bool:
        """Return true if the params require plotting."""
        return self.plot_all or self.plot != ""

    def plot_test(self) -> Callable[[str], bool]:
        """Return a function to test if a gene should be plotted."""
        targets = set(self.plot.split(","))

        def should_plot(gene: str) -> bool:
            if self.plot_all:
                return True
            return gene in targets

        return should_plot
