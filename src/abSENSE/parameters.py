"""Wrapper for parameters from click."""
from __future__ import annotations

from dataclasses import dataclass
from typing import TextIO


@dataclass
class AbsenseParameters:
    """Collection of parameters from the command line."""

    distances: TextIO
    bitscores: TextIO
    e_value: float
    include_only: str | None
    gene_lengths: TextIO | None
    db_lengths: TextIO | None
    predict_all: bool
    plot_all: bool
    out_dir: str | None
    start_time: str

    @property
    def output_directory(self) -> str:
        """Get the output directory supplied or the default."""
        if self.out_dir is None:
            return f"abSENSE_results_{self.start_time}"

        return self.out_dir

    @property
    def default_gene_length(self) -> float:
        """Default gene length."""
        return 400

    @property
    def default_db_length(self) -> float:
        """Default database length."""
        return 8_000_000
