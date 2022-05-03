from dataclasses import dataclass
from datetime import datetime
from typing import Optional, TextIO


@dataclass
class AbsenseParameters:
    """Collection of parameters from the command line."""
    distances: TextIO
    bitscores: TextIO
    e_value: float
    include_only: Optional[str]
    gene_lengths: Optional[TextIO]
    db_lengths: Optional[TextIO]
    predict_all: bool
    out_dir: Optional[str]
    start_time: str

    @property
    def output_directory(self) -> str:
        """Get the output directory supplied or the default."""
        if self.out_dir is None:
            return f"abSENSE_results_{self.start_time}"

        return self.out_dir

    @property
    def default_gene_length(self) -> float:
        return 400

    @property
    def default_db_length(self) -> float:
        return 8_000_000
