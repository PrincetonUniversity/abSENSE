"""Result objects to represent various outcomes from analysis."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pandas as pd

from abSENSE.recorder import FileRecorder


class FitResult:
    """Base class for fitting results."""

    def __init__(self, gene: str):
        self.gene = gene
        self.status = ""

    def __eq__(self, other: object) -> bool:
        """Test if this fitting result is the same as other."""
        return self.__dict__ == other.__dict__

    def record_to(self, recorder: FileRecorder) -> None:
        """Record this result to the FileRecorder.

        Args:
            recorder: a FileRecorder
        """
        recorder.write_gene(self.gene)


class ErrorResult(FitResult):
    """An analysis result when an analysis error occurs."""

    def __init__(self, gene: str, predictions: list[float] | None = None):
        super().__init__(gene)
        self.predictions = predictions
        self.status = "Analysis Error"

    def record_to(self, recorder: FileRecorder) -> None:
        super().record_to(recorder)
        recorder.analysis_error(predictions=self.predictions)


class NotEnoughDataResult(FitResult):
    """Analysis result when not enough data is present."""

    def __init__(self, gene: str):
        super().__init__(gene)
        self.status = "Not Enough Data"

    def record_to(self, recorder: FileRecorder) -> None:
        super().record_to(recorder)
        recorder.not_enough_data()


class SampledResult(FitResult):
    """Analysis result when the fit is successful."""

    def __init__(
        self,
        gene: str,
        result: pd.DataFrame,
        a_fit: float,
        b_fit: float,
        bit_threshold: pd.DataFrame,
        correlation: float,
        covariance: npt.NDArray[np.float64],
    ):
        super().__init__(gene)
        self.result = result
        self.a_fit = a_fit
        self.b_fit = b_fit
        self.bit_threshold = bit_threshold
        self.correlation = correlation
        self.covariance = covariance

    def record_to(self, recorder: FileRecorder) -> None:
        super().record_to(recorder)
        recorder.plot(result=self)

        recorder.write_params(self.a_fit, self.b_fit)

        for _, row in self.result.round(2).iterrows():
            recorder.write_result(
                prediction=row.prediction,
                high=row.high_interval,
                low=row.low_interval,
                pval=row.p_values,
                realscore=row.score,
                is_considered=row.in_fit,
                is_ambiguous=row.ambiguous,
            )

        recorder.finalize_row()
