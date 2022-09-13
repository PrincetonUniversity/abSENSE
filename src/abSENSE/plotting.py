"""Module for producing plots for abSENSE results."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.patches import Patch

from abSENSE.constants import GREY, ORANGE
from abSENSE.utilities import exponential, find_confidence_interval, sample_parameters


class FitPlot:
    """Wrapper of matplotlib to specialize for abSENSEE."""

    def __init__(self) -> None:
        """Initialize new plot."""
        self.figure, self.axes = plt.subplots(figsize=(11.5, 6.0), tight_layout=True)
        self.interval = 0.01

    def title(self, gene: str, a_fit: float, b_fit: float, correlation: float) -> None:
        """Place title on plot.

        Args:
            gene: gene name
            a_fit: a fit value in exponential function
            b_fit: b fit value in exponential function
            correlation: goodness of fit
        """
        self.axes.set_title(
            (
                f"Gene: {gene}\n"
                f"a = {round(a_fit, 1)}, b = {round(b_fit, 2)}\n"
                f"$r^2$ = {round(correlation**2, 2)}"
            ),
            color="black",
            fontsize=13,
            fontweight="bold",
        )

    def scores(
        self,
        distance: list[float],
        score: list[float],
        xerr: list[float],
    ) -> None:
        """Plot score values.

        Args:
            distance: x values
            score: y values
            xerr: uncertainty in x values, None to skip plotting errorbar
        """
        # passing in na values raises a runtime error for errorbar
        self.axes.errorbar(
            distance,
            np.zeros_like(distance),
            xerr=xerr,
            ecolor="red",
            fmt="none",
            elinewidth=4,
        )
        self.axes.scatter(
            x=distance,
            y=score,
            s=40,
            c="black",
            label="Bitscores of detected orthologs used in fit",
        )

    def fit(
        self,
        distance: npt.NDArray[np.float64],
        a_fit: float,
        b_fit: float,
        covariance: npt.NDArray[np.float64],
        bit_threshold: float,
    ) -> None:
        """Plot fit and confidence interval.

        Args:
            distance: target x values
            a_fit: a fit value of exponential
            b_fit: b fit value of exponential
            covariance: covariance matrix of a/b fit
            bit_threshold: bitscore threshold of dataset
        """
        distance, high, low = self._interpolate(distance, a_fit, b_fit, covariance)

        prediction = exponential(distance, a_fit, b_fit)
        self.axes.plot(distance, prediction, color="red", label="Predicted bitscore")
        self.axes.plot(distance, high, color="black")
        self.axes.plot(distance, low, color="black")

        self.axes.fill_between(
            distance,
            high,
            low,
            facecolor="blue",
            alpha=0.2,
            label="99% confidence interval",
        )

        self.interval = (max(distance) - min(distance)) / 100
        self.axes.set_xlim([-self.interval, max(distance) + self.interval])
        self.axes.set_ylim([0, max(prediction) * 1.1])
        self.axes.axhline(
            y=bit_threshold,
            linestyle="dashed",
            c="black",
            label="Detectability threshold",
        )

    @staticmethod
    def _interpolate(
        distance: npt.NDArray[np.float64],
        a_fit: float,
        b_fit: float,
        covariance: npt.NDArray[np.float64],
    ) -> tuple[npt.NDArray[np.float64], pd.Series[float], pd.Series[float]]:
        """Interpolate fits for smooth plotting.

        Args:
            distance: target x values
            a_fit: a fit value of exponential
            b_fit: b fit value of exponential
            covariance: covariance matrix of a/b fit

        Returns:
            new distance with 100 points
            high values of estimate
            low values of estimate
        """
        distance = np.linspace(distance.min(), distance.max(), num=100, endpoint=True)
        random = np.random.default_rng()
        result = find_confidence_interval(
            random,
            distance.reshape(-1, 1),
            sample_parameters(random, a_fit, b_fit, covariance),
        )

        return distance, result.high_interval, result.low_interval

    def set_axes(self, distance: list[float], species: pd.Index) -> None:
        """Set axes labels and ticks.

        Args:
            distance: target x values
            species: species names
        """
        self.axes.set_ylabel(
            "Bitscore",
            fontsize=13,
            labelpad=10,
        )

        self.axes.set_xlabel(
            "Evolutionary distance from focal species",
            fontsize=13,
            labelpad=10,
        )

        self.axes.spines["right"].set_visible(False)
        self.axes.spines["top"].set_visible(False)
        self.axes.tick_params(axis="x", width=2, length=7, direction="inout")
        self.axes.set_xticks(distance, species, fontsize=10, rotation=90)

        self.axes.tick_params(axis="y", labelsize=10)

    def highlight_not_in_fit(self, distance: list[float], in_fit: list[bool]) -> None:
        """Highlight species not included in fit.

        Args:
            distance: target x values
            in_fit: specify if each species was included in fit
        """
        for i, is_in_fit in enumerate(in_fit):
            if not is_in_fit:
                self.axes.axvspan(
                    distance[i] - self.interval,
                    distance[i] + self.interval,
                    facecolor=ORANGE,
                    alpha=0.3,
                    capstyle="round",
                )
                self.axes.get_xticklabels()[i].set_color(ORANGE)
                self.axes.get_xticklabels()[i].set_weight("bold")

    def highlight_ambiguous(self, ambiguous: list[bool]) -> None:
        """Mark species whose homology was ambiguous.

        Args:
            ambiguous: list of bools specifying if the species was ambiguous
        """
        for i, is_ambiguous in enumerate(ambiguous):
            if is_ambiguous:
                self.axes.get_xticklabels()[i].set_color(GREY)

    def show_label(self, any_not_in_fit: bool, any_ambiguous: bool) -> None:
        """Produce legend for the plot.

        Args:
            any_not_in_fit: true if any species were not present in fit
            any_ambiguous: true if any species were ambiguous
        """
        handles, labels = self.axes.get_legend_handles_labels()

        if any_ambiguous:
            handles.append(Patch(facecolor=GREY, alpha=0.3))
            labels.append("Homolog detected, but orthology ambiguous")

        if any_not_in_fit:
            handles.append(Patch(facecolor=ORANGE, alpha=0.3))
            labels.append("No homolog detected!")

        self.axes.legend(handles, labels, fontsize=9)

    def save(self, filename: str, format: str = "svg") -> None:
        """Save this figure to file.

        Args:
            filename: the file to save to
        """
        self.figure.savefig(filename, format=format)

    @staticmethod
    def show() -> None:
        """Show the plot."""
        plt.show()

    def generate_plot(self, result: SampledResult):
        self.title(result.gene, result.a_fit, result.b_fit, result.correlation)

        fit = result.result

        self.scores(fit.distance, fit.score, fit.dist_stdev)

        self.fit(
            fit.distance,
            a_fit=result.a_fit,
            b_fit=result.b_fit,
            covariance=result.covariance,
            bit_threshold=result.bit_threshold,
        )

        self.set_axes(fit.distance, fit.index)
        # test if the species is also ambiguous, that takes precedence
        self.highlight_not_in_fit(fit.distance, (fit.in_fit | fit.ambiguous))
        self.highlight_ambiguous(fit.ambiguous)
        self.show_label(
            any_not_in_fit=not (fit.in_fit | fit.ambiguous).all(),
            any_ambiguous=fit.ambiguous.any(),
        )
