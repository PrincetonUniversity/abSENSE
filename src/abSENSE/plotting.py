import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import scipy.interpolate

from abSENSE.constants import exponential, ORANGE, GREY, sample_parameters, find_confidence_interval


class FitPlot():
    def __init__(self):
        self.figure, self.axes = plt.subplots(figsize=(11.5, 6.0), tight_layout=True)
        self.interval = 0.01

    def title(self, gene, a_fit, b_fit, correlation):
        self.axes.set_title(
            (
                f'Gene: {gene}\n'
                f'a = {round(a_fit, 1)}, b = {round(b_fit, 2)}\n'
                f'$r^2$ = {round(correlation**2, 2)}'
            ),
            color="black",
            fontsize=13,
            fontweight="bold",
        )

    def scores(self, distance, score):
        self.axes.scatter(
            x=distance,
            y=score,
            s=40,
            c="black",
            label="Bitscores of detected orthologs used in fit",
        )

    def fit(self, distance, a_fit, b_fit, covariance, bit_threshold):
        distance, high, low = self._interpolate(distance, a_fit, b_fit, covariance)

        prediction = exponential(distance, a_fit, b_fit)
        self.axes.plot(distance, prediction, color='red', label='Predicted bitscore')
        self.axes.plot(distance, high, color='black')
        self.axes.plot(distance, low, color='black')

        self.axes.fill_between(
            distance,
            high,
            low,
            facecolor='blue',
            alpha=0.2,
            label='99% confidence interval',
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

    def _interpolate(self, distance, a_fit, b_fit, covariance):
        distance = np.linspace(distance.min(), distance.max(), num=100, endpoint=True)
        random = np.random.default_rng()
        result = find_confidence_interval(
            random,
            distance.reshape(-1, 1),
            sample_parameters(random, a_fit, b_fit, covariance),
        )

        return distance, result.high_interval, result.low_interval

    def set_axes(self, distance, species):
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

        self.axes.spines['right'].set_visible(False)
        self.axes.spines['top'].set_visible(False)
        self.axes.tick_params(axis="x", width=2, length=7, direction="inout")
        self.axes.set_xticks(distance, species, fontsize=10, rotation=90)

        self.axes.tick_params(axis='y', labelsize=10)

    def highlight_not_in_fit(self, distance, in_fit):
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

    def highlight_ambiguous(self, ambiguous):
        for i, is_ambiguous in enumerate(ambiguous):
            if is_ambiguous:
                self.axes.get_xticklabels()[i].set_color(GREY)

    def show_label(self, any_not_in_fit, any_ambiguous):
        handles, labels = self.axes.get_legend_handles_labels()

        if any_ambiguous:
            handles.append(Patch(facecolor=GREY, alpha=0.3))
            labels.append('Homolog detected, but orthology ambiguous')

        if any_not_in_fit:
            handles.append(Patch(facecolor=ORANGE, alpha=0.3))
            labels.append('No homolog detected!')

        self.axes.legend(handles, labels, fontsize=9)

    def save(self, filename):
        self.figure.savefig(filename, format='svg')

    def show(self):
        plt.show()
