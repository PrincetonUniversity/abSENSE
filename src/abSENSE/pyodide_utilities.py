"""Helper methods for performing analyses with pyodide."""

from io import StringIO
from typing import TextIO
import pandas as pd

from abSENSE.analyzer import AbsenseAnalyzer
from abSENSE.parameters import AbsenseParameters
from abSENSE.plotting import FitPlot
from abSENSE.results import SampledResult

def get_plots_from_text(
        data: str,
        e_value: float,
        gene_length: float,
        db_length: float,
    ):
    """Given text input, yield plots of each fit."""
    df = pd.read_csv(StringIO(data), sep=r'[\t|,]', engine='python')

    distances = StringIO()
    df.iloc[:, 0:2].to_csv(distances, sep='\t', header=False, index=False)
    distances.seek(0)

    bitscores = StringIO()
    bits = df.iloc[:, 2:].T
    bits.columns=df.iloc[:, 0]
    bits.index.name='Gene'
    bits.to_csv(bitscores, sep='\t', header=True, index=True, na_rep='N/A')
    bitscores.seek(0)

    return get_plots_from_files(
        bitscores=bitscores,
        distances=distances,
        e_value=e_value,
        gene_length=gene_length,
        db_length=db_length,
    )

def get_plots_from_files(
        bitscores: TextIO,
        distances: TextIO,
        e_value: float,
        gene_length: float,
        db_length: float,
    ):
    """Given file inputs, yield plots of each fit."""
    params = AbsenseParameters(
        bitscores=bitscores,
        distances=distances,
        e_value=e_value,
        default_gene_length=gene_length,
        default_db_length=db_length,
        include_only=None,
        out_dir='pyodide',
        plot='',
        plot_all=True,
        predict_all=True,
        start_time='',
        validate=False,
        # these are expected to be files
        gene_lengths=None,
        db_lengths=None,
    )
    analyzer = AbsenseAnalyzer(params)
    for result in analyzer.fit_genes():
        if isinstance(result, SampledResult):
            plot = FitPlot()
            plot.generate_plot(result)
            yield plot
