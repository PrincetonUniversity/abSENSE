"""Helper methods for performing analyses with pyodide."""

import base64
from io import BytesIO, StringIO
from typing import TextIO
import pandas as pd

from abSENSE.analyzer import AbsenseAnalyzer
from abSENSE.parameters import AbsenseParameters
from abSENSE.plotting import FitPlot
from abSENSE.results import SampledResult, ErrorResult, NotEnoughDataResult

def get_analyzer_from_text(
        data: str,
        e_value: float,
        gene_length: float,
        db_length: float,
    ):
    """Given text input, return analyzer."""
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

    return get_analyzer_from_files(
        bitscores=bitscores,
        distances=distances,
        e_value=e_value,
        gene_length=gene_length,
        db_length=db_length,
    )

def get_analyzer_from_files(
        bitscores: TextIO,
        distances: TextIO,
        e_value: float,
        gene_length: float,
        db_length: float,
    ):
    """Given file inputs, return analyzer."""
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
    return AbsenseAnalyzer(params)

def get_plot_and_table(gene: str, analyzer: AbsenseAnalyzer):
    bitscore = analyzer.bitscores.loc[gene]
    result = analyzer.fit_gene(gene, bitscore)
    plot = FitPlot()
    missing_homologs = []
    if isinstance(result, SampledResult):
        plot.generate_plot(result)
        for species, row in result.result[result.result.score == 0].iterrows():
            missing_homologs.append(
                {
                    'species': species,
                    'p_value': f"{row.p_values:0.2f}",
                    'bitscore': f"{row.prediction:0.2f}",
                    'interval': f"{row.low_interval:0.2f} - {row.high_interval:0.2f}",
                }
            )
    else:
        if isinstance(result, ErrorResult):
            plot.write_error(result.gene, "Analysis error.")
        elif isinstance(result, NotEnoughDataResult):
            plot.write_error(result.gene, "Not enough data.")

    return plot, missing_homologs


def generate_png_string(plot: FitPlot):
    buf = BytesIO()
    plot.save(buf, format='png')
    buf.seek(0)
    base64_str = base64.b64encode(buf.read()).decode('UTF-8')
    return f"data:image/png;base64,{base64_str}"
