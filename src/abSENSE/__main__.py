from __future__ import annotations

from datetime import datetime

import click
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import chi2

from abSENSE.analyzer import AbsenseAnalyzer
from abSENSE.parameters import AbsenseParameters
from abSENSE.recorder import FileRecorder


@click.command(help="A method to interpred undetected homologs")
@click.option('--distances', type=click.File('rt'), required=True,
              help="tsv file containing pairwise evolutionary distances "
              "between focal species and each of the other species.")
@click.option('--bitscores', type=click.File('rt'), required=True,
              help="tsv file containing bitscores between focal "
              "species gene and orthologs in other species.")
@click.option('--e-value', default=0.001,
              help="E-value threshold. Default 0.001.")
@click.option('--include-only', type=str,
              help="Comma separated list of species to include in analysis.")
@click.option('--gene-lengths', type=click.File('rt'),
              help="tsv file containing amino acid lengths of all genes in "
              "bitscores. Used to accurately calculate E-value threshold. "
              "Default is 400aa for all genes. "
              "Only large deviations will qualitatively affect results.")
@click.option('--db-lengths', type=click.File('rt'),
              help="tsv file containing amino acid size of each species' database."
              "Used to accurately calculate E-value threshold. "
              "Default is 400aa/gene * 20,000 genes = 8000000 for all species, "
              "intended to be the size of an average proteome. "
              "Only large deviations will significantly affect results.",
              )
@click.option('--predict-all', is_flag=True,
              help="If set, predicts bitscores and P(detectable) of homologs "
              "in all species, including those in which homologs were "
              "detected. By default, will make predictions only for homologs "
              "that seem to be absent.")
@click.option('--out-dir', type=click.Path(file_okay=False),
              help="Name of output directory. "
              "Default is date and time when analysis was run.")
def main(**args):
    """CLI wrapper."""
    now = datetime.now()
    start_time = now.strftime("%m.%d.%Y_%H.%M")

    params = AbsenseParameters(**args, start_time=start_time)

    perform_analysis(params)


def perform_analysis(params: AbsenseParameters):
    """Wrapper for unit testing."""
    analyzer = AbsenseAnalyzer(params)
    total_genes = analyzer.total_genes()

    with FileRecorder(params, analyzer.species).open() as recorder:
        recorder.write_info(params)
        recorder.write_headers()

        print("Running!")

        for i, result in enumerate(analyzer.fit_genes()):
            status = result.status()
            if status != '':
                status = '\t-- ' + status
            print(f"gene {i+1} out of {total_genes}: {result.gene}{status}")

            result.record_to(recorder)


if __name__ == '__main__':
    main()
