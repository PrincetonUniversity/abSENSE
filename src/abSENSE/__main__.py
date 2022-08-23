"""Main CLI entry point."""
from __future__ import annotations

from datetime import datetime
from typing import Any

import click

from abSENSE.analyzer import AbsenseAnalyzer
from abSENSE.parameters import AbsenseParameters
from abSENSE.recorder import Recorder


@click.command(help="A method to interpret undetected homologs")
@click.option(
    "--distances",
    type=click.File("rt"),
    required=True,
    help="tsv file containing pairwise evolutionary distances "
    "between focal species and each of the other species.",
)
@click.option(
    "--bitscores",
    type=click.File("rt"),
    required=True,
    help="tsv file containing bitscores between focal "
    "species gene and orthologs in other species.",
)
@click.option("--e-value", default=0.001, help="E-value threshold. Default 0.001.")
@click.option(
    "--include-only",
    type=str,
    help="Comma separated list of species to include in analysis.",
)
@click.option(
    "--gene-lengths",
    type=click.File("rt"),
    help="tsv file containing amino acid lengths of all genes in "
    "bitscores. Used to accurately calculate E-value threshold. "
    "Default is 400aa for all genes. "
    "Only large deviations will qualitatively affect results.",
)
@click.option(
    "--db-lengths",
    type=click.File("rt"),
    help="tsv file containing amino acid size of each species' database."
    "Used to accurately calculate E-value threshold. "
    "Default is 400aa/gene * 20,000 genes = 8000000 for all species, "
    "intended to be the size of an average proteome. "
    "Only large deviations will significantly affect results.",
)
@click.option(
    "--predict-all",
    is_flag=True,
    help="If set, predicts bitscores and P(detectable) of homologs "
    "in all species, including those in which homologs were "
    "detected. By default, will make predictions only for homologs "
    "that seem to be absent.",
)
@click.option(
    "--plot-all",
    is_flag=True,
    help="If set, will plot fitting results for each gene present in "
    "the bitscores file. By default, will not produce plots.",
)
@click.option(
    "--plot",
    default="",
    help="A comma separated list of genes to plot",
)
@click.option(
    "--validate",
    is_flag=True,
    help="If set, will remove each gene in the bitscore file and tabulate "
    "statistics on how well the current distances fit the data.",
)
@click.option(
    "--out-dir",
    type=click.Path(file_okay=False),
    help="Name of output directory. " "Default is date and time when analysis was run.",
)
def main(**args: Any) -> None:
    """CLI wrapper."""
    now = datetime.now()
    start_time = now.strftime("%m.%d.%Y_%H.%M")

    params = AbsenseParameters(**args, start_time=start_time)

    perform_analysis(params)


def perform_analysis(params: AbsenseParameters) -> None:
    """Wrapper for unit testing."""
    analyzer = AbsenseAnalyzer(params)
    total_genes = analyzer.total_genes()

    with Recorder.build_recorder(params, analyzer.species).open() as recorder:
        recorder.write_info(params)
        recorder.write_headers()

        click.echo("Running!")

        for i, result in enumerate(analyzer.fit_genes()):
            status = result.status
            if status != "":
                status = "\t-- " + status
            click.echo(f"gene {i+1} out of {total_genes}: {result.gene}{status}")

            result.record_to(recorder)


if __name__ == "__main__":
    main()
