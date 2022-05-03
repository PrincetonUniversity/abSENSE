from __future__ import annotations

import glob
import inspect
import math
import os
import random
import re
import sys
import warnings
from datetime import datetime

import click
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import chi2

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
    distancefile = np.transpose(np.genfromtxt(params.distances, dtype=str, delimiter="\t"))
    bitscores = np.genfromtxt(params.bitscores, dtype=str, delimiter="\t")

    speciesorder = distancefile[0]
    rawdistances = distancefile[1].astype(float)
    genelist = bitscores[1:, 0]  # skip header

    if params.gene_lengths is None:
        genelengthfilefound = False
        defgenelen = params.default_gene_length
    else:
        genelengths = np.genfromtxt(params.gene_lengths, dtype=str, delimiter="\t")
        genelengthfilefound = True

    if params.db_lengths is None:
        speciesdblengthfilefound = False
        defdbsize = params.default_db_length
        speciesdblengths = np.transpose(
            np.vstack((speciesorder, [float(defdbsize)] * len(speciesorder)))
        )
    else:
        speciesdblengths = np.genfromtxt(params.db_lengths, dtype=str, delimiter="\t")
        speciesdblengthfilefound = True


    # Determine order of species in the bitscore file, in case they're not the same as in the distance file, from which they're extracted
    # Also determine the locations in the ordering to be used of the species to omit from curve fit, if any
    ordervec = (
        []
    )  # first index is location of first species in speciesorder in the file; and so on
    if params.include_only is None:
        pred_specs = []
        pred_spec_locs = []
    else:
        pred_specs = re.split(",", params.include_only)
        pred_spec_locs = []

    for i, species in enumerate(speciesorder):
        found = False
        for j, score in enumerate(bitscores[0]):
            if species in score:
                ordervec.append(j)
                found = True
                if species in pred_specs:
                    pred_spec_locs.append(i)
        if not found:
            raise ValueError(
                "One or more species names in header of bitscore file do not "
                "match species names in header of distance file! "
                f"The first I encountered was {species}. Quitting. \n"
            )

    invordervec = (
        []
    )  # first index is location of first species in file in speciesorder; and so on
    for score in bitscores[0]:
        for j, species in enumerate(speciesorder):
            if score in species:
                invordervec.append(j)


    # Find species db sizes in the right order, from either file or manual input (done above)
    speciestotallengths = []
    for species in speciesorder:
        found = False
        for species_db in speciesdblengths:
            if species in species_db[0]:
                speciestotallengths.append(float(species_db[1]))
                found = True
        if not found:
            if speciesdblengthfilefound:
                raise ValueError(
                    "One or more species names in your database size file do not "
                    "match species names in distance file! The first I "
                    f"encountered was {species}. Quitting. \n"
                )

    with FileRecorder(
            params,
            [speciesorder[i] for i in invordervec],
    ).open() as recorder:

        # Ignore warning that sometimes happen as a result of stochastic 
        # sampling but that doesn't affect overall computation
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")

        recorder.write_info(params, pred_specs)
        recorder.write_headers()

        print("Running!")

        for i, gene in enumerate(genelist):
            # report current position and gene name
            print(f"gene {i+1} out of {len(genelist)}: {gene}")

            recorder.write_gene(gene)

            # make new arrays to put truncated (not using omitted species) scores, distances
            genebitscores = []
            truncdistances = []

            # make new arrays to indicate which species/distances are ambiguous
            # orthology but have a homolog; don't predict these later
            ambigdists = []

            # if gene length input file given, look for length of current gene
            # if not given, assume default value (specified above)
            lengthfound = False
            if genelengthfilefound:
                for z in range(1, len(genelengths)):  # skip header
                    if gene in genelengths[z]:
                        seqlen = float(genelengths[z][1])
                        lengthfound = True
                        break
                if not lengthfound:
                    raise ValueError(
                        f"Gene {gene} not found in specified " "gene length file! Quitting \n"
                    )
            else:
                seqlen = float(defgenelen)

            # put scores for current gene in bitscore file in right order
            # i + 1 because header skipped in gene list formation, so one behind now
            orderedscores = [bitscores[i+1][index] for index in ordervec]

            # append score of species and corresponding distance of species to gene-specific distance, score vectors if:
            # score isn't 0 (can't distinguish absence from loss from bad data etc)
            # score isn't 'N/A' or some other string (ie indicator of unclear orthology or otherwise absent, or generally unclear what it is)
            # score isn't from species that is to be excluded from fit
            for k, score in enumerate(orderedscores):
                if len(pred_spec_locs) > 0:
                    if k in pred_spec_locs and score != "0" and isfloat(score):
                        genebitscores.append(float(score))
                        truncdistances.append(rawdistances[k])
                    elif score == "N/A":
                        ambigdists.append(rawdistances[k])
                else:
                    if score != "0" and isfloat(score):
                        genebitscores.append(float(score))
                        truncdistances.append(rawdistances[k])
                    elif score == "N/A":
                        ambigdists.append(rawdistances[k])

            if len(truncdistances) <= 2:
                recorder.not_enough_data()
                continue

            try:
                (a, b), covar = curve_fit(
                    func,
                    truncdistances,
                    genebitscores,
                    bounds=((-np.inf, 0), (np.inf, np.inf)),
                )

            except RuntimeError:
                recorder.analysis_error()
                continue

            parout = parameter_CI_find(a, b, covar)
            if parout == "failed":
                recorder.analysis_error(predictions=[func(dist, a, b) for dist in rawdistances])
                continue

            testavals, testbvals = parout
            for j in range(0, len(rawdistances)):
                bitthresh = -1 * math.log(
                    params.e_value / (seqlen * speciestotallengths[invordervec[j]]), 2
                )

                prediction = round(func(rawdistances[invordervec[j]], a, b), 2)
                lowprediction, highprediction, pval = PI_find(
                    testavals, testbvals, rawdistances[invordervec[j]], bitthresh
                )

                is_considered = rawdistances[invordervec[j]] in truncdistances
                realscore = 0
                if is_considered:
                    realscore = genebitscores[
                            truncdistances.index(rawdistances[invordervec[j]])
                        ]

                recorder.write_result(
                    prediction=prediction,
                    high=highprediction,
                    low=lowprediction,
                    pval=pval,
                    realscore=realscore,
                    is_considered=is_considered,
                    is_ambiguous=rawdistances[invordervec[j]] in ambigdists,
                )

            recorder.write_params(a, b)
            recorder.finalize_row()


# curve to fit
def func(x, a, b):
    return a * np.exp(-b * x)


def isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def parameter_CI_find(mla, mlb, covar):
    """
    function to, where possible, use maximum likelihood estimates of a and b
    parameter plus estimated covariance matrix to directly sample from the
    probability distribution of a and b (assume Gaussian with mean of max
    likelihood estimates and given covariance structure)
    """
    if True not in np.isinf(covar):

        testavals = []
        testbvals = []
        for i in range(0, 200):
            a, b = np.random.multivariate_normal([mla, mlb], covar)
            testavals.append(a)
            testbvals.append(b)

        return testavals, testbvals
    else:
        return "failed"


def PI_find(testavals, testbvals, currx, bitthresh):
    """Gives an empirical estimate of the prediction interval.
    function to take each of the sampled a, b values and use them to sample
    directly from the distribution of scores taking into account the Gaussian
    noise (a function of distance, a, b)."""

    # sample from score distribution: Gaussian with mean a, b and noise determined by distance (currx), a, b
    PIsamples = []
    for test_a_val, test_b_val in zip(testavals, testbvals):
        detval = func(currx, test_a_val, test_b_val)
        estnoise = np.sqrt(
            test_a_val
            * (1 - math.exp(-1 * test_b_val * currx))
            * (math.exp(-1 * test_b_val * currx))
        )
        if estnoise > 0:
            for _ in range(0, 200):
                PIsamples.append(detval + np.random.normal(0, estnoise))
        else:
            PIsamples.append(detval)
    # compute mean of sample
    mean = np.mean(PIsamples)
    # compute std dev of sample
    std = np.std(PIsamples)

    # calculate this analytically from std estimate
    pval = stats.norm.cdf(bitthresh, mean, std)

    # calculate 99% CI
    (lowint, highint) = stats.norm.interval(0.99, mean, std)

    return lowint, highint, pval


if __name__ == '__main__':
    main()
