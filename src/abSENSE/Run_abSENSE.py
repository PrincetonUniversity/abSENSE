from __future__ import annotations

import argparse
import glob
import inspect
import math
import os
import random
import re
import sys
import warnings
from datetime import datetime

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import chi2

from recorder import File_Recorder


def main():
    now = datetime.now()
    starttime = now.strftime("%m.%d.%Y_%H.%M")

    parser = argparse.ArgumentParser(description="abSENSE arguments:")

    parser.add_argument(
        "--distfile",
        type=str,
        required=True,
        help="Required. Name of file containing pairwise evolutionary distances between focal species and each of the other species",
    )
    parser.add_argument(
        "--scorefile",
        type=str,
        required=True,
        help="Required. Name of file containing bitscores between focal species gene and orthologs in other species",
    )
    parser.add_argument(
        "--Eval",
        default=0.001,
        type=float,
        help="Optional. E-value threshold. Scientific notation (e.g. 10E-5) accepted. Default 0.001.",
    )
    parser.add_argument(
        "--includeonly",
        type=str,
        help="Optional. Species whose orthologs' bitscores will be included in fit; all others will be omitted. Default is all species. Format as species names, exactly as in input files, separated by commas (no spaces).",
    )
    parser.add_argument(
        "--genelenfile",
        type=str,
        help="Optional. File containing lengths (aa) of all genes to be analyzed. Used to accurately calculate E-value threshold. Default is 400aa for all genes. Only large deviations will qualitatively affect results.",
    )
    parser.add_argument(
        "--dblenfile",
        type=str,
        help="Optional. File containing size (aa) of databases on which the anticipated homology searches will be performed. Species-specific. Used to accurately calculate E-value threshold. Default is 400aa/gene * 20,000 genes for each species, intended to be the size of an average proteome. Only large deviations will significantly affect results.",
    )
    parser.add_argument(
        "--predall",
        type=bool,
        default=False,
        help="Optional. True: Predicts bitscores and P(detectable) of homologs in all species, including those in which homologs were actually detected. Default is False: only make predictions for homologs that seem to be absent.",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="abSENSE_results_" + starttime,
        help="Optional. Name of directory for output data. Default is date and time when analysis was run.",
    )
    args = parser.parse_args()

    distancefilecheck = glob.glob(args.distfile)
    if len(distancefilecheck) == 0:
        sys.exit(
            "Distance file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n"
        )
    else:
        distancefile = np.transpose(np.genfromtxt(args.distfile, dtype=str, delimiter="\t"))

    scorefilecheck = glob.glob(args.scorefile)
    if len(scorefilecheck) == 0:
        sys.exit(
            "Bitscore file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n"
        )
    else:
        bitscores = np.genfromtxt(args.scorefile, dtype=str, delimiter="\t")

    speciesorder = distancefile[0]
    rawdistances = distancefile[1].astype(float)
    genelist = bitscores[1:, 0]  # skip header

    ethresh = args.Eval

    if args.genelenfile is None:
        genelengthfilefound = False
        defgenelen = float(400)
    else:
        genelengthfilefound = True

    if args.dblenfile is None:
        speciesdblengthfilefound = False
        defdbsize = float(8000000)
        speciesdblengths = np.transpose(
            np.vstack((speciesorder, [float(defdbsize)] * len(speciesorder)))
        )
    else:
        speciesdblengthfilefound = True

    if genelengthfilefound:
        genelengthfilecheck = glob.glob(args.genelenfile)
        if len(genelengthfilecheck) != 1:
            sys.exit(
                "Gene length file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n"
            )
        else:
            genelengths = np.genfromtxt(args.genelenfile, dtype=str, delimiter="\t")

    if speciesdblengthfilefound:
        speciesdblengthfilecheck = glob.glob(args.dblenfile)
        if len(speciesdblengthfilecheck) != 1:
            sys.exit(
                "Species database size file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n"
            )
        else:
            speciesdblengths = np.genfromtxt(args.dblenfile, dtype=str, delimiter="\t")

    # Determine order of species in the bitscore file, in case they're not the same as in the distance file, from which they're extracted
    # Also determine the locations in the ordering to be used of the species to omit from curve fit, if any
    ordervec = (
        []
    )  # first index is location of first species in speciesorder in the file; and so on
    if args.includeonly is None:
        pred_specs = []
        pred_spec_locs = []
    else:
        pred_specs = re.split(",", args.includeonly)
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
            sys.exit(
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
                sys.exit(
                    "One or more species names in your database size file do not "
                    "match species names in distance file! The first I "
                    f"encountered was {species}. Quitting. \n"
                )

    with File_Recorder(
            args.out,
            [speciesorder[i] for i in invordervec],
            args.predall,
    ).open() as recorder:

        # Ignore warning that sometimes happen as a result of stochastic 
        # sampling but that doesn't affect overall computation
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")

        recorder.write_info(args, defgenelen, defdbsize, pred_specs)
        recorder.write_headers()

        print("Running!")

        for i, gene in enumerate(genelist):
            # report current position and gene name
            print(f"gene {i} out of {len(genelist)}: {gene}")

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
                    sys.exit(
                        f"Gene {gene} not found in specified " "gene length file! Quitting \n"
                    )
            else:
                seqlen = float(defgenelen)

            # put scores for current gene in bitscore file in right order
            orderedscores = []
            for index in ordervec:  # ordervec starts at 1
                # i + 1 because header skipped in gene list formation, so one behind now
                orderedscores.append(bitscores[i + 1][index])

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
                    ethresh / (seqlen * speciestotallengths[invordervec[j]]), 2
                )

                prediction = round(func(rawdistances[invordervec[j]], a, b), 2)
                lowprediction, highprediction, pval = PI_find(
                    testavals, testbvals, rawdistances[invordervec[j]], bitthresh
                )

                is_truncated = rawdistances[invordervec[j]] in truncdistances
                realscore = 0
                if is_truncated:
                    realscore = genebitscores[
                            truncdistances.index(rawdistances[invordervec[j]])
                        ]

                recorder.write_result(
                    prediction=prediction,
                    high=highprediction,
                    low=lowprediction,
                    pval=pval,
                    realscore=realscore,
                    is_truncated=is_truncated,
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
