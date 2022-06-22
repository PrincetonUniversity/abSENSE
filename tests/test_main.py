from __future__ import annotations

from io import StringIO

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import abSENSE.__main__ as main


def test_defaults(faked_recorder_non_closing, default_params):
    _, files = faked_recorder_non_closing
    with pytest.warns(UserWarning) as record:
        main.perform_analysis(default_params)

    assert len(record) == 1
    assert record[0].message.args[0] == (
        "Only one estimate of distances provided, "
        "additional values can improve result accuracy"
    )

    files["./parameters.tsv"].seek(0)
    result = pd.read_csv(files["./parameters.tsv"], delimiter="\t", comment="#")
    expected = pd.read_csv(
        StringIO(
            "# the best-fit parameters (performed using only bitscores from species not omitted from the fit; see run info file) for a and b for each gene\n"
            "Gene\ta\tb\n"
            "NP_010181.2\t2359.910340705106\t0.676126045968598\n"
            "NP_009362.1\t1053.5985693169105\t0.4073235374987012\n"
            "NP_014555.1\t732.8505436789137\t0.5595110997139245\n"
            "NP_116682.3\t352.34734049379813\t14.460626046962682\n"
            "NP_011284.1\t303.0830341558144\t5.143646534570729\n"
            "NP_011878.1\t589.7752939170018\t2.729326889253197\n"
            "NP_013320.1\t1552.4815439696647\t5.120090782287544\n"
            "NP_014160.2\t968.5487558750334\t0.9669086263255712\n"
            "NP_014890.1\t383.70820281955355\t4.700780867702744\n"
        ),
        delimiter="\t",
        comment="#",
    )
    assert_frame_equal(result, expected)

    assert files["./predicted_bitscores.tsv"].getvalue() == (
        "# maximum likelihood bitscore predictions for each tested gene in each species\n"
        "Gene\tS_cer\tS_par\tS_mik\tS_kud\tS_bay\tS_castellii\tK_waltii\tA_gossyppi\tK_lactis\tA_nidulans\tS_pombe\tY_lipolytica\n"
        "NP_010181.2\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\n"
        "NP_009362.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\n"
        "NP_014555.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\n"
        "NP_116682.3\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t73.91\t1.85\t0.28\t0.2\t0.11\t0.0\t0.0\t0.0\n"
        "NP_011284.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t46.85\t23.88\t21.11\t17.27\t2.91\t2.64\t2.24\n"
        "NP_011878.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t218.98\t153.16\t143.44\t128.96\t50.16\t47.62\t43.64\n"
        "NP_013320.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t15.24\t13.83\t11.74\n"
        "NP_014160.2\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t600.73\tOrtholog_detected\t565.23\tOrtholog_detected\t397.15\tOrtholog_detected\n"
        "NP_014890.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t37.63\t33.61\t27.98\t5.5\t5.03\t4.33\n"
    )

    assert files["./failure_probabilities.tsv"].getvalue() == (
        "# the probability of a homolog being undetected at the specified significance threshold (see run info file) in each tested gene in each species\n"
        "Gene\tS_cer\tS_par\tS_mik\tS_kud\tS_bay\tS_castellii\tK_waltii\tA_gossyppi\tK_lactis\tA_nidulans\tS_pombe\tY_lipolytica\n"
        "NP_010181.2\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\n"
        "NP_009362.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\n"
        "NP_014555.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\tOrtholog_detected\tHomolog_detected(orthology_ambiguous)\n"
        "NP_116682.3\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t0.0\t1.0\t1.0\t1.0\t1.0\t1.0\t1.0\t1.0\n"
        "NP_011284.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t0.29\t0.98\t0.99\t1.0\t1.0\t1.0\t1.0\n"
        "NP_011878.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t0.0\t0.0\t0.0\t0.0\t0.23\t0.28\t0.39\n"
        "NP_013320.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t1.0\t1.0\t1.0\n"
        "NP_014160.2\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t0.0\tOrtholog_detected\t0.0\tOrtholog_detected\t0.0\tOrtholog_detected\n"
        "NP_014890.1\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\tOrtholog_detected\t0.55\t0.65\t0.8\t1.0\t1.0\t1.0\n"
    )

    assert files["./run_info.txt"].getvalue() == (
        "abSENSE analysis run on NOW\n"
        "Input bitscore file: bitscores\n"
        "Input distance file: distances\n"
        "Gene length (for all genes): 400 (default)\n"
        "Database length (for all species): 8000000 (default)\n"
        "Species used in fit: All (default)\n"
        "E-value threshold: 0.001\n"
    )
