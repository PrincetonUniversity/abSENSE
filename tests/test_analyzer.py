from __future__ import annotations

from io import StringIO

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from abSENSE.analyzer import AbsenseAnalyzer
from abSENSE.exceptions import MissingGeneException, MissingSpeciesException
from abSENSE.recorder import Recorder
from abSENSE.results import (
    ErrorResult,
    NotEnoughDataResult,
    SampledResult,
    ValidationResult,
)
from abSENSE.utilities import find_confidence_interval, sample_parameters


def test_init_defaults(default_params):
    with pytest.warns(UserWarning) as record:
        analyzer = AbsenseAnalyzer(default_params)

    assert analyzer.odr_model is None

    assert len(record) == 1
    assert record[0].message.args[0] == (
        "Only one estimate of distances provided, "
        "additional values can improve result accuracy"
    )

    assert (analyzer.gene_lengths.index == analyzer.genes).all()
    assert (analyzer.gene_lengths["length"] == default_params.default_gene_length).all()

    assert (analyzer.db_lengths.index == analyzer.species).all()
    assert (analyzer.db_lengths["length"] == default_params.default_db_length).all()


def test_init_defaults_zero_variance(default_params):
    default_params.distances = StringIO(
        "S_cer\t0\t0\t0\t0\n"
        "S_par\t0.051\t0.051\t0.051\t0.051\n"  # zero variance
        "S_mik\t0.088\t0.098\t0.078\t0.098\n"
        "S_kud\t0.104\t0.114\t0.094\t0.114\n"
        "S_bay\t0.108\t0.118\t0.098\t0.107\n"
        "S_castellii\t0.363\t0.353\t0.373\t0.366\n"
        "K_waltii\t0.494\t0.484\t0.504\t0.496\n"
        "A_gossyppi\t0.518\t0.508\t0.528\t0.520\n"
        "K_lactis\t0.557\t0.547\t0.567\t0.559\n"
        "A_nidulans\t0.903\t0.893\t0.913\t0.905\n"
        "S_pombe\t0.922\t0.912\t0.932\t0.924\n"
        "Y_lipolytica\t0.954\t0.944\t0.964\t0.956\n"
    )
    with pytest.warns(UserWarning) as record:
        analyzer = AbsenseAnalyzer(default_params)

    assert analyzer.odr_model is None

    assert len(record) == 1
    assert record[0].message.args[0] == (
        "One or more estimates have 0 " "variance, skipping ODR"
    )

    assert (analyzer.gene_lengths.index == analyzer.genes).all()
    assert (analyzer.gene_lengths["length"] == default_params.default_gene_length).all()

    assert (analyzer.db_lengths.index == analyzer.species).all()
    assert (analyzer.db_lengths["length"] == default_params.default_db_length).all()


@pytest.mark.filterwarnings("ignore: Only one estimate of distances")
def test_init_with_include_only(default_params):
    # include only (or the distance index on default)
    # determines the order of bitscore columns
    default_params.include_only = "S_cer,S_mik,S_par"
    analyzer = AbsenseAnalyzer(default_params)
    assert (analyzer.species == "S_cer,S_mik,S_par".split(",")).all()
    assert (analyzer.bitscores.columns == "S_cer,S_mik,S_par".split(",")).all()


@pytest.mark.filterwarnings("ignore: Only one estimate of distances")
def test_init_with_include_only_missing(default_params):
    default_params.include_only = "S_cer,S_par,S_mik,NOTHERE"
    with pytest.raises(MissingSpeciesException) as error:
        AbsenseAnalyzer(default_params)
    assert str(error.value) == (
        "Unable to find all requested species in bitscores. " "Missing: ['NOTHERE']"
    )


@pytest.mark.filterwarnings("ignore: Only one estimate of distances")
def test_init_with_missing_db_species(default_params):
    default_params.include_only = "S_cer,S_par,S_mik,S_bay"
    default_params.db_lengths = StringIO(
        "#Species\tDBsize\n" "S_par\t2606124\n" "S_mik\t2603233\n" "S_cer\t2615464\n"
    )
    with pytest.raises(MissingSpeciesException) as error:
        AbsenseAnalyzer(default_params)
    assert str(error.value) == (
        "Unable to find all requested species in database lengths. "
        "Missing: ['S_bay']"
    )


@pytest.mark.filterwarnings("ignore: Only one estimate of distances")
def test_init_with_missing_gene_lengths(default_params):
    default_params.gene_lengths = StringIO(
        "#GeneID\tGeneLength\n"
        "NP_001018029.1\t66\n"
        "NP_001018030.1\t360\n"
        "NP_001018031.2\t113\n"
    )
    with pytest.raises(MissingGeneException) as error:
        AbsenseAnalyzer(default_params)
    assert str(error.value).startswith(
        "Unable to find all requested genes in gene lengths. "
        "Missing: ['NP_010181.2', "
    )


@pytest.fixture()
def analyzer_with_files(fungi_database_lengths, fungi_gene_lengths, default_params):
    default_params.gene_lengths = StringIO(fungi_gene_lengths)
    default_params.db_lenghts = StringIO(fungi_database_lengths)
    with pytest.warns(UserWarning):
        return AbsenseAnalyzer(default_params)


@pytest.fixture()
def analyzer_with_replicates(
    fungi_database_lengths,
    fungi_gene_lengths,
    quick_distances_replicates,
    default_params,
):
    default_params.distances = StringIO(quick_distances_replicates)
    default_params.distances.name = "distances"
    default_params.gene_lengths = StringIO(fungi_gene_lengths)
    default_params.db_lenghts = StringIO(fungi_database_lengths)
    return AbsenseAnalyzer(default_params)


def test_fit_gene_not_enough(analyzer_with_files):
    result = analyzer_with_files.fit_gene(
        "NP_001018030.1",
        pd.Series(
            [0, 1, 2] + [np.nan] * 9,
            index=analyzer_with_files.species,
        ),
    )
    assert result == NotEnoughDataResult("NP_001018030.1")


def test_fit_gene_analysis_error(analyzer_with_files, mocker):
    mocker.patch("scipy.optimize.curve_fit", side_effect=RuntimeError)
    result = analyzer_with_files.fit_gene(
        "NP_001018030.1",
        pd.Series(
            [3, 2, 1, 0] + [np.nan] * 8,
            index=analyzer_with_files.species,
        ),
    )
    assert result == ErrorResult("NP_001018030.1")


def test_fit_gene_infinite_covariance(analyzer_with_files, mocker):
    mocker.patch(
        "scipy.optimize.curve_fit", return_value=((1, 1), [[np.inf, 0], [0, 1]])
    )
    result = analyzer_with_files.fit_gene(
        "NP_001018030.1",
        pd.Series(
            [3, 2, 1, 0] + [np.nan] * 8,
            index=analyzer_with_files.species,
        ),
    )
    assert isinstance(result, ErrorResult)
    assert result.gene == "NP_001018030.1"
    assert_series_equal(
        result.predictions, np.exp(-1 * analyzer_with_files.distances)["distance"]
    )


def test_fit_normal(analyzer_with_files):
    result = next(analyzer_with_files.fit_genes())
    assert isinstance(result, SampledResult)
    assert result.gene == "NP_010181.2"
    assert result.a_fit == pytest.approx(2359.91, abs=1e-2)
    assert result.b_fit == pytest.approx(0.67612, abs=1e-4)
    assert result.correlation == pytest.approx(-0.976, abs=1e-2)
    assert result.bit_threshold == pytest.approx(39.55, abs=1e-2)


def test_fit_normal_replicates(analyzer_with_replicates):
    assert analyzer_with_replicates.odr_model is not None
    result = next(analyzer_with_replicates.fit_genes())
    assert isinstance(result, SampledResult)
    assert result.gene == "NP_010181.2"
    assert result.a_fit == pytest.approx(2424.33, abs=1e-2)
    assert result.b_fit == pytest.approx(0.74619, abs=1e-4)
    assert result.correlation == pytest.approx(-0.976, abs=1e-2)
    assert result.bit_threshold == pytest.approx(39.55, abs=1e-2)


def test_sample_parameters(analyzer_with_files):
    result = sample_parameters(analyzer_with_files.random, 1, 1, [[1, 0], [0, 1]])
    a_vals = result[:, 0]
    b_vals = result[:, 1]
    assert np.mean(a_vals) == pytest.approx(1, abs=1e-1)
    assert np.mean(b_vals) == pytest.approx(1, abs=1e-1)
    assert result.shape == (200, 2)


def test_bit_threshold(analyzer_with_files):
    old_result = -1 * np.log2(
        analyzer_with_files.e_value / (100 * analyzer_with_files.db_lengths)
    )
    old_result = old_result.rename(columns={"length": "bit_threshold"})
    result = analyzer_with_files.bit_threshold(100)
    assert_frame_equal(result, old_result)

    old_result = -1 * np.log2(
        analyzer_with_files.e_value / (234 * analyzer_with_files.db_lengths)
    )
    old_result = old_result.rename(columns={"length": "bit_threshold"})
    result = analyzer_with_files.bit_threshold(234)
    assert_frame_equal(result, old_result)


def test_find_confidence_interval(analyzer_with_files):
    bit_threshold = pd.DataFrame(
        {"bit_threshold": 40.1},
        index=analyzer_with_files.species,
    )
    result = find_confidence_interval(
        analyzer_with_files.random,
        np.expand_dims(analyzer_with_files.distances["distance"].to_numpy(), axis=1),
        np.array(
            [
                [40.1, 1],
                [40.2, 2],
                [40.3, 3],
                [40.4, -1],
            ]
        ),
        bit_threshold,
        analyzer_with_files.species,
    )

    cer = result.iloc[0, :]
    assert cer["p_values"] == pytest.approx(0.089856, abs=1e-4)
    assert cer["low_interval"] == pytest.approx(39.962014, abs=1e-4)
    assert cer["high_interval"] == pytest.approx(40.537986, abs=1e-4)


def test_validation(mocker, default_params, analyzer_with_files):
    analyzer_with_files.validate = True
    default_params.validate = True
    default_params.gene_lengths.name = "GENE"

    mocker.patch("abSENSE.recorder.os.makedirs")
    default_params.out_dir = "."
    recorder = Recorder.build_recorder(default_params, analyzer_with_files.species)
    # replace files with stringIO for easier testing
    files = {}

    def new_stringIO(*args, **kwargs):
        result = StringIO()
        result.close = lambda: None
        result.name = args[0]
        files[args[0]] = result
        return result

    mocker.patch("abSENSE.recorder.open", side_effect=new_stringIO)

    with recorder.open() as rec:
        rec.write_info(default_params)
        for result in analyzer_with_files.fit_genes():
            assert isinstance(result, ValidationResult)
            result.record_to(rec)

    assert files["./run_info.txt"].getvalue().split("\n")[-2] == (
        "Run in validation mode"
    )

    files["./validation.tsv"].seek(0)
    result = pd.read_csv(files["./validation.tsv"], sep="\t", comment="#", index_col=0)
    assert (result["mre"].isna()).sum() == 1  # just focal species
    assert pd.isna(result.iloc[0, -1])  # at first index, last column
