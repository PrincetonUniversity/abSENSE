from abSENSE.analyzer import AbsenseAnalyzer
from abSENSE.exceptions import MissingGeneException, MissingSpeciesException
from abSENSE.results import SampledResult, ErrorResult, NotEnoughDataResult
from abSENSE.utilities import sample_parameters, find_confidence_interval
from io import StringIO
import pytest
import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal as afe, assert_series_equal as ase


def test_init_defaults(default_params):
    analyzer = AbsenseAnalyzer(default_params)

    assert (analyzer.gene_lengths.index == analyzer.genes).all()
    assert (analyzer.gene_lengths['length'] == default_params.default_gene_length).all()

    assert (analyzer.db_lengths.index == analyzer.species).all()
    assert (analyzer.db_lengths['length'] == default_params.default_db_length).all()


def test_init_with_include_only(default_params):
    # include only (or the distance index on default)
    # determines the order of bitscore columns
    default_params.include_only = 'S_cer,S_mik,S_par'
    analyzer = AbsenseAnalyzer(default_params)
    assert (analyzer.species == 'S_cer,S_mik,S_par'.split(',')).all()
    assert (analyzer.bitscores.columns == 'S_cer,S_mik,S_par'.split(',')).all()


def test_init_with_include_only_missing(default_params):
    default_params.include_only = 'S_cer,S_par,S_mik,NOTHERE'
    with pytest.raises(MissingSpeciesException) as error:
        AbsenseAnalyzer(default_params)
    assert str(error.value) == ('Unable to find all requested species in bitscores. '
                                "Missing: ['NOTHERE']")


def test_init_with_missing_db_species(default_params):
    default_params.include_only = 'S_cer,S_par,S_mik,S_bay'
    default_params.db_lengths = StringIO(
        '#Species\tDBsize\n'
        'S_par\t2606124\n'
        'S_mik\t2603233\n'
        'S_cer\t2615464\n'
    )
    with pytest.raises(MissingSpeciesException) as error:
        AbsenseAnalyzer(default_params)
    assert str(error.value) == ('Unable to find all requested species in database lengths. '
                                "Missing: ['S_bay']")

def test_init_with_missing_gene_lengths(default_params):
    default_params.gene_lengths = StringIO(
        '#GeneID\tGeneLength\n'
        'NP_001018029.1\t66\n'
        'NP_001018030.1\t360\n'
        'NP_001018031.2\t113\n'
    )
    with pytest.raises(MissingGeneException) as error:
        AbsenseAnalyzer(default_params)
    assert str(error.value).startswith('Unable to find all requested genes in gene lengths. '
                                       "Missing: ['NP_010181.2', ")


@pytest.fixture()
def analyzer_with_files(fungi_database_lengths, fungi_gene_lengths, default_params):
    default_params.gene_lengths = StringIO(fungi_gene_lengths)
    default_params.db_lenghts = StringIO(fungi_database_lengths)
    return AbsenseAnalyzer(default_params)


def test_fit_genes(analyzer_with_files):
    analyzer_with_files.fit_genes()


def test_fit_gene_not_enough(analyzer_with_files):
    result = analyzer_with_files.fit_gene('NP_001018030.1', pd.Series(
        [0, 1, 2] + [np.nan]*9,
        index=analyzer_with_files.species,
    ))
    assert result == NotEnoughDataResult('NP_001018030.1')


def test_fit_gene_analysis_error(analyzer_with_files, mocker):
    mocker.patch('scipy.optimize.curve_fit', side_effect=RuntimeError)
    result = analyzer_with_files.fit_gene('NP_001018030.1', pd.Series(
        [3, 2, 1, 0] + [np.nan]*8,
        index=analyzer_with_files.species,
    ))
    assert result == ErrorResult('NP_001018030.1')



def test_fit_gene_infinite_covariance(analyzer_with_files, mocker):
    mocker.patch('scipy.optimize.curve_fit',
                 return_value=((1, 1), [[np.inf, 0], [0, 1]]))
    result = analyzer_with_files.fit_gene('NP_001018030.1', pd.Series(
        [3, 2, 1, 0] + [np.nan]*8,
        index=analyzer_with_files.species,
    ))
    assert isinstance(result, ErrorResult)
    assert result.gene == 'NP_001018030.1'
    ase(result.predictions, np.exp(-1 * analyzer_with_files.distances)['distance'])


def test_fit_normal(analyzer_with_files):
    result = next(analyzer_with_files.fit_genes())
    assert isinstance(result, SampledResult)
    assert result.gene == 'NP_010181.2'
    assert result.a_fit == pytest.approx(2359.91, abs=1e-2)
    assert result.b_fit == pytest.approx(0.67612, abs=1e-4)
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
    old_result = -1 * np.log2(analyzer_with_files.e_value / (100 * analyzer_with_files.db_lengths))
    old_result = old_result.rename(columns={'length': 'bit_threshold'})
    result = analyzer_with_files.bit_threshold(100)
    afe(result, old_result)

    old_result = -1 * np.log2(analyzer_with_files.e_value / (234 * analyzer_with_files.db_lengths))
    old_result = old_result.rename(columns={'length': 'bit_threshold'})
    result = analyzer_with_files.bit_threshold(234)
    afe(result, old_result)


def test_find_confidence_interval(analyzer_with_files):
    bit_threshold = pd.DataFrame(
        {'bit_threshold': 40.1},
        index=analyzer_with_files.species,
    )
    result = find_confidence_interval(
        analyzer_with_files.random,
        analyzer_with_files.distances.to_numpy(),
        np.array([
            [40.1, 1],
            [40.2, 2],
            [40.3, 3],
            [40.4, -1],
        ]),
        bit_threshold,
        analyzer_with_files.species,
    )

    cer = result.iloc[0, :]
    assert cer['p_values'] == pytest.approx(0.089856, abs=1e-4)
    assert cer['low_interval'] == pytest.approx(39.962014, abs=1e-4)
    assert cer['high_interval'] == pytest.approx(40.537986, abs=1e-4)
