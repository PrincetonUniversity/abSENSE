import pytest
from io import StringIO

from abSENSE.recorder import FileRecorder
from abSENSE.parameters import AbsenseParameters


@pytest.fixture()
def faked_recorder(mocker, default_params):
    mocker.patch('abSENSE.recorder.os.makedirs')
    default_params.out_dir = '.'
    default_params.predict_all = True
    recorder = FileRecorder(default_params, ['s1', 's2', 's3'])
    # replace files with stringIO for easier testing
    files = {}
    def new_stringIO(*args, **kwargs):
        result = StringIO()
        files[args[0]] = result
        return result
    mocker.patch('abSENSE.recorder.open', side_effect=new_stringIO)
    return recorder, files


@pytest.fixture()
def quick_bitscores():
    return (
        'Gene\tS_cer\tS_par\tS_mik\tS_bay\tS_kud\tS_castellii\tK_waltii\tK_lactis\tA_gossyppi\tY_lipolytica\tA_nidulans\tS_pombe\n'
        'NP_010181.2\t2284.0\t2254.0\t2234.0\t2200.0\t2197.0\t1895.0\t1767.0\t1691.0\t1736.0\t1297.0\t1099.0\t1194.0\n'
        'NP_009362.1\t1027.0\t1016.0\t1018.0\t1008.0\t1010.0\t909.0\t903.0\t894.0\t875.0\t714.0\t691.0\t684.0\n'
        'NP_014555.1\t712.0\t698.0\t696.0\t684.0\t688.0\tN/A\t617.0\t611.0\tN/A\tN/A\tN/A\t344.0\n'
        'NP_116682.3\t352.0\t171.0\t92.8\t0\t82.0\t0\t0\t0\t0\t0\t0\t0\n'
        'NP_011284.1\t308.0\t225.0\t179.0\t191.0\t178.0\t0\t0\t0\t0\t0\t0\t0\n'
        'NP_011878.1\t588.0\t522.0\t446.0\t439.0\t455.0\t0\t0\t0\t0\t0\t0\t0\n'
        'NP_013320.1\t1491.0\t1205.0\t1013.0\t910.0\t996.0\t141.0\t101.0\t102.0\t63.2\t0\t0\t0\n'
        'NP_014160.2\t941.0\t927.0\t889.0\t885.0\t893.0\t638.0\t0\t0\t657.0\t375.0\t378.0\t0\n'
        'NP_014890.1\t394.0\t264.0\t254.0\t230.0\t271.0\t58.9\t0\t0\t0\t0\t0\t0\n'
    )


@pytest.fixture()
def quick_distances():
    return (
        'S_cer\t0\n'
        'S_par\t0.051\n'
        'S_mik\t0.088\n'
        'S_kud\t0.104\n'
        'S_bay\t0.108\n'
        'S_castellii\t0.363\n'
        'K_waltii\t0.494\n'
        'A_gossyppi\t0.518\n'
        'K_lactis\t0.557\n'
        'A_nidulans\t0.903\n'
        'S_pombe\t0.922\n'
        'Y_lipolytica\t0.954\n'
    )


@pytest.fixture
def default_params(quick_bitscores, quick_distances):
    return AbsenseParameters(
        distances=StringIO(quick_distances),
        bitscores=StringIO(quick_bitscores),
        e_value=0.001,
        predict_all=False,
        include_only=None,
        gene_lengths=None,
        db_lengths=None,
        out_dir=None,
        start_time='NOW',
    )
