from __future__ import annotations

from io import StringIO

import pytest

from abSENSE.parameters import AbsenseParameters
from abSENSE.recorder import Recorder


@pytest.fixture()
def faked_recorder(mocker, default_params):
    mocker.patch("abSENSE.recorder.os.makedirs")
    default_params.out_dir = "."
    default_params.predict_all = True
    recorder = Recorder.build_recorder(default_params, ["s1", "s2", "s3"])
    # replace files with stringIO for easier testing
    files = {}

    def new_stringIO(*args, **kwargs):
        result = StringIO()
        result.name = args[0]
        files[args[0]] = result
        return result

    mocker.patch("abSENSE.recorder.open", side_effect=new_stringIO)
    return recorder, files


@pytest.fixture()
def faked_recorder_non_closing(mocker, default_params):
    mocker.patch("abSENSE.recorder.os.makedirs")
    default_params.out_dir = "."
    recorder = Recorder.build_recorder(default_params, ["s1", "s2", "s3"])
    # replace files with stringIO for easier testing
    files = {}

    def new_stringIO(*args, **kwargs):
        result = StringIO()
        result.close = lambda: None
        result.name = args[0]
        files[args[0]] = result
        return result

    mocker.patch("abSENSE.recorder.open", side_effect=new_stringIO)
    return recorder, files


@pytest.fixture()
def quick_bitscores():
    return (
        "Gene\tS_cer\tS_par\tS_mik\tS_bay\tS_kud\tS_castellii\tK_waltii\tK_lactis\tA_gossyppi\tY_lipolytica\tA_nidulans\tS_pombe\n"
        "NP_010181.2\t2284.0\t2254.0\t2234.0\t2200.0\t2197.0\t1895.0\t1767.0\t1691.0\t1736.0\t1297.0\t1099.0\t1194.0\n"
        "NP_009362.1\t1027.0\t1016.0\t1018.0\t1008.0\t1010.0\t909.0\t903.0\t894.0\t875.0\t714.0\t691.0\t684.0\n"
        "NP_014555.1\t712.0\t698.0\t696.0\t684.0\t688.0\tN/A\t617.0\t611.0\tN/A\tN/A\tN/A\t344.0\n"
        "NP_116682.3\t352.0\t171.0\t92.8\t0\t82.0\t0\t0\t0\t0\t0\t0\t0\n"
        "NP_011284.1\t308.0\t225.0\t179.0\t191.0\t178.0\t0\t0\t0\t0\t0\t0\t0\n"
        "NP_011878.1\t588.0\t522.0\t446.0\t439.0\t455.0\t0\t0\t0\t0\t0\t0\t0\n"
        "NP_013320.1\t1491.0\t1205.0\t1013.0\t910.0\t996.0\t141.0\t101.0\t102.0\t63.2\t0\t0\t0\n"
        "NP_014160.2\t941.0\t927.0\t889.0\t885.0\t893.0\t638.0\t0\t0\t657.0\t375.0\t378.0\t0\n"
        "NP_014890.1\t394.0\t264.0\t254.0\t230.0\t271.0\t58.9\t0\t0\t0\t0\t0\t0\n"
    )


@pytest.fixture()
def quick_distances():
    return (
        "#Species\tDistance\n"
        "S_cer\t0\n"
        "S_par\t0.051\n"
        "S_mik\t0.088\n"
        "S_kud\t0.104\n"
        "S_bay\t0.108\n"
        "S_castellii\t0.363\n"
        "K_waltii\t0.494\n"
        "A_gossyppi\t0.518\n"
        "K_lactis\t0.557\n"
        "A_nidulans\t0.903\n"
        "S_pombe\t0.922\n"
        "Y_lipolytica\t0.954\n"
    )


@pytest.fixture()
def quick_distances_replicates():
    return (
        "S_cer\t0\t0\t0\t0\n"
        "S_par\t0.051\t0.061\t0.041\t0.052\n"
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


@pytest.fixture()
def fungi_database_lengths():
    return (
        "#Species\tDBsize\n"
        "A_nidulans\t5073867\n"
        "Y_lipolytica\t3139837\n"
        "A_gossyppi\t2335136\n"
        "K_lactis\t2462503\n"
        "S_pombe\t2381887\n"
        "S_castellii\t2768407\n"
        "S_mik\t2612415\n"
        "S_par\t2606124\n"
        "S_kud\t2603233\n"
        "S_cer\t2615464\n"
        "S_bay\t2595487\n"
        "K_waltii\t2424992\n"
    )


@pytest.fixture()
def fungi_gene_lengths():
    return (
        "#GeneID\tGeneLength\n"
        "NP_001018029.1\t66\n"
        "NP_001018030.1\t360\n"
        "NP_001018031.2\t113\n"
        "NP_001018032.1\t52\n"
        "NP_001018033.3\t92\n"
        "NP_001027023.1\t163\n"
        "NP_001027534.1\t73\n"
        "NP_001032571.3\t67\n"
        "NP_001032572.1\t78\n"
        "NP_001032573.1\t79\n"
        "NP_010181.2\t101\n"
        "NP_009362.1\t102\n"
        "NP_014555.1\t103\n"
        "NP_116682.3\t104\n"
        "NP_011284.1\t105\n"
        "NP_011878.1\t106\n"
        "NP_013320.1\t107\n"
        "NP_014160.2\t108\n"
        "NP_014890.1\t109\n"
    )


@pytest.fixture
def default_params(quick_bitscores, quick_distances):
    distances = StringIO(quick_distances)
    distances.name = "distances"
    bitscores = StringIO(quick_bitscores)
    bitscores.name = "bitscores"
    return AbsenseParameters(
        distances=distances,
        bitscores=bitscores,
        e_value=0.001,
        predict_all=False,
        plot_all=False,
        validate=False,
        plot="",
        include_only=None,
        gene_lengths=None,
        db_lengths=None,
        out_dir=None,
        start_time="NOW",
    )
