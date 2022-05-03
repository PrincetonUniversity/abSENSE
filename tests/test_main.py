from abSENSE.parameters import AbsenseParameters
import abSENSE.__main__ as main
from io import StringIO
import pytest


def test_defaults(faked_recorder, default_params):
    recorder, files = faked_recorder
    main.perform_analysis(default_params)
