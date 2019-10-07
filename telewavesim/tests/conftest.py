import pytest
from telewavesim import conf as cf
from telewavesim import utils as ut
from pkg_resources import resource_filename


@pytest.fixture(scope="module")
def load_params():
    modfile = resource_filename('telewavesim',
                                'examples/models/model_Audet2016.txt')
    ut.read_model(modfile)

    cf.nt = 2000
    cf.dt = 0.05
    cf.slow = 0.04
    cf.baz = 0.
    cf.dp = 1000.
    cf.wvtype = 'P'

    ut.model2for()
    ut.wave2for()
    ut.obs2for()

    return
