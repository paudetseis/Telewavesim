import pytest
from telewavesim import conf as cf
from telewavesim import utils as ut

@pytest.fixture(scope="module")
def load_params():

    try:
        modfile='tests/test_model_Audet2016.txt'
        ut.read_model(modfile)
    except:
        modfile='test_model_Audet2016.txt'
        ut.read_model(modfile)

    cf.nt = 2000
    cf.dt = 0.05
    cf.slow = 0.04
    cf.baz = 0.
    cf.dp = 1000.
    cf.wvtype='P'

    ut.model2for()
    ut.wave2for()
    ut.obs2for()

    return