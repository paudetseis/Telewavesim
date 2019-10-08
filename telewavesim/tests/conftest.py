import pytest
from telewavesim import conf as cf
from telewavesim import utils as ut
from pkg_resources import resource_filename


@pytest.fixture(scope="module")
def load_params():
    modfile = resource_filename('telewavesim',
                                'examples/models/model_Audet2016.txt')
    model = ut.read_model(modfile)

    cf.nt = 2000
    cf.wvtype = 'P'
    dt = 0.05
    slow = 0.04
    baz = 0.
    dp = 1000
    c = 1.5
    rhof = 1027

    model2for(model)
    wave2for(dt, slow, baz)
    obs2for(dp, c, rhof)

    return
