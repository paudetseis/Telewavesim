import telewavesim
import numpy as np
from telewavesim import conf as cf
from telewavesim import utils as ut
from telewavesim.rmat_f import conf as cf_f

from telewavesim.tests.conftest import load_params

def test_model2for_conf(load_params):

    error_msg = "Failed! Python model conf values are {} and Fortran conf values are {}".format(
        [cf.a, cf.rho, cf.thickn],[cf_f.a[:,:,:,:,0:cf.nlay], cf_f.rho[0:cf.nlay], cf_f.thickn[0:cf.nlay]])
    assert np.all(cf.a==cf_f.a[:,:,:,:,0:cf.nlay]), error_msg
    assert np.all(cf.rho==cf_f.rho[0:cf.nlay]), error_msg
    assert np.all(cf.thickn==cf_f.thickn[0:cf.nlay]), error_msg

def test_wave2for_conf(load_params):

    error_msg = "Failed! Python wave conf values are {} and Fortran conf values are {}".format(
        [cf.dt, cf.slow, cf.baz],[cf_f.dt, cf_f.slow, cf_f.baz])
    assert cf.dt==cf_f.dt, error_msg
    assert cf.slow==cf_f.slow, error_msg
    assert cf.baz==cf_f.baz, error_msg

def test_obs2for_conf(load_params):

    error_msg = "Failed! Python obs conf values are {} and Fortran conf values are {}".format(
        [cf.dp, cf.c, cf.rhof],[cf_f.dp, cf_f.c, cf_f.rhof])
    assert cf.dp==cf_f.dp, error_msg
    assert cf.c==cf_f.c, error_msg
    assert cf.rhof==cf_f.rhof, error_msg

def test_tensor_shape():
    a = 0.
    b = 0.
    tr = 0.
    pl = 0.
    ani = 0.
    typ = 'atg'
    cc_iso = ut.set_iso_tensor(a, b)
    cc_tri = ut.set_tri_tensor(a, b, tr, pl, ani)
    cc_ani, rho = ut.set_aniso_tensor(tr, pl, typ=typ)
    cc = np.zeros((3, 3, 3, 3))

    error_msg = "Failed! shape of stiffness matrices {}".format(cc_iso.shape,cc_tri.shape,cc_ani.shape)

    assert cc_iso.shape==cc.shape, error_msg
    assert cc_tri.shape==cc.shape, error_msg
    assert cc_ani.shape==cc.shape, error_msg

def test_rotate():
    a = 6.
    b = 3.6
    tr = 0.
    pl = 0.
    ani = 5.
    cc_tri = ut.set_tri_tensor(a, b, tr, pl, ani)
    aa_tri = ut.rot_tensor(cc_tri,2.*np.pi,2.*np.pi,2.*np.pi)

    assert np.allclose(cc_tri, aa_tri)
