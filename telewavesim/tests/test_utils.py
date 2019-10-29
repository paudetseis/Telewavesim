import numpy as np
from telewavesim import utils as ut


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

    error_msg = "Failed! shape of stiffness matrices {}".format(
        cc_iso.shape, cc_tri.shape, cc_ani.shape)

    assert cc_iso.shape == cc.shape, error_msg
    assert cc_tri.shape == cc.shape, error_msg
    assert cc_ani.shape == cc.shape, error_msg


def test_rotate():
    a = 6.
    b = 3.6
    tr = 0.
    pl = 0.
    ani = 5.
    cc_tri = ut.set_tri_tensor(a, b, tr, pl, ani)
    aa_tri = ut.rot_tensor(cc_tri, 2.*np.pi, 2.*np.pi, 2.*np.pi)

    assert np.allclose(cc_tri, aa_tri)
