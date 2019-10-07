from telewavesim import utils as ut
from telewavesim.rmat_f import plane as pw_f
from telewavesim import conf as cf
import numpy as np
import pyfftw

from telewavesim.tests.conftest import load_params

def test_plane_obs(load_params):

    yx, yy, yz = pw_f.plane_obs(cf.nt, cf.nlay, np.array(cf.wvtype, dtype='c'))
    ux = np.real(pyfftw.interfaces.numpy_fft.fft(yx))
    uy = np.real(pyfftw.interfaces.numpy_fft.fft(yy))
    uz = np.real(pyfftw.interfaces.numpy_fft.fft(yz))

    # seismogram should be maximized on vertical component
    assert np.max(np.abs(uz)) > np.max(np.abs(ux)) > np.max(np.abs(uy)), \
        'Failed! Energy is not maximized on vertical component'

    # tangential component should all be close to zero
    assert np.allclose(uy, np.zeros(len(uy))), 'non-zero values in uy'

def test_plane_land(load_params):

    yx, yy, yz = pw_f.plane_land(cf.nt, cf.nlay, np.array(cf.wvtype, dtype='c'))
    ux = np.real(pyfftw.interfaces.numpy_fft.fft(yx))
    uy = np.real(pyfftw.interfaces.numpy_fft.fft(yy))
    uz = np.real(pyfftw.interfaces.numpy_fft.fft(yz))

    # seismogram should be maximized on vertical component
    assert np.max(np.abs(uz)) > np.max(np.abs(ux)) > np.max(np.abs(uy)), \
        'Failed! Energy is not maximized on vertical component'

    # tangential component should all be close to zero
    assert np.allclose(uy, np.zeros(len(uy))), 'non-zero values in uy'

    trxyz = ut.get_trxyz(ux, uy, uz)
    tfs = ut.tf_from_xyz(trxyz)

    nt = tfs[0].stats.npts

    # zero-lag should be maximized on radial component
    assert tfs[0].data[int(nt/2)] > tfs[1].data[int(nt/2)], \
        'Failed! Zero-lag is not maximized on radial component'

