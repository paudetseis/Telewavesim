from telewavesim import utils as ut
from telewavesim.rmat_f import plane as pw_f
from pkg_resources import resource_filename
import numpy as np
from numpy.fft import fft


def test_plane_obs():
    modfile = resource_filename('telewavesim',
                                'examples/models/model_Audet2016.txt')
    model = ut.read_model(modfile)
    npts = 2000
    wvtype = 'P'
    dt = 0.05
    slow = 0.04
    baz = 0.
    dp = 1000
    c = 1.5
    rhof = 1027

    ut.model2for(model)
    ut.wave2for(dt, slow, baz)
    ut.obs2for(dp, c, rhof)
    yx, yy, yz = pw_f.plane_obs(npts, model.nlay, np.array(wvtype, dtype='c'))
    ux = np.real(fft(yx))
    uy = np.real(fft(yy))
    uz = np.real(fft(yz))
    # seismogram should be maximized on vertical component
    assert np.max(np.abs(uz)) > np.max(np.abs(ux)) > np.max(np.abs(uy)), \
        'Failed! Energy is not maximized on vertical component'

    # tangential component should all be close to zero
    assert np.allclose(uy, np.zeros(len(uy))), 'non-zero values in uy'


def test_plane_land():
    modfile = resource_filename('telewavesim',
                                'examples/models/model_Audet2016.txt')
    model = ut.read_model(modfile)
    npts = 2000
    wvtype = 'P'
    dt = 0.05
    slow = 0.04
    baz = 0.
    ut.model2for(model)
    ut.wave2for(dt, slow, baz)
    yx, yy, yz = pw_f.plane_land(npts, model.nlay, np.array(wvtype, dtype='c'))
    ux = np.real(fft(yx))
    uy = np.real(fft(yy))
    uz = np.real(fft(yz))

    # seismogram should be maximized on vertical component
    assert np.max(np.abs(uz)) > np.max(np.abs(ux)) > np.max(np.abs(uy)), \
        'Failed! Energy is not maximized on vertical component'

    # tangential component should all be close to zero
    assert np.allclose(uy, np.zeros(len(uy))), 'non-zero values in uy'

    trxyz = ut.get_trxyz(ux, uy, uz, npts, dt, slow, baz, wvtype)
    tfs = ut.tf_from_xyz(trxyz)

    nt = tfs[0].stats.npts

    assert nt == npts

    # zero-lag should be maximized on radial component
    assert tfs[0].data[int(nt/2)] > tfs[1].data[int(nt/2)], \
        'Failed! Zero-lag is not maximized on radial component'
