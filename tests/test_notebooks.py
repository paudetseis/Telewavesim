# Copyright 2019 Pascal Audet, Tom Eulenfeld

# This file is part of Telewavesim.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from pkg_resources import resource_filename
import tempfile
from os.path import join


def test_Porter2011():
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    from obspy.core import Stream
    from telewavesim import utils as ut
    from telewavesim import wiggle as wg

    modfile = resource_filename('telewavesim',
                                'examples/models/model_Porter2011.txt')
    wvtype = 'P'
    npts = 3000  # Number of samples
    dt = 0.01  # Sample distance in seconds
    slow = 0.06  # Horizontal slowness (or ray parameter) in s/km
    baz = np.arange(0., 360., 10.)
    model = ut.read_model(modfile)
    trR = Stream()
    trT = Stream()
    # Loop over range of data
    for bb in baz:
        # Calculate the plane waves seismograms
        trxyz = ut.run_plane(model, slow, npts, dt, bb, wvtype=wvtype,
                             obs=False)
        # Then the transfer functions in Z-R-T coordinate system
        tfs = ut.tf_from_xyz(trxyz, pvh=False)
        # Append to streams
        trR.append(tfs[0])
        trT.append(tfs[1])
    # Set frequency corners in Hz
    f1 = 0.01
    f2 = 1.0
    # Filter to get wave-like traces
    trR.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
    trT.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
    # Stack over all traces
    trR_stack, trT_stack = ut.stack_all(trR, trT, pws=True)
    # Plot as wiggles
    with tempfile.TemporaryDirectory() as tempdir:
        wg.rf_wiggles_baz(trR, trT, trR_stack, trT_stack, 'test', btyp='baz',
                          scale=1.e3, tmin=-5., tmax=8., save=True,
                          ftitle=join(tempdir, 'porter2011.png'),
                          wvtype='P')


def test_Audet2016():
    import matplotlib
    matplotlib.use('Agg')
    from obspy.core import Stream
    from obspy.signal.rotate import rotate_ne_rt
    from telewavesim import utils as ut
    from telewavesim import wiggle as wg
    modfile = resource_filename('telewavesim',
                                'examples/models/model_Audet2016.txt')
    wvtype = 'P'
    npts = 3000  # Number of samples
    dt = 0.01  # Sample distance in seconds
    dp = 2000.  # Deployment depth below sea level in meters
    c = 1.500    # P-wave velocity in salt water (km/s)
    rhof = 1027.  # Density of salt water (kg/m^3)
    slow = 0.06  # Horizontal slowness (or ray parameter) in s/km
    # Back-azimuth direction in degrees
    # (has no influence if model is isotropic)
    baz = 0.
    model = ut.read_model(modfile)
    assert list(model.rho) == [2800.0, 2800.0, 3200.0]
    t1 = ut.calc_ttime(model, slow, wvtype=wvtype)
    assert round(t1, 1) == 1.1
    trxyz = ut.run_plane(model, slow, npts, dt, baz=baz, wvtype=wvtype,
                         obs=True, dp=dp, c=c, rhof=rhof)
    tfs = ut.tf_from_xyz(trxyz, pvh=False)
    ntr = trxyz[0]  # North component
    etr = trxyz[1]  # East component
    ztr = trxyz[2]  # Vertical component
    rtr = ntr.copy()  # Radial component
    ttr = etr.copy()  # Transverse component
    rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)
    strf = Stream(traces=[tfs[0], ztr, rtr])
    # Set frequency corners in Hz
    f1 = 0.1
    f2 = 1.0
    # Plot as wiggles
    with tempfile.TemporaryDirectory() as tempdir:
        wg.pw_wiggles_Audet2016(strf, t1=t1, tmax=10., f1=f1, f2=f2,
                                ftitle=join(tempdir, 'audet2016'),
                                scale=2.e-7, save=True)


def test_SKS():
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    from obspy.core import Stream
    from obspy.signal.rotate import rotate_ne_rt
    from telewavesim import utils as ut
    from telewavesim import wiggle as wg

    modfile = resource_filename('telewavesim',
                                'examples/models/model_SKS.txt')
    wvtype = 'SV'
    npts = 3000  # Number of samples
    dt = 0.05  # Sample distance in seconds
    slow = 0.04  # Horizontal slowness (or ray parameter) in s/km
    baz = np.arange(0., 360., 10.)
    model = ut.read_model(modfile)
    t1 = ut.calc_ttime(model, slow, wvtype=wvtype)
    assert round(t1, 1) == 21.6
    trR = Stream()
    trT = Stream()
    # Loop over range of data
    for bb in baz:
        # Calculate the plane wave seismograms
        trxyz = ut.run_plane(model, slow, npts, dt, bb, wvtype=wvtype)
        # Extract East, North and Vertical
        ntr = trxyz[0]
        etr = trxyz[1]
        ztr = trxyz[2]
        # Copy to radial and transverse
        rtr = ntr.copy()
        ttr = etr.copy()
        # Rotate to radial and transverse
        rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, bb)
        # Append to streams
        trR.append(rtr)
        trT.append(ttr)

    # Set frequency corners in Hz
    f1 = 0.01
    f2 = 0.2
    # Filter to get wave-like traces
    trR.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
    trT.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
    # Plot as wiggles
    with tempfile.TemporaryDirectory() as tempdir:
        wg.pw_wiggles_baz(trR, trT, 'test', btyp='baz', scale=0.05,
                          t1=t1, tmin=0., tmax=40, save=True,
                          ftitle=join(tempdir, 'sks'),
                          wvtype='SV')
