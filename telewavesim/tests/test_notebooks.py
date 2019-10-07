# Copyright 2019 Pascal Audet, Tom Eulenfeld

# This file is part of Telewavesim.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

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
    import numpy as np
    from obspy.core import Stream
    from telewavesim import utils as ut
    from telewavesim import wiggle as wg
    from telewavesim import conf as cf

    modfile = resource_filename('telewavesim',
                                'examples/models/model_Porter2011.txt')
    cf.wvtype = 'P'
    cf.nt = 3000  # Number of samples
    cf.dt = 0.01  # Sample distance in seconds
    cf.slow = 0.06  # Horizontal slowness (or ray parameter) in s/km
    baz = np.arange(0., 360., 10.)
    ut.read_model(modfile)
    trR = Stream()
    trT = Stream()
    # Loop over range of data
    for bb in baz:
        # Pass baz global variable
        cf.baz = bb
        # Calculate the plane waves seismograms
        trxyz = ut.run_plane()
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
                          ftitle=join(tempdir, 'porter2011.png'))


def test_Autdet2016():
    from obspy.core import Stream
    from obspy.signal.rotate import rotate_ne_rt
    from telewavesim import utils as ut
    from telewavesim import wiggle as wg
    from telewavesim import conf as cf
    modfile = resource_filename('telewavesim',
                                'examples/models/model_Audet2016.txt')
    cf.wvtype = 'P'
    cf.nt = 3000  # Number of samples
    cf.dt = 0.01  # Sample distance in seconds
    cf.dp = 2000.  # Deployment depth below sea level in meters
    cf.c = 1500.    # P-wave velocity in salt water (m/s)
    cf.rhof = 1027.  # Density of salt water (kg/m^3)
    cf.slow = 0.06  # Horizontal slowness (or ray parameter) in s/km
    # Back-azimuth direction in degrees (has no influence if model is isotropic)
    cf.baz = 0.
    ut.read_model(modfile)
    assert cf.rho == [2800.0, 2800.0, 3200.0]
    t1 = ut.calc_ttime(cf.slow)
    assert round(t1, 1) == 1.1
    trxyz = ut.run_plane(obs=True)
    tfs = ut.tf_from_xyz(trxyz, pvh=False)
    ntr = trxyz[0]  # North component
    etr = trxyz[1]  # East component
    ztr = trxyz[2]  # Vertical component
    rtr = ntr.copy()  # Radial component
    ttr = etr.copy()  # Transverse component
    rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, cf.baz)
    strf = Stream(traces=[tfs[0], ztr, rtr])
    # Set frequency corners in Hz
    f1 = 0.1
    f2 = 1.0
    # Plot as wiggles
    with tempfile.TemporaryDirectory() as tempdir:
        wg.pw_wiggles_Audet2016(strf, t1=t1, tmax=10., f1=f1, f2=f2,
                                scale=2.e-7, save=True,
                                ftitle=join(tempdir, 'porter2011.png'))


def test_SKS():
    import numpy as np
    from obspy.core import Stream
    from obspy.signal.rotate import rotate_ne_rt
    from telewavesim import utils as ut
    from telewavesim import wiggle as wg
    from telewavesim import conf as cf

    modfile = resource_filename('telewavesim',
                                'examples/models/model_SKS.txt')
    cf.wvtype = 'SV'
    cf.nt = 3000  # Number of samples
    cf.dt = 0.05  # Sample distance in seconds
    cf.slow = 0.04  # Horizontal slowness (or ray parameter) in s/km
    baz = np.arange(0., 360., 10.)
    ut.read_model(modfile)
    t1 = ut.calc_ttime(cf.slow)
    assert round(t1, 1) == 21.6
    trR = Stream()
    trT = Stream()
    # Loop over range of data
    for bb in baz:
        # Pass baz global variable
        cf.baz = bb
        # Calculate the plane wave seismograms
        trxyz = ut.run_plane()
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
                          ftitle=join(tempdir, 'sks.png'))
