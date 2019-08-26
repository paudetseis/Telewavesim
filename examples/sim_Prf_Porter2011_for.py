'''
PROGRAM sim_P.py

Generates synthetic P waves.

Pascal Audet, 2019

'''

# Import modules and functions
import numpy as np
from obspy.core import Stream, Trace
from obspy.signal.rotate import rotate_ne_rt
from telewavesim import green as gp
from telewavesim import utils as ut
from telewavesim import wiggle as wg
from telewavesim import conf as cf
from telewavesim.rmat_f import green as gp_f
from telewavesim.rmat_f import conf as cf_f

modfile = 'model_Porter2011.txt'
cf.wvtype = 'P'

"""
Set parameters to use in generation of receiver functions

"""
# Frequencies
f1 = 0.01
f2 = 1.

# Sampling values
cf.nt = 3000 
cf.dt = 0.01

# Wave parameters
slow = 0.06
baz = np.arange(0., 360., 10.)

# Set parameters from model file
ut.read_model(modfile)

# Pass global variables to fortran
ut.model2for()

# Calculate time shift (for plotting Green's functions)
t1 = ut.calc_ttime(slow)
print('Predicted arrival time from model: {0:4.1f} sec'.format(t1))

# Initialize streams
trR_f = Stream(); trT_f = Stream()

# Loop over range of data - here simple loop over backazimuths
for bb in baz:

    # Pass baz and slow to global variables
    cf.baz = bb
    cf.slow = slow
    ut.wave2for()

    # Calculate the Green's functions
    ux, uy, uz = gp_f.green_land(cf.nt,cf.nlay,np.array(cf.wvtype, dtype='c'))

    # Store in traces
    tux = Trace(data=ux)
    tuy = Trace(data=uy)
    tuz = Trace(data=uz)

    # Update trace header
    tux = ut.update_stats(tux, cf.nt, cf.dt, cf.slow, cf.baz)
    tuy = ut.update_stats(tuy, cf.nt, cf.dt, cf.slow, cf.baz)
    tuz = ut.update_stats(tuz, cf.nt, cf.dt, cf.slow, cf.baz)

    # Append to stream
    trxyz_f = Stream(traces=[tux, tuy, tuz])

    # Then the transfer functions
    tfs_f = gp.tf_from_xyz(trxyz_f,pvh=False)

    # Append to streams
    trR_f.append(tfs_f[0]); trT_f.append(tfs_f[1])

# Filter to get wave-like traces
trR_f.filter('bandpass',freqmin=f1, freqmax=f2, corners=2, zerophase=True)
trT_f.filter('bandpass',freqmin=f1, freqmax=f2, corners=2, zerophase=True)

# Stack over all traces
trR_f_stack, trT_f_stack = ut.stack_all(trR_f, trT_f, pws=True)

# Plot as wiggles
wg.rf_wiggles_baz(trR_f, trT_f, trR_f[0], trT_f[0], 'test', btyp='baz', scale=1.e3, tmin=-5., tmax=8., save=False, title='_fortran')
    
