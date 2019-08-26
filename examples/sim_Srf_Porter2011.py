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
cf.wvtype = 'SV'

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

# Calculate time shift (for plotting Green's functions - only isotropic cases)
t1 = ut.calc_ttime(slow)
print('Predicted arrival time from model: {0:4.1f} sec'.format(t1))

# Initialize streams
trR = Stream(); trT = Stream()

# Loop over range of data - here simple loop over backazimuths
for bb in baz:

    # Pass baz and slow to global variables
    cf.baz = bb
    cf.slow = slow

    # Calculate the Green's functions
    trxyz = gp.green_land()

    # Then the transfer functions
    tfs = gp.tf_from_xyz(trxyz,pvh=False)

    # Append to streams
    trR.append(tfs[0]); trT.append(tfs[1])

# Filter to get wave-like traces
trR.filter('bandpass',freqmin=f1, freqmax=f2, corners=2, zerophase=True)
trT.filter('bandpass',freqmin=f1, freqmax=f2, corners=2, zerophase=True)

# Stack over all traces
trR_stack, trT_stack = ut.stack_all(trR, trT, pws=True)

# Plot as wiggles
wg.rf_wiggles_baz(trR, trT, trR_stack, trT_stack, 'test', btyp='baz', scale=1.e3, tmin=-5., tmax=8., save=False, title='_python')
    
