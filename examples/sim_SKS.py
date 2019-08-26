'''
PROGRAM sim_SKS.py

Generates synthetic SKS waves.

Pascal Audet, 2019

'''

# Import modules and functions
import numpy as np
from obspy.core import Stream
from obspy.signal.rotate import rotate_ne_rt
from telewavesim import green as gp
from telewavesim import utils as ut
from telewavesim import wiggle as wg
from telewavesim import conf as cf

modfile = 'model_SKS.txt'
cf.wvtype = 'SV'

"""
Set parameters of time series to use in generation of receiver functions

"""
# Frequencies
f1 = 0.01
f2 = 0.2

# Sampling values
cf.nt = 4000 
cf.dt = 0.2

# Ray parameters
slow = 0.04
baz = np.arange(0., 360., 10.)

# Set parameters from model file
ut.read_model(modfile)

# Calculate time shift (for plotting Green's functions - only isotropic cases)
t1 = ut.calc_ttime(slow)
print('Predicted arrival time from model: {0:4.1f} sec'.format(t1))

# Initialize streams
trR = Stream()
trT = Stream()

# Loop over range of data - here simple look over backazimuths
for bb in baz:

    # Pass baz and slow to global variables
    cf.baz = bb
    cf.slow = slow

    # Calculate the Green's functions
    trxyz = gp.green_land()

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
    
# Filter to get wave-like traces
trR.filter('bandpass',freqmin=f1, freqmax=f2, corners=4, zerophase=True)
trT.filter('bandpass',freqmin=f1, freqmax=f2, corners=4, zerophase=True)

# Stack over all traces
trR_stack, trT_stack = ut.stack_all(trR, trT, pws=True)

# Plot as wiggles
wg.gf_wiggles_baz(trR, trT, 'test', btyp='baz', scale=0.025, t1=t1, tmin=0., tmax=40, save=False, title='_synt')
    
