'''
PROGRAM sim_obs_Audet2016.py

Generates synthetic receiver functions for an ocean-bottom 
seismograph (OBS) station. The reference model is described 
in Audet, GRI, 2016 (see also the file model_Audet2016.txt 
for details).

This version uses the numpy implementation of the software.

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

modfile = 'model_Audet2016.txt'
cf.wvtype = 'P'

"""
Set parameters of time series to use in generation of receiver functions

"""
# Frequencies
f1 = 0.01
f2 = 1.

# Sampling values
cf.nt = 3000
cf.dt = 0.01

# OBS parameters
cf.dp = 2000.

# Parameters of fluid layer
cf.c = 1500.    # P-wave velocity in salt water
cf.rhof = 1027. # Density of salt water

# Ray parameters
slow = 0.06
baz = 0.

# Set model parameters from file
ut.read_model(modfile)

# Calculate time shift (for plotting Green's functions - only isotropic cases)
t1 = ut.calc_ttime(slow)
print('Predicted arrival time from model: {0:4.1f} sec'.format(t1))

# Pass baz and slow to global variables
cf.baz = baz
cf.slow = slow

# Calculate the Green's functions
trxyz = gp.green_obs()

# Then the transfer functions
tfs = gp.tf_from_xyz(trxyz, pvh=False)

# Extract East, North and Vertical
ntr = trxyz[0]
etr = trxyz[1]
ztr = trxyz[2]

# Copy to radial and transverse
rtr = ntr.copy()
ttr = etr.copy()

# Rotate to radial and transverse
rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)

# Store to Stream
strf = Stream(traces=[tfs[0],ztr,rtr])

# Plot as wiggles
wg.gf_wiggles_Audet2016(strf, t1=t1, tmax=10., f1=f1, f2=f2, scale=2.e-7, save=False, title='')
    
