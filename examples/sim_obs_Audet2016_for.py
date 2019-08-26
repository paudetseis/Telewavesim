'''
PROGRAM sim_obs_Audet2016_for.py

Generates synthetic receiver functions for an ocean-bottom 
seismograph (OBS) station. The reference model is described 
in Audet, GRI, 2016 (see also the file model_Audet2016.txt 
for details).

This version uses the fortran implementation of the software,
wrapped using numpy f2py.

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

# Pass global variables to fortran
ut.obs2for()

# Ray parameters
slow = 0.06
baz = 0.

# Set model parameters from file
ut.read_model(modfile)

# Pass global variables to fortran
ut.model2for()

# Calculate time shift (for plotting Green's functions - only isotropic cases)
t1 = ut.calc_ttime(slow)
print('Predicted arrival time from model: {0:4.1f} sec'.format(t1))

# Pass baz and slow to global variables
cf.baz = baz
cf.slow = slow
ut.wave2for()

# Calculate the Green's functions
ux, uy, uz = gp_f.green_obs(cf.nt,cf.nlay,np.array(cf.wvtype, dtype='c'))

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
tfs = gp.tf_from_xyz(trxyz_f, pvh=False)

# Extract East, North and Vertical
ntr = trxyz_f[0]
etr = trxyz_f[1]
ztr = trxyz_f[2]

# Copy to radial and transverse
rtr = ntr.copy()
ttr = etr.copy()

# Rotate to radial and transverse
rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)

# Store to Stream
strf = Stream(traces=[tfs[0],ztr,rtr])

# Plot as wiggles
wg.gf_wiggles_Audet2016(strf, t1=t1, tmax=10., f1=f1, f2=f2, scale=2.e-7, save=False, title='')
    
