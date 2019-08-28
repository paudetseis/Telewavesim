"""
Configuration script to set up global variables

"""

import numpy as np

# wavetype
global wvtype
wvtype = None

# wavefield variables
global slow, baz, nt, dt
slow = None; baz = None
nt = None; dt = None

# model variables
global a, rho, isoflg, thickn, nlay, aa, bb, ani, tr, pl
a = None; rho = None; isoflg = None; thickn = None

# obs variables
global dp, c, rhof
dp = None; c = 1500.; rhof = 1027.

# Global variables for rmatrix
global evecs, evals
global Tui, Tdi, Rui, Rdi
