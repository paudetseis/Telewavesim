# Copyright 2019 Pascal Audet

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

"""

Configuration script to set up global variables

Variables are:
    ``wavetype``:
        - wvtype (str): Incident wavetype (``'P'``, ``'SV'``, ``'SH'``, ``'Si'``)

    ``wavefield parameters``:
        - slow (float): Slowness (s/km)
        - baz (float): Back-azimuth (degree)
        - nt (int): Number of samples in time series
        - dt (float): Sampling rate (Hz)

    ``model parameters``:
        - a (np.ndarray): Elastic thickness (shape ``(3, 3, 3, 3, nlay)``)
        - rho (np.ndarray): Density (kg/m^3) (shape ``(nlay)``)
        - isoflg (list of str): Flags for type of layer material (dimension ``nlay``)
        - thickn (np.ndarray): Thickness of layers (km) (shape ``(nlay)``)
        - nlay (int): Number of layers
        - aa (np.ndarray): P-wave velocity (km/s) (shape ``(nlay)``)
        - bb (np.ndarray): S-wave velocity (km/s) (shape ``(nlay)``)
        - tr (np.ndarray): Trend of symmetry axis (degree) (shape ``(nlay)``)
        - pl (np.ndarray): Plunge of symmetry axis (degree) (shape ``(nlay)``)

    ``obs parameters``:
        - dp (float, optional): Deployment depth below sea level (km)
        - c (float): P-wave velocity of salt water (default = ``1050.0`` km/s)
        - rhof (float): Density of salt water (default = ``1027.0`` kg/m^3)

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
