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

'''
Compute R/T coefficient for solid-fluid interface.

Original code by M. G. Bostock.

See paper by Bostock & Trehu, BSSA (2012) for details.

'''

import numpy as np
from numpy.linalg import solve
from telewavesim import rmatrix as rm


def rtfluid(p1, p2, c, rhof, a, b, rhos):
    """ 
    Function to compute reflection and transmission coefficients for 
    a fluid-solid boundary

    Args:
        p1 (float): x-component of horizontal slowness (s/km)
        p2 (float): y-component of horizontal slowness (s/km)
        c (float): P-wave velocity of salt water (km/s)
        rhof (float): Density of salt water (kg/m^3)
        a (float): P-wave velocity of topmost solid layer (km/s)
        b (float): S-wave velocity of topmost solid layer (km/s)
        rhos (float): Density of topmost solid layer (kg/m^3)

    Returns:
        (tuple): Tuple containing R/T coefficients f
            * TdPP (float):
            * TdSP (float):
            * RdPP (float): 
            * TuPP (float):
            * TuPS (float): 
            * RuPP (float): 
            * RuSP (float): 
            * RuPS (float): 
            * RuSS (float):

    """

    pp = p1*p1 + p2*p2

    # Fluid properties.
    qf = np.sqrt(1./c**2 - pp)
    Ff = np.array([[qf*c, -qf*c],[rhof*c, rhof*c]])

    # Solid properties (call isotroc for consistency with 
    # Rmatrix on conventions of signs etc.)
    q, N = rm.isotroc(a, b, rhos, p1, p2)
    Fs = np.array([[N[0,0],N[0,1],N[0,3],N[0,4]],
        [N[2,0],N[2,1],N[2,3],N[2,4]],
        [N[3,0],N[3,1],N[3,3],N[3,4]],
        [N[5,0],N[5,1],N[5,3],N[5,4]]])
    Fl = np.array([Fs[1,:],Fs[3,:]])

    # Define matrix L that encapsulates continuity conditions (continuity of
    # vertical displacement and vertical traction) as well as 
    # vanishing horizontal traction at fluid solid interface.
    Fl = solve(Ff,Fl)
    L = np.vstack((Fl,Fs[2,:]))

    # P-incidence from above.
    dum = solve([[L[0,0],L[0,1],0],
        [L[1,0],L[1,1],-1],
        [L[2,0],L[2,1],0]],[1,0,0])

    TdPP = dum[0]
    TdSP = dum[1]
    RdPP = dum[2]

    # P-incidence from below.
    dum = solve([[0,L[0,0],L[0,1]],
        [-1,L[1,0],L[1,1]],
        [0,L[2,0],L[2,1]]],[-L[0,2],-L[1,2],-L[2,2]])

    TuPP = dum[0]
    RuPP = dum[1]
    RuSP = dum[2]

    # S-incidence from below.
    dum = solve([[0.,L[0,0],L[0,1]],
        [-1.,L[1,0],L[1,1]],
        [0.,L[2,0],L[2,1]]],[-L[0,3],-L[1,3],-L[2,3]])

    TuPS = dum[0]
    RuPS = dum[1]
    RuSS = dum[2]

    return TdPP, TdSP, RdPP, TuPP, TuPS, RuPP, RuSP, RuPS, RuSS



