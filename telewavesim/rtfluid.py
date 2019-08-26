'''
MODULE rt_fluid.py

Function translated from M. Bostock's Matlab code.

See paper by Bostock & Trehu, BSSA (2012) for details.

'''

import numpy as np
from numpy.linalg import solve
from telewavesim import rmatrix as rm


def rtfluid(p1, p2, c, rhof, a, b, rhos):
    """ 
    Function to compute reflection and transmission coefficients for 
    a fluid-solid boundary

    p is slowness, c is fluid velocity, rhof is fluid density, (a, b) are P and
    S velocities of solid and rhos is solid density

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



