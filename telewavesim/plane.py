'''

Calculate plane wave seismograms for stacks of layers using \
the propagator matrix approach of Kennett (1983), implemented by \
Colin J. Thomson (1995).

Original code by C. Thomson.

'''

import numpy as np
import pyfftw
from numpy.linalg import solve, inv
from telewavesim import rmatrix as rm
from telewavesim import rtfluid as rtf
from telewavesim import conf as cf
from telewavesim import utils as ut
from obspy.core import Trace, Stream


def plane_land():
    """
    Function to generate plane wave seismograms from a stack of layers
    for land surface stations at a single slowness. Also handles anisotropy. 
    All model and time series properties are passed through the 
    configuration module ``conf``. Handles any number of layers.

    Returns:
        (tuple): Tuple containing:
            * ux (np.ndarray): x-component displacement seismogram
            * uy (np.ndarray): y-component displacement seismogram
            * uz (np.ndarray): z-conponent displacement seismogram

    """

    # Local versions of global parameters
    dt = cf.dt
    n2 = cf.nt
    slow = cf.slow
    baz = cf.baz
    nlay = cf.nlay

    # Other local parameters
    ifree = 1
    iter = -1
    omg = 1.+0.001j
    om0 = 1.+0.j

    # Frequency and time axes
    omega = np.arange(0,n2/2+1)*omg*2.*np.pi/(n2*dt)
    time = np.arange(0,n2*dt,dt)

    # Slowness vector
    psi = (baz - 180.)*np.pi/180.
    p1 = slow*1.e-3*np.cos(psi)
    p2 = slow*1.e-3*np.sin(psi)

    # Produce R/T matrices for interfaces in solid
    Tus, Rus, Tds, Rds = rm.rmatrix(p1,p2,om0,nlay,ifree,iter)

    # Get partitions of matrix D (Kennett, p.214)
    md = cf.evecs[0:3,0:3,0]
    mu = cf.evecs[0:3,3:6,0]
    nd = cf.evecs[3:6,0:3,0]
    nu = cf.evecs[3:6,3:6,0]
        
    # Free surface matrix. Add +0. to remove negative sign of -0.
    Ruf0 = -np.dot(np.linalg.inv(nd),nu) + 0. 
    
    # Initialize upgoing wavefield
    wup = np.zeros((3,len(omega)),dtype=complex)
    wuv = np.zeros((3,len(omega)),dtype=complex)
    wuh = np.zeros((3,len(omega)),dtype=complex)

    # Initialize displacement vector
    y = np.zeros((3,len(omega)),dtype=complex)

    a = pyfftw.empty_aligned(len(omega), dtype='float')

    # Setup loop for rmatrix
    for iw in range(0,int(n2/2+1)):

        # R/T matrices for solid
        Tus, Rus, Tds, Rds = rm.rmatrix(p1,p2,omega[iw],nlay,ifree,iw)
    
        # Upgoing wavefield at free surface
        wup[0:3,iw] = np.dot(np.dot(inv(np.identity(3) - \
                np.dot(Rds,Ruf0)),Tus),np.array([1,0,0]))
        wuv[0:3,iw] = np.dot(np.dot(inv(np.identity(3) - \
                np.dot(Rds,Ruf0)),Tus),np.array([0,1,0]))
        wuh[0:3,iw] = np.dot(np.dot(inv(np.identity(3) - \
                np.dot(Rds,Ruf0)),Tus),np.array([0,0,1]))

        # Displacement vector (S wave is defined as isotropic source)
        if cf.wvtype=='P':
            y[0:3,iw] = np.dot(mu + np.dot(md,Ruf0),wup[0:3,iw])
        elif cf.wvtype=='Si':
            y[0:3,iw] = np.dot(mu + np.dot(md,Ruf0),0.5*(wuv[0:3,iw]+wuh[0:3,iw]))
        elif cf.wvtype=='SV':
            y[0:3,iw] = np.dot(mu + np.dot(md,Ruf0),wuv[0:3,iw])
        elif cf.wvtype=='SH':
            y[0:3,iw] = np.dot(mu + np.dot(md,Ruf0),wuh[0:3,iw])
        else:
            print('Error - wave type can only be P, Si, SV or SH')
            return
    
    # Get displacements in time domain
    nend = int(n2/2)
    ux = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[0,0:nend+1],\
            np.conj(np.flipud(y[0,1:nend])))))).transpose()
    uy = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[1,0:nend+1],\
            np.conj(np.flipud(y[1,1:nend])))))).transpose()
    uz = -np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[2,0:nend+1],\
            np.conj(np.flipud(y[2,1:nend])))))).transpose()

    # Return displacement stream
    return ux, uy, uz


def plane_obs():
    """
    Function to generate plane wave seismograms from a stack of layers
    for ocean-bottom stations at a single slowness. Also handles anisotropy. 
    All model and time series properties are passed through the 
    configuration module ``conf``. Handles any number of layers.

    Returns:
        (tuple): Tuple containing:
            * ux (np.ndarray): x-component displacement seismogram
            * uy (np.ndarray): y-component displacement seismogram
            * uz (np.ndarray): z-conponent displacement seismogram

    """

    # Local versions of global parameters
    dt = cf.dt
    nn = cf.nt
    slow = cf.slow
    baz = cf.baz
    nlay = cf.nlay
    h = cf.dp
    a0 = np.sqrt(cf.a[2,2,2,2,0])
    b0 = np.sqrt(cf.a[1,2,1,2,0])
    r0 = cf.rho[0]

    # Other local parameters
    ifree = 1
    iter = -1
    omg = 1.+0.001j
    om0 = 1.+0.j
    n2 = int(nn/2)

    # Frequency and time axes
    omega = np.arange(0,n2+1)*omg*2.*np.pi/(nn*dt)
    time = np.arange(0,nn*dt,dt)

    # Slowness vector
    psi = (baz - 180.)*np.pi/180.
    p1 = slow*1.e-3*np.cos(psi)
    p2 = slow*1.e-3*np.sin(psi)
    pp = p1*p1 + p2*p2

    # Properties of fluid layer above
    c = cf.c
    rhof = cf.rhof
    qf = np.sqrt(1/c**2 - pp)
    
    # Reflection matrix for water column
    rf = np.array([[-1.,0.,0.],[0.,0.,0.],[0.,0.,0.]])

    # Wave type
    if cf.wvtype=='P':
        const = np.array([0.,0.,0.,1.,0.,0.])
    elif cf.wvtype=='Si':
        const = np.array([0.,0.,0.,0.,0.5,0.5])
    elif cf.wvtype=='SV':
        const = np.array([0.,0.,0.,0.,1.,0.])
    elif cf.wvtype=='SH':
        const = np.array([0.,0.,0.,0.,0.,1.])
    else: 
        print('Error - wave type can only be P, Si, SV or SH')
        return

    # Produce R/T matrices for interfaces in solid
    Tus, Rus, Tds, Rds = rm.rmatrix(p1,p2,om0,nlay,ifree,iter)
    md = cf.evecs[0:3,0:3,0]
    mu = cf.evecs[0:3,3:6,0]
    nd = cf.evecs[3:6,0:3,0]
    nu = cf.evecs[3:6,3:6,0]
    F1 = np.vstack((np.hstack((md,mu)),np.hstack((nd,nu))))
    
    # Produce R/T matrices for ocean bottom
    TdPP, TdSP, RdPP, TuPP, TuPS, RuPP, RuSP, RuPS, RuSS = \
            rtf.rtfluid(p1, p2, c, rhof, a0, b0, r0)
    tu = np.array([[TuPP,TuPS,0.],[0.,0.,0.],[0.,0.,0.]])
    td = np.array([[TdPP,0.,0.],[TdSP,0.,0.],[0.,0.,0.]])
    ru = np.array([[RuPP,RuPS,0.],[RuSP,RuSS,0.],[0.,0.,1.]])
    rd = np.array([[RdPP,0.,0.],[0.,0.,0.],[0.,0.,0.]])

    # Initialize wavefield
    wave = np.zeros(6,dtype=complex)
    
    # Initialize displacement vector
    y = np.zeros((6,len(omega)),dtype=complex)
    
    a = pyfftw.empty_aligned(len(omega), dtype='float')

    # Setup loop for rmatrix
    for iw in range(0,int(n2+1)):

        # R/T matrices for solid
        Tus, Rus, Tds, Rds = rm.rmatrix(p1,p2,omega[iw],nlay,ifree,iw)
    
        # Phase income through ocean layer
        eu = np.exp(np.array([[1j*omega[iw]*qf*h,0.,0.],[0.,0.,0.],[0.,0.,0.]]))
        ed = eu
    
        # Upgoing wavefield on fluid side of OB
        Tup = np.dot(tu,solve(np.identity(3) - np.dot(Rds,ru),Tus))
        Rdp = rd + np.dot(tu,np.dot(Rds,(solve(np.identity(3) - 
            np.dot(ru,Rds),td))))
    
        # Upgoing wavefield on solid side of OB using wave propgator cast 
        # in terms of interfacial R/T coefficients.
        Ruh = ru + np.dot(td,np.dot(ed,np.dot(rf,np.dot(eu, 
            solve(np.identity(3) - 
                np.dot(rd,np.dot(ed,np.dot(rf,eu))),tu)))))

        wave[3:6] = np.dot(solve(np.identity(3) - 
            np.dot(Rds,Ruh),Tus),const[3:6])
        wave[0:3] = np.dot(Ruh,wave[3:6])

        # Displacement vector
        y[:,iw] = np.dot(F1,wave)

    # Get displacements in time domain
    ux = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[0,0:n2+1],\
            np.conj(np.flipud(y[0,1:n2])))))).transpose()
    uy = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[1,0:n2+1],\
            np.conj(np.flipud(y[1,1:n2])))))).transpose()
    uz = -np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[2,0:n2+1],\
            np.conj(np.flipud(y[2,1:n2])))))).transpose()

    # Return displacement stream
    return ux, uy, uz