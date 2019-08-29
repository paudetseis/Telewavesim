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
from obspy.signal.rotate import rotate_ne_rt


def plane_land():
    """
    Function to generate plane wave seismograms from a stack of layers
    for land surface stations. Also handles anisotropy. All
    model and time series properties are passed through the 
    configuration module ``conf``. Handles any number of layers.

    Returns:
        trxyz (obspy.stream): ObsPy ``Stream`` containing displacement traces for
                given model and slowness vector

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
    ux_time = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[0,0:nend+1],\
            np.conj(np.flipud(y[0,1:nend])))))).transpose()
    uy_time = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[1,0:nend+1],\
            np.conj(np.flipud(y[1,1:nend])))))).transpose()
    uz_time = -np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[2,0:nend+1],\
            np.conj(np.flipud(y[2,1:nend])))))).transpose()

    # Store in traces
    tux = Trace(data=ux_time)
    tuy = Trace(data=uy_time)
    tuz = Trace(data=uz_time)

    # Update trace header
    tux = ut.update_stats(tux, n2, dt, slow, baz)
    tuy = ut.update_stats(tuy, n2, dt, slow, baz)
    tuz = ut.update_stats(tuz, n2, dt, slow, baz)

    # Append to stream
    trxyz = Stream(traces=[tux, tuy, tuz])

    # Return displacement stream
    return trxyz


def plane_obs():
    """
    Function to generate plane wave seismograms from a stack of layers
    for ocean-bottom stations. Also handles anisotropy. All
    model and time series properties are passed through the 
    configuration module ``conf``. Handles any number of layers.

    Returns:
        trxyz (obspy.stream): ObsPy ``Stream`` containing displacement traces for
                given model and slowness vector

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
    ux_time = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[0,0:n2+1],\
            np.conj(np.flipud(y[0,1:n2])))))).transpose()
    uy_time = np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[1,0:n2+1],\
            np.conj(np.flipud(y[1,1:n2])))))).transpose()
    uz_time = -np.real(pyfftw.interfaces.numpy_fft.fft(np.hstack((y[2,0:n2+1],\
            np.conj(np.flipud(y[2,1:n2])))))).transpose()

    # Store in traces
    tux = Trace(data=ux_time)
    tuy = Trace(data=uy_time)
    tuz = Trace(data=uz_time)

    # Update trace header
    tux = ut.update_stats(tux, nn, dt, slow, baz)
    tuy = ut.update_stats(tuy, nn, dt, slow, baz)
    tuz = ut.update_stats(tuz, nn, dt, slow, baz)

    # Append to stream
    trxyz = Stream(traces=[tux, tuy, tuz])

    # Return displacement stream
    return trxyz
   

def tf_from_xyz(trxyz,pvh=False):
    """
    Function to generate transfer functions from displacement traces. 

    Args:
        trxyz (obspy.stream): Obspy ``Stream`` object in cartesian coordinate system

    Returns:
        (obspy.stream): tfs: Stream containing Radial and Transverse transfer functions

    """

    # Extract East, North and Vertical
    ntr = trxyz[0]
    etr = trxyz[1]
    ztr = trxyz[2]
    baz = cf.baz

    # Copy to radial and transverse
    rtr = ntr.copy()
    ttr = etr.copy()

    # Rotate to radial and transverse
    rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)
    a = pyfftw.empty_aligned(len(rtr.data), dtype='float')
    # print(rtr.data, ttr.data)

    if pvh:
        vp = np.sqrt(cf.a[2,2,2,2,0])/1.e3
        vs = np.sqrt(cf.a[1,2,1,2,0])/1.e3
        trP, trV, trH = tru.rotate_zrt_pvh(ztr, rtr, ttr, vp=vp, vs=vs)
        
        tfr = trV.copy(); tfr.data = np.zeros(len(tfr.data))
        tft = trH.copy(); tft.data = np.zeros(len(tft.data))
        ftfv = pyfftw.interfaces.numpy_fft.fft(trV.data)
        ftfh = pyfftw.interfaces.numpy_fft.fft(trH.data)
        ftfp = pyfftw.interfaces.numpy_fft.fft(trP.data)

        if cf.wvtype=='P':
            # Transfer function
            tfr.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(ftfv,ftfp))))
            tft.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(ftfh,ftfp))))
        elif cf.wvtype=='Si': 
            tfr.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfp,ftfv))))
            tft.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfp,ftfh))))
        elif cf.wvtype=='SV':
            tfr.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfp,ftfv))))
        elif cf.wvtype=='SH':
            tft.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfp,ftfh))))
    else:
        tfr = rtr.copy(); tfr.data = np.zeros(len(tfr.data))
        tft = ttr.copy(); tft.data = np.zeros(len(tft.data))
        ftfr = pyfftw.interfaces.numpy_fft.fft(rtr.data)
        ftft = pyfftw.interfaces.numpy_fft.fft(ttr.data)
        ftfz = pyfftw.interfaces.numpy_fft.fft(ztr.data)

        if cf.wvtype=='P':
            # Transfer function
            tfr.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(ftfr,ftfz))))
            tft.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(ftft,ftfz))))
        elif cf.wvtype=='Si':
            tfr.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfz,ftfr))))
            tft.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfz,ftft))))
        elif cf.wvtype=='SV':
            tfr.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfz,ftfr))))
        elif cf.wvtype=='SH':
            tft.data = np.fft.fftshift(np.real(pyfftw.interfaces.numpy_fft.ifft(np.divide(-ftfz,ftft))))

    # Store in stream
    tfs = Stream(traces=[tfr, tft])

    # Return stream
    return tfs


