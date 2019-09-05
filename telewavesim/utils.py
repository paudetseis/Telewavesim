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

Utility functions to interact with ``telewavesim`` modules.

'''
import sys
import itertools
import numpy as np
import pyfftw
from scipy.signal import hilbert
from obspy.core import Trace, Stream
from obspy.signal.rotate import rotate_ne_rt
from telewavesim import conf as cf 
from telewavesim import elast as es
from telewavesim.rmat_f import conf as cf_f
from telewavesim.rmat_f import plane as pw_f


def set_iso_tensor(a, b):
    """
    Function to generate tensor for isotropic material. 

    Args:
        a (float): P-wave velocity (km/s)
        b (float): S-wave velocity (km/s)

    Returns:
        (np.ndarray): cc: Elastic stiffness matrix in Voigt components (GPa) \
        (shape ``(3, 3, 3, 3)``)

    """

    a = a*1.e3
    b = b*1.e3
    C = es.iso_tensor(a, b)

    # Convert Voigt to full tensor
    cc = voigt2cc(C)

    return cc


def set_tri_tensor(a, b, tr, pl, ani):
    """
    Function to generate tensor for transverse isotropy. The
    tensor is rotated using the trend and plunge of the symmetry
    axis.

    Args:
        a (float): P-wave velocity (km/s)
        b (float): S-wave velocity (km/s)
        tr (float): Trend angle of symmetry axis (degree)
        pl (float): Plunge angle of symmetry axis (degree)
        ani (float): Percent anisotropy 

    Returns:
        (np.ndarray): cc: Elastic stiffness matrix in Voigt components (GPa) \
        (shape ``(3, 3, 3, 3)``)


    """

    # Trend and plunge of symmetry axis
    tr = -tr*np.pi/180.
    pl = (90. - pl)*np.pi/180.

    # Percent anisotropy
    da = (a*1.e3)*ani/100.
    db = (b*1.e3)*ani/100.
    
    # Set up matrix elements
    AA = (a*1.e3 - da/2.)**2
    CC = (a*1.e3 + da/2.)**2
    LL = (b*1.e3 + db/2.)**2
    NN = (b*1.e3 - db/2.)**2 
    AC = (a*1.e3)**2
    FF = -LL + np.sqrt((2.*AC)**2 - 2.*AC*(AA + CC + 2.*LL) + 
            (AA + LL)*(CC + LL))
    # eta = FF/(AA - 2.*LL)

    # Get tensor with horizontal axis
    cc = es.tri_tensor(AA, CC, FF, LL, NN)

    # Rotate tensor using trend and plunge
    cc = rot_tensor(cc, pl, tr, 0.)

    # Return tensor
    return cc


def set_aniso_tensor(tr, pl, typ='atg'):
    """
    Function to generate tensor for anisotropic minerals. The \
    tensor is rotated using the trend and plunge of the symmetry \
    axis.

    Args:
        tr (float): Trend angle of symmetry axis (degree)
        pl (float): Plunge angle of symmetry axis (degree)
        type (str, optional): Type of elastic material

    Returns:
        (tuple): Tuple containing:
            * cc (np.ndarray): Elastic stiffness matrix in Voigt components (GPa)\
            (shape ``(3, 3, 3, 3)``)
            * rho (float): Density (kg/m^3)

    """

    # Trend and plunge of symmetry axis
    tr = -tr*np.pi/180.
    pl = (90. - pl)*np.pi/180.

    # Get tensor with horizontal axis

    # Minerals
    if typ=='atg':
        C, rho = es.antigorite()
    elif typ=='bt':
        C, rho = es.biotite()
    elif typ=='cpx':
        C, rho = es.clinopyroxene_92()
    elif typ=='dol':
        C, rho = es.dolomite()
    elif typ=='ep':
        C, rho = es.epidote()
    elif typ=='grt':
        C, rho = es.garnet()
    elif typ=='gln':
        C, rho = es.glaucophane()
    elif typ=='hbl':
        C, rho = es.hornblende()
    elif typ=='jade':
        C, rho = es.jadeite()
    elif typ=='lws':
        C, rho = es.lawsonite()
    elif typ=='lz':
        C, rho = es.lizardite()
    elif typ=='ms':
        C, rho = es.muscovite()
    elif typ=='ol':
        C, rho = es.olivine()
    elif typ=='opx':
        C, rho = es.orthopyroxene()
    elif typ=='plag':
        C, rho = es.plagioclase_06()
    elif typ=='qtz':
        C, rho = es.quartz()
    elif typ=='zo':
        C, rho = es.zoisite()

    # Rocks
    elif typ=='BS_f':
        C, rho = es.blueschist_felsic()
    elif typ=='BS_m':
        C, rho = es.blueschist_mafic()
    elif typ=='EC_f':
        C, rho = es.eclogite_foliated()
    elif typ=='EC_m':
        C, rho = es.eclogite_massive()
    elif typ=='HB':
        C, rho = es.harzburgite()
    elif typ=='SP_37':
        C, rho = es.serpentinite_37()
    elif typ=='SP_80':
        C, rho = es.serpentinite_80()
    elif typ=='LHZ':
        C, rho = es.lherzolite()

    else:
        print('type of mineral/rock not implemented')
        return

    # Convert Voigt to full tensor
    cc = voigt2cc(C)*1.e9/rho

    # Rotate tensor using trend and plunge
    cc = rot_tensor(cc, pl, tr, 0.)

    # Return tensor
    return cc, rho


def full_3x3_to_Voigt_6_index(i, j):
    """
    Conversion of tensor to Voigt notation for indices

    """
    if i == j:
        return i
    return 6-i-j

def voigt2cc(C):
    """
    Convert the Voigt representation of the stiffness matrix to the full
    3x3x3x3 representation.

    Args:
        C (np.ndarray): Stiffness matrix (shape ``(6, 6)``)

    Returns:
        (np.ndarray): cc: Elastic tensor (shape ``(3, 3, 3, 3)``)

    """
    
    C = np.asarray(C)
    cc = np.zeros((3,3,3,3), dtype=float)
    for i, j, k, l in itertools.product(range(3), range(3), range(3), range(3)):
        Voigt_i = full_3x3_to_Voigt_6_index(i, j)
        Voigt_j = full_3x3_to_Voigt_6_index(k, l)
        cc[i, j, k, l] = C[Voigt_i, Voigt_j]
    return cc


def cc2voigt(cc):
    """
    Convert from the full 3x3x3x3 representation of the stiffness matrix
    to the representation in Voigt notation. 

    Args:
        cc (np.ndarray): Elastic tensor (shape ``(3, 3, 3, 3)``)

    Returns:
        (np.ndarray): C: Stiffness matrix (shape ``(6, 6)``)

    """

    Voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

    tol = 1e-3

    cc = np.asarray(cc)
    C = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            C[i,j] = cc[k,l,m,n]

    return C


def VRH_average(C):
    """
    Performs a Voigt-Reuss-Hill average of the anisotropic
    stifness matrix to the bulk modulus K and the shear modulus
    G.

    Args:
        C (np.ndarray): Stiffness matrix (shape ``(6, 6)``)

    Returns:     
        (tuple): Tuple containing: 
            * Kvoigt (float): Voigt average bulk modulus (GPa)
            * Gvoigt (float): Voigt average shear modulus (GPa)
            * Kreuss (float): Reuss average bulk modulus (GPa)
            * Greuss (float): Reuss average shear modulus (GPa)
            * Kvrh (float):   Voigt-Reuss-Hill average bulk modulus (GPa)
            * Gvrh (float):   Voigt-Reuss-Hill average shear modulus (GPa)

    """

    # Compliance matrix
    S = np.linalg.inv(C)

    # Voigt averaging
    Kvoigt = (C[0,0] + C[1,1] + C[2,2] + 2.*C[0,1] + 2.*C[0,2] + 2.*C[1,2])/9.
    Gvoigt = (C[0,0] + C[1,1] + C[2,2] - C[0,1] - C[0,2] - C[1,2] + 3.*C[3,3] + \
        3.*C[4,4] + 3.*C[5,5])/15.

    # Reuss averaging
    Kreuss = 1./(S[0,0] + S[1,1] + S[2,2] + 2.*S[0,1] + 2.*S[0,2] + 2.*S[1,2])
    Greuss = 15./(4.*S[0,0] + 4.*S[1,1] + 4.*S[2,2] - 4.*S[0,1] - 4.*S[0,2] - \
        4.*S[1,2] + 3.*S[3,3] + 3.*S[4,4] + 3.*S[5,5])

    # Voigt-Reuss-Hill average
    Kvrh = (Kvoigt + Kreuss)/2.
    Gvrh = (Gvoigt + Greuss)/2.

    return Kvoigt, Gvoigt, Kreuss, Greuss, Kvrh, Gvrh


def mod2vel(K,G,rho):
    """

    Calculates the isotropic P and S wave velocities from given
    bulk (K) and shear (G) moduli and density (rho) in kg/m^3

    Args:
        K (float): Bulk modulus (GPa)
        G (float): Shear modulus (GPa)
        rho (float): Density (kg/m^3)
    Returns:
        (tuple): tuple containing:

            * Vp (float): P-wave velocity (km/s)
            * Vs (float): S-wave velocity (km/s)
    """

    Vp = np.sqrt((K + 4.*G/3.)/rho)
    Vs = np.sqrt(G/rho)

    return Vp, Vs


def rot_tensor(a,alpha,beta,gam):
    """

    Performs a rotation of the tensor cc (c_ijkl) about three angles (alpha, 
    beta, gamma)

    Args:
        a (np.ndarray): Elastic tensor with shape ``(3, 3, 3, 3)``
        alpha (float): Angle in radians
        beta (float): Angle in radians
        gam (float): Angle in radians

    Returns:
        (np.ndarray): aa: Rotated tensor with shape ``(3, 3, 3, 3)``

    .. note::

        The three angles (``alpha``, ``beta``, ``gam``) correspond to rotation about the 
        x_2, x_3, x_1 axes. Note that the sequence of the rotation is important: 
        (AB ~= BA). In this case we rotate about x_2 first, x_3 second and x_1 third. 
        
        For trend and plunge of symmetry axis (e.g., tri_tensor): 
            
            ``alpha`` = plunge

            ``beta`` = trend

    """

    rot = np.zeros((3,3))
    aa = np.zeros((3,3,3,3))

    rot[0,0] = np.cos(alpha)*np.cos(beta)
    rot[0,1] = np.sin(beta)
    rot[0,2] = np.sin(alpha)*np.cos(beta)

    rot[1,0] = -np.cos(gam)*np.sin(beta)*np.cos(alpha) - \
            np.sin(gam)*np.sin(alpha)
    rot[1,1] = np.cos(gam)*np.cos(beta)
    rot[1,2] = -np.cos(gam)*np.sin(beta)*np.sin(alpha) + \
            np.sin(gam)*np.cos(alpha)
    rot[2,0] = np.sin(gam)*np.sin(beta)*np.cos(alpha) - \
            np.cos(gam)*np.sin(alpha)
    rot[2,1] = -np.sin(gam)*np.cos(beta)
    rot[2,2] = np.sin(gam)*np.sin(beta)*np.sin(alpha) + \
            np.cos(gam)*np.cos(alpha)

    #
    #  c_ijkl ---> c_mnrs
    #
    for m in range(3):
        for n in range(3):
            for r in range(3):
                for s in range(3):
                    asum=0.0
                    for i in range(3):
                        for j in range(3):
                            for k in range(3):
                                for l in range(3):
                                    rr = rot[m,i]*rot[n,j]*rot[r,k]*rot[s,l]
                                    asum = asum + rr*a[i,j,k,l]
                    aa[m,n,r,s] = asum

    return aa


def rotate_zrt_pvh(trZ, trR, trT, vp=6., vs=3.5):
    """
    Rotates traces from `Z-R-T` orientation to `P-SV-SH` wave mode.

    Args:
        trZ (obspy.trace): Vertical component
        trR (obspy.trace): Radial component
        trT (obspy.trace): Transverse component
        vp (float, optional): P-wave velocity used for rotation
        vs (float, optional): S-wave velocity used for rotation

    Returns:
        (tuple): tuple containing:

            * trP (obspy.trace): Compressional (P) wave mode
            * trV (obspy.trace): Vertically polarized shear (SV) wave mode
            * trH (obspy.trace): Horizontally polarized shear (SH) wave mode

    """
    # Copy traces
    trP = trZ.copy()
    trV = trR.copy()
    trH = trT.copy()

    # Vertical slownesses
    qp = np.sqrt(1/vp/vp - cf.slow*cf.slow) 
    qs = np.sqrt(1/vs/vs - cf.slow*cf.slow) 

    # Elements of rotation matrix
    m11 = cf.slow*vs*vs/vp
    m12 = -(1 - 2*vs*vs*cf.slow*cf.slow)/(2*vp*qp)
    m21 = (1 - 2*vs*vs*cf.slow*cf.slow)/(2*vs*qs)
    m22 = cf.slow*vs

    # Rotation matrix
    rot = np.array([[-m11, m12], [-m21, m22]])
    
    # Vector of Radial and Vertical
    r_z = np.array([trR.data,trZ.data])
                
    # Rotation
    vec = np.dot(rot, r_z)
                
    # Extract P and SV components
    trP.data = vec[0,:]
    trV.data = vec[1,:]
    trH.data = -trT.data/2.
                
    return trP, trV, trH


def stack_all(st1, st2, pws=False):
    """
    Stacks all traces in two ``Stream`` objects.

    Args:
        st1 (obspy.stream): Stream 1
        st2 (obspy.stream,): Stream 2
        pws (bool, optional): Enables Phase-Weighted Stacking

    Returns:
        (tuple): tuple containing:

            * stack1 (obspy.trace): Stacked trace for Stream 1
            * stack2 (obspy.trace): Stacked trace for Stream 2

    """

    print()
    print('Stacking ALL traces in streams')

    # Copy stats from stream
    str_stats = st1[0].stats

    # Initialize arrays
    tmp1 = np.zeros(len(st1[0].data))
    tmp2 = np.zeros(len(st2[0].data))
    weight1 = np.zeros(len(st1[0].data), dtype=complex)
    weight2 = np.zeros(len(st2[0].data), dtype=complex)

    # Stack all traces
    for tr in st1:
        tmp1 += tr.data
        hilb1 = hilbert(tr.data)
        phase1 = np.arctan2(hilb1.imag, hilb1.real)
        weight1 += np.exp(1j*phase1)

    for tr in st2:
        tmp2 += tr.data
        hilb2 = hilbert(tr.data)
        phase2 = np.arctan2(hilb2.imag, hilb2.real)
        weight2 += np.exp(1j*phase2)

    # Normalize
    tmp1 = tmp1/np.float(len(st1))
    tmp2 = tmp2/np.float(len(st2))

    # Phase-weighting
    if pws:
        weight1 = weight1/np.float(len(st1))
        weight2 = weight2/np.float(len(st2))
        weight1 = np.real(abs(weight1))
        weight2 = np.real(abs(weight2))
    else:
        weight1 = np.ones(len(st1[0].data))
        weight2 = np.ones(len(st1[0].data))

    # Put back into traces
    stack1 = Trace(data=weight1*tmp1,header=str_stats)
    stack2 = Trace(data=weight2*tmp2,header=str_stats)

    return stack1, stack2


def calc_ttime(slow):
    """
    Calculates total propagation time through model.

    Args:
        slow (float): Slowness value (s/km)

    Returns:
        (float): t1: Time in seconds
    """


    t1 = 0.

    for i in range(cf.nlay-1):
        if cf.isoflg[i] == 'iso':    
            a0 = cf.a[2,2,2,2,i]
            b0 = cf.a[1,2,1,2,i]
        else:
            cc = cc2voigt(cf.a[:,:,:,:,i])
            rho = cf.rho[i]
            K1,G1,K2,G2,K,G = VRH_average(cc*rho)
            a0, b0 = mod2vel(K,G,rho)
            a0 = a0**2
            b0 = b0**2
        if cf.wvtype=='P':
            t1 += cf.thickn[i]*np.sqrt(1./a0 - (slow*1.e-3)**2)
        elif cf.wvtype=='Si' or cf.wvtype=='SV' or cf.wvtype=='SH':
            t1 += cf.thickn[i]*np.sqrt(1./b0 - (slow*1.e-3)**2)

    return t1


def read_model(modfile):
    """
    Reads model parameters from file that are passed 
        through the configuration module ``conf``. 

    Returns: 
        None: Parameters are now global variables shared
                between all other modules

    """

    h = []; r = []; a = []; b = []; fl = []; ani = []; tr = []; pl = []

    # Read file line by line and populate lists
    try:
        open(modfile)
    except:
        raise(Exception('model file cannot be opened: ',modfile))

    with open(modfile) as fileobj:
        for line in fileobj:
            if not line.rstrip().startswith('#'):
                model = line.rstrip().split()
                h.append(np.float64(model[0])*1.e3)
                r.append(np.float64(model[1]))
                a.append(np.float64(model[2]))
                b.append(np.float64(model[3]))
                fl.append(model[4])
                ani.append(np.float64(model[5]))
                tr.append(np.float64(model[6]))
                pl.append(np.float64(model[7]))

    # Pass configuration parameters
    cf.nlay = len(h)
    cf.thickn = h
    cf.rho = r
    cf.isoflg = fl

    cf.a = np.zeros((3,3,3,3,cf.nlay))
    cf.evecs = np.zeros((6,6,cf.nlay),dtype=complex)
    cf.evals = np.zeros((6,cf.nlay),dtype=complex)
    cf.Tui = np.zeros((3,3,cf.nlay),dtype=complex)
    cf.Rui = np.zeros((3,3,cf.nlay),dtype=complex)
    cf.Tdi = np.zeros((3,3,cf.nlay),dtype=complex)
    cf.Rdi = np.zeros((3,3,cf.nlay),dtype=complex)

    mins = ['atg', 'bt', 'cpx', 'dol', 'ep', 'grt', 'gln', 'hbl', 'jade',\
            'lws', 'lz', 'ms', 'ol', 'opx', 'plag', 'qtz', 'zo']

    rocks = ['BS_f', 'BS_m', 'EC_f', 'EC_m', 'HB', 'LHZ', 'SP_37', 'SP_80']

    for j in range(cf.nlay):
        if fl[j]=='iso':
            cc = set_iso_tensor(a[j],b[j])
            cf.a[:,:,:,:,j] = cc
        elif fl[j]=='tri':
            cc = set_tri_tensor(a[j],b[j],tr[j],pl[j],ani[j])
            cf.a[:,:,:,:,j] = cc
        elif fl[j] in mins or fl[j] in rocks:
            cc, rho = set_aniso_tensor(tr[j],pl[j],typ=fl[j])
            cf.a[:,:,:,:,j] = cc
            cf.rho[j] = rho
        else:
            print('\nFlag not defined: use either "iso", "tri" or one among\n')
            print(mins,rocks)
            print()
            raise(Exception())

    return


def check_cf(obs=False):
    """
    Checks whether or not all required global variables are set and throws an Exception if not.

    Args:
        obs (bool, optional): Whether the analysis is done for an OBS case or not.

    :raises ExceptionError: Throws ExceptionError if not all variables are set.
    """
    lst = [cf.a, cf.rho, cf.thickn, cf.isoflg, cf.dt, cf.nt, cf.slow, cf.baz]
    check = [f is None for f in lst]
    if sum(check)/len(check)>0.:
        raise Exception("global variables not all set. Set all of the following variables through the conf module: 'a', 'rho', 'thickn', 'isoflg', 'dt', 'nt', 'slow', 'baz'")

    if obs:
        lst = [cf.dp, cf.c, cf.rhof]
        check = [f is None for f in lst]
        if sum(check)/len(check)>0.:
            raise Exception("global variables not all set for OBS case. Set all of the following variables through the conf module: 'dp', 'c', 'rhof'")


def model2for():
    """
    Passes global model variables to Fortran ``conf`` module.

    Returns:
        None

    Variables to pass are ``a``, ``rho``, ``thickn``, ``isoflg``
    """

    nlaymx = cf_f.nlaymx
    cf_f.a = np.zeros((3,3,3,3,nlaymx))
    cf_f.rho = np.zeros((nlaymx))
    cf_f.thickn = np.zeros((nlaymx))
    cf_f.isoflg = np.zeros((nlaymx), dtype='int')

    for i in range(cf.nlay):
        cf_f.a[:,:,:,:,i] = cf.a[:,:,:,:,i]
        cf_f.rho[i] = cf.rho[i]
        cf_f.thickn[i] = cf.thickn[i]
        if cf.isoflg[i]=='iso':
            cf_f.isoflg[i] = 1


def wave2for():
    """
    Passes global wavefield variables to Fortran ``conf`` module.

    Returns:
        None

    Variables to pass are ``dt``, ``slow``, ``baz``
    """

    cf_f.dt = cf.dt
    cf_f.slow = cf.slow
    cf_f.baz = cf.baz


def obs2for():
    """
    Passes global OBS-related variables to Fortran ``conf`` module.

    Returns:
        None

    Variables to pass are ``dp``, ``c``, ``rhof``
    """
    cf_f.dp = cf.dp
    cf_f.c = cf.c
    cf_f.rhof = cf.rhof


def run_plane(obs=False):
    """
    Function to run the ``plane`` module and return 3-component seismograms as an ``obspy``
    ``Stream`` object. 

    Args:
        fortran (book, option): Whether or not the Fortran modules are used
        obs (bool, optional): Whether or not the analysis is done for an OBS stations

    Returns:
        (obspy.stream): trxyz: Stream containing 3-component displacement seismograms

    """

    # Check if all variables are set. If not, throw an Exception and stop
    check_cf(obs)

    # Pass  variables to Fortran conf
    model2for()
    wave2for()

    # Run the ``plane`` module depending on land or OBS case.
    if obs:

        # If OBS, then further pass OBS-related paramters to Fortran conf
        obs2for()

        # Get the Fourier transform of seismograms for ``obs``case
        yx, yy, yz = pw_f.plane_obs(cf.nt,cf.nlay,np.array(cf.wvtype, dtype='c'))

    else:

        # Get the Fourier transform of seismograms for ``land`` case
        yx, yy, yz = pw_f.plane_land(cf.nt,cf.nlay,np.array(cf.wvtype, dtype='c'))

    # Transfer displacement seismograms to an ``obspy`` ``Stream`` object.
    trxyz = get_trxyz(yx, yy, yz)

    return trxyz


def get_trxyz(yx, yy, yz):
    """
    Function to store displacement seismograms into ``obspy`` ``Trace`` obsjects and 
    then an ``obspy`` ``Stream`` object. 

    Args:
        ux (np.ndarray): x-component displacement seismogram
        uy (np.ndarray): y-component displacement seismogram
        uz (np.ndarray): z-component displacement seismogram

    Returns:
        (obspy.stream): trxyz: Stream containing 3-component displacement seismograms

    """

    # Get displacements in time domain
    ux = np.real(pyfftw.interfaces.numpy_fft.fft(yx))
    uy = np.real(pyfftw.interfaces.numpy_fft.fft(yy))
    uz = -np.real(pyfftw.interfaces.numpy_fft.fft(yz))

    # Store in traces
    tux = Trace(data=ux)
    tuy = Trace(data=uy)
    tuz = Trace(data=uz)

    # Update trace header
    tux = update_stats(tux, cf.nt, cf.dt, cf.slow, cf.baz)
    tuy = update_stats(tuy, cf.nt, cf.dt, cf.slow, cf.baz)
    tuz = update_stats(tuz, cf.nt, cf.dt, cf.slow, cf.baz)

    # Append to stream
    trxyz = Stream(traces=[tux, tuy, tuz])

    return trxyz


def tf_from_xyz(trxyz, pvh=False):
    """
    Function to generate transfer functions from displacement traces. 

    Args:
        trxyz (obspy.stream): Obspy ``Stream`` object in cartesian coordinate system
        pvh (bool, optional): Whether to rotate from Z-R-T coordinate system to P-SV-SH wave mode

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
        trP, trV, trH = rotate_zrt_pvh(ztr, rtr, ttr, vp=vp, vs=vs)
        
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


def update_stats(tr, nt, dt, slow, baz):
    """
    Updates the ``stats`` doctionary from an obspy ``Trace`` object.

    Args:
        tr (obspy.trace): Trace object to update
        nt (int): Number of samples
        dt (float): Sampling rate
        slow (float): Slowness value (s/km)
        baz (float): Back-azimuth value (degree)

    Returns:
        (obspy.trace): tr: Trace with updated stats
    """

    tr.stats.delta = dt
    tr.stats.slow = slow
    tr.stats.baz = baz

    return tr