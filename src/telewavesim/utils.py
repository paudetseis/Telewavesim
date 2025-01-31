# Copyright 2019 Pascal Audet, Tom Eulenfeld

# This file is part of Telewavesim.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

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
import itertools
import numpy as np
from numpy import sin, cos
from numpy.fft import fft, fftshift, ifft
from scipy.signal import hilbert
from obspy.core import Trace, Stream
from obspy.signal.rotate import rotate_ne_rt
from telewavesim import elast as es
from telewavesim.rmat_f import conf as cf_f
from telewavesim.rmat_f import plane as pw_f


MINERALS = ['atg', 'bt', 'cpx', 'dol', 'ep', 'grt', 'gln', 'hbl', 'jade',
            'lws', 'lz', 'ms', 'ol', 'opx', 'plag', 'qtz', 'zo']

ROCKS = ['BS_f', 'BS_m', 'EC_f', 'EC_m', 'HB', 'LHZ', 'SP_37', 'SP_80']


def set_iso_tensor(a, b):
    """
    Function to generate tensor for isotropic material.

    Args:
        a (float): P-wave velocity (km/s)
        b (float): S-wave velocity (km/s)

    Returns:
        (np.ndarray): cc: Elastic tensor (GPa /density) \
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
        (np.ndarray): cc: Elastic tensor (GPa / kg/m^3) \
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
    cc = np.zeros((3, 3, 3, 3))

    cc[0, 0, 0, 0] = AA
    cc[1, 1, 1, 1] = AA
    cc[2, 2, 2, 2] = CC

    cc[0, 0, 1, 1] = (AA - 2.*NN)
    cc[1, 1, 0, 0] = (AA - 2.*NN)

    cc[0, 0, 2, 2] = FF
    cc[2, 2, 0, 0] = FF

    cc[1, 1, 2, 2] = FF
    cc[2, 2, 1, 1] = FF

    cc[1, 2, 1, 2] = LL
    cc[2, 1, 2, 1] = LL
    cc[2, 1, 1, 2] = LL
    cc[1, 2, 2, 1] = LL

    cc[2, 0, 2, 0] = LL
    cc[0, 2, 0, 2] = LL
    cc[2, 0, 0, 2] = LL
    cc[0, 2, 2, 0] = LL

    cc[0, 1, 0, 1] = NN
    cc[1, 0, 1, 0] = NN
    cc[0, 1, 1, 0] = NN
    cc[1, 0, 0, 1] = NN

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
            * cc (np.ndarray): Elastic tensor (GPa / kg/m^3)\
            (shape ``(3, 3, 3, 3)``)
            * rho (float): Density (kg/m^3)

    """

    # Trend and plunge of symmetry axis
    tr = -tr*np.pi/180.
    pl = (90. - pl)*np.pi/180.

    # Get tensor with horizontal axis

    # Minerals
    if typ == 'atg':
        C, rho = es.antigorite()
    elif typ == 'bt':
        C, rho = es.biotite()
    elif typ == 'cpx':
        C, rho = es.clinopyroxene_92()
    elif typ == 'dol':
        C, rho = es.dolomite()
    elif typ == 'ep':
        C, rho = es.epidote()
    elif typ == 'grt':
        C, rho = es.garnet()
    elif typ == 'gln':
        C, rho = es.glaucophane()
    elif typ == 'hbl':
        C, rho = es.hornblende()
    elif typ == 'jade':
        C, rho = es.jadeite()
    elif typ == 'lws':
        C, rho = es.lawsonite()
    elif typ == 'lz':
        C, rho = es.lizardite()
    elif typ == 'ms':
        C, rho = es.muscovite()
    elif typ == 'ol':
        C, rho = es.olivine()
    elif typ == 'opx':
        C, rho = es.orthopyroxene()
    elif typ == 'plag':
        C, rho = es.plagioclase_06()
    elif typ == 'qtz':
        C, rho = es.quartz()
    elif typ == 'zo':
        C, rho = es.zoisite()

    # Rocks
    elif typ == 'BS_f':
        C, rho = es.blueschist_felsic()
    elif typ == 'BS_m':
        C, rho = es.blueschist_mafic()
    elif typ == 'EC_f':
        C, rho = es.eclogite_foliated()
    elif typ == 'EC_m':
        C, rho = es.eclogite_massive()
    elif typ == 'HB':
        C, rho = es.harzburgite()
    elif typ == 'SP_37':
        C, rho = es.serpentinite_37()
    elif typ == 'SP_80':
        C, rho = es.serpentinite_80()
    elif typ == 'LHZ':
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
    3x3x3x3 tensor representation.

    Args:
        C (np.ndarray): Stiffness matrix (shape ``(6, 6)``)

    Returns:
        (np.ndarray): cc: Elastic tensor (shape ``(3, 3, 3, 3)``)

    """

    C = np.asarray(C)
    cc = np.zeros((3, 3, 3, 3), dtype=float)
    for i, j, k, l in itertools.product(range(3), repeat=4):
        Voigt_i = full_3x3_to_Voigt_6_index(i, j)
        Voigt_j = full_3x3_to_Voigt_6_index(k, l)
        cc[i, j, k, l] = C[Voigt_i, Voigt_j]
    return cc


def cc2voigt(cc):
    """
    Convert from the full 3x3x3x3 tensor representation
    to the Voigt notation of the stiffness matrix.

    Args:
        cc (np.ndarray): Elastic tensor (shape ``(3, 3, 3, 3)``)

    Returns:
        (np.ndarray): C: Stiffness matrix (shape ``(6, 6)``)

    """

    Voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

    cc = np.asarray(cc)
    C = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            C[i, j] = cc[k, l, m, n]

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

    Example
    -------
    >>> from telewavesim import utils
    >>> cc, rho = utils.set_aniso_tensor(0., 0., typ='atg')
    >>> C = utils.cc2voigt(cc)
    >>> utils.VRH_average(C*rho)
    (75655555555.555557, 48113333333.333336, 61245706544.967415,
    28835098086.844658, 68450631050.26149, 38474215710.088997)

    """

    # Compliance matrix
    S = np.linalg.inv(C)

    # Voigt averaging
    Kvoigt = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * C[0, 1] + 2 * C[0, 2] +
              2 * C[1, 2]) / 9
    Gvoigt = (C[0, 0] + C[1, 1] + C[2, 2] - C[0, 1] - C[0, 2] - C[1, 2] +
              3 * C[3, 3] + 3 * C[4, 4] + 3 * C[5, 5]) / 15

    # Reuss averaging
    Kreuss = 1 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * S[0, 1] + 2 * S[0, 2] +
                  2 * S[1, 2])
    Greuss = 15 / (4 * S[0, 0] + 4 * S[1, 1] + 4 * S[2, 2] - 4 * S[0, 1] -
                   4 * S[0, 2] - 4 * S[1, 2] + 3 * S[3, 3] + 3 * S[4, 4] +
                   3 * S[5, 5])

    # Voigt-Reuss-Hill average
    Kvrh = (Kvoigt + Kreuss) / 2
    Gvrh = (Gvoigt + Greuss) / 2

    return Kvoigt, Gvoigt, Kreuss, Greuss, Kvrh, Gvrh


def mod2vel(K, G, rho):
    """

    Calculates the isotropic P and S wave velocities from given
    bulk (K) and shear (G) moduli and density (rho) in kg/m^3

    Args:
        K (float): Bulk modulus (GPa)
        G (float): Shear modulus (GPa)
        rho (float): Density (kg/m^3)
    Returns:
        (tuple): tuple containing:

            * Vp (float): P-wave velocity (m/s)
            * Vs (float): S-wave velocity (m/s)

    Example
    -------
    >>> from telewavesim import utils
    >>> cc, rho = utils.set_aniso_tensor(0., 0., typ='atg')
    >>> C = utils.cc2voigt(cc)
    >>> K, G = utils.VRH_average(C*rho)[4:6]
    >>> utils.mod2vel(K, G, rho)
    (6760.617471753726, 3832.0771334254896)

    """

    Vp = np.sqrt((K + 4.*G/3.)/rho)
    Vs = np.sqrt(G/rho)

    return Vp, Vs


def rot_tensor(a, alpha, beta, gam):
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

        The three angles (``alpha``, ``beta``, ``gam``) correspond to rotation
        about the x_2, x_3, x_1 axes. Note that the sequence of the rotation
        is important: (AB ~= BA). In this case we rotate about x_2 first,
        x_3 second and x_1 third.

        For trend and plunge of symmetry axis (e.g., tri_tensor):

            ``alpha`` = plunge

            ``beta`` = trend

    """

    rot = np.zeros((3, 3))
    aa = np.zeros((3, 3, 3, 3))

    rot[0, 0] = cos(alpha) * cos(beta)
    rot[0, 1] = sin(beta)
    rot[0, 2] = sin(alpha) * cos(beta)

    rot[1, 0] = -cos(gam) * sin(beta) * cos(alpha) - sin(gam) * sin(alpha)
    rot[1, 1] = +cos(gam) * cos(beta)
    rot[1, 2] = -cos(gam) * sin(beta) * sin(alpha) + sin(gam) * cos(alpha)
    rot[2, 0] = +sin(gam) * sin(beta) * cos(alpha) - cos(gam) * sin(alpha)
    rot[2, 1] = -sin(gam) * cos(beta)
    rot[2, 2] = +sin(gam) * sin(beta) * sin(alpha) + cos(gam) * cos(alpha)

    #
    #  c_ijkl ---> c_mnrs
    #
    for m, n, r, s in itertools.product(range(3), repeat=4):
        asum = 0.0
        for i, j, k, l in itertools.product(range(3), repeat=4):
            rr = rot[m, i] * rot[n, j] * rot[r, k] * rot[s, l]
            asum = asum + rr * a[i, j, k, l]
        aa[m, n, r, s] = asum

    return aa


def rotate_zrt_pvh(trZ, trR, trT, slow, vp=None, vs=None):
    """
    Rotates traces from `Z-R-T` orientation to `P-SV-SH` wave mode.

    Args:
        trZ (obspy.trace): Vertical component
        trR (obspy.trace): Radial component
        trT (obspy.trace): Transverse component
        slow (float): slowness of wave
        vp (float, optional): P-wave velocity used for rotation
        vs (float, optional): S-wave velocity used for rotation

    Returns:
        (tuple): tuple containing:

            * trP (obspy.trace): Compressional (P) wave mode
            * trV (obspy.trace): Vertically polarized shear (SV) wave mode
            * trH (obspy.trace): Horizontally polarized shear (SH) wave mode

    """
    if vp is None:
        vp = 6.0
    if vs is None:
        vs = 3.5
    # Copy traces
    trP = trZ.copy()
    trV = trR.copy()
    trH = trT.copy()

    # Vertical slownesses
    qp = np.sqrt(1/vp**2 - slow**2)
    qs = np.sqrt(1/vs**2 - slow**2)

    # Elements of rotation matrix
    m11 = slow*vs*vs/vp
    m12 = -(1 - 2*vs*vs*slow*slow)/(2*vp*qp)
    m21 = (1 - 2*vs*vs*slow*slow)/(2*vs*qs)
    m22 = slow*vs

    # Rotation matrix
    rot = np.array([[-m11, m12], [-m21, m22]])

    # Vector of Radial and Vertical
    r_z = np.array([trR.data, trZ.data])

    # Rotation
    vec = np.dot(rot, r_z)

    # Extract P and SV components
    trP.data = vec[0, :]
    trV.data = vec[1, :]
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
    stack1 = Trace(data=weight1*tmp1, header=str_stats)
    stack2 = Trace(data=weight2*tmp2, header=str_stats)

    return stack1, stack2


def calc_ttime(model, slow, wvtype='P'):
    """
    Calculates total propagation time through model given the corresponding
    ``'P'`` or ``'S'`` wave type. The bottom layer is irrelevant in this
    calculation. All ``'S'`` wave types will return the same predicted time.
    This function is useful mainly for plotting purposes. For example, to show
    the first phase arrival at a time of 0, the traces can be shifted by the
    total propagation time through the model.

    Args:
        model (Model): Model object
        slow (float): Slowness value (s/km)
        wvtype (str): Incident wavetype (``'P'``, ``'SV'``, ``'SH'``, ``'Si'``)

    Returns:
        (float): t1: Time in seconds

    Example
    -------
    >>> from telewavesim import utils
    >>> # Define two-layer model with identical material
    >>> model = utils.Model([10, 0], None, 0, 0, 'atg', 0, 0, 0)
    >>> # Only topmost layer is useful for travel time calculation
    >>> wvtype = 'P'
    >>> slow = 0.06     # s/km
    >>> utils.calc_ttime(model, slow, wvtype)
    1.3519981570791182

    """

    t1 = 0.

    for i in range(model.nlay-1):
        if model.isoflg[i] == 'iso':
            a0 = model.a[2, 2, 2, 2, i]
            b0 = model.a[1, 2, 1, 2, i]
        else:
            cc = cc2voigt(model.a[:, :, :, :, i])
            rho = model.rho[i]
            K1, G1, K2, G2, K, G = VRH_average(cc*rho)
            a0, b0 = mod2vel(K, G, rho)
            a0 = a0**2
            b0 = b0**2
        if wvtype == 'P':
            t1 += 1000*model.thickn[i]*np.sqrt(1./a0 - (slow*1.e-3)**2)
        elif wvtype == 'Si' or wvtype == 'SV' or wvtype == 'SH':
            t1 += 1000*model.thickn[i]*np.sqrt(1./b0 - (slow*1.e-3)**2)
        else:
            raise ValueError('Invalid wave type')

    return t1


class Model(object):
    """
    ``model parameters``:
        - thickn (np.ndarray): Thickness of layers (km) (shape ``(nlay)``)
        - rho (np.ndarray): Density (kg/m^3) (shape ``(nlay)``)
        - vp (np.ndarray): P-wave velocity (km/s) (shape ``(nlay)``)
        - vs (np.ndarray): S-wave velocity (km/s) (shape ``(nlay)``)
        - isoflg (list of str, optional, defaut: ``'iso'``):
            Flags for type of layer material (dimension ``nlay``)
        - ani (np.ndarray, optional): Anisotropy (percent) (shape ``(nlay)``)
        - tr (np.ndarray, optional):
            Trend of symmetry axis (degree) (shape ``(nlay)``)
        - pl (np.ndarray, optional):
            Plunge of symmetry axis (degree) (shape ``(nlay)``)

        - nlay (int): Number of layers
        - a (np.ndarray): Elastic thickness (shape ``(3, 3, 3, 3, nlay)``)
    """

    def __init__(self, thickn, rho, vp, vs, isoflg='iso',
                 ani=None, tr=None, pl=None):
        def _get_val(v):
            return (np.array([v] * self.nlay if isinstance(v, (int, float))
                             else v) if v is not None else None)
        self.nlay = len(thickn)
        self.thickn = np.array(thickn)
        self.rho = np.array(rho) if rho is not None else [None] * self.nlay
        self.vp = np.array(vp)
        self.vs = np.array(vs)
        self.isoflg = (list(isoflg) if not isinstance(isoflg, str)
                       else [isoflg] * self.nlay)
        self.ani = _get_val(ani)
        self.tr = _get_val(tr)
        self.pl = _get_val(pl)
        self.update_tensor()

    def update_tensor(self):
        """
        Update the elastic thickness tensor ``a``.

        Needs to be called when model parameters change.
        """
        self.nlay = len(self.thickn)
        self.a = np.zeros((3, 3, 3, 3, self.nlay))

        for j in range(self.nlay):
            if self.isoflg[j] == 'iso':
                cc = set_iso_tensor(self.vp[j], self.vs[j])
                self.a[:, :, :, :, j] = cc
            elif self.isoflg[j] == 'tri':
                cc = set_tri_tensor(self.vp[j], self.vs[j],
                                    self.tr[j], self.pl[j], self.ani[j])
                self.a[:, :, :, :, j] = cc
            elif self.isoflg[j] in MINERALS or self.isoflg[j] in ROCKS:
                cc, rho = set_aniso_tensor(self.tr[j], self.pl[j],
                                           typ=self.isoflg[j])
                self.a[:, :, :, :, j] = cc
                self.rho[j] = rho
            else:
                msg = ('\nFlag not defined: use either "iso", "tri" or one '
                       'among\n%s\n%s\n')
                raise ValueError(msg % (MINERALS, ROCKS))


def read_model(modfile, encoding=None):
    """
    Reads model parameters from file and returns a Model object.

    Returns:
        Model object
    """
    values = np.genfromtxt(modfile, dtype=None, encoding=encoding)
    return Model(*zip(*values))


def model2for(model):
    """
    Passes global model variables to Fortran ``conf`` module.

    Returns:
        None

    Variables to pass are ``a``, ``rho``, ``thickn``, ``isoflg``
    """

    nlaymx = cf_f.nlaymx
    cf_f.a = np.zeros((3, 3, 3, 3, nlaymx))
    cf_f.rho = np.zeros((nlaymx))
    cf_f.thickn = np.zeros((nlaymx))
    cf_f.isoflg = np.zeros((nlaymx), dtype='int')

    for i in range(model.nlay):
        cf_f.a[:, :, :, :, i] = model.a[:, :, :, :, i]
        cf_f.rho[i] = model.rho[i]
        cf_f.thickn[i] = 1000. * model.thickn[i]
        if model.isoflg[i] == 'iso':
            cf_f.isoflg[i] = 1
        else:
            cf_f.isoflg[i] = 0


def wave2for(dt, slow, baz):
    """
    Passes global wavefield variables to Fortran ``conf`` module.

    Returns:
        None

    Variables to pass are ``dt``, ``slow``, ``baz``
    """

    cf_f.dt = dt
    cf_f.slow = slow
    cf_f.baz = baz


def obs2for(dp, c, rhof):
    """
    Passes global OBS-related variables to Fortran ``conf`` module.

    Returns:
        None

    Variables to pass are ``dp``, ``c``, ``rhof``
    """
    cf_f.dp = dp
    cf_f.c = 1000 * c
    cf_f.rhof = rhof


def run_plane(model, slow, npts, dt, baz=0, wvtype='P',
              obs=False, dp=100., c=1.5, rhof=1027):
    """
    Function to run the ``plane`` module and return 3-component seismograms as
    an ``obspy.Stream`` object. This function builds the seismic response
    spectrum by frequency using the matrix propagation approach.
    Required arguments are the seismic model, slowness value, and the sampling
    properties (maximum number of samples and
    the sampling distance in seconds).

    By default, the function uses a back-azimuth value of 0 degree,
    which is suitable for events coming from the North pole or isotropic
    seismic velocity models (i.e., those that do not vary with direction of
    incoming waves).
    For anisotropic velocity models, users need to specify the back-azimuth
    value in degrees. Furthermore, the default type of the incoming
    teleseismic body wave is ``'P'`` for compressional wave. Other options are
    ``'SV'``, ``'SH'``, or ``'Si'`` for vertically-polarized shear wave,
    horizontally-polarized shear wave or isotropic shear wave, respectively.
    Wave modes cannot be mixed.

    Finally, it is possible to simulate the seismic response for ocean-bottom
    seismic (OBS) stations using the flag ``obs=True``. If the flag is set to
    ``True``, the user can specify the water depth below sea level
    (in meters, positive value) as well as the properties of the sea water
    (defaults are acoustic wavespeed of 1.5 km/s and density of 1027 kg/m^3).

    Args:
        model (Model):
            Instance of the ``Model`` class that contains the physical
            properties of subsurface layers.
        slow (float): Slowness (s/km)
        baz (float): Back-azimuth (degree)
        npts (int): Number of samples in time series
        dt (float): Sampling distance (s)
        baz (float, optional): Back-azimuth (degree)
        wvtype (str, optional, default: ``'P'``):
            Incident wavetype (``'P'``, ``'SV'``, ``'SH'``, ``'Si'``)
        obs (bool, optional):
            Whether or not the analysis is done for an OBS stations
        dp (float, optional): Deployment depth below sea level (m)
        c (float, optional):
            P-wave velocity of salt water (default = ``1.5`` km/s)
        rhof (float, optional):
            Density of salt water (default = ``1027.0`` kg/m^3)


    Returns:
        (obspy.stream):
            trxyz: Stream containing 3-component displacement seismograms


    Example
    -------

    Basic example:

    >>> from telewavesim import utils
    >>> # Define three-layer model with isotropic crust and antigorite upper mantle layer over isotropic half-space
    >>> model = utils.Model([20, 10, 0], [2800., None, 3300.], [4.6, 0, 6.0], [2.6, 0, 3.6], ['iso', 'atg', 'iso'], [0, 0, 0], [0, 0, 0], [0, 0, 0])
    >>> slow = 0.06     # s/km
    >>> npts = 1500
    >>> dt = 0.025      # s
    >>> st = utils.run_plane(model, slow, npts, dt)

    >>> type(st)
    <class 'obspy.core.stream.Stream'>
    >>> print(st)
    3 Trace(s) in Stream:
    ...N | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    ...E | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    ...Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    >>> st.plot(size=(600, 450))

    .. figure:: ../telewavesim/examples/picture/Figure_land.png
       :align: center

    OBS station:

    >>> from telewavesim import utils
    >>> # Define two-layer model with foliated eclogitic crust over isotropic half-space
    >>> model = utils.Model([20, 0], [None, 3300.], [0, 6.0], [0, 3.6], ['EC_f', 'iso'], [0, 0], [0, 0], [0, 0])
    >>> slow = 0.06     # s/km
    >>> npts = 3000
    >>> dt = 0.01      # s
    >>> wvtype = 'SV'
    >>> baz = 45.
    >>> dp = 1000.
    >>> st = utils.run_plane(model, slow, npts, dt, baz=baz, wvtype=wvtype, obs=True, dp=dp)
    >>> st.plot(size=(600, 450))

    .. figure:: ../telewavesim/examples/picture/Figure_obs.png
       :align: center

    """

    # Pass  variables to Fortran conf
    model2for(model)
    wave2for(dt, slow, baz)

    # Run the ``plane`` module depending on land or OBS case.
    if obs:

        # If OBS, then further pass OBS-related paramters to Fortran conf
        obs2for(dp, c, rhof)

        # Get the Fourier transform of seismograms for ``obs``case
        yx, yy, yz = pw_f.plane_obs(
            npts, model.nlay, np.array(wvtype, dtype='c'))

    else:

        # Get the Fourier transform of seismograms for ``land`` case
        yx, yy, yz = pw_f.plane_land(
            npts, model.nlay, np.array(wvtype, dtype='c'))

    # Transfer displacement seismograms to an ``obspy`` ``Stream`` object.
    trxyz = get_trxyz(yx, yy, yz, npts, dt, slow, baz, wvtype)

    return trxyz


def get_trxyz(yx, yy, yz, npts, dt, slow, baz, wvtype):
    """
    Function to store displacement seismograms into ``obspy.Trace`` objects
    and then an ``obspy`` ``Stream`` object.

    Args:
        ux (np.ndarray): x-component displacement seismogram
        uy (np.ndarray): y-component displacement seismogram
        uz (np.ndarray): z-component displacement seismogram

    Returns:
        (obspy.stream): trxyz: Stream containing 3-component displacement
          seismograms

    """

    # Get displacements in time domain
    ux = np.real(fft(yx))
    uy = np.real(fft(yy))
    uz = -np.real(fft(yz))

    # Store in traces
    tux = Trace(data=ux)
    tuy = Trace(data=uy)
    tuz = Trace(data=uz)

    # Update trace header
    tux = update_stats(tux, dt, slow, baz, wvtype, 'N')
    tuy = update_stats(tuy, dt, slow, baz, wvtype, 'E')
    tuz = update_stats(tuz, dt, slow, baz, wvtype, 'Z')

    # Append to stream
    trxyz = Stream(traces=[tux, tuy, tuz])

    return trxyz


def tf_from_xyz(trxyz, pvh=False, vp=None, vs=None):
    """
    Function to generate transfer functions from displacement traces.

    Args:
        trxyz (obspy.stream):
            Obspy ``Stream`` object in cartesian coordinate system
        pvh (bool, optional):
            Whether to rotate from Z-R-T coordinate system to P-SV-SH wave mode
        vp (float, optional):
            Vp velocity at surface for rotation to P-SV-SH system
        vs (float, optional):
            Vs velocity at surface for rotation to P-SV-SH system

    Returns:
        (obspy.stream):
            tfs: Stream containing Radial and Transverse transfer functions

    """

    # Extract East, North and Vertical
    ntr = trxyz[0]
    etr = trxyz[1]
    ztr = trxyz[2]
    baz = ntr.stats.baz
    slow = ntr.stats.slow
    wvtype = ntr.stats.wvtype

    # Copy to radial and transverse
    rtr = ntr.copy()
    ttr = etr.copy()

    # Rotate to radial and transverse
    rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)

    if pvh:
        trP, trV, trH = rotate_zrt_pvh(ztr, rtr, ttr, slow, vp=vp, vs=vs)

        tfr = trV.copy()
        tfr.data = np.zeros(len(tfr.data))
        tft = trH.copy()
        tft.data = np.zeros(len(tft.data))
        ftfv = fft(trV.data)
        ftfh = fft(trH.data)
        ftfp = fft(trP.data)

        if wvtype == 'P':
            # Transfer function
            tfr.data = fftshift(np.real(ifft(np.divide(ftfv, ftfp))))
            tft.data = fftshift(np.real(ifft(np.divide(ftfh, ftfp))))
        elif wvtype == 'Si':
            tfr.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfv))))
            tft.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfh))))
        elif wvtype == 'SV':
            tfr.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfv))))
        elif wvtype == 'SH':
            tft.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfh))))
    else:
        tfr = rtr.copy()
        tfr.data = np.zeros(len(tfr.data))
        tft = ttr.copy()
        tft.data = np.zeros(len(tft.data))
        ftfr = fft(rtr.data)
        ftft = fft(ttr.data)
        ftfz = fft(ztr.data)

        if wvtype == 'P':
            # Transfer function
            tfr.data = fftshift(np.real(ifft(np.divide(ftfr, ftfz))))
            tft.data = fftshift(np.real(ifft(np.divide(ftft, ftfz))))
        elif wvtype == 'Si':
            tfr.data = fftshift(np.real(ifft(np.divide(-ftfz, ftfr))))
            tft.data = fftshift(np.real(ifft(np.divide(-ftfz, ftft))))
        elif wvtype == 'SV':
            tfr.data = fftshift(np.real(ifft(np.divide(-ftfz, ftfr))))
        elif wvtype == 'SH':
            tft.data = fftshift(np.real(ifft(np.divide(-ftfz, ftft))))

    # Store in stream
    tfs = Stream(traces=[tfr, tft])

    # Return stream
    return tfs


def update_stats(tr, dt, slow, baz, wvtype, cha):
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
    tr.stats.wvtype = wvtype
    tr.stats.channel = cha

    return tr
