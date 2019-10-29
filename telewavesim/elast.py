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

Functions to define elastic stiffness matrices. Most references for elastic
constants are included in the paper by
`Brownlee et al. (2017) <https://doi.org/10.1002/2017TC004625>`_.

'''

import numpy as np

def iso_tensor(a, b):
    """
    Elastic constants (GPa /density) of isotropic material in Voigt notation

    Args:
        a (float): P-wave velocity (km/s)
        b (float): S-wave velocity (km/s)

    Returns:
        (np.ndarray): C: Elastic stiffness matrix (shape ``(6, 6)``)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.iso_tensor(6.e3, 3.6e3)
    array([[36000000., 10080000., 10080000.,        0.,        0.,        0.],
           [10080000., 36000000., 10080000.,        0.,        0.,        0.],
           [10080000., 10080000., 36000000.,        0.,        0.,        0.],
           [       0.,        0.,        0., 12960000.,        0.,        0.],
           [       0.,        0.,        0.,        0., 12960000.,        0.],
           [       0.,        0.,        0.,        0.,        0., 12960000.]])

    """

    # Velocity to Lame parameters
    mu = b*b
    lam = a*a - 2.*mu

    # Components of Voigt matrix
    Cii = 2.*mu + lam
    Cij = lam
    Cjj = mu

    C = np.zeros((6,6), dtype=float)
    C[0,0] = Cii;       C[0,1] = Cij;       C[0,2] = Cij;       C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = Cii;       C[1,2] = Cij;       C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = Cii;       C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = Cjj;       C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = Cjj;       C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = Cjj

    return C


def antigorite():
    """
    Elastic constants of antigorite mineral (GPa) from
    Bezacier, EPSL, 2010, in Voigt notation.

    - Abbreviation: ``'atg'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2620 kg/m^3)

    Example
    -------

    >>> from telewavesim import elast
    >>> elast.antigorite()[0]
    array([[208.4,  66.2,  15.9,   0. ,   2.4,   0. ],
           [ 66.2, 201.6,   5. ,   0. ,  -4.4,   0. ],
           [ 15.9,   5. ,  96.7,   0. ,   2.5,   0. ],
           [  0. ,   0. ,   0. ,  17.4,   0. , -13.1],
           [  2.4,  -4.4,   2.5,   0. ,  18.3,   0. ],
           [  0. ,   0. ,   0. , -13.1,   0. ,  65. ]])
    >>> elast.antigorite()[1]
    2620.0

    """

    rho = 2620.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 208.4;     C[0,1] = 66.2;      C[0,2] = 15.9;      C[0,3] = 0.;        C[0,4] = 2.4;       C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 201.6;     C[1,2] = 5.;        C[1,3] = 0.;        C[1,4] = -4.4;      C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 96.7;      C[2,3] = 0.;        C[2,4] = 2.5;       C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 17.4;      C[3,4] = 0.;        C[3,5] = -13.1
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 18.3;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 65.

    return C, rho


def biotite():
    """
    Elastic constants of biotite mineral (GPa) from
    Aleksandrov and Ryzhova (1986), in Voigt notation

    - Abbreviation: ``'bt'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2800 kg/m^3)

    Example
    -------

    >>> from telewavesim import elast
    >>> elast.biotite()[0]
    array([[186. ,  32.4,  11.6,   0. ,   0. ,   0. ],
           [ 32.4, 186. ,  11.6,   0. ,   0. ,   0. ],
           [ 11.6,  11.6,  54. ,   0. ,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   5.8,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   5.8,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   0. ,  76.8]])
    >>> elast.biotite()[1]
    2800.0

    """

    rho = 2800.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 186.;      C[0,1] = 32.4;      C[0,2] = 11.6;      C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 186.;      C[1,2] = 11.6;      C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 54.;       C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 5.8;       C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 5.8;       C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 76.8

    return C, rho


def blueschist_felsic():
    """
    Elastic constants of Felsic Blueschist (GPa) from
    Cao et al. (2013), in Voigt notation

    - Abbreviation: ``'BS_f'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2970 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.blueschist_felsic()[0]
    array([[ 1.4985e+02,  3.8700e+01,  3.2590e+01, -1.5000e-01, -1.0000e+00,
            -1.9000e-01],
           [ 3.8700e+01,  1.6355e+02,  3.0030e+01,  1.0500e+00, -1.8100e+00,
            -1.7800e+00],
           [ 3.2590e+01,  3.0030e+01,  1.2162e+02,  2.2000e-01, -9.5000e-01,
            -1.3000e-01],
           [-1.5000e-01,  1.0500e+00,  2.2000e-01,  4.8030e+01, -6.3000e-01,
            -1.1400e+00],
           [-1.0000e+00, -1.8100e+00, -9.5000e-01, -6.3000e-01,  4.8620e+01,
            -1.0000e-02],
           [-1.9000e-01, -1.7800e+00, -1.3000e-01, -1.1400e+00, -1.0000e-02,
             5.8420e+01]])
    >>> elast.blueschist_felsic()[1]
    2970.0

    """

    rho = 2970.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 149.85;    C[0,1] = 38.7;      C[0,2] = 32.59;     C[0,3] = -0.15;     C[0,4] = -1.;       C[0,5] = -0.19
    C[1,0] = C[0,1];    C[1,1] = 163.55;    C[1,2] = 30.03;     C[1,3] = 1.05;      C[1,4] = -1.81;     C[1,5] = -1.78
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 121.62;    C[2,3] = 0.22;      C[2,4] = -0.95;     C[2,5] = -0.13
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 48.03;     C[3,4] = -0.63;     C[3,5] = -1.14
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 48.62;     C[4,5] = -0.01
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 58.42

    return C, rho


def blueschist_mafic():
    """
    Elastic constants of Mafic Blueschist (GPa) from
    Cao et al. (2013), in Voigt notation

    - Abbreviation: ``'BS_m'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3190 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.blueschist_mafic()[0]
    array([[ 1.9079e+02,  6.2280e+01,  5.2940e+01, -4.4000e-01,  4.6800e+00,
             6.0000e-01],
           [ 6.2280e+01,  2.1838e+02,  5.3100e+01, -8.7000e-01,  1.5700e+00,
             2.8000e-01],
           [ 5.2940e+01,  5.3100e+01,  1.5804e+02, -4.4000e-01,  2.6600e+00,
            -3.5000e-01],
           [-4.4000e-01, -8.7000e-01, -4.4000e-01,  6.0860e+01, -2.9000e-01,
             1.8600e+00],
           [ 4.6800e+00,  1.5700e+00,  2.6600e+00, -2.9000e-01,  5.8940e+01,
            -2.0000e-01],
           [ 6.0000e-01,  2.8000e-01, -3.5000e-01,  1.8600e+00, -2.0000e-01,
             6.9630e+01]])
    >>> elast.blueschist_mafic()[1]
    3190.0

    """

    rho = 3190.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 190.79;    C[0,1] = 62.28;     C[0,2] = 52.94;     C[0,3] = -0.44;     C[0,4] = 4.68;      C[0,5] = 0.6
    C[1,0] = C[0,1];    C[1,1] = 218.38;    C[1,2] = 53.1;      C[1,3] = -0.87;     C[1,4] = 1.57;      C[1,5] = 0.28
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 158.04;    C[2,3] = -0.44;     C[2,4] = 2.66;      C[2,5] = -0.35
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 60.86;     C[3,4] = -0.29;     C[3,5] = 1.86
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 58.94;     C[4,5] = -0.2
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 69.63

    return C, rho


def clinopyroxene_92():
    """
    Elastic constants of clinopyroxene mineral (GPa) from
    Baghat et al. (1992), in Voigt notation

    - Abbreviation: ``'cpx'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3327 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.clinopyroxene_92()[0]
    array([[257.3,  85.9,  76.2,   0. ,   7.1,   0. ],
           [ 85.9, 216.2,  71.8,   0. ,  13.3,   0. ],
           [ 76.2,  71.8, 260.2,   0. ,  33.7,   0. ],
           [  0. ,   0. ,   0. ,  80.2,   0. ,  10.2],
           [  7.1,  13.3,  33.7,   0. ,  70.6,   0. ],
           [  0. ,   0. ,   0. ,  10.2,   0. ,  85.8]])
    >>> elast.clinopyroxene_92()[1]
    3327.0

    """

    rho = 3327.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 257.3;     C[0,1] = 85.9;      C[0,2] = 76.2;      C[0,3] = 0.;        C[0,4] = 7.1;       C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 216.2;     C[1,2] = 71.8;      C[1,3] = 0.;        C[1,4] = 13.3;      C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 260.2;     C[2,3] = 0.;        C[2,4] = 33.7;      C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 80.2;      C[3,4] = 0.;        C[3,5] = 10.2
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 70.6;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 85.8

    return C, rho


def clinopyroxene_98():
    """
    Elastic constants of clinopyroxene mineral (GPa) from
    Collins and Brown (1998), in Voigt notation

    - Abbreviation: ``None``: not currently used

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3190 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.clinopyroxene_98()[0]
    array([[237.8,  83.5,  80. ,   0. ,   9. ,   0. ],
           [ 83.5, 183.6,  59.9,   0. ,   9.5,   0. ],
           [ 80. ,  59.9, 229.5,   0. ,  48.1,   0. ],
           [  0. ,   0. ,   0. ,  76.5,   0. ,   8.4],
           [  9. ,   9.5,  48.1,   0. ,  73. ,   0. ],
           [  0. ,   0. ,   0. ,   8.4,   0. ,  81.6]])
    >>> elast.clinopyroxene_98()[1]
    3190.0

    """

    rho = 3190.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 237.8;     C[0,1] = 83.5;      C[0,2] = 80.;       C[0,3] = 0.;        C[0,4] = 9.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 183.6;     C[1,2] = 59.9;      C[1,3] = 0.;        C[1,4] = 9.5;       C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 229.5;     C[2,3] = 0.;        C[2,4] = 48.1;      C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 76.5;      C[3,4] = 0.;        C[3,5] = 8.4
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 73.;       C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 81.6

    return C, rho


def dolomite():
    """
    Elastic constants of dolomite mineral (GPa) from
    Humbert and Plicque (1972), in Voigt notation

    - Abbreviation: ``'dol'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2840 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.dolomite()[0]
    array([[205. ,  71. ,  57.4, -19.5,  13.7,   0. ],
           [ 71. , 205. ,  57.4,  19.5, -13.7,   0. ],
           [ 57.4,  57.4, 113. ,   0. ,   0. ,   0. ],
           [-19.5,  19.5,   0. ,  39.8,   0. , -13.7],
           [ 13.7, -13.7,   0. ,   0. ,  39.8, -19.5],
           [  0. ,   0. ,   0. , -13.7, -19.5,  67. ]])
    >>> elast.dolomite()[1]
    2840.0
    """

    rho = 2840.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 205.;      C[0,1] = 71.;       C[0,2] = 57.4;      C[0,3] = -19.5;     C[0,4] = 13.7;      C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 205.;      C[1,2] = 57.4;      C[1,3] = 19.5;      C[1,4] = -13.7;     C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 113.;      C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 39.8;      C[3,4] = 0.;        C[3,5] = -13.7
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 39.8;      C[4,5] = -19.5
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 67.

    return C, rho


def eclogite_foliated():
    """
    Elastic constants of Foliated Eclogite rock (GPa) from
    Cao et al. (2013), in Voigt notation

    - Abbreviation: ``'EC_f'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3300 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.eclogite_foliated()[0]
    array([[ 2.0345e+02,  6.7760e+01,  6.4470e+01,  8.0000e-02,  1.9000e+00,
            -4.0000e-01],
           [ 6.7760e+01,  2.2058e+02,  6.3650e+01,  4.6000e-01,  5.9000e-01,
             6.0000e-02],
           [ 6.4470e+01,  6.3650e+01,  1.8975e+02,  1.3000e-01,  9.5000e-01,
            -2.0000e-01],
           [ 8.0000e-02,  4.6000e-01,  1.3000e-01,  6.6320e+01, -2.7000e-01,
             7.3000e-01],
           [ 1.9000e+00,  5.9000e-01,  9.5000e-01, -2.7000e-01,  6.5770e+01,
            -2.0000e-02],
           [-4.0000e-01,  6.0000e-02, -2.0000e-01,  7.3000e-01, -2.0000e-02,
             7.0750e+01]])
    >>> elast.eclogite_foliated()[1]
    3300.0

    """

    rho = 3300.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 203.45;    C[0,1] = 67.76;     C[0,2] = 64.47;     C[0,3] = 0.08;      C[0,4] = 1.9;       C[0,5] = -0.4
    C[1,0] = C[0,1];    C[1,1] = 220.58;    C[1,2] = 63.65;     C[1,3] = 0.46;      C[1,4] = 0.59;      C[1,5] = 0.06
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 189.75;    C[2,3] = 0.13;      C[2,4] = 0.95;      C[2,5] = -0.2
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 66.32;     C[3,4] = -0.27;     C[3,5] = 0.73
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 65.77;     C[4,5] = -0.02
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 70.75

    return C, rho


def eclogite_massive():
    """
    Elastic constants of Massive Eclogite rock (GPa) from
    Cao et al. (2013), in Voigt notation

    - Abbreviation: ``'EC_m'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3490 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.eclogite_massive()[0]
    array([[ 2.3885e+02,  8.2010e+01,  8.1440e+01,  3.0000e-01, -2.0000e-02,
             5.0000e-01],
           [ 8.2010e+01,  2.4212e+02,  8.1110e+01, -6.6000e-01,  3.3000e-01,
             1.2000e-01],
           [ 8.1440e+01,  8.1110e+01,  2.3557e+02, -2.8000e-01,  2.2000e-01,
             3.1000e-01],
           [ 3.0000e-01, -6.6000e-01, -2.8000e-01,  7.8720e+01,  2.7000e-01,
             0.0000e+00],
           [-2.0000e-02,  3.3000e-01,  2.2000e-01,  2.7000e-01,  7.8370e+01,
             2.5000e-01],
           [ 5.0000e-01,  1.2000e-01,  3.1000e-01,  0.0000e+00,  2.5000e-01,
             7.7910e+01]])
    >>> elast.eclogite_massive()[1]
    3490.0

    """

    rho = 3490.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 238.85;    C[0,1] = 82.01;     C[0,2] = 81.44;     C[0,3] = 0.3;       C[0,4] = -0.02;     C[0,5] = 0.5
    C[1,0] = C[0,1];    C[1,1] = 242.12;    C[1,2] = 81.11;     C[1,3] = -0.66;     C[1,4] = 0.33;      C[1,5] = 0.12
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 235.57;    C[2,3] = -0.28;     C[2,4] = 0.22;      C[2,5] = 0.31
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 78.72;     C[3,4] = 0.27;      C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 78.37;     C[4,5] = 0.25
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 77.91

    return C, rho


def epidote():
    """
    Elastic constants of epidote mineral (GPa) from
    Aleksandrakov et al. (1974), in Voigt notation

    - Abbreviation: ``'ep'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3465 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.epidote()[0]
    array([[211.5,  65.6,  43.2,   0. ,  -6.5,   0. ],
           [ 65.6, 239. ,  43.6,   0. , -10.4,   0. ],
           [ 43.2,  43.6, 202.1,   0. , -20. ,   0. ],
           [  0. ,   0. ,   0. ,  39.1,   0. ,  -2.3],
           [ -6.5, -10.4, -20. ,   0. ,  43.4,   0. ],
           [  0. ,   0. ,   0. ,  -2.3,   0. ,  79.5]])
    >>> elast.epidote()[1]
    3465.0

    """

    rho = 3465.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 211.5;     C[0,1] = 65.6;      C[0,2] = 43.2;      C[0,3] = 0.;        C[0,4] = -6.5;      C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 239.;      C[1,2] = 43.6;      C[1,3] = 0.;        C[1,4] = -10.4;     C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 202.1;     C[2,3] = 0.;        C[2,4] = -20.;      C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 39.1;      C[3,4] = 0.;        C[3,5] = -2.3
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 43.4;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 79.5

    return C, rho


def garnet():
    """
    Elastic constants of garnet mineral (GPa) from
    Babuska et al. (1978), in Voigt notation

    - Abbreviation: ``'grt'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3660 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.garnet()[0]
    array([[306.2, 112.5, 112.5,   0. ,   0. ,   0. ],
           [112.5, 306.2, 112.5,   0. ,   0. ,   0. ],
           [112.5, 112.5, 306.2,   0. ,   0. ,   0. ],
           [  0. ,   0. ,   0. ,  92.7,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,  92.7,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   0. ,  92.7]])
    >>> elast.garnet()[1]
    3660.0

    """

    rho = 3660.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 306.2;     C[0,1] = 112.5;     C[0,2] = 112.5;     C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 306.2;     C[1,2] = 112.5;     C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 306.2;     C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 92.7;      C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 92.7;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 92.7

    return C, rho


def glaucophane():
    """
    Elastic constants of glaucophane mineral (GPa) from
    Bezacier et al. (2010), in Voigt notation

    - Abbreviation: ``'gln'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3070 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.glaucophane()[0]
    array([[122.3 ,  45.7 ,  37.2 ,   0.  ,   2.3 ,   0.  ],
           [ 45.7 , 231.5 ,  74.9 ,   0.  ,  -4.8 ,   0.  ],
           [ 37.2 ,  74.9 , 254.6 ,   0.  ,  -2.37,   0.  ],
           [  0.  ,   0.  ,   0.  ,  79.6 ,   0.  ,   8.9 ],
           [  2.3 ,  -4.8 ,  -2.37,   0.  ,  52.8 ,   0.  ],
           [  0.  ,   0.  ,   0.  ,   8.9 ,   0.  ,  51.2 ]])
    >>> elast.glaucophane()[1]
    3070.0

    """

    rho = 3070.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 122.3;     C[0,1] = 45.7;      C[0,2] = 37.2;      C[0,3] = 0.;        C[0,4] = 2.3;       C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 231.5;     C[1,2] = 74.9;      C[1,3] = 0.;        C[1,4] = -4.8;      C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 254.6;     C[2,3] = 0.;        C[2,4] = -2.37;     C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 79.6;      C[3,4] = 0.;        C[3,5] = 8.9
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 52.8;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 51.2

    return C, rho


def harzburgite():
    """
    Elastic constants of harzburgite rock (GPa) from
    Covey-Crump et al. (2003), in Voigt notation

    - Abbreviation: ``'HB'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3200 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.harzburgite()[0]
    array([[226.5 ,  75.34,  74.73,  -0.27,  -2.  ,   1.85],
           [ 75.34, 242.8 ,  73.68,  -3.6 ,  -1.91,   4.14],
           [ 74.73,  73.68, 230.  ,  -4.36,  -4.27,  -0.27],
           [ -0.27,  -3.6 ,  -4.36,  80.75,   1.81,  -2.19],
           [ -2.  ,  -1.91,  -4.27,   1.81,  76.94,  -1.88],
           [  1.85,   4.14,  -0.27,  -2.19,  -1.88,  79.15]])
    >>> elast.harzburgite()[1]
    3200.0

    """

    rho = 3200.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 226.5;     C[0,1] = 75.34;     C[0,2] = 74.73;     C[0,3] = -0.27;     C[0,4] = -2.00;     C[0,5] = 1.85
    C[1,0] = C[0,1];    C[1,1] = 242.8;     C[1,2] = 73.68;     C[1,3] = -3.6;      C[1,4] = -1.91;     C[1,5] = 4.14
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 230.;      C[2,3] = -4.36;     C[2,4] = -4.27;     C[2,5] = -0.27
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 80.75;     C[3,4] = 1.81;      C[3,5] = -2.19
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 76.94;     C[4,5] = -1.88
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 79.15

    return C, rho


def hornblende():
    """
    Elastic constants of hornblende mineral (GPa) from
    Aleksandrov and Ryzhova (1986), in Voigt notation

    - Abbreviation: ``'hbl'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3200 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.hornblende()[0]
    array([[116. ,  49.9,  61.4,   0. ,   4.3,   0. ],
           [ 49.9, 159.7,  65.5,   0. ,  -2.5,   0. ],
           [ 61.4,  65.5, 191.6,   0. ,  10. ,   0. ],
           [  0. ,   0. ,   0. ,  57.4,   0. ,  -6.2],
           [  4.3,  -2.5,  10. ,   0. ,  31.8,   0. ],
           [  0. ,   0. ,   0. ,  -6.2,   0. ,  36.8]])
    >>> elast.hornblende()[1]
    3200.0

    """

    rho = 3200.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 116.;      C[0,1] = 49.9;      C[0,2] = 61.4;      C[0,3] = 0.;        C[0,4] = 4.3;       C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 159.7;     C[1,2] = 65.5;      C[1,3] = 0.;        C[1,4] = -2.5;      C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 191.6;     C[2,3] = 0.;        C[2,4] = 10.;       C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 57.4;      C[3,4] = 0.;        C[3,5] = -6.2
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 31.8;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 36.8

    return C, rho


def jadeite():
    """
    Elastic constants of jadeite mineral (GPa) from
    Kandelin and Weiner (1988), in Voigt notation

    - Abbreviation: ``'jade'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3330 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.jadeite()[0]
    array([[274.,  94.,  71.,   0.,   4.,   0.],
           [ 94., 253.,  82.,   0.,  14.,   0.],
           [ 71.,  82., 282.,   0.,  28.,   0.],
           [  0.,   0.,   0.,  88.,   0.,  13.],
           [  4.,  14.,  28.,   0.,  65.,   0.],
           [  0.,   0.,   0.,  13.,   0.,  94.]])
    >>> elast.jadeite()[1]
    3330.0

    """

    rho = 3330.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 274.;      C[0,1] = 94.;       C[0,2] = 71.;       C[0,3] = 0.;        C[0,4] = 4.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 253.;      C[1,2] = 82.;       C[1,3] = 0.;        C[1,4] = 14.;       C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 282.;      C[2,3] = 0.;        C[2,4] = 28.;       C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 88.;       C[3,4] = 0.;        C[3,5] = 13.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 65.;       C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 94.

    return C, rho


def lawsonite():
    """
    Elastic constants of jadeite mineral (GPa) from
    Kandelin and Weiner (1988), in Voigt notation

    - Abbreviation: ``'lws'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3090 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.lawsonite()[0]
    array([[214.,  69.,  82.,   0.,   0.,   0.],
           [ 69., 226.,  65.,   0.,   0.,   0.],
           [ 82.,  65., 259.,   0.,   0.,   0.],
           [  0.,   0.,   0.,  60.,   0.,   0.],
           [  0.,   0.,   0.,   0.,  65.,   0.],
           [  0.,   0.,   0.,   0.,   0.,  17.]])
    >>> elast.lawsonite()[1]
    3090.0

    """

    rho = 3090.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 214.;      C[0,1] = 69.;       C[0,2] = 82.;       C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 226.;      C[1,2] = 65.;       C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 259.;      C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 60.;       C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 65.;       C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 17.

    return C, rho


def lherzolite():
    """
    Elastic constants of lherzolite rock (GPa) from
    Peselnick et al. (1974), in Voigt notation

    - Abbreviation: ``'LHZ'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3270 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.lherzolite()[0]
    array([[ 1.8740e+02,  6.3710e+01,  6.3870e+01,  7.8000e-01,  2.0200e+00,
            -3.2000e+00],
           [ 6.3710e+01,  2.1125e+02,  6.4500e+01, -3.0700e+00,  8.7000e-01,
            -5.7800e+00],
           [ 6.3870e+01,  6.4500e+01,  1.9000e+02,  3.8000e-01,  2.3800e+00,
            -1.2000e-01],
           [ 7.8000e-01, -3.0700e+00,  3.8000e-01,  6.7900e+01, -2.1200e+00,
             1.6000e+00],
           [ 2.0200e+00,  8.7000e-01,  2.3800e+00, -2.1200e+00,  6.3120e+01,
            -5.5000e-01],
           [-3.2000e+00, -5.7800e+00, -1.2000e-01,  1.6000e+00, -5.5000e-01,
             6.6830e+01]])
    >>> elast.lherzolite()[1]
    3270.0

    """

    rho = 3270.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 187.4;     C[0,1] = 63.71;     C[0,2] = 63.87;     C[0,3] = 0.78;      C[0,4] = 2.02;      C[0,5] = -3.2
    C[1,0] = C[0,1];    C[1,1] = 211.25;    C[1,2] = 64.5;      C[1,3] = -3.07;     C[1,4] = 0.87;      C[1,5] = -5.78
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 190.;      C[2,3] = 0.38;      C[2,4] = 2.38;      C[2,5] = -0.12
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 67.9;      C[3,4] = -2.12;     C[3,5] = 1.6
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 63.12;     C[4,5] = -0.55
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 66.83

    return C, rho


def lizardite_atom():
    """
    Elastic constants of lizardite mineral (GPa) from
    Auzende et al., Phys. Chem. Min. 2006 from atomistic calculations.

    - Abbreviation: ``None``: Not currently used.

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2515 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.lizardite_atom()[0]
    array([[ 2.29080e+02,  8.90440e+01,  1.35580e+01, -1.00000e-04,
             4.60250e+00,  1.00000e-04],
           [ 8.90440e+01,  2.29080e+02,  1.35570e+01, -1.00000e-04,
            -4.60160e+00,  1.00000e-04],
           [ 1.35580e+01,  1.35570e+01,  4.58380e+01, -1.00000e-04,
             1.50000e-03,  1.00000e-04],
           [-1.00000e-04, -1.00000e-04, -1.00000e-04,  1.27650e+01,
            -1.00000e-04, -4.45980e+00],
           [ 4.60250e+00, -4.60160e+00,  1.50000e-03, -1.00000e-04,
             1.27740e+01,  1.00000e-04],
           [ 1.00000e-04,  1.00000e-04,  1.00000e-04, -4.45980e+00,
             1.00000e-04,  7.00166e+01]])
    >>> elast.lizardite_atom()[1]
    2515.5

    """

    rho = 2515.5

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 229.08;    C[0,1] = 89.044;    C[0,2] = 13.558;    C[0,3] = -0.0001;   C[0,4] = 4.6025;    C[0,5] = 0.0001
    C[1,0] = C[0,1];    C[1,1] = 229.08;    C[1,2] = 13.557;    C[1,3] = -0.0001;   C[1,4] = -4.6016;   C[1,5] = 0.0001
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 45.838;    C[2,3] = -0.0001;   C[2,4] = 0.0015;    C[2,5] = 0.0001
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 12.765;    C[3,4] = -0.0001;   C[3,5] = -4.4598
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 12.774;    C[4,5] = 0.0001
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 70.0166

    return C, rho


def lizardite():
    """
    Elastic constants of lizardite mineral (GPa) from
    Reynard, GRL, 2007 from Density Functional Theory.

    - Abbreviation: ``'lz'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2610 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.lizardite()[0]
    array([[245. ,  50. ,  31. ,   0. ,   0. ,   0. ],
           [ 50. , 245. ,  31. ,   0. ,   0. ,   0. ],
           [ 31. ,  31. ,  23. ,   0. ,   0. ,   0. ],
           [  0. ,   0. ,   0. ,  11.6,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,  11.6,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   0. ,  97.5]])
    >>> elast.lizardite()[1]
    2610.0

    """

    rho = 2610.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 245.;      C[0,1] = 50.;       C[0,2] = 31.;       C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 245.;      C[1,2] = 31.;       C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 23.;       C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 11.6;      C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 11.6;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 97.5

    return C, rho


def muscovite():
    """
    Elastic constants of muscovite mineral (GPa) from
    Vaughan and Guggenheim (1986), in Voigt notation

    - Abbreviation: ``'ms'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2834 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.muscovite()[0]
    array([[181. ,  48.8,  25.6,   0. , -14.2,   0. ],
           [ 48.8, 178.4,  21.2,   0. ,   1.1,   0. ],
           [ 25.6,  21.2,  58.6,   0. ,   1. ,   0. ],
           [  0. ,   0. ,   0. ,  16.5,   0. ,  -5.2],
           [-14.2,   1.1,   1. ,   0. ,  19.5,   0. ],
           [  0. ,   0. ,   0. ,  -5.2,   0. ,  72. ]])
    >>> elast.muscovite()[1]
    2834.0

    """

    rho = 2834.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 181.;      C[0,1] = 48.8;      C[0,2] = 25.6;      C[0,3] = 0.;        C[0,4] = -14.2;     C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 178.4;     C[1,2] = 21.2;      C[1,3] = 0.;        C[1,4] = 1.1;       C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 58.6;      C[2,3] = 0.;        C[2,4] = 1.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 16.5;      C[3,4] = 0.;        C[3,5] = -5.2
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 19.5;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 72.

    return C, rho


def olivine():
    """
    Elastic constants of olivine mineral (GPa) from
    Abrahamson et al. (1997), in Voigt notation

    - Abbreviation: ``'ol'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3355 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.olivine()[0]
    array([[320.5 ,  68.15,  71.6 ,   0.  ,   0.  ,   0.  ],
           [ 68.15, 196.5 ,  76.8 ,   0.  ,   0.  ,   0.  ],
           [ 71.6 ,  76.8 , 233.5 ,   0.  ,   0.  ,   0.  ],
           [  0.  ,   0.  ,   0.  ,  64.  ,   0.  ,   0.  ],
           [  0.  ,   0.  ,   0.  ,   0.  ,  77.  ,   0.  ],
           [  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,  78.7 ]])
    >>> elast.olivine()[1]
    3355.0

    """

    rho = 3355.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 320.5;     C[0,1] = 68.15;     C[0,2] = 71.6;      C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 196.5;     C[1,2] = 76.8;      C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 233.5;     C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 64.;       C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 77.;       C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 78.7

    return C, rho


def orthopyroxene():
    """
    Elastic constants of orthopyroxene mineral (GPa) from
    Chai et al. (1977), in Voigt notation

    - Abbreviation: ``'opx'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3304 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.orthopyroxene()[0]
    array([[236.9,  79.6,  63.2,   0. ,   0. ,   0. ],
           [ 79.6, 180.5,  56.8,   0. ,   0. ,   0. ],
           [ 63.2,  56.8, 230.4,   0. ,   0. ,   0. ],
           [  0. ,   0. ,   0. ,  84.3,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,  79.4,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   0. ,  80.1]])
    >>> elast.orthopyroxene()[1]
    3304.0

    """

    rho = 3304.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 236.9;     C[0,1] = 79.6;      C[0,2] = 63.2;      C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 180.5;     C[1,2] = 56.8;      C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 230.4;     C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 84.3;      C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 79.4;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 80.1

    return C, rho


def plagioclase_64():
    """
    Elastic constants of plagioclase mineral (GPa) from
    Ryzhova (1964), in Voigt notation

    - Abbreviation: ``None``: Not currently used.

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2700 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.plagioclase_64()[0]
    array([[ 81.8,  39.3,  40.7,   0. ,  -9. ,   0. ],
           [ 39.3, 145. ,  34.1,   0. ,  -7.9,   0. ],
           [ 40.7,  34.1, 133. ,   0. , -18.5,   0. ],
           [  0. ,   0. ,   0. ,  17.7,   0. ,  -0.8],
           [ -9. ,  -7.9, -18.5,   0. ,  31.2,   0. ],
           [  0. ,   0. ,   0. ,  -0.8,   0. ,  33.3]])
    >>> elast.plagioclase_64()[1]
    2700.0

    """

    rho = 2700.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 81.8;      C[0,1] = 39.3;      C[0,2] = 40.7;      C[0,3] = 0.;        C[0,4] = -9.;       C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 145.;      C[1,2] = 34.1;      C[1,3] = 0.;        C[1,4] = -7.9;      C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 133.;      C[2,3] = 0.;        C[2,4] = -18.5;     C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 17.7;      C[3,4] = 0.;        C[3,5] = -0.8
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 31.2;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 33.3

    return C, rho


def plagioclase_06():
    """
    Elastic constants of plagioclase mineral (GPa) from
    Brown et al. (2006), in Voigt notation

    - Abbreviation: ``'plag'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2700 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.plagioclase_06()[0]
    array([[ 69.9  ,  33.24 ,  31.56 ,   5.28 ,  -2.46 ,  -0.72 ],
           [ 33.24 , 183.28 ,   7.53 ,   5.31 ,  -7.6  ,  -0.423],
           [ 31.56 ,   7.53 , 175.65 , -17.48 ,   5.86 , -11.29 ],
           [  5.28 ,   5.31 , -17.48 ,  26.93 ,  -3.94 ,  -6.56 ],
           [ -2.46 ,  -7.6  ,   5.86 ,  -3.94 ,  26.91 ,   0.98 ],
           [ -0.72 ,  -0.423, -11.29 ,  -6.56 ,   0.98 ,  33.39 ]])
    >>> elast.plagioclase_06()[1]
    2700.0

    """

    rho = 2700.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 69.9;      C[0,1] = 33.24;     C[0,2] = 31.56;     C[0,3] = 5.28;      C[0,4] = -2.46;     C[0,5] = -0.72
    C[1,0] = C[0,1];    C[1,1] = 183.28;    C[1,2] = 7.53;      C[1,3] = 5.31;      C[1,4] = -7.6;      C[1,5] = -0.423
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 175.65;    C[2,3] = -17.48;    C[2,4] = 5.86;      C[2,5] = -11.29
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 26.93;     C[3,4] = -3.94;     C[3,5] = -6.56
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 26.91;     C[4,5] = 0.98
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 33.39

    return C, rho


def quartz():
    """
    Elastic constants of quartz mineral (GPa) from
    Lakshanov et al. (2007), in Voigt notation

    - Abbreviation: ``'qtz'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2649 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.quartz()[0]
    array([[ 86.9,   7.6,  12. ,  17.8,   0. ,   0. ],
           [  7.6,  86.9,  12. , -17.8,   0. ,   0. ],
           [ 12. ,  12. , 106.4,   0. ,   0. ,   0. ],
           [ 17.8, -17.8,   0. ,  59.5,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,  59.5, -17.8],
           [  0. ,   0. ,   0. ,   0. , -17.8,  39.6]])
    >>> elast.quartz()[1]
    2649.0

    """

    rho = 2649.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 86.9;      C[0,1] = 7.6;       C[0,2] = 12.;       C[0,3] = 17.8;      C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 86.9;      C[1,2] = 12.;       C[1,3] = -17.8;     C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 106.4;     C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 59.5;      C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 59.5;      C[4,5] = -17.8
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 39.6

    return C, rho


def serpentinite_37():
    """
    Elastic constants of serpentinite rock sample HPS-M (GPa) from
    Watanabe et al., 2011, in Voigt notation.

    - Mineralogy:
        ``Ol`` (57.7%), ``Atg`` (36.9%), ``Trm`` (4.5%), ``Mgt`` (1.1%)

    - Abbreviation: ``'SP_37'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3000 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.serpentinite_37()[0]
    array([[ 2.0552e+02,  6.6360e+01,  6.2290e+01, -1.0000e-01, -1.4800e+00,
             3.8600e+00],
           [ 6.6360e+01,  1.9579e+02,  6.2530e+01, -3.7000e-01,  2.0000e-01,
             1.5400e+00],
           [ 6.2290e+01,  6.2530e+01,  1.9330e+02, -1.7800e+00, -2.4000e-01,
             8.3000e-01],
           [-1.0000e-01, -3.7000e-01, -1.7800e+00,  6.6170e+01,  1.4700e+00,
            -5.7000e-01],
           [-1.4800e+00,  2.0000e-01, -2.4000e-01,  1.4700e+00,  6.4700e+01,
            -8.4000e-01],
           [ 3.8600e+00,  1.5400e+00,  8.3000e-01, -5.7000e-01, -8.4000e-01,
             6.7830e+01]])
    >>> elast.serpentinite_37()[1]
    3000.0

    """

    rho = 3000.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 205.52;    C[0,1] = 66.36;     C[0,2] = 62.29;     C[0,3] = -0.1;      C[0,4] = -1.48;    C[0,5] = 3.86
    C[1,0] = C[0,1];    C[1,1] = 195.79;    C[1,2] = 62.53;     C[1,3] = -0.37;     C[1,4] = 0.2;      C[1,5] = 1.54
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 193.30;    C[2,3] = -1.78;     C[2,4] = -0.24;    C[2,5] = 0.83
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 66.17;     C[3,4] = 1.47;     C[3,5] = -0.57
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 64.70;    C[4,5] = -0.84
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];   C[5,5] = 67.83

    return C, rho


def serpentinite_80():
    """
    Elastic constants of serpentinite rock sample HKB-B (GPa) from
    Watanabe et al., 2011, in Voigt notation

    - Mineralogy: ``Ol`` (12.0%), ``Atg`` (80.2%), ``Mgt`` (7.8%)

    - Abbreviation: ``'SP_80'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (2800 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.serpentinite_80()[0]
    array([[ 1.9225e+02,  4.9350e+01,  4.1700e+01, -4.5500e+00,  8.0400e+00,
             9.7800e+00],
           [ 4.9350e+01,  1.5690e+02,  4.2360e+01, -6.9100e+00,  7.1000e-01,
             1.8400e+00],
           [ 4.1700e+01,  4.2360e+01,  1.4162e+02, -4.2800e+00,  1.1100e+00,
             1.9000e-01],
           [-4.5500e+00, -6.9100e+00, -4.2800e+00,  5.3480e+01,  1.0000e-02,
            -6.0000e-02],
           [ 8.0400e+00,  7.1000e-01,  1.1100e+00,  1.0000e-02,  5.1910e+01,
            -3.7200e+00],
           [ 9.7800e+00,  1.8400e+00,  1.9000e-01, -6.0000e-02, -3.7200e+00,
             5.9130e+01]])
    >>> elast.serpentinite_80()[1]
    2800.0

    """

    rho = 2800.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 192.25;    C[0,1] = 49.35;     C[0,2] = 41.7;      C[0,3] = -4.55;     C[0,4] = 8.04;      C[0,5] = 9.78
    C[1,0] = C[0,1];    C[1,1] = 156.9;     C[1,2] = 42.36;     C[1,3] = -6.91;     C[1,4] = 0.71;      C[1,5] = 1.84
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 141.62;    C[2,3] = -4.28;     C[2,4] = 1.11;      C[2,5] = 0.19
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 53.48;     C[3,4] = 0.01;      C[3,5] = -0.06
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 51.91;     C[4,5] = -3.72
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 59.13

    return C, rho


def zoisite():
    """
    Elastic constants of zoisite mineral (GPa) from
    Mao et al. (2007), in Voigt notation

    - Abbreviation: ``'zo'``

    Returns:
        (tuple): tuple containing:
            * C (np.ndarray): Elastic stiffness matrix (shape ``(6, 6)``)
            * rho (float): Density (3343 kg/m^3)

    Example
    -------
    >>> from telewavesim import elast
    >>> elast.zoisite()[0]
    array([[279.8,  94.7,  88.7,   0. ,   0. ,   0. ],
           [ 94.7, 249.2,  27.5,   0. ,   0. ,   0. ],
           [ 88.7,  27.5, 209.4,   0. ,   0. ,   0. ],
           [  0. ,   0. ,   0. ,  51.8,   0. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,  81.4,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   0. ,  66.3]])
    >>> elast.zoisite()[1]
    3343.0

    """

    rho = 3343.

    C = np.zeros((6,6), dtype=float)
    C[0,0] = 279.8;     C[0,1] = 94.7;      C[0,2] = 88.7;      C[0,3] = 0.;        C[0,4] = 0.;        C[0,5] = 0.
    C[1,0] = C[0,1];    C[1,1] = 249.2;     C[1,2] = 27.5;      C[1,3] = 0.;        C[1,4] = 0.;        C[1,5] = 0.
    C[2,0] = C[0,2];    C[2,1] = C[1,2];    C[2,2] = 209.4;     C[2,3] = 0.;        C[2,4] = 0.;        C[2,5] = 0.
    C[3,0] = C[0,3];    C[3,1] = C[1,3];    C[3,2] = C[2,3];    C[3,3] = 51.8;      C[3,4] = 0.;        C[3,5] = 0.
    C[4,0] = C[0,4];    C[4,1] = C[1,4];    C[4,2] = C[2,4];    C[4,3] = C[3,4];    C[4,4] = 81.4;      C[4,5] = 0.
    C[5,0] = C[0,5];    C[5,1] = C[1,5];    C[5,2] = C[2,5];    C[5,3] = C[3,5];    C[5,4] = C[4,5];    C[5,5] = 66.3

    return C, rho

