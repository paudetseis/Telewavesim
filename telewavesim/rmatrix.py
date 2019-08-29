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
Compute R/T matrices for stacks of generally anisotropic layers.

Original code by C. J. Thomson

'''

import numpy as np
from numpy.linalg import solve, inv, eig
from telewavesim import conf as cf


def rmatrix(px,py,omega,nlay,ifree,iter):
    """
    Computes Reflection and Transmission matrices for upgoing waves
    through a stack of layers. Elastic medium can be isotropic or
    anisotropic.

    Args:
        px (float): x-component of horizontal slowness (s/km)
        py (float): y-component of horizontal slowness (s/km)
        omega (float): Angular frequency (rad/s)
        nlay (int): Number of layers
        ifree (int): Flag for free-surface calculation (always == 1)
        iter (int): Iteration

    Returns:
        (tuple): Tuple containing R/T matrices for up- or doing-going wavefield
            * Tus (np.ndarray): shape ``(3, 3)``
            * Rus (np.ndarray): shape ``(3, 3)``
            * Tds (np.ndarray): shape ``(3, 3)``
            * Rds (np.ndarray): shape ``(3. 3)``

    """

    # Definitions.
    Tds = np.zeros((3,3))
    Rds = np.zeros((3,3))
    Tus = np.zeros((3,3))
    Rus = np.zeros((3,3))
 
    # Start at bottom layer.
    ilay=nlay-1

    # Only compute e-vecs, etc, if on first pass for this slowness
    # as they are frequency independent.
    if iter == -1:
        if cf.isoflg[ilay] == 'iso':
            vs = np.sqrt(cf.a[1,2,1,2,ilay])
            vp = np.sqrt(cf.a[2,2,2,2,ilay])
            eval1, evec1 = isotroc(vp,vs,cf.rho[ilay],px,py)
            cf.evals[:,ilay] = eval1  
            cf.evecs[:,:,ilay] = evec1
        else:
            c = cf.a[:,:,:,:,ilay]*cf.rho[ilay]
            eval1, evec1 = anisotroc(c,cf.rho[ilay],px,py)
            cf.evals[:,ilay] = eval1  
            cf.evecs[:,:,ilay] = evec1

    # Invert eigenvalue matrix.
    eveci1 = xeveci(nlay,ilay,cf.evecs[:,:,ilay])

    # Next layer and loop point
    for ilay in list(reversed(range(0,nlay-1))):

        # Check if on first pass for this slowness.
        if iter == -1:
            if cf.isoflg[ilay] == 'iso':
                vs = np.sqrt(cf.a[1,2,1,2,ilay])
                vp = np.sqrt(cf.a[2,2,2,2,ilay])
                eval2, evec2 = isotroc(vp,vs,cf.rho[ilay],px,py)
                cf.evals[:,ilay] = eval2
                cf.evecs[:,:,ilay] = evec2
            else:
                c = cf.a[:,:,:,:,ilay]*cf.rho[ilay]
                eval2, evec2 = anisotroc(c,cf.rho[ilay],px,py)
                cf.evals[:,ilay] = eval2
                cf.evecs[:,:,ilay] = evec2

            # Scattering matrix for interface.
            Q = np.dot(eveci1,evec2)

            # R/T matrices. 
            Tu = inv(Q[3:6,3:6])
            Ru = np.dot(Q[0:3,3:6],Tu) 
            Rd = np.dot(-Tu,Q[3:6,0:3])
            Td = Q[0:3,0:3] - np.dot(np.dot(Q[0:3,3:6],Tu),Q[3:6,0:3])

            # Store on first pass.
            cf.Tui[:,:,ilay] = Tu
            cf.Rui[:,:,ilay] = Ru
            cf.Rdi[:,:,ilay] = Rd
            cf.Tdi[:,:,ilay] = Td

        else:

            # Retrieve on subsequent passes.
            Tu = cf.Tui[:,:,ilay]
            Ru = cf.Rui[:,:,ilay]
            Rd = cf.Rdi[:,:,ilay]
            Td = cf.Tdi[:,:,ilay]

        # Addition rule for previous and present R/T matrices
        Tus, Rus, Tds, Rds = addit(nlay,ilay,Tus,Rus,Tds,Rds,Tu,Ru,Td,Rd)

        # Layer up to next interface.
        if ilay > 0:
            Tus, Tds, Rds = layer(omega,cf.thickn[ilay],cf.evals[:,ilay],Tus,Tds,Rds)

            if iter == -1:
                eval1 = eval2
                evec1 = evec2
                eveci1 = xeveci(nlay,ilay,evec1)

    # Free surface at top of first layer -- so include phase shift from
    # first interface up to free surface. Note this does not however include
    # free surface reflection coefficients.
    if ifree == 1:

        # Tus, Tds, Rds = layer(omega,cf.thickn[ilay],cf.evals[:,ilay],Tus,Tds,Rds)
        Tus, Tds, Rds = layer(omega,cf.thickn[0],cf.evals[:,0],Tus,Tds,Rds)

    return Tus, Rus, Tds, Rds


def layer(omega,xh,eval,Tu,Td,Rd):
    """
    Includes the phase shift across a layer (Rmatrix theory notes, 
    section 4, equations (4.6)).

    Args:
        omega (float): Angular frequency (rad/s)
        xh (float): Thickness of layer (km)
        eval (float): Eigen-values for layer (shape ``(6)``)
        Tu (np.ndarray): shape ``(3, 3)``
        Td (np.ndarray): shape ``(3, 3)``
        Rd (np.ndarray): shape ``(3, 3)``

    Returns:
        (tuple): Tuple containing:
            Tu (np.ndarray): shape ``(3, 3)``
            Td (np.ndarray): shape ``(3, 3)``
            Rd (np.ndarray): shape ``(3, 3)``

    """

    # Elements of layer matrix for downgoing waves.
    Ed = np.diag(np.exp(1j*omega*eval[0:3]*xh))

    # Elements of layer matrix for downgoing waves.
    Eui = np.diag(np.exp(-1j*omega*eval[3:6]*xh))

    # Update R/T coefficients.
    Td = np.dot(Td,Ed)
    Rd = np.dot(np.dot(Eui,Rd),Ed)
    Tu = np.dot(Eui,Tu)

    return Tu, Td, Rd


def addit(nlay,ilay,Tus,Rus,Tds,Rds,Tu,Ru,Td,Rd):
    """
    Addition/recursion formulas (Rmatrix theory notes, section 4,
    equations (4.8)).

    Args:
        nlay (int): Number of layers
        ilay (int): Index of layer
        Tus (np.ndarray): shape ``(3, 3)`` 
        Rus (np.ndarray): shape ``(3, 3)``
        Tds (np.ndarray): shape ``(3, 3)``
        Rds (np.ndarray): shape ``(3, 3)``
        Tu (np.ndarray): shape ``(3, 3)``
        Ru (np.ndarray): shape ``(3, 3)``
        Td (np.ndarray): shape ``(3, 3)``
        Rd (np.ndarray): shape ``(3, 3)``


    Returns:
        (tuple): Tuple containing:
            Tu (np.ndarray): shape ``(3, 3)``
            Td (np.ndarray): shape ``(3, 3)``
            Rd (np.ndarray): shape ``(3, 3)``

    """

    Ximat = np.identity(3)

    if ilay == nlay-2: 
        Tus = Tu
        Rus = Ru
        Tds = Td
        Rds = Rd
          
        return Tus, Rus, Tds, Rds

    else:
        Rvrb = Ximat - np.dot(Rds,Ru)
        Rvrbi = inv(Rvrb)
        Tusw = np.dot(Tu,solve(Rvrbi,Tus))
        Rdsw = Rd + np.dot(np.dot(Tu,solve(Rvrbi,Rds)),Td)
        Rusw = Rus + np.dot(np.dot(Tds,Ru),solve(Rvrbi,Tus))
        Tdsw = np.dot(np.dot(Tds,(Ximat + np.dot(Ru,solve(Rvrbi,Rds)))),Td)

        Tus = Tusw
        Rus = Rusw
        Tds = Tdsw
        Rds = Rdsw

    return Tus, Rus, Tds, Rds


def xeveci(nlay,ilay,evec):
    """
    Inverts eigen-vector matrix.

    Args:
        nlay (int): Number of layers
        ilay (int): Index of layer
        evec (np.ndarray): Eigen-vectors for layer (shape ``(6, 6)``)

    Returns:
        (np.ndarray): eveci: Inverse of eigen-vector matrix (shape ``(6, 6)``)

    """

    # Create N'*J1. Note just need simple transpose here NOT conjugate.
    eveci_1 = np.hstack((evec[3:6,0:3].transpose(),evec[0:3,0:3].transpose()))
    eveci_2 = np.hstack((evec[3:6,3:6].transpose(),evec[0:3,3:6].transpose()))
    eveci = np.vstack((eveci_1,eveci_2))

    # Determine scaling factors.
    wrk = np.dot(eveci,evec)
 
    # Normalize by 1/diag(wrk).
    for i in range(6):
        eveci[i,:] = eveci[i,:]/wrk[i,i]

    return eveci


def isotroc(vp, vs, rho, p1, p2):
    """ 
    Analytic construction of the fundamental matrix for isotropic
    media from Fryer and Frazer (1987, eq (4.16)) but modified for
    a different Fourier transform sign convention.

    Args:
        vp (float): P-wave velocity (km/s)
        vs (float): S-wave velocity (km/s)
        p1 (float): x-component of horizontal slowness (s/km)
        p2 (float): y-component of horizontal slowness (s/km)

    Returns:
        (tuple): Tuple containing:
            * q (np.ndarray): Eigen-value vector (shape ``(6)``)
            * N (np.ndarray): Eigen-vector matrix (shape ``(6, 6)``)

    """

    # Correction for normal incidence
    if p1 == 0. and p2 == 0.:
        p1 = 1.e-150
    
    mu = rho*vs*vs
    pp = p1*p1 + p2*p2
    qdp = np.sqrt(1./(vp*vp) - pp)
    qds = np.sqrt(1./(vs*vs) - pp)
    qup = -qdp
    qus = -qds

    # Sort eigenvalues and eigenvectors (don't forget normal vs
    # conjugate transpose).
    q = np.array([qdp, qds, qds, qup, qus, qus])

    N50 = rho-2.*mu*pp
    N = np.transpose(np.array(
        [[p1, p2, qdp, 2.*mu*p1*qdp, 2.*mu*p2*qdp, N50],
        [p1, p2, -pp/qds, p1*N50/qds, p2*N50/qds, -2.*mu*pp],
        [-p2, p1, 0., -p2*qds*mu, p1*qds*mu, 0.],
        [p1, p2, qup, 2.*mu*p1*qup, 2.*mu*p2*qup, N50],
        [p1, p2, -pp/qus, p1*N50/qus, p2*N50/qus, -2.*mu*pp],
        [-p2, p1, 0., -p2*qus*mu, p1*qus*mu, 0.]])) + 0.

    # Normalize vectors wrt displacement magnitude.
    for j in range(6):
        norm = np.sqrt(np.dot(N[0:3,j].transpose(),np.conj(N[0:3,j])))
        N[:,j] = N[:,j]/norm

    return q, N


def anisotroc(c,rho,p1,p2):
    """
    Function ANISOTROC takes stiffness tensor c and horizontal components
    of slowness and produces fundamental eigenvector matrix N and diagonal
    matrix of eigenvalues q. The quantities are ordered as qP, qS1, qS2.

    Args:
        c (float): Elastic tensor (GPa) (shape ``(3, 3, 3, 3)``)
        rho (float): Density (kg/m^3) 
        p1 (float): x-component of horizontal slowness (s/km)
        p2 (float): y-component of horizontal slowness (s/km)

    Returns:
        (tuple): Tuple containing:
            * q (np.ndarray): Eigen-value vector (shape ``(6)``)
            * N (np.ndarray): Eigen-vector matrix (shape ``(6, 6)``)

    """

    # Define (Cij)kl. (Cij)kl=c(k,i,l,j) 
    C11 = np.reshape(c[:,0,:,0],(3,3),order='F')
    C12 = np.reshape(c[:,0,:,1],(3,3),order='F')
    C13 = np.reshape(c[:,0,:,2],(3,3),order='F')
    C21 = np.reshape(c[:,1,:,0],(3,3),order='F')
    C22 = np.reshape(c[:,1,:,1],(3,3),order='F')
    C23 = np.reshape(c[:,1,:,2],(3,3),order='F')
    C31 = np.reshape(c[:,2,:,0],(3,3),order='F')
    C32 = np.reshape(c[:,2,:,1],(3,3),order='F')
    C33 = np.reshape(c[:,2,:,2],(3,3),order='F')

    # Construct system matrix A.
    iC33 = inv(C33)
    T = np.dot(-p1*C13-p2*C23,iC33)
    S = rho*np.identity(3) - (p1*p1*(C11 - np.dot(C13,np.dot(iC33,C31))) \
            + p2*p1*(C12 - np.dot(C13,np.dot(iC33,C32))) \
            + p1*p2*(C21 - np.dot(C23,np.dot(iC33,C31))) \
            + p2*p2*(C22 - np.dot(C23,np.dot(iC33,C32))))
    A = np.vstack((np.hstack((T.transpose(),iC33)),np.hstack((S,T))))

    # Diagonalize.
    eva, eve = eig(A)

    # Here you have to be careful that you have properly
    # sorted things. If all eigenvalues are real then a 
    # simple sort will suffice, such that positive ones
    # represent upgoing waves (elements 4:6) while negative 
    # values are downgoing (elements 1:3). If there are 
    # complex values then you must take care to sort correctly.
    # It doesn t really matter how you sort as long as first
    # 3 columns are downgoing/decaying and second 3 are upgoing/
    # growing (in fact can be very difficult to sort 2 qS waves).
    #q0 = np.diag(eva)
    q0 = eva

    N = np.zeros((6,6))
    if np.isreal(q0.all()):
        x = np.sort(q0)
        j = np.argsort(q0)
        q = np.array([x[3],x[4],x[5],x[2],x[1],x[0]])
        N = eve[:,[j[3],j[4],j[5],j[2],j[1],j[0]]]

    else:
        # Never had to deal with this case yet - UNVERIFIED
        imp = np.where(np.imag(q0)>0)[0]
        #imp = find(imag(q0) > 0)
        imn = np.where(np.imag(q0)<0)[0]
        #imn = find(imag(q0) < 0)
        irp = np.where(np.imag(q0)==0 and np.real(q0)>0)
        #irp = find(imag(q0) == 0 & real(q0) > 0)
        irn = np.where(np.imag(q0)==0 and np.real(q0)<0)
        #irn = find(imag(q0) == 0 & real(q0) < 0)
        q = q0[np.concatenate((imp,irp,imn,irn))]
        #q = q0([imp;irp;imn;irn])
        N = eve[:,[np.concatenate((imp,irp,imn,irn))]]
    
    # Normalize vectors wrt displacement magnitude.
    for j in range(6):
        norm = np.sqrt(np.dot(N[0:3,j].transpose(),N[0:3,j]))
        N[:,j] = N[:,j]/norm


    return q, N


