! Copyright 2019 Pascal Audet

! This file is part of Telewavesim.

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!===========================================================================
!
! MODULE conf
!
! Configuration module that contains global variables used in rmat and plane 
! modules to interface with the Python codes.
!
!===========================================================================

      MODULE conf

      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793d0
      INTEGER, PARAMETER :: nlaymx = 30
!
! Model parameters
!
      DOUBLE PRECISION :: a(3,3,3,3,nlaymx), thickn(nlaymx)
      DOUBLE PRECISION :: rho(nlaymx), dp, c, rhof
      INTEGER :: isoflg(nlaymx)
!
! Wavefield parameters
!
      DOUBLE PRECISION :: dt, slow, baz
!
! Eigen values and eigen vectors
!
      DOUBLE COMPLEX :: evals(6,nlaymx), evecs(6,6,nlaymx)
!
! R/T matrices
!
      DOUBLE COMPLEX :: tui(3,3,nlaymx), rui(3,3,nlaymx), &
                        tdi(3,3,nlaymx), rdi(3,3,nlaymx)

      END MODULE conf


!===========================================================================
!
! MODULE rmat
!
! Contains subroutines rmatrix and rtfluid that compute R/T matrices
! for stacks of generally anisotropic layers.
!
!===========================================================================

      MODULE rmat

      CONTAINS

!---------------------------------------------------------------------------
! Subroutine rmatrix
!
! Subroutine to compute Reflection and Transmission matrices for upgoing 
! waves through a stack of layers. Elastic medium can be isotropic or
! anisotropic.
!
! Finds e-vals and e-vecs for an isotropic layer.
! Based on modified (corrected) version of Sean Guest s routine isoeig.
!
! Program to compute r/t matrices for a stack of homogeneous
! layers. Uses reflection matrix algorithm of B. L. N. Kennett.
! General anisotropy permitted.
!
! C. J. Thomson
!
! Original software available here:
! http://sw3d.cz/software/sw3dcd20/rmatrix/rmatrix.htm
!---------------------------------------------------------------------------
      SUBROUTINE rmatrix(px, py, omega, nlay, tus, rus, tds, rds, iter) 

      USE conf

      IMPLICIT NONE

      DOUBLE PRECISION :: px, py
      DOUBLE COMPLEX :: omega
      INTEGER :: nlay, iter
      DOUBLE PRECISION :: aa(3,3,3,3)
!
! Eigenvalue, eigenvector and inverse eigenvector arrays
!
      INTEGER :: ilay
      DOUBLE COMPLEX :: eval1(6), evec1(6,6), eveci1(6,6)
      DOUBLE COMPLEX :: eval2(6), evec2(6,6)
!
! Partitions of wave propagator
!
      DOUBLE COMPLEX :: qq(6,6)
!
! r/t coefficient matrices - single interface and stack of layers
!
      DOUBLE COMPLEX :: tu(3,3), ru(3,3), td(3,3), rd(3,3)
      DOUBLE COMPLEX :: tus(3,3), rus(3,3), tds(3,3), rds(3,3)
!
! Workspace
!
      DOUBLE PRECISION :: r0, r1, vp, vs
      DOUBLE COMPLEX :: d0
!
! Python bindings
!
!f2py DOUBLE PRECISION, intent(in) :: px, py
!f2py DOUBLE COMPLEX, intent(in) :: omega
!f2py INTEGER, intent(in) :: nlay
!f2py DOUBLE COMPLEX, intent(out) :: tus, rus, tds,rds
!f2py INTEGER, intent(in) :: iter

!
! Assign output arrays
!
        r0 = DBLE(0.d0)
        r1 = DBLE(1.d0)
        d0 = DCMPLX(r0,r0)
        tus = d0
        rus = d0
        tds = d0
        rds = d0
!
! Start at bottom layer
!
        ilay = nlay
!
! Only compute e-vecs, etc, if on first pass for this slowness
!     as they are frequency independent
!
        IF (iter.eq.0) THEN
          IF (isoflg(ilay).eq.1) THEN
            vs = SQRT(a(2,3,2,3,ilay))
            vp = SQRT(a(3,3,3,3,ilay))
            CALL isotroc(vp,vs,rho(ilay),px,py,eval1,evec1)
            evals(:,ilay) = eval1
            evecs(:,:,ilay) = evec1
          ELSE
            aa = a(:,:,:,:,ilay)*rho(ilay)
            CALL anisotroc(aa,rho(ilay),px,py,eval1,evec1)
            evals(:,ilay) = eval1
            evecs(:,:,ilay) = evec1
          END IF
        END IF
!
! Invert eigenvalue matrix
!
        CALL xeveci(evecs(:,:,ilay),eveci1)
!
! Next layer and loop point
!
        DO ilay = nlay-1, 1, -1
!
! Get eigenvalues and eigenvectors
!
          IF (iter.eq.0) THEN
            IF (isoflg(ilay).eq.1) THEN
              vs = SQRT(a(2,3,2,3,ilay))
              vp = SQRT(a(3,3,3,3,ilay))
              CALL isotroc(vp,vs,rho(ilay),px,py,eval2,evec2)
              evals(:,ilay) = eval2
              evecs(:,:,ilay) = evec2
            ELSE
              aa = a(:,:,:,:,ilay)*rho(ilay)
              CALL anisotroc(aa,rho(ilay),px,py,eval2,evec2)
              evals(:,ilay) = eval2
              evecs(:,:,ilay) = evec2
            END IF
!
! Scattering matrix for interface
!
            qq = MATMUL(eveci1,evec2)
!
! r/t matrices for interface from scattering matrix
!
            CALL xinv3(qq(4:6,4:6),tu)
            ru = MATMUL(qq(1:3,4:6),tu)
            rd = MATMUL(-tu,qq(4:6,1:3))
            td = qq(1:3,1:3) - MATMUL(MATMUL(qq(1:3,4:6),tu), &
                qq(4:6,1:3))
!
! Store on first pass.
!
            Tui(:,:,ilay) = Tu
            Rui(:,:,ilay) = Ru
            Rdi(:,:,ilay) = Rd
            Tdi(:,:,ilay) = Td

          ELSE
!
! Retrieve on subsequent passes.
!
            Tu = Tui(:,:,ilay)
            Ru = Rui(:,:,ilay)
            Rd = Rdi(:,:,ilay)
            Td = Tdi(:,:,ilay)

          END IF
!
! Addition rule for previous and present r/t matrices
!
          CALL addit(nlay,ilay,tus,rus,tds,rds,tu,ru,td,rd)
!
! Layer up to next interface
!
          IF (ilay.gt.1) THEN

            CALL layer(omega,thickn(ilay),evals(:,ilay), &
                       tus,tds,rds)

            IF (iter.eq.0) THEN
!
! Transfer data
!              
              eval1 = eval2
              evec1 = evec2
!
! Invert eigenvalue matrix
!
              CALL xeveci(evec1,eveci1)

            END IF
          END IF
        END DO
!
! Free surface at top of first layer -- so include phase shift from
!     first interface up to free surface
!
        CALL layer(omega,thickn(1),evals(:,1),tus,tds,rds)

        RETURN

      END SUBROUTINE rmatrix


!---------------------------------------------------------------------------
! Subroutine rtfluid
!
! Subroutine to compute reflection and transmission coefficients for 
! a fluid-solid boundary
!
! p1, p2 are horizontal slowness values, c is fluid velocity, rhof is fluid 
! density, (a, b) are P and S velocities of solid and rhos is solid density
!
! See paper by Bostock & Trehu, BSSA (2012) for details.
!---------------------------------------------------------------------------
      SUBROUTINE rtfluid(p1, p2, a0, b0, rhos, RTmat)

      USE conf

      IMPLICIT NONE

      INTEGER :: info, ipiv(2), ipiv3(3)
      DOUBLE PRECISION :: p1, p2, a0, b0, rhos
!
! r/t coefficients
!
      DOUBLE PRECISION :: tdpp, tdsp, rdpp, tupp, tups, rupp, rusp, rups, russ
      DOUBLE PRECISION :: RTmat(9)
!
! Eigenvalue, eigenvector and F matrices
!
      DOUBLE COMPLEX :: q(6), N(6,6)
      DOUBLE PRECISION :: Fs(4,4), Fl(2,4), L(3,4)
      DOUBLE PRECISION :: Ff(2,2)
!
! Workspace
!
      DOUBLE PRECISION :: pp, qf, dum(3), dum33(3,3)
      DOUBLE PRECISION :: r0, r1
!
! Python bindings
!
!f2py DOUBLE PRECISION, intent(in) :: p1, py, a, b, rhos
!f2py DOUBLE PRECISION, intent(out) :: RTmat(9)

        r0 = DBLE(0.d0)
        r1 = DBLE(1.d0)

        pp = p1*p1 + p2*p2
!
! Fluid properties.
!
        qf = SQRT(1.d0/(c*c) - pp)
        Ff = TRANSPOSE(RESHAPE((/qf*c, -qf*c, rhof*c, rhof*c/),(/2,2/)))
!
! Solid properties (call isotroc for consistency with 
! Rmatrix on conventions of signs etc.)
!
        CALL isotroc(a0, b0, rhos, p1, p2, q, N)
        Fs = TRANSPOSE(REAL(RESHAPE((/N(1,1), N(1,2), N(1,4), N(1,5), &
                      N(3,1), N(3,2), N(3,4), N(3,5), &
                      N(4,1), N(4,2), N(4,4), N(4,5), &
                      N(6,1), N(6,2), N(6,4), N(6,5)/), (/4,4/))))
        Fl(1,:) = Fs(2,:)
        Fl(2,:) = Fs(4,:) 
!
! Define matrix L that encapsulates continuity conditions (continuity of
! vertical displacement and vertical traction) as well as 
! vanishing horizontal traction at fluid solid interface.
!
        CALL dgesv(2,4,Ff,2,ipiv,Fl,2,info)
        L(1:2,:) = Fl(1:2,:)
        L(3,:) = Fs(3,:)

        ! P-incidence from above.
        dum33 = TRANSPOSE(RESHAPE((/L(1,1), L(1,2), r0, &
                          L(2,1), L(2,2), -r1, &
                          L(3,1), L(3,2), r0/), (/3,3/)))
        dum = (/r1, r0, r0/)
        CALL dgesv(3,1,dum33,3,ipiv3,dum,3,info)
                             
        tdpp = dum(1)
        tdsp = dum(2)
        rdpp = dum(3)

        ! P-incidence from below.
        dum33 = TRANSPOSE(RESHAPE((/r0, L(1,1), L(1,2), &
                          -r1, L(2,1), L(2,2), &
                          r0, L(3,1), L(3,2)/), (/3,3/)))
        dum = (/-L(1,3), -L(2,3), -L(3,3)/)
        CALL dgesv(3,1,dum33,3,ipiv3,dum,3,info)

        tupp = dum(1)
        rupp = dum(2)
        rusp = dum(3)

        ! S-incidence from below.
        dum33 = TRANSPOSE(RESHAPE((/r0, L(1,1), L(1,2), &
                          -r1, L(2,1), L(2,2), &
                          r0, L(3,1), L(3,2)/), (/3,3/)))
        dum = (/-L(1,4), -L(2,4), -L(3,4)/)
        CALL dgesv(3,1,dum33,3,ipiv3,dum,3,info)

        tups = dum(1)
        rups = dum(2)
        russ = dum(3)

        RTmat = (/tdpp, tdsp, rdpp, tupp, tups, rupp, rusp, rups, russ/)

        RETURN

      END SUBROUTINE rtfluid

      END MODULE rmat


!===========================================================================
!
! MODULE plane
!
! Contains subroutines plane_land and plane_obs that runs rmat subroutines
! to obtain 3-component seismograms.
!
!===========================================================================

      MODULE plane

      CONTAINS

!---------------------------------------------------------------------------
! Subroutine plane_land
!
! Subroutine to generate plane wave seismograms from a stack of layers
! for land surface stations. Also handles anisotropy. All
! model and time series properties are passed through the 
! configuration module 'conf'. 
!
! Returns displacement traces for given model and slowness vector.
!---------------------------------------------------------------------------
      SUBROUTINE plane_land(nt, nlay, wvtype, yx, yy, yz) 

      USE conf
      USE rmat

      IMPLICIT NONE

      DOUBLE PRECISION :: p1, p2, psi
      DOUBLE COMPLEX :: omg, om0
      INTEGER :: nt, nlay, iter
      DOUBLE PRECISION :: time(nt)
      DOUBLE COMPLEX :: omega(int(nt/2)+1)
      DOUBLE COMPLEX :: tus(3,3), rus(3,3), tds(3,3), rds(3,3)
      DOUBLE COMPLEX :: md(3,3), mu(3,3), nd(3,3), nu(3,3), Ruf0(3,3)
      DOUBLE COMPLEX :: ndi(3,3)
      DOUBLE COMPLEX :: wup(3,int(nt/2)+1), wuv(3,int(nt/2)+1), wuh(3,int(nt/2)+1)
      DOUBLE COMPLEX :: y(3,nt), yx(nt), yy(nt), yz(nt)
      CHARACTER(len=2) :: wvtype

      INTEGER :: i, iw, n2
      DOUBLE PRECISION :: eye(3,3),up(3),uv(3),uh(3)
      DOUBLE COMPLEX :: tmp(3,3), d0
      DOUBLE PRECISION :: r0, r1
!
! Python bindings
!
!f2py DOUBLE COMPLEX, intent(out) :: yx, yy, yz
!f2py INTEGER, intent(in) :: nt, nlay
!f2py intent(in) :: wvtype

        r0 = DBLE(0.d0)
        r1 = DBLE(1.d0)
        d0 = DCMPLX(r0,r0)
        eye = RESHAPE( (/r1, r0, r0, r0, r1, r0, r0, r0, r1/), (/3,3/))
        up = (/r1, r0, r0/)
        uv = (/r0, r1, r0/)
        uh = (/r0, r0, r1/)
        n2 = int(nt/2)
        iter = 0
!
! local parameters
!
        omg = DCMPLX(r1, 0.001d0)
        om0 = DCMPLX(r1, r0)
!
! Frequency and time axes
!
        omega = (/(i, i = 0, n2, 1)/)*omg*2.*pi/(nt*dt)
        time = (/(i, i = 0, nt-1, 1)/)*dt
!
! Slowness vector
!
        psi = (baz - 180.d0)*pi/180.d0
        p1 = slow*1.d-3*dcos(psi)
        p2 = slow*1.d-3*dsin(psi)
!
! Produce R/T matrices for interfaces in solid
!
        CALL rmatrix(p1,p2,om0,nlay,tus,rus,tds,rds,iter)
!
! Get partitions of matrix D (Kennett, p.214)
!
        md = evecs(1:3,1:3,1)
        mu = evecs(1:3,4:6,1)
        nd = evecs(4:6,1:3,1)
        nu = evecs(4:6,4:6,1)
!       
! Free surface matrix. Add +0. to remove negative sign of -0.
! 
        CALL xinv3(nd, ndi)
        Ruf0 = -MATMUL(ndi,nu) + d0
!       
! Initialize upgoing wavefield
!
        wup = d0
        wuv = d0
        wuh = d0
!
! Initialize displacement vector
!
        y = d0
!
! Set up loop for rmatrix
!
        DO iw = 1, n2+1
!
! R/T matrices for solid
!
          CALL rmatrix(p1,p2,omega(iw),nlay,tus,rus,tds,rds,iw)
!
! Upgoing wavefield at free surface
!
          CALL xinv3(eye - MATMUL(rds,ruf0),tmp)
          wup(1:3,iw) = MATMUL(MATMUL(tmp,tus),up)
          wuv(1:3,iw) = MATMUL(MATMUL(tmp,tus),uv)
          wuh(1:3,iw) = MATMUL(MATMUL(tmp,tus),uh)
!
! Displacement vector (Si wave is defined as isotropic source)
!
          IF (wvtype(1:1)=='P') THEN
              y(1:3,iw) = MATMUL(mu + MATMUL(md,Ruf0),wup(1:3,iw))
          ELSEIF (wvtype(1:2)=='Si') THEN
              y(1:3,iw) = MATMUL(mu + MATMUL(md,Ruf0),0.5*(wuv(1:3,iw) &
                          + wuh(1:3,iw)))
          ELSE IF (wvtype(1:2)=='SV') THEN
              y(1:3,iw) = MATMUL(mu + MATMUL(md,Ruf0),wuv(1:3,iw))
          ELSE IF (wvtype(1:2)=='SH') THEN
              y(1:3,iw) = MATMUL(mu + MATMUL(md,Ruf0),wuh(1:3,iw))
          ELSE
              PRINT*,'Error - wave type can only be P, Si, SV or SH'
              RETURN
          END IF
        END DO
!
! Get displacements in time domain
!
        DO iw = n2+2, nt
          y(:,iw) = DCONJG(y(:,nt-iw+2))
        END DO
        yx = y(1,:)
        yy = y(2,:)
        yz = y(3,:)

        RETURN

      END SUBROUTINE plane_land


!---------------------------------------------------------------------------
! Subroutine plane_obs
!
! Subroutine to generate plane wave seismograms from a stack of layers
! for ocean-bottom stations. Also handles anisotropy. All
! model and time series properties are passed through the 
! configuration module 'conf'. 
!
! Returns displacement traces for given model and slowness vector.
!---------------------------------------------------------------------------
      SUBROUTINE plane_obs(nt, nlay, wvtype, yx, yy, yz)

      USE conf
      USE rmat

      IMPLICIT NONE

      DOUBLE PRECISION :: p1, p2, psi
      DOUBLE COMPLEX :: omg, om0
      INTEGER :: nt, nlay, iter
      DOUBLE PRECISION :: time(nt)
      DOUBLE COMPLEX :: omega(int(nt/2)+1)
      DOUBLE COMPLEX :: tus(3,3), rus(3,3), tds(3,3), rds(3,3)
      CHARACTER(len=2) :: wvtype
!
! Workspace
!
      INTEGER :: i, iw, n2, ipiv3(3), info
      DOUBLE PRECISION :: up(3), uv(3), uh(3), rf(3,3), cons(6)
      DOUBLE PRECISION :: a0, b0, rho0, h, qf, pp
      DOUBLE COMPLEX :: md(3,3), mu(3,3), nd(3,3), nu(3,3)
      DOUBLE COMPLEX :: tup(3,3), rdp(3,3), ruh(3,3), eu(3,3), ed(3,3)
      DOUBLE COMPLEX :: tu(3,3), ru(3,3), td(3,3), rd(3,3)
      DOUBLE COMPLEX :: wave(6)
      DOUBLE COMPLEX :: dum3_1(3,3), dum3_2(3,3)
!
! R/T matrices for fluid-solid
!
      DOUBLE PRECISION :: tdpp, tdsp, rdpp, tupp, tups, rupp, rusp, rups
      DOUBLE PRECISION :: russ
      DOUBLE PRECISION :: RTmat(9)
      DOUBLE COMPLEX :: F1(6,6), eye(3,3)
      DOUBLE COMPLEX :: y(6,nt), yx(nt), yy(nt), yz(nt)
!
! Workspace
!
      DOUBLE COMPLEX :: d0, d1, di
      DOUBLE PRECISION :: r0, r1
!
! Python bindings
!
!f2py DOUBLE COMPLEX, intent(out) :: yx, yy, yz
!f2py INTEGER, intent(in) :: nt, nlay
!f2py intent(in) :: wvtype

        r0 = DBLE(0.d0)
        r1 = DBLE(1.d0)
        d0 = DCMPLX(r0, r0)
        d1 = DCMPLX(r1, r0)
        di = DCMPLX(r0, r1)

        eye = RESHAPE((/d1, d0, d0, d0, d1, d0, d0, d0, d1/),(/3,3/))
        up = (/r1, r0, r0/)
        uv = (/r0, r1, r0/)
        uh = (/r0, r0, r1/)
        h = dp
        a0 = SQRT(a(3,3,3,3,1))
        b0 = SQRT(a(2,3,2,3,1))
        rho0 = rho(1)
        n2 = int(nt/2)
        iter = 0
!
! local parameters
!
        omg = DCMPLX(r1, 0.001d0)
        om0 = DCMPLX(r1, r0)
!
! Frequency and time axes
!
        omega = (/(i, i = 0, n2, 1)/)*omg*2.*pi/(nt*dt)
        time = (/(i, i = 0, nt-1, 1)/)*dt
!
! Slowness vector
!
        psi = (baz - 180.d0)*pi/180.d0
        p1 = slow*1.d-3*dcos(psi)
        p2 = slow*1.d-3*dsin(psi)
        pp = p1*p1 + p2*p2
!
! Vertical slowness
!
        qf = SQRT(r1/(c*c) - pp)
!
! Reflection matrix for water column
!
        rf = RESHAPE((/-r1, r0, r0, r0, r0, r0, r0, r0, r0/),(/3,3/))
!
! Wave type
!
        IF (wvtype(1:1)=='P') THEN
          cons = (/r0, r0, r0, r1, r0, r0/)
        ELSE IF (wvtype(1:2)=='Si') THEN
          cons = (/r0, r0, r0, r0, 0.5d0*r1, 0.5d0*r1/)
        ELSE IF (wvtype(1:2)=='SV') THEN
          cons = (/r0, r0, r0, r0, r1, r0/)
        ELSE IF (wvtype(1:2)=='SH') THEN
          cons = (/r0, r0, r0, r0, r0, r1/)
        ELSE 
          PRINT*,'Error - wave type can only be P, Si, SV or SH'
          RETURN
        END IF
!
! Produce R/T matrices for interfaces in solid
!

        CALL rmatrix(p1,p2,om0,nlay,tus,rus,tds,rds,iter)
!
! Get partitions of matrix D (Kennett, p.214)
!
        md = evecs(1:3,1:3,1)
        mu = evecs(1:3,4:6,1)
        nd = evecs(4:6,1:3,1)
        nu = evecs(4:6,4:6,1)
!
! Matrix F1
!
        F1(1:3,1:3) = md
        F1(1:3,4:6) = mu
        F1(4:6,1:3) = nd
        F1(4:6,4:6) = nu
!        
! Produce R/T matrices for ocean bottom
!
        CALL rtfluid(p1,p2,a0,b0,rho0,RTmat)
        tdpp = RTmat(1); tdsp = RTmat(2); rdpp = RTmat(3)
        tupp = RTmat(4); tups = RTmat(5); rupp = RTmat(6)
        rusp = RTmat(7); rups = RTmat(8); russ = RTmat(9)

        tu = DCMPLX(TRANSPOSE(RESHAPE((/tupp, tups, r0, r0, r0, r0, &
             r0, r0, r0/),(/3,3/))))
        td = DCMPLX(TRANSPOSE(RESHAPE((/tdpp, r0, r0, tdsp, r0, r0, &
             r0, r0, r0/),(/3,3/))))
        ru = DCMPLX(TRANSPOSE(RESHAPE((/rupp, rups, r0, rusp, russ, &
             r0, r0, r0, r1/),(/3,3/))))
        rd = DCMPLX(TRANSPOSE(RESHAPE((/rdpp, r0, r0, r0, r0, r0, r0, & 
             r0, r0/),(/3,3/))))
!        
! Assign upgoing and downgoing wavefields
!
        wave = d0
!        
! Assign displacement vector
!
        y = d0
!
! Set up loop for rmatrix
!
        DO iw = 1,n2+1
!
! R/T matrices for solid
!
          CALL rmatrix(p1,p2,omega(iw),nlay,tus,rus,tds,rds,iw)
!        
! Phase income through ocean layer
!
          eu = cdexp(RESHAPE((/di*omega(iw)*qf*h,d0,d0,&
                    d0,d0,d0,d0,d0,d0/),(/3,3/)))
          ed = eu
!       
! Upgoing wavefield on fluid side of OB
!
          dum3_1 = eye - MATMUL(rds,ru)
          dum3_2 = Tus
          CALL zgesv(3,3,dum3_1,3,ipiv3,dum3_2,3,info)
          tup = MATMUL(tu,dum3_2)

          dum3_1 = eye - MATMUL(ru,rds)
          dum3_2 = td
          CALL zgesv(3,3,dum3_1,3,ipiv3,dum3_2,3,info)
          rdp = rd + MATMUL(tu, MATMUL(rds, dum3_2))
!        
! Upgoing wavefield on solid side of OB using wave propgator cast 
! in terms of interfacial R/T coefficients.
!
          dum3_1 = eye - MATMUL(rd,MATMUL(ed,MATMUL(rf,eu)))
          dum3_2 = tu
          CALL zgesv(3,3,dum3_1,3,ipiv3,dum3_2,3,info)
          ruh = ru + MATMUL(td,MATMUL(ed,MATMUL(rf,MATMUL(eu,dum3_2))))

          dum3_1 = eye - MATMUL(rds,ruh)
          dum3_2 = Tus
          CALL zgesv(3,3,dum3_1,3,ipiv3,dum3_2,3,info)

          wave(4:6) = MATMUL(dum3_2,cons(4:6))
          wave(1:3) = MATMUL(ruh,wave(4:6))

          y(:,iw) = MATMUL(F1,wave)

        END DO

!
! Get displacements in time domain
!
        DO iw = n2+2, nt
          y(:,iw) = DCONJG(y(:,nt-iw+2))
        END DO
        yx = y(1,:)
        yy = y(2,:)
        yz = y(3,:)

        RETURN

      END SUBROUTINE plane_obs

      END MODULE plane
