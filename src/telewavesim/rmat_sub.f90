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

!----------------------------------------------------------
! Subroutine layer
!
! Copyright (c) 1996 C. J. Thomson.
! All rights reserved by the author.
!
! Includes the phase shIFt across a layer (theory notes, 
! section 4, equations (4.6)).
!----------------------------------------------------------
      SUBROUTINE layer(omega, xh, eval, tu, td, rd)
      IMPLICIT NONE

      DOUBLE COMPLEX :: omega
      DOUBLE COMPLEX :: xi
      DOUBLE COMPLEX :: eval(6)
      DOUBLE COMPLEX :: tu(3,3), td(3,3), rd(3,3)
      DOUBLE COMPLEX :: ed(3,3), eui(3,3)
      DOUBLE PRECISION :: xh

        xi = DCMPLX(0.d0, 1.d0)
        ed = DCMPLX(0.d0, 0.d0)
        eui = DCMPLX(0.d0, 0.d0)
!
! Elements of layer matrix for downgoing waves
!
        ed(1,1) = CDEXP(xi*omega*eval(1)*xh)
        ed(2,2) = CDEXP(xi*omega*eval(2)*xh)
        ed(3,3) = CDEXP(xi*omega*eval(3)*xh)
!
! Elements of inverse layer matrix for upgoing waves
!
        eui(1,1) = CDEXP(-xi*omega*eval(4)*xh)
        eui(2,2) = CDEXP(-xi*omega*eval(5)*xh)
        eui(3,3) = CDEXP(-xi*omega*eval(6)*xh)
!
! Propagate matrices
!
        td = MATMUL(td,ed)
        rd = MATMUL(MATMUL(eui,rd),ed)
        tu = MATMUL(eui,tu)

        RETURN
      END SUBROUTINE layer


!----------------------------------------------------------
! Subroutine addit
!
! Copyright (c) 1996 C. J. Thomson.
! All rights reserved by the author.
!
! Addition/recursion formulas (see theory notes, section 4,
! equations (4.8)).
!----------------------------------------------------------
      SUBROUTINE addit(nlay, ilay, tus, rus, tds, rds, &
                       tu, ru, td, rd)
      IMPLICIT NONE

      DOUBLE COMPLEX :: tus(3,3), rus(3,3), tds(3,3), rds(3,3)
      DOUBLE COMPLEX :: tu(3,3), ru(3,3), td(3,3), rd(3,3)
      DOUBLE COMPLEX :: tusw(3,3), rusw(3,3), tdsw(3,3), rdsw(3,3)
      DOUBLE COMPLEX :: ximat(3,3), reverb(3,3), reverbi(3,3)
      INTEGER :: ilay, nlay
      DOUBLE COMPLEX :: d0, d1
      DOUBLE PRECISION :: r0, r1

        r0 = DBLE(0.d0)
        r1 = DBLE(1.d0)
        d0 = DCMPLX(r0,r0)
        d1 = DCMPLX(r1,r0)

        Ximat = RESHAPE((/d1, d0, d0, d0, d1, d0, d0, d0, d1/), (/3,3/))

        IF (ilay.eq.nlay-1) THEN
!
! First interface -- trivial
!
          tus = tu
          rus = ru
          tds = td
          rds = rd
        ELSE
!
! Symbolic notation...    
!
!     reverb = ximat - rds * ru
          reverb = ximat - MATMUL(rds,ru)
          CALL xinv3(reverb,reverbi)
!     tusw = tu * reverbi * tus
          tusw = MATMUL(tu,MATMUL(reverb,tus))
!     rdsw = rd + tu * reverbi * rds * td
          rdsw = rd + MATMUL(MATMUL(tu,MATMUL(reverb,rds)),td)
!     rusw = rus + tds * ru * reverbi * tus
          rusw = rus + MATMUL(MATMUL(tds,ru),MATMUL(reverb,tus))
!     tdsw = tds * [ ximat + ru*reverbi*rds ] * td
          tdsw = MATMUL(MATMUL(tds,(ximat + &
            MATMUL(ru,MATMUL(reverb,rds)))),td)
!
! copy over
!
          tus = tusw
          rus = rusw
          tds = tdsw
          rds = rdsw
        END IF

        RETURN
      END SUBROUTINE addit


!----------------------------------------------------------
! Subroutine xeveci
!
! Copyright (c) 1996 C. J. Thomson.
! All rights reserved by the author.
!
! Inverts eigen-vector matrix.
!----------------------------------------------------------
      SUBROUTINE xeveci(evec, eveci)
      IMPLICIT NONE

      DOUBLE COMPLEX :: evec(6,6), eveci(6,6), wrk(6,6)
      INTEGER :: i, j

        DO i = 1, 3
          DO j = 1, 3
            eveci(i,j)=evec(j+3,i)
            eveci(i,j+3)=evec(j,i)
            eveci(i+3,j)=evec(j+3,i+3)
            eveci(i+3,j+3)=evec(j,i+3)
          END DO
        END DO
        wrk = MATMUL(eveci,evec)
        DO i = 1, 6
          DO j = 1, 6
            eveci(i,j) = eveci(i,j)/wrk(i,i)
          END DO
        END DO

        RETURN
      END SUBROUTINE xeveci


!----------------------------------------------------------
! Subroutine isotroc
!
! Calculates eigen values and eigen vectors of propagator 
! matrix where the elasticity tensor is isotropic.
!----------------------------------------------------------
      SUBROUTINE isotroc(vp, vs, dens, p11, p22, q, N)
      IMPLICIT NONE

      DOUBLE PRECISION :: vp, vs, dens, p1, p2
      DOUBLE COMPLEX :: p11, p22
      DOUBLE COMPLEX :: q(6), N(6,6)
      DOUBLE PRECISION :: mu, pp, qdp, qds, qup, qus, N50

        ! IF (p11.eq.0..and.p22.eq.0.) THEN
        !   pp=1.d-20
        ! END IF

        p1 = REAL(p11)
        p2 = REAL(p22)
        mu = dens*vs*vs
        pp = p1*p1 + p2*p2
        qdp = SQRT(1./(vp*vp) - pp)
        qds = SQRT(1./(vs*vs) - pp)
        qup = -qdp
        qus = -qds
!
! Sort eigenvalues and eigenvectors (don't forget normal vs
! conjugate transpose).
!
        q = (/qdp, qds, qds, qup, qus, qus/)

        N50 = dens-2.*mu*pp
        N = DCMPLX(RESHAPE((/p1, p2, qdp, 2.*mu*p1*qdp, 2.*mu*p2*qdp, &
            N50, p1, p2, -pp/qds, p1*N50/qds, p2*N50/qds, -2.*mu*pp, &
            -p2, p1, 0.d0, -p2*qds*mu, p1*qds*mu, 0.d0, &
            p1, p2, qup, 2.*mu*p1*qup, 2.*mu*p2*qup, N50, &
            p1, p2, -pp/qus, p1*N50/qus, p2*N50/qus, -2.*mu*pp, &
            -p2, p1, 0.d0, -p2*qus*mu, p1*qus*mu, 0.d0/), &
            (/6,6/))) 
!
! Normalize wrt displacements
!
        CALL cnorm(6,N)

        RETURN  

      END SUBROUTINE isotroc


!----------------------------------------------------------
! Subroutine anisotroc
!
! Calculates eigen values and eigen vectors of propagator 
! matrix where the elasticity tensor is anisotropic.
!----------------------------------------------------------
      SUBROUTINE anisotroc(c, dens, p11, p22, q, N)
      IMPLICIT NONE

      DOUBLE PRECISION :: dens, c(3,3,3,3), p1, p2
      DOUBLE COMPLEX :: q(6), N(6,6)
      DOUBLE COMPLEX :: p11, p22
      DOUBLE PRECISION :: C11(3,3), C12(3,3), C13(3,3)
      DOUBLE PRECISION :: C21(3,3), C22(3,3), C23(3,3)
      DOUBLE PRECISION :: C31(3,3), C32(3,3), C33(3,3)
      DOUBLE PRECISION :: A(6,6), T(3,3), S(3,3), iC33(3,3)
      DOUBLE PRECISION :: eye(3,3)
!
! Lapack variable for matrix diagonalization
!
      INTEGER :: lwork, ldvl, ldvr, info, i, j
      DOUBLE PRECISION :: vl(6,6), vr(6,6), wr(6), wi(6)
      INTEGER, PARAMETER :: lwmax = 1000
      DOUBLE PRECISION :: work(lwmax)

        eye = RESHAPE( (/1.d0, 0.d0, 0.d0, &
                         0.d0, 1.d0, 0.d0, &
                         0.d0, 0.d0, 1.d0 /), (/3,3/))
        p1 = REAL(p11)
        p2 = REAL(p22)
!
! Get Woodhouse matrices
!        
        C11(1:3,1:3) = c(1:3,1,1:3,1)
        C12(1:3,1:3) = c(1:3,1,1:3,2)
        C13(1:3,1:3) = c(1:3,1,1:3,3)
        C21(1:3,1:3) = c(1:3,2,1:3,1)
        C22(1:3,1:3) = c(1:3,2,1:3,2)
        C23(1:3,1:3) = c(1:3,2,1:3,3)
        C31(1:3,1:3) = c(1:3,3,1:3,1)
        C32(1:3,1:3) = c(1:3,3,1:3,2)
        C33(1:3,1:3) = c(1:3,3,1:3,3)
!
! Construct system matrix A
!
        CALL rxinv3(C33,iC33)
        T = -p1*C13 - p2*C23
        T = MATMUL(T,iC33)

        S = dens*eye - (p1*p1*(C11-MATMUL(C13,MATMUL(iC33,C31))) + &
                        p2*p1*(C12-MATMUL(C13,MATMUL(iC33,C32))) + &
                        p1*p2*(C21-MATMUL(C23,MATMUL(iC33,C31))) + &
                        p2*p2*(C22-MATMUL(C23,MATMUL(iC33,C32))))

        A(1:3,1:3) = TRANSPOSE(T)
        A(1:3,4:6) = iC33
        A(4:6,1:3) = S
        A(4:6,4:6) = T
!
! Query the optimal workspace
!
        ldvl = 6; ldvr = 6
        lwork = -1
        CALL dgeev('N','V',6,A,6,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
        lwork = MIN(lwmax, int(work(1)))
!
! Solve eigenproblem
!
        CALL dgeev('N','V',6,A,6,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
!
! Assign complex eigenvalues and eigenvectors 
! Note:right-handed vectors, to match numpy.linalg.eig 
!
        q = DCMPLX(wr,wi)
        DO i=1,6
          j=1
          DO WHILE (j.le.6)
            IF (wi(j).eq.0.d0) THEN
              N(i,j) = DCMPLX(vr(i,j),0.d0)
              j = j + 1
            ELSE
              N(i,j) = DCMPLX(vr(i,j),vr(i,j+1))
              N(i,j+1) = DCMPLX(vr(i,j),-vr(i,j+1))
              j = j + 2
            END IF
          END DO
        END DO
!
! Sort eigenvalues and eigenvectors and normalize wrt displacements
!
        CALL sort_evec(q,N,6)
        CALL cnorm(6,N)

        RETURN
      END SUBROUTINE anisotroc


!----------------------------------------------------------
! Subroutine xinv3
!
! Copyright (c) 1996 C. J. Thomson.
! All rights reserved by the author.
!
! Complex matrix inversion ainv=inv(a)
! dimension 3.
!----------------------------------------------------------
      SUBROUTINE xinv3(a, ainv)
      IMPLICIT NONE

      DOUBLE COMPLEX :: a(3,3)
      DOUBLE COMPLEX :: ainv(3,3)
      DOUBLE COMPLEX :: co(3,3), deta, detb, detc
      INTEGER :: i, j
!
! Inversion of 3x3 matrix
!
        co(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))
        co(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
        co(1,3) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))

        co(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
        co(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))
        co(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))

        co(3,1) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))
        co(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
        co(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))

        deta = a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
        detb = a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
        detc = a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

        DO i = 1,3
          DO j = 1,3
            ainv(i,j) = co(j,i)/deta
          END DO
        END DO

      END SUBROUTINE xinv3


!----------------------------------------------------------
! Subroutine rxinv3
!
! Copyright (c) 1996 C. J. Thomson.
! All rights reserved by the author.
!
! Real matrix inversion ainv=inv(a)
! dimension 3.
!----------------------------------------------------------
      SUBROUTINE rxinv3(a, ainv)
      IMPLICIT NONE

      DOUBLE PRECISION :: a(3,3)
      DOUBLE PRECISION :: ainv(3,3)
      DOUBLE PRECISION :: co(3,3), deta, detb, detc
      INTEGER :: i, j
!
! Inversion of 3x3 matrix
!
        co(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))
        co(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
        co(1,3) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))

        co(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
        co(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))
        co(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))

        co(3,1) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))
        co(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
        co(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))

        deta = a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
        detb = a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
        detc = a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

        DO i = 1,3
          DO j = 1,3
            ainv(i,j) = co(j,i)/deta
          END DO
        END DO

      END SUBROUTINE rxinv3

!----------------------------------------------------------
! Subroutine xmult3
!
! Copyright (c) 1996 C. J. Thomson.
! All rights reserved by the author.
!
! Complex matrix multiplication a*b=c
! dimension 3 ... allows overwriting.
!----------------------------------------------------------
      SUBROUTINE xmult3(a, b, c)
      IMPLICIT NONE

      DOUBLE COMPLEX :: a(3,3), b(3,3), c(3,3), wrk(3,3)
      INTEGER :: i, j, k

        DO i = 1,3
          DO j = 1,3
            wrk(i,j) = DCMPLX(0.d0, 0.d0)
            DO k = 1,3
              wrk(i,j) = wrk(i,j) + a(i,k)*b(k,j)
            END DO
          END DO
        END DO

        DO i = 1,3
          DO j = 1,3
            c(i,j) = wrk(i,j)
          END DO
        END DO

        RETURN

      END SUBROUTINE xmult3


!----------------------------------------------------------
! Subroutine sort_evec
!
! Sorts eigenvalues and eigenvectors in the same manner 
! as the MATLAB anisotroc.m routine. If all evals are real, 
! the order will be like [1 2 3 -1 -2 3]; otherwise, 
! evanescent waves go first, e.g. [2i,1+2i, 1 -2i, -1-2i, -1].
!----------------------------------------------------------
        SUBROUTINE sort_evec(q, evec, n)
        IMPLICIT NONE

        INTEGER :: n, i, j, indx(n)
        INTEGER :: imagpos(n), imagneg(n)
        INTEGER :: realpos(n), realneg(n)
        INTEGER :: nrp, nrn, nip, nin
        DOUBLE PRECISION :: evalr(n), evali(n), ztol
        DOUBLE PRECISION :: evecr(n,n), eveci(n,n)
        DOUBLE COMPLEX :: q(n), evec(n,n)

          ztol = 1.d-7
          nrp = 0; realpos = 0
          nrn = 0; realneg = 0
          nip = 0; imagpos = 0
          nin = 0; imagneg = 0
!
! Separate into realN and imaginary
!
          evalr = REAL(q)
          evali = AIMAG(q)
          evecr = REAL(evec)
          eveci = AIMAG(evec)
!
! Divide eigenvalues up into real positive, real negative, 
! complex positive, complex negative
!
          DO i = 1, n
            ! is it real?
            IF (abs(evali(i)/evalr(i)) .lt. ztol) THEN
                DO j = 1, n
                  eveci(j,i) = 0.d0
                END DO
                IF (evalr(i).ge.0.d0) THEN
                    nrp = nrp + 1
                    realpos(nrp) = i
                ELSE
                    nrn = nrn + 1
                    realneg(nrn) = i
                END IF
            ELSE
                IF (evali(i).ge.0.d0) THEN
                    nip = nip + 1
                    imagpos(nip) = i
                ELSE
                    nin = nin + 1
                    imagneg(nin) = i
                END IF
            END IF
          END DO
!
! Sort sub-groups
!
          CALL sort_r(evalr,realpos,n,nrp)
          CALL sort_r(evalr,realneg,n,nrn)
          CALL sort_r(evali,imagpos,n,nip)
          CALL sort_r(evali,imagneg,n,nin)
!
! Assemble sub-ranges
!
          DO i = 1, nip
            indx(i) = imagpos(i)
          END DO
          DO i = 1, nrp
            indx(i+nip) = realpos(i)
          END DO
          DO i = 1, nin
            indx(i+nip+nrp) = imagneg(nin-i+1)
          END DO
          DO i = 1, nrn
            indx(i+nip+nrp+nin) = realneg(nrn-i+1)
          END DO

          CALL reorder_evec(evalr,evali,evecr,eveci,indx,q,evec,n)

          RETURN

        END SUBROUTINE sort_evec


!----------------------------------------------------------
! Subroutine reorder_evec
!
! Reorders eigenvectors based on ordered indices.
!----------------------------------------------------------
      SUBROUTINE reorder_evec(evalr, evali, evecr, eveci, indx, &
                 eval, evec, n)
      IMPLICIT NONE

      INTEGER :: n, indx(n), j, k
      DOUBLE PRECISION :: evalr(n), evali(n), evecr(n,n), eveci(n,n)
      DOUBLE COMPLEX :: eval(n), evec(n,n)
        
        DO j = 1, n
          eval(j) = DCMPLX(evalr(indx(j)),evali(indx(j)))
          DO k = 1, n
            evec(k,j) = DCMPLX(evecr(k,indx(j)),eveci(k,indx(j)))
          END DO
        END DO

        RETURN
        
      END SUBROUTINE reorder_evec


!----------------------------------------------------------
! Subroutine sort_r
!
! Sort an array of n real-valued elements, from 
! smallest (most negative) to largest (most positive). 
! Actual elements aren't moved -- indx returns new positions. 
! Uses a basic algorithm (insertion sort) since we're sorting 
! few elements. m is the number of elements in a; n is the 
! number of elements in indx (n <= m).
!
! Indx should be pre-initialized with initial locations.
!----------------------------------------------------------
        SUBROUTINE sort_r(a, indx, m, n)
        IMPLICIT NONE

        INTEGER :: m, n
        INTEGER :: indx(n), j, pos, i, k
        DOUBLE PRECISION :: a(m), val
        
          IF (n .le. 1) THEN
            RETURN
          END IF
!
! Loop over k to ensure recursivity
! Without it, values ordered 3, 2, 1 
! would be re-ordered 3, 2, 1 -> 2, 3, 1 -> 2, 1, 3
!
          DO k = 1, n - 1
            DO j = 2, n
              i = indx(j)
              val = a(i)          
              pos = j
              IF (pos .gt. 1) THEN
                IF (val .lt. a(indx(pos-1))) THEN
                  indx(pos) = indx(pos-1)
                  pos = pos-1
                END IF
              END IF
              indx(pos) = i
            END DO
          END DO

          RETURN

        END SUBROUTINE sort_r
        

!----------------------------------------------------------
! Subroutine cnorm
!
! Calculates norm of complex eigen-vector with respect
! to first 3 vectors (i.e., displacements in the case of
! rmatrix).
!----------------------------------------------------------
        SUBROUTINE cnorm(nm, N)
        IMPLICIT NONE

        INTEGER :: nm, i
        DOUBLE COMPLEX :: N(nm,nm)
        DOUBLE PRECISION :: norm

          DO i = 1, nm
           norm = DSQRT(abs(SUM(N(1:3,i)*DCONJG(N(1:3,i)))))
           N(1:6,i) = N(1:6,i)/norm
          END DO

          RETURN

        END SUBROUTINE cnorm
