
! Code converted using TO_F90 by Alan Miller
! Date: 2012-03-16  Time: 11:07:48

!> \file
!! General linear algebra routines.
!!
!! ***** Collection of utility routines from V. Blobel *****
!!
!!     V. Blobel, Univ. Hamburg
!!     Numerical subprograms used in MP-II: matrix equations,
!!        and matrix products, double precision
!!
!!     Solution by inversion
!!        SQMINV
!!        SQMINL  for LARGE matrix, use OpenMP (CHK)
!!
!!     Solution by diagonalization
!!        DEVROT, DEVPRT, DEFSOL,DEVINV
!!
!!     Solution by Cholesky decomposition of symmetric matrix
!!        CHOLDC
!!
!!     Solution by Cholesky decomposition of variable-band matrix
!!        VABDEC
!!
!!     Solution by Cholesky decomposition of bordered band matrix
!!        SQMIBB  (CHK)
!!
!!     Matrix/vector products
!!        DBDOT     dot vector product
!!        DBAXPY    multiplication and addition
!!        DBSVX     symmetric matrix vector
!!        DBSVX     LARGE symmetric matrix vector (CHK)
!!        DBGAX     general matrix vector
!!        DBAVAT    AVAT product
!!        DBMPRV    print parameter and matrix
!!        DBPRV     print matrix  (CHK)
!!
!!     Chi^2 cut values
!!        CHINDL
!!
!!     Accurate summation (moved to pede.f90)
!!        ADDSUM
!!
!!     Sorting
!!        HEAPF    heap sort reals direct
!!        SORT1K   sort 1-dim key-array (CHK)
!!        SORT2K   sort 2-dim key-array
!!

!----------------------------------------------------------------------
!> Matrix inversion and solution.
!!
!! Obtain solution of a system of linear equations with symmetric
!! matrix (V * X = B) and the inverse.
!!
!! Method of solution is by elimination selecting the pivot on the
!! diagonal each stage. The rank of the matrix is returned in NRANK.
!! For NRANK ne N, all remaining rows and cols of the resulting
!! matrix V and the corresponding elements of B are set to zero.
!!
!! \param [in,out] V symmetric N-by-N matrix in symmetric storage mode
!!                   (V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, ...),
!!                   replaced by inverse matrix
!! \param [in,out] B N-vector, replaced by solution vector
!! \param [in]     N size of V, B
!! \param [out]    NRANK rank of matrix V
!! \param [out]    DIAG  double precision scratch array
!! \param [out]    NEXT  integer aux array

SUBROUTINE sqminv(v,b,n,nrank,diag,next)   ! matrix inversion
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ij
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jl
    INTEGER :: k
    INTEGER :: kk
    INTEGER :: l
    INTEGER :: last
    INTEGER :: lk
    INTEGER :: next0

    DOUBLE PRECISION, INTENT(IN OUT)         :: v(*)
    DOUBLE PRECISION, INTENT(OUT)            :: b(n)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(OUT)                     :: nrank
    DOUBLE PRECISION, INTENT(OUT)            :: diag(n)
    INTEGER, INTENT(OUT)                     :: next(n)
    DOUBLE PRECISION :: vkk
    DOUBLE PRECISION :: vjk

    DOUBLE PRECISION, PARAMETER :: eps=1.0D-10
    !     ...
    next0=1
    l=1
    DO i=1,n
        next(i)=i+1                ! set "next" pointer
        diag(i)=ABS(v((i*i+i)/2))  ! save abs of diagonal elements
    END DO
    next(n)=-1                  ! end flag

    nrank=0
    DO i=1,n                    ! start of loop
        k  =0
        vkk=0.0D0
  
        j=next0
        last=0
05      IF(j > 0) THEN
            jj=(j*j+j)/2
            IF(ABS(v(jj)) > MAX(ABS(vkk),eps*diag(j))) THEN
                vkk=v(jj)
                k=j
                l=last
            END IF
            last=j
            j=next(last)
            GO TO 05
        END IF
  
        IF(k /= 0) THEN            ! pivot found
            kk=(k*k+k)/2
            IF(l == 0) THEN
                next0=next(k)
            ELSE
                next(l)=next(k)
            END IF
            next(k)=0               ! index is used, reset
            nrank=nrank+1           ! increase rank and ...
            vkk    =1.0/vkk
            v(kk)  =-vkk
            b(k)   =b(k)*vkk
            jk     =kk-k
            jl     =0
            DO j=1,n                ! elimination
                IF(j == k) THEN
                    jk=kk
                    jl=jl+j
                ELSE
                    IF(j < k) THEN
                        jk=jk+1
                    ELSE
                        jk=jk+j-1
                    END IF
                    vjk  =v(jk)
                    v(jk)=vkk*vjk
                    b(j) =b(j)-b(k)*vjk
                    lk   =kk-k
                    DO l=1,j
                        jl=jl+1
                        IF(l == k) THEN
                            lk=kk
                        ELSE
                            IF(l < k) THEN
                                lk=lk+1
                            ELSE
                                lk=lk+l-1
                            END IF
                            v(jl)=v(jl)-v(lk)*vjk
                        END IF
                    END DO
                END IF
            END DO
        ELSE
            DO k=1,n
                IF(next(k) /= 0) THEN
                    b(k)=0.0D0       ! clear vector element
                    DO j=1,k
                        IF(next(j) /= 0) v((k*k-k)/2+j)=0.0D0  ! clear matrix row/col
                    END DO
                END IF
            END DO
            GO TO 10
        END IF
    END DO             ! end of loop
    10   DO ij=1,(n*n+n)/2
        v(ij)=-v(ij)      ! finally reverse sign of all matrix elements
    END DO
END SUBROUTINE sqminv

!> Matrix inversion for LARGE matrices.
!!
!! Like SQMINV, additional parallelization with OpenMP.
!!
!! \param [in,out] V symmetric N-by-N matrix in symmetric storage mode
!!                   (V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, ...),
!!                   replaced by inverse matrix
!! \param [in,out] B N-vector, replaced by solution vector
!! \param [in]     N size of V, B
!! \param [out]    NRANK rank of matrix V
!! \param [out]    DIAG  double precision scratch array
!! \param [out]    NEXT  integer aux array

SUBROUTINE sqminl(v,b,n,nrank,diag,next)   !
    USE mpdef

    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    INTEGER :: last
    INTEGER :: next0

    DOUBLE PRECISION, INTENT(IN OUT)         :: v(*)
    DOUBLE PRECISION, INTENT(OUT)            :: b(n)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(OUT)                     :: nrank
    DOUBLE PRECISION, INTENT(OUT)            :: diag(n)
    INTEGER, INTENT(OUT)                     :: next(n)
    INTEGER*8 :: i8
    INTEGER*8 :: j8
    INTEGER*8 :: jj
    INTEGER*8 :: k8
    INTEGER*8 :: kk
    INTEGER*8 :: kkmk
    INTEGER*8 :: jk
    INTEGER*8 :: jl
    INTEGER*8 :: llk
    INTEGER*8 :: ljl
    
    DOUBLE PRECISION :: vkk
    DOUBLE PRECISION :: vjk

    DOUBLE PRECISION, PARAMETER :: eps=1.0D-10
    !     ...
    next0=1
    l=1
    DO i=1,n
        i8=int8(i)
        next(i)=i+1                ! set "next" pointer
        diag(i)=ABS(v((i8*i8+i8)/2))  ! save abs of diagonal elements
    END DO
    next(n)=-1                  ! end flag

    nrank=0
    DO i=1,n                    ! start of loop
        k  =0
        vkk=0.0D0
        j=next0
        last=0
05      IF(j > 0) THEN
            j8=int8(j)
            jj=(j8*j8+j8)/2
            IF(ABS(v(jj)) > MAX(ABS(vkk),eps*diag(j))) THEN
                vkk=v(jj)
                k=j
                l=last
            END IF
            last=j
            j=next(last)
            GO TO 05
        END IF
  
        IF(k /= 0) THEN            ! pivot found
            k8=int8(k)
            kk=(k8*k8+k8)/2
            kkmk=kk-k8
            IF(l == 0) THEN
                next0=next(k)
            ELSE
                next(l)=next(k)
            END IF
            next(k)=0               ! index is used, reset
            nrank=nrank+1           ! increase rank and ...
            vkk    =1.0/vkk
            v(kk)  =-vkk
            b(k)   =b(k)*vkk
            ! elimination
            jk     =kkmk
            DO j=1,n
                IF(j == k) THEN
                    jk=kk
                ELSE
                    IF(j < k) THEN
                        jk=jk+1
                    ELSE
                        jk=jk+int8(j)-1
                    END IF
                    v(jk)=v(jk)*vkk
                END IF
            END DO
            ! parallelize row loop
            ! slot of 128 'J' for next idle thread
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(JL,JK,L,LJL,LLK,VJK,J8) &
            !$OMP SCHEDULE(DYNAMIC,128)
            DO j=n,1,-1
                j8=int8(j)
                jl=j8*(j8-1)/2
                IF(j /= k) THEN
                    IF(j < k) THEN
                        jk=kkmk+j8
                    ELSE
                        jk=k8+jl
                    END IF
                    vjk  =v(jk)/vkk
                    b(j) =b(j)-b(k)*vjk
                    ljl=jl
                    llk=kkmk
                    DO l=1,MIN(j,k-1)
                        ljl=ljl+1
                        llk=llk+1
                        v(ljl)=v(ljl)-v(llk)*vjk
                    END DO
                    ljl=ljl+1
                    llk=kk
                    DO l=k+1,j
                        ljl=ljl+1
                        llk=llk+l-1
                        v(ljl)=v(ljl)-v(llk)*vjk
                    END DO
                END IF
            END DO
        !$OMP END PARALLEL DO
        ELSE
            DO k=1,n
                k8=int8(k)
                kk=(k8*k8-k8)/2
                IF(next(k) /= 0) THEN
                    b(k)=0.0D0       ! clear vector element
                    DO j=1,k
                        IF(next(j) /= 0) v(kk+int8(j))=0.0D0  ! clear matrix row/col
                    END DO
                END IF
            END DO
            GO TO 10
        END IF
    END DO             ! end of loop
    10   DO jj=1,(int8(n)*int8(n)+int8(n))/2
        v(jj)=-v(jj)      ! finally reverse sign of all matrix elements
    END DO
END SUBROUTINE sqminl

!> Diagonalization.
!!
!!  Determination  of  eigenvalues  and  eigenvectors of
!!  symmetric matrix V by  Householder method
!!
!! \param [in]  n      size of matrix
!! \param [out] diag   diagonal elements
!! \param [out] u      transformation matrix
!! \param [in]  v      symmetric matrix, unchanged
!! \param [out] work   work array
!! \param [out] iwork  work array

SUBROUTINE devrot(n,diag,u,v,work,iwork)   ! diagonalization
    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(OUT)            :: diag(n)
    DOUBLE PRECISION, INTENT(OUT)            :: u(n,n)
    DOUBLE PRECISION, INTENT(IN)             :: v(*)
    DOUBLE PRECISION, INTENT(OUT)            :: work(n)
    INTEGER, INTENT(OUT)                     :: iwork(n)


    INTEGER, PARAMETER :: itmax=30
    DOUBLE PRECISION, PARAMETER :: tol=1.0D-16
    DOUBLE PRECISION, PARAMETER :: eps=1.0D-16

    DOUBLE PRECISION :: f
    DOUBLE PRECISION :: g
    DOUBLE PRECISION :: h
    DOUBLE PRECISION :: sh
    DOUBLE PRECISION :: hh
    DOUBLE PRECISION :: b
    DOUBLE PRECISION :: p
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: s
    DOUBLE PRECISION :: c
    DOUBLE PRECISION :: workd

    INTEGER :: ij
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    INTEGER :: m
    INTEGER :: ll
    !     ...
    !     1. part: symmetric matrix V reduced to tridiagonal from
    ij=0
    DO i=1,n
        DO j=1,i
            ij=ij+1
            u(i,j)=v(ij)  ! copy half of symmetric matirx
        END DO
    END DO

    DO i=n,2,-1
        l=i-2
        f=u(i,i-1)
        g=0.0D0
        IF(l /= 0) THEN
            DO k=1,l
                IF(ABS(u(i,k)) > tol) g=g+u(i,k)*u(i,k)
            END DO
            h=g+f*f
        END IF
        IF(g < tol) THEN   ! G too small
            work(i)=f           ! skip transformation
            h   =0.0D0
        ELSE
            l=l+1
            sh=SQRT(h)
            IF(f >= 0.0D0) sh=-sh
            g=sh
            work(i)=sh
            h=h-f*g
            u(i,i-1)=f-g
            f=0.0D0
            DO j=1,l
                u(j,i)=u(i,j)/h
                g=0.0D0
                !          form element of a u
                DO k=1,j
                    IF(ABS(u(j,k)) > tol.AND.ABS(u(i,k)) > tol) THEN
                        g=g+u(j,k)*u(i,k)
                    END IF
                END DO
                DO k=j+1,l
                    IF(ABS(u(k,j)) > tol.AND.ABS(u(i,k)) > tol) THEN
                        g=g+u(k,j)*u(i,k)
                    END IF
                END DO
                work(j)=g/h
                f=f+g*u(j,i)
            END DO
            !         form k
            hh=f/(h+h)
            !         form reduced a
            DO j=1,l
                f=u(i,j)
                work(j)=work(j)-hh*f
                g=work(j)
                DO k=1,j
                    u(j,k)=u(j,k)-f*work(k)-g*u(i,k)
                END DO
            END DO
        END IF
        diag(i)=h
    END DO

    diag(1)=0.0D0
    work(1)=0.0D0

    !     accumulation of transformation matrices
    DO i=1,n
        IF(diag(i) /= 0.0) THEN
            DO j=1,i-1
                g=0.0D0
                DO k=1,i-1
                    g=g+u(i,k)*u(k,j)
                END DO
                DO k=1,i-1
                    u(k,j)=u(k,j)-g*u(k,i)
                END DO
            END DO
        END IF
        diag(i)=u(i,i)
        u(i,i)=1.0D0
        DO j=1,i-1
            u(i,j)=0.0D0
            u(j,i)=0.0D0
        END DO
    END DO

    !     2. part: diagonalization of tridiagonal matrix
    DO i=2,n
        work(i-1)=work(i)
    END DO
    work(n)=0.0D0
    b=0.0D0
    f=0.0D0

    DO l=1,n
        j=0
        h=eps*(ABS(diag(l))+ABS(work(l)))
        IF(b < h) b=h
        DO m=l,n
            IF(ABS(work(m)) <= b) GO TO 10 ! look for small sub-diagonal element
        END DO
        m=l
10      IF(m == l) GO TO 30
        !      next iteration
20      IF(j == itmax) THEN
            WRITE(*,*) 'DEVROT: Iteration limit reached'
            STOP
        END IF
        j=j+1
        g=diag(l)
        p=(diag(l+1)-g)/(2.0D0*work(l))
        r=SQRT(1.0D0+p*p)
        diag(l)=work(l)
        IF(p < 0.0D0) diag(l)=diag(l)/(p-r)
        IF(p >= 0.0D0) diag(l)=diag(l)/(p+r)
        h=g-diag(l)
        DO i=l+1,n
            diag(i)=diag(i)-h
        END DO
        f=f+h
        !      QL transformation
        p=diag(m)
        c=1.0D0
        s=0.0D0
        DO i=m-1,l,-1          ! reverse loop
            g=c*work(i)
            h=c*p
            IF(ABS(p) >= ABS(work(i))) THEN
                c=work(i)/p
                r=SQRT(1.0D0+c*c)
                work(i+1)=s*p*r
                s=c/r
                c=1.0D0/r
            ELSE
                c=p/work(i)
                r=SQRT(1.0D0+c*c)
                work(i+1)=s*work(i)*r
                s=1.0D0/r
                c=c/r
            END IF
            p=c*diag(i)-s*g
            diag(i+1)=h+s*(c*g+s*diag(i))
            !       form vector
            DO k=1,n
                h=u(k,i+1)
                u(k,i+1)=s*u(k,i)+c*h
                u(k,i)=c*u(k,i)-s*h
            END DO
        END DO
        work(l)=s*p
        diag(l)=c*p
        IF(ABS(work(l)) > b) GO TO 20 ! next iteration
30      diag(l)=diag(l)+f
    END DO
    DO i=1,n
        iwork(i)=i
    END DO

    m=1
40  m=1+3*m    ! determine initial increment
    IF(m <= n) GO TO 40
50  m=m/3
    DO j=1,n-m ! sort with increment M
        l=j
60      IF(diag(iwork(l+m)) > diag(iwork(l))) THEN ! compare
            ll=iwork(l+m)      ! exchange the two index values
            iwork(l+m)=iwork(l)
            iwork(l)=ll
            l=l-m
            IF(l > 0) GO TO 60
        END IF
    END DO
    IF(m > 1) GO TO 50

    DO i=1,n
        IF(iwork(i) /= i) THEN
            !         move vector from position I to the work area
            workd=diag(i)
            DO l=1,n
                work(l)=u(l,i)
            END DO
            k=i
70          j=k
            k=iwork(j)
            iwork(j)=j
            IF(k /= i) THEN
                !            move vector from position K to the (free) position J
                diag(j)=diag(k)
                DO l=1,n
                    u(l,j)=u(l,k)
                END DO
                GO TO 70
            END IF
            !         move vector from the work area to position J
            diag(j)=workd
            DO l=1,n
                u(l,j)=work(l)
            END DO
        END IF
    END DO
END SUBROUTINE devrot

!> Calculate significances.
SUBROUTINE devsig(n,diag,u,b,coef)
    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN)             :: diag(n)
    DOUBLE PRECISION, INTENT(IN)             :: u(n,n)
    DOUBLE PRECISION, INTENT(IN)             :: b(n)
    DOUBLE PRECISION, INTENT(OUT)            :: coef(n)
    INTEGER :: i
    INTEGER :: j
    DOUBLE PRECISION :: dsum
    !     ...
    DO i=1,n
        coef(i)=0.0D0
        IF(diag(i) > 0.0D0) THEN
            dsum=0.0D0
            DO j=1,n
                dsum=dsum+u(j,i)*b(j)
            END DO
            coef(i)=ABS(dsum)/SQRT(diag(i))
        END IF
    END DO
END SUBROUTINE devsig


!> Solution by diagonalization.
!!
!! Solution of matrix equation   V * X = B after diagonalization of V.
!!
!! \param [in]   N       size of matrix
!! \param [in]   DIAG    diagonal elements
!! \param [in]   U       transformation matrix
!! \param [in]   B       r.h.s. of matrix equation (unchanged)
!! \param [out]  X       solution vector
!! \param [out]  WORK    work array

SUBROUTINE devsol(n,diag,u,b,x,work)
    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN)             :: diag(n)
    DOUBLE PRECISION, INTENT(IN)             :: u(n,n)
    DOUBLE PRECISION, INTENT(IN)             :: b(n)
    DOUBLE PRECISION, INTENT(OUT)            :: x(n)
    DOUBLE PRECISION, INTENT(OUT)            :: work(n)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: jj
    DOUBLE PRECISION :: s
    !     ...
    DO j=1,n
        s=0.0D0
        work(j)=0.0D0
        IF(diag(j) /= 0.0D0) THEN
            DO i=1,n
                !           j-th eigenvector is U(.,J)
                s=s+u(i,j)*b(i)
            END DO
            work(j)=s/diag(j)
        END IF
    END DO

    DO j=1,n
        s=0.0D0
        DO jj=1,n
            s=s+u(j,jj)*work(jj)
        END DO
        x(j)=s
    END DO
!      WRITE(*,*) 'DEVSOL'
!      WRITE(*,*) 'X ',X
END SUBROUTINE devsol

!> Inversion by diagonalization.
!! Get inverse matrix V from DIAG and U.
!!
!! \param [in]  N     size of matrix
!! \param [in]  DIAG  diagonal elements
!! \param [in]  U     transformation matrix
!! \param [out] V     smmmetric matrix

SUBROUTINE devinv(n,diag,u,v)

    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ij
    INTEGER :: j
    INTEGER :: k

    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN)             :: diag(n)
    DOUBLE PRECISION, INTENT(IN)             :: u(n,n)
    DOUBLE PRECISION, INTENT(OUT)            :: v(*)
    DOUBLE PRECISION:: dsum
    !     ...
    ij=0
    DO i=1,n
        DO j=1,i
            ij=ij+1
            dsum=0.0D0
            DO k=1,n
                IF(diag(k) /= 0.0D0) THEN
                    dsum=dsum+u(i,k)*u(j,k)/diag(k)
                END IF
            END DO
            v(ij)=dsum
        END DO
    END DO
END SUBROUTINE devinv


!> Cholesky decomposition.
!!
!! Cholesky decomposition of the matrix G:      G =  L  D  L^T
!!
!!  - G = symmetric matrix, in symmetric storage mode
!!
!!  - L = unit triangular matrix (1's on diagonal)
!!
!!  - D = diagonal matrix (elements store on diagonal of L)
!!
!! The sqrts of the usual Cholesky decomposition are avoided by D.
!! Matrices L and D are stored in the place of matrix G; after the
!! decomposition, the solution of matrix equations and the computation
!! of the inverse of the (original) matrix G are done by CHOLSL and CHOLIN.
!!
!! \param [in,out] g symmetric matrix, replaced by D,L
!! \param [in]     n size of matrix
!!
SUBROUTINE choldc(g,n)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ii
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: k
    INTEGER :: kk
    
    DOUBLE PRECISION, INTENT(IN OUT)         :: g(*)
    INTEGER, INTENT(IN)                      :: n

    DOUBLE PRECISION :: ratio
    !     ...
    ii=0
    DO i=1,n
        ii=ii+i
        IF(g(ii) /= 0.0) g(ii)=1.0/g(ii)  ! (I,I) div !
        jj=ii
        DO j=i+1,n
            ratio=g(i+jj)*g(ii)              ! (I,J) (I,I)
            kk=jj
            DO k=j,n
                g(kk+j)=g(kk+j)-g(kk+i)*ratio   ! (K,J) (K,I)
                kk=kk+k
            END DO ! K
            g(i+jj)=ratio                    ! (I,J)
            jj=jj+j
        END DO ! J
    END DO ! I
    RETURN
END SUBROUTINE choldc

!> Solution after decomposition.
!!
!! The matrix equation  G X = B is solved for X, where the matrix
!! G in the argument is already decomposed by CHOLDC. The vector B
!! is called X in the argument and the content is replaced by the
!! resulting vector X.
!!
!! \param [in]      g  decomposed symmetric matrix
!! \param [in,out]  x  r.h.s vector B, replaced by solution vector X
!! \param [in]      n  size of matrix
!!
SUBROUTINE cholsl(g,x,n)
    IMPLICIT NONE
    DOUBLE PRECISION :: dsum
    INTEGER :: i
    INTEGER :: ii
    INTEGER :: k
    INTEGER :: kk
    
    DOUBLE PRECISION, INTENT(IN)            :: g(*)
    DOUBLE PRECISION, INTENT(IN OUT)        :: x(n)
    INTEGER, INTENT(IN)                     :: n

    ii=0
    DO i=1,n
        dsum=x(i)
        DO k=1,i-1
            dsum=dsum-g(k+ii)*x(k)             ! (K,I)
        END DO
        x(i)=dsum
        ii=ii+i
    END DO
    DO i=n,1,-1
        dsum=x(i)*g(ii)                    ! (I,I)
        kk=ii
        DO k=i+1,n
            dsum=dsum-g(kk+i)*x(k)             ! (K,I)
            kk=kk+k
        END DO
        x(i)=dsum
        ii=ii-i
    END DO
    RETURN
END SUBROUTINE cholsl

!> Inversion after decomposition.
!!
!! The inverse of the (original) matrix G is computed and stored
!! in symmetric storage mode in matrix V. Arrays G and V must be
!! different arrays.
!!
!! \param [in]      g  decomposed symmetric matrix
!! \param [in,out]  v  inverse matrix
!! \param [in]      n  size of matrix
!!
SUBROUTINE cholin(g,v,n)
    IMPLICIT NONE
    DOUBLE PRECISION :: dsum
    INTEGER :: i
    INTEGER :: ii
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    INTEGER :: m

    DOUBLE PRECISION, INTENT(IN)            :: g(*)
    DOUBLE PRECISION, INTENT( OUT)          :: v(*)
    INTEGER, INTENT(IN)                     :: n

    ii=(n*n-n)/2
    DO i=n,1,-1
        dsum=g(ii+i)                       ! (I,I)
        DO j=i,1,-1
            DO k=j+1,n
                l=MIN(i,k)
                m=MAX(i,k)
                dsum=dsum-g(j+(k*k-k)/2)*v(l+(m*m-m)/2) ! (J,K) (I,K)
            END DO
            v(ii+j)=dsum                      ! (I,J)
            dsum=0.0D0
        END DO
        ii=ii-i+1
    END DO
END SUBROUTINE cholin

!     variable band matrix operations ----------------------------------

!> Variable band matrix decomposition.
!!
!! Decomposition: A = L D L^T
!!
!! Variable-band matrix row Doolittle decomposition.
!! A variable-band NxN symmetric matrix, also called skyline, is stored
!! row by row in the array VAL(.). For each row every coefficient
!! between the first non-zero element in the row and the diagonal is
!! stored.
!! The pointer array ILPTR(N) contains the indices in VAL(.) of the
!! diagonal elements. ILPTR(1) is always 1, and ILPTR(N) is equal
!! to the total number of coefficients stored, called the profile.
!! The form of a variable-band matrix is preserved in the L D L^T
!! decomposition no fill-in is created ahead in any row or ahead of the
!! first entry in any column, but existing zero-values will become
!! non-zero. The decomposition is done "in-place".
!!
!! \param [in]      n      size of matrix
!! \param [in,out]  val    variable-band matrix, replaced by D,L
!! \param [in]      ilptr  pointer array

SUBROUTINE vabdec(n,val,ilptr)
    IMPLICIT NONE

    INTEGER :: i
    INTEGER :: in
    INTEGER :: j
    INTEGER :: k
    INTEGER :: kj
    INTEGER :: mj
    INTEGER :: mk
    DOUBLE PRECISION :: sn
    DOUBLE PRECISION :: beta
    DOUBLE PRECISION :: delta
    DOUBLE PRECISION :: theta

    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN OUT)         :: val(*)
    INTEGER, INTENT(IN)                      :: ilptr(n)
    
    DOUBLE PRECISION :: dgamma
    DOUBLE PRECISION :: xi
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: prd
    DOUBLE PRECISION :: valkj

    DOUBLE PRECISION, PARAMETER :: one=1.0D0
    DOUBLE PRECISION, PARAMETER :: two=2.0D0
    !     DATA EPS/0.0D0/
    DATA eps/2.22044605E-16/
    SAVE eps
    !     ...
    IF(eps == 0.0D0) THEN
        eps    = two**(-12)
10      eps    = eps/two
        prd=one
        prd=prd+one*eps ! CALL dbaxpy(1,one,eps,prd)
        IF(prd > one) GO TO 10
        eps=eps*two                      ! EPS is machine presision
        WRITE(*,*) 'Machine precision is ',eps
    END IF

    WRITE(*,*) 'Variable band matrix Cholesky decomposition'

    dgamma=0.0D0
    i=1
    DO j=1,ilptr(n) ! loop thrugh all matrix elements
        IF(ilptr(i) == j) THEN ! diagonal element
            IF(val(j) <= 0.0D0) GO TO 01   ! exit loop for negative diag
            dgamma=MAX(dgamma,ABS(val(j)))  ! max diagonal element
            i=i+1
        END IF
    END DO
    i=n+1
01  in=i-1      ! IN positive diagonal elements
    WRITE(*,*) '  ',in,' positive diagonal elements'
    xi=0.0D0
    i=1
    DO j=1,ilptr(in)          ! loop for positive diagonal elements
        ! through all matrix elements
        IF(ilptr(i) == j) THEN   ! diagonal element
            i=i+1
        ELSE
            xi=MAX(xi,ABS(val(j))) ! Xi = abs(max) off-diagonal element
        END IF
    END DO

    delta=eps*MAX(1.0D0,dgamma+xi)
    sn=1.0D0
    IF(n > 1) sn=1.0D0/SQRT(dfloat(n*n-1))
    beta=SQRT(MAX(eps,dgamma,xi*sn))           ! beta
    WRITE(*,*) '   DELTA and BETA ',delta,beta

    DO k=2,n
        mk=k-ilptr(k)+ilptr(k-1)+1
  
        theta=0.0D0
  
        DO j=mk,k
            mj=j-ilptr(j)+ilptr(j-1)+1
            kj=ilptr(k)-k+j        ! index kj
    
            DO i=MAX(mj,mk),j-1
                val(kj)=val(kj) &     ! L_kj := L_kj - L_ki D_ii L_ji
                    -val(ilptr(k)-k+i)*val(ilptr(i))*val(ilptr(j)-j+i)
      
            END DO !
    
            theta=MAX(theta,ABS(val(kj)))  ! maximum value of row
    
            IF(j /= k) THEN
                IF(val(ilptr(j)) /= 0.0D0) THEN
                    val(kj)=val(kj)/val(ilptr(j))
                ELSE
                    val(kj)=0.0D0
                END IF
            END IF                                ! L_kj := L_kj/D_jj ! D_kk
    
            IF(j == k) THEN
                valkj=val(kj)
                IF(k <= in) THEN
                    val(kj)=MAX(ABS(val(kj)),(theta/beta)**2,delta)
                    IF(valkj /= val(kj)) THEN
                        WRITE(*,*) '   Index K=',k
                        WRITE(*,*) '   ',valkj,val(kj), (theta/beta)**2,delta,theta
                    END IF
                END IF
            END IF
        END DO ! J
  
    END DO ! K

    DO k=1,n
        IF(val(ilptr(k)) /= 0.0D0) val(ilptr(k))=1.0D0/val(ilptr(k))
    END DO

    RETURN
END SUBROUTINE vabdec

!>  Variable band matrix print minimum and maximum.
!!
!! \param [in]      n      size of matrix
!! \param [in,out]  val    variable-band matrix, replaced by D,L
!! \param [in]      ilptr  pointer array

SUBROUTINE vabmmm(n,val,ilptr)
    IMPLICIT NONE
    INTEGER :: k
    INTEGER :: kp
    INTEGER :: kr
    INTEGER :: ks
    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN OUT)         :: val(*)
    INTEGER, INTENT(IN)                      :: ilptr(n)
    kr=1
    ks=1
    kp=1
    DO k=1,n
        IF(val(ilptr(k)) > val(ilptr(ks))) ks=k
        IF(val(ilptr(k)) < val(ilptr(kr))) kr=k
        IF(val(ilptr(k)) > 0.0.AND.val(ilptr(k)) < val(ilptr(kp))) kp=k
    END DO
    WRITE(*,*) '   Index value ',ks,val(ilptr(ks))
    WRITE(*,*) '   Index value ',kp,val(ilptr(kp))
    WRITE(*,*) '   Index value ',kr,val(ilptr(kr))

    RETURN
END SUBROUTINE vabmmm

!> Variable band matrix solution.
!!
!! The matrix equation  A X = B  is solved. The matrix is assumed to
!! decomposed before using VABDEC. The array X(N) contains on entry
!! the right-hand-side B(N); at return it contains the solution.
!!
!! \param [in]      n      size of matrix
!! \param [in,out]  val    decomposed variable-band matrix
!! \param [in]      ilptr  pointer array
!! \param [in,out]  x      r.h.s vector B, replaced by solution vector X

SUBROUTINE vabslv(n,val,ilptr,x)
    IMPLICIT NONE
    INTEGER :: j
    INTEGER :: k
    INTEGER :: mk

    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN OUT)         :: val(*)
    INTEGER, INTENT(IN)                      :: ilptr(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: x(n)
    !     ...
    DO k=1,n                  ! forward loop
        mk=k-ilptr(k)+ilptr(k-1)+1
        DO j=mk,k-1
            x(k)=x(k)-val(ilptr(k)-k+j)*x(j)  ! X_k := X_k - L_kj B_j
        END DO
    END DO ! K

    DO k=1,n                  ! divide by diagonal elements
        x(k)=x(k)*val(ilptr(k))            ! X_k := X_k*D_kk
    END DO

    DO k=n,1,-1               ! backward loop
        mk=k-ilptr(k)+ilptr(k-1)+1
        DO j=mk,k-1
            x(j)=x(j)-val(ilptr(k)-k+j)*x(k)  ! X_j := X_j - L_kj X_k
        END DO
    END DO ! K
END SUBROUTINE vabslv

!     matrix/vector products -------------------------------------------

!> Dot product.
!!
!! \param [in] n   vector size
!! \param [in] dx  vector
!! \param [in] dy  vector
!! \return dot product dx*dy

DOUBLE PRECISION FUNCTION dbdot(n,dx,dy)
    IMPLICIT NONE
    INTEGER:: i
    DOUBLE PRECISION :: dtemp

    INTEGER, INTENT(IN)          :: n
    DOUBLE PRECISION, INTENT(IN) :: dx(*)
    DOUBLE PRECISION, INTENT(IN) :: dy(*)
    !     ...
    dtemp=0.0D0
    DO i = 1,MOD(n,5)
        dtemp=dtemp+dx(i)*dy(i)
    END DO
    DO i =MOD(n,5)+1,n,5
        dtemp=dtemp+dx(i)*dy(i)+dx(i+1)*dy(i+1)+dx(i+2)*dy(i+2)  &
            +dx(i+3)*dy(i+3)+dx(i+4)*dy(i+4)
    END DO
    dbdot=dtemp
END FUNCTION dbdot

!> Multiply, addition.
!!
!! Constant times vector added to a vector: DY:=DY+DA*DX
!!
!! \param [in]     n   vector size
!! \param [in]     dx  vector
!! \param [in,out] dy  vector
!! \param [in]     da  scalar

SUBROUTINE dbaxpy(n,da,dx,dy)
    IMPLICIT NONE
    INTEGER:: i

    INTEGER, INTENT(IN)              :: n
    DOUBLE PRECISION, INTENT(IN)     :: dx(*)
    DOUBLE PRECISION, INTENT(IN OUT) :: dy(*)
    DOUBLE PRECISION, INTENT(IN)     :: da
    !     ...
    DO i=1,MOD(n,4)
        dy(i)=dy(i)+da*dx(i)
    END DO
    DO i=MOD(n,4)+1,n,4
        dy(i  )=dy(i  )+da*dx(i  )
        dy(i+1)=dy(i+1)+da*dx(i+1)
        dy(i+2)=dy(i+2)+da*dx(i+2)
        dy(i+3)=dy(i+3)+da*dx(i+3)
    END DO
END SUBROUTINE dbaxpy

!> Product symmetric matrix, vector.
!!
!! Multiply symmetric N-by-N matrix and N-vector.
!!
!! \param[in]  v  symmetric matrix
!! \param[in]  a  vector
!! \param[out] b  vector B = V * A
!! \param[in]  n  size of matrix

SUBROUTINE dbsvx(v,a,b,n)                  !
    IMPLICIT NONE
    INTEGER:: i
    INTEGER:: ij
    INTEGER:: ijs
    INTEGER:: j

    !         B   :=    V   *    A
    !         N        N*N       N

    INTEGER, INTENT(IN)           :: n
    DOUBLE PRECISION, INTENT(IN)  :: v(*)
    DOUBLE PRECISION, INTENT(IN)  :: a(*)
    DOUBLE PRECISION, INTENT(OUT) :: b(*)

    DOUBLE PRECISION:: dsum
    ijs=1
    DO i=1,n
        dsum=0.0
        ij=ijs
        DO j=1,n
            dsum=dsum+v(ij)*a(j)
            IF(j < i) THEN
                ij=ij+1
            ELSE
                ij=ij+j
            END IF
        END DO
        b(i)=dsum
        ijs=ijs+i
    END DO
END SUBROUTINE dbsvx

!> Product LARGE symmetric matrix, vector.
!!
!! Multiply LARGE symmetric N-by-N matrix and N-vector:
!!
!! \param[in]  v  symmetric matrix
!! \param[in]  a  vector
!! \param[out] b  product vector B = V * A
!! \param[in]  n  size of matrix

SUBROUTINE dbsvxl(v,a,b,n)                  ! LARGE symm. matrix, vector
    IMPLICIT NONE
    INTEGER:: i
    INTEGER:: j

    !         B   :=    V   *    A
    !         N        N*N       N

    INTEGER, INTENT(IN)           :: n
    DOUBLE PRECISION, INTENT(IN)  :: v(*)
    DOUBLE PRECISION, INTENT(IN)  :: a(*)
    DOUBLE PRECISION, INTENT(OUT) :: b(*)

    DOUBLE PRECISION:: dsum
    INTEGER*8 :: ij
    INTEGER*8 :: ijs
    ijs=1
    DO i=1,n
        dsum=0.0
        ij=ijs
        DO j=1,n
            dsum=dsum+v(ij)*a(j)
            IF(j < i) THEN
                ij=ij+1
            ELSE
                ij=ij+int8(j)
            END IF
        END DO
        b(i)=dsum
        ijs=ijs+int8(i)
    END DO
END SUBROUTINE dbsvxl

!> Multiply general M-by-N matrix A and N-vector X.
!!
!! \param [in]  A general M-by-N matrix (A11 A12 ... A1N  A21 A22 ...)
!! \param [in]  X N vector
!! \param [out] Y = M vector
!! \param [in]  M rows of A
!! \param [in]  N columns of A

SUBROUTINE dbgax(a,x,y,m,n)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ij
    INTEGER :: j

    DOUBLE PRECISION, INTENT(IN)             :: a(*)
    DOUBLE PRECISION, INTENT(IN)             :: x(*)
    DOUBLE PRECISION, INTENT(OUT)            :: y(*)
    INTEGER, INTENT(IN)                      :: m
    INTEGER, INTENT(IN)                      :: n

    !     ...
    ij=0
    DO i=1,m
        y(i)=0.0D0
        DO j=1,n
            ij=ij+1
            y(i)=y(i)+a(ij)*x(j)
        END DO
    END DO
END SUBROUTINE dbgax

!> A V AT product (similarity).
!!
!! Multiply symmetric N-by-N matrix from the left with general M-by-N
!! matrix and from the right with the transposed of the same general
!! matrix to form symmetric M-by-M matrix (used for error propagation).
!!
!! \param [in]     V  symmetric N-by-N matrix
!! \param [in]     A  general M-by-N matrix
!! \param [in,out] W  symmetric M-by-M matrix
!! \param [in]     MS rows of A (-rows: don't reset W)
!! \param [in]     N  columns of A
!!
SUBROUTINE dbavat(v,a,w,n,ms)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ij
    INTEGER :: ijs
    INTEGER :: il
    INTEGER :: j
    INTEGER :: jk
    INTEGER :: k
    INTEGER :: l
    INTEGER :: lk
    INTEGER :: lkl
    INTEGER :: m

    DOUBLE PRECISION, INTENT(IN)             :: v(*)
    DOUBLE PRECISION, INTENT(IN)             :: a(*)
    DOUBLE PRECISION, INTENT(INOUT)          :: w(*)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(IN)                      :: ms

    DOUBLE PRECISION :: cik
    !     ...
    m=ms
    IF (m > 0) THEN
        DO i=1,(m*m+m)/2
            w(i)=0.0D0             ! reset output matrix
        END DO
    ELSE
        m=-m
    END IF

    il=-n
    ijs=0
    DO i=1,m                 ! do I
        ijs=ijs+i-1             !
        il=il+n                 !
        lkl=0                   !
        DO k=1,n                !   do K
            cik=0.0D0              !
            lkl=lkl+k-1            !
            lk=lkl                 !
            DO l=1,k               !     do L
                lk=lk+1               !     .
                cik=cik+a(il+l)*v(lk) !     .
            END DO                 !     end do L
            DO l=k+1,n             !     do L
                lk=lk+l-1             !     .
                cik=cik+a(il+l)*v(lk) !     .
            END DO                 !     end do L
            jk=k                   !
            ij=ijs                 !
            DO j=1,i               !     do J
                ij=ij+1               !     .
                w(ij)=w(ij)+cik*a(jk) !     .
                jk=jk+n               !     .
            END DO                 !     end do J
        END DO                  !   end do K
    END DO                   ! end do I
END SUBROUTINE dbavat

!> Print symmetric matrix, vector.
!!
!! Prints the n-vector X and the symmetric N-by-N  covariance  matrix
!! V, the latter as a correlation matrix.
!!
!! \param[in]  lun  unit number
!! \param[in]  x    vector
!! \param[in]  v    symmetric matrix
!! \param[in]  n    size of matrix, vector

SUBROUTINE dbmprv(lun,x,v,n)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ii
    INTEGER :: ij
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: l
    INTEGER :: m
    INTEGER :: mc(15)
    REAL :: pd
    REAL :: rho
    REAL :: err

    INTEGER, INTENT(IN)          :: lun
    DOUBLE PRECISION, INTENT(IN) :: x(*)
    DOUBLE PRECISION, INTENT(IN) :: v(*)
    INTEGER, INTENT(IN)          :: n

    WRITE(lun,103)
    WRITE(lun,101)
    ii=0
    DO i=1,n
        ij=ii
        ii=ii+i
        ERR=0.0
        IF(v(ii) > 0.0) ERR=SQRT(REAL(v(ii)))
        l=0
        jj=0
        DO j=1,i
            jj=jj+j
            ij=ij+1
            rho=0.0
            pd=REAL(v(ii)*v(jj))
            IF(pd > 0.0) rho=REAL(v(ij))/SQRT(pd)
            l=l+1
            mc(l)=INT(100.0*ABS(rho)+0.5)
            IF(rho < 0.0) mc(l)=-mc(l)
            IF(j == i.OR.l == 15) THEN
                IF(j <= 15) THEN
                    IF(j == i) THEN
                        WRITE(lun,102) i,x(i),ERR,(mc(m),m=1,l-1)
                    ELSE
                        WRITE(lun,102) i,x(i),ERR,(mc(m),m=1,l)
                    END IF
                ELSE
                    IF(j == i) THEN
                        WRITE(lun,103) (mc(m),m=1,l-1)
                    ELSE
                        WRITE(lun,103) (mc(m),m=1,l)
                    END IF
                    l=0
                END IF
            END IF
        END DO
    END DO
    WRITE(lun,104)
    !  100 RETURN
    RETURN
101 FORMAT(9X,'Param',7X,'error',7X,'correlation coefficients'/)
102 FORMAT(1X,i5,2G12.4,1X,15I5)
103 FORMAT(31X,15I5)
104 FORMAT(33X,'(correlation coefficients in percent)')
END SUBROUTINE dbmprv

!> Print symmetric matrix.
!!
!! Prints the symmetric N-by-N matrix V.
!!
!! \param[in]  lun  unit number
!! \param[in]  v    symmetric matrix
!! \param[in]  n    size of matrix,

SUBROUTINE dbprv(lun,v,n)
    IMPLICIT NONE
    INTEGER, PARAMETER :: istp=6
    INTEGER :: i
    INTEGER :: ip
    INTEGER :: ipe
    INTEGER :: ipn
    INTEGER :: ips
    INTEGER :: k

    INTEGER, INTENT(IN)          :: lun
    DOUBLE PRECISION, INTENT(IN) :: v(*)
    INTEGER, INTENT(IN)          :: n

    WRITE(lun,101)

    DO i=1,n
        ips=(i*i-i)/2
        ipe=ips+i
        ip =ips
100 CONTINUE
    ipn=ip+istp
    WRITE(lun,102), i, ip+1-ips, (v(k),k=ip+1,MIN(ipn,ipe))
    IF (ipn < ipe) THEN
        ip=ipn
        GO TO 100
    END IF
END DO
RETURN
101 FORMAT(1X,'--- DBPRV -----------------------------------')
102 FORMAT(1X,2I3,6G12.4)
END SUBROUTINE dbprv

!     sort -------------------------------------------------------------

!> Heap sort direct (real).
!!
!! Real keys A(*), sorted at return.
!!
!! \param[in,out] a array of keys
!! \param[in]     n number of keys

SUBROUTINE heapf(a,n)
    IMPLICIT NONE
    !
    INTEGER :: i
    INTEGER :: j
    INTEGER :: l
    INTEGER :: r
    REAL :: at    ! pivot key value

    REAL, INTENT(IN OUT) :: a(*)
    INTEGER, INTENT(IN)  :: n
    !     ...
    IF(n <= 1) RETURN
    l=n/2+1
    r=n
10  IF(l > 1) THEN
        l=l-1
        at  =a(l)
    ELSE
        at  =a(r)
        a(r)=a(1)
        r=r-1
        IF(r == 1) THEN
            a(1)=at
            RETURN
        END IF
    END IF
    i=l
    j=l+l
20  IF(j <= r) THEN
        IF(j < r) THEN
            IF(a(j) < a(j+1)) j=j+1
        END IF
        IF(at < a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
        ELSE
            j=r+1
        END IF
        GO TO 20
    END IF
    a(i)=at
    GO TO 10
END SUBROUTINE heapf

!> Quick sort 1.
!!
!! Quick sort of A(1,N) integer.
!!
!! \param[in,out] a vector of integers, sorted at return
!! \param[in]     n size of vector

SUBROUTINE sort1k(a,n)
    IMPLICIT NONE
    INTEGER :: nlev          ! stack size
    PARAMETER (nlev=2*32) ! ... for N = 2**32 = 4.3 10**9
    INTEGER :: i
    INTEGER :: j
    INTEGER :: l
    INTEGER :: r
    INTEGER :: lev
    INTEGER :: lr(nlev)
    INTEGER :: lrh
    INTEGER :: maxlev
    INTEGER :: a1    ! pivot key
    INTEGER :: at    ! pivot key

    INTEGER, INTENT(IN OUT) :: a(*)
    INTEGER, INTENT(IN)     :: n
    !     ...
    IF (n <= 0) RETURN
    maxlev=0
    lev=0
    l=1
    r=n
10  IF(r-l == 1) THEN     ! sort two elements L and R
        IF(a(l) > a(r)) THEN
            at=a(l)       ! exchange L <-> R
            a(l)=a(r)
            a(r)=at
        END IF
        r=l
    END IF
    IF(r == l) THEN
        IF(lev <= 0) THEN
            !            WRITE(*,*) 'SORT1K (quicksort): maxlevel used/available =',
            !     +                 MAXLEV,'/64'
            RETURN
        END IF
        lev=lev-2
        l=lr(lev+1)
        r=lr(lev+2)
    ELSE
        !        LRH=(L+R)/2
        lrh=(l/2)+(r/2)          ! avoid bit overflow
        IF(MOD(l,2) == 1.AND.MOD(r,2) == 1) lrh=lrh+1
        a1=a(lrh)      ! middle
        i=l-1            ! find limits [J,I] with [L,R]
        j=r+1
20      i=i+1
        IF(a(i) < a1) GO TO 20
30      j=j-1
        IF(a(j) > a1) GO TO 30
        IF(i <= j) THEN
            at=a(i)     ! exchange I <-> J
            a(i)=a(j)
            a(j)=at
            GO TO 20
        END IF
        IF(lev+2 > nlev) STOP 'SORT1K (quicksort): stack overflow'
        IF(r-i < j-l) THEN
            lr(lev+1)=l
            lr(lev+2)=j
            l=i
        ELSE
            lr(lev+1)=i
            lr(lev+2)=r
            r=j
        END IF
        lev=lev+2
        maxlev=MAX(maxlev,lev)
    END IF
    GO TO 10
END SUBROUTINE sort1k

!> Quick sort 2.
!!
!! Quick sort of A(2,N) integer.
!!
!! \param[in,out] a vector (pair) of integers, sorted at return
!! \param[in]     n size of vector

SUBROUTINE sort2k(a,n)
    IMPLICIT NONE
    INTEGER:: nlev          ! stack size
    PARAMETER (nlev=2*32) ! ... for N = 2**32 = 4.3 10**9
    INTEGER:: i
    INTEGER::j
    INTEGER::l
    INTEGER::r
    INTEGER::lev
    INTEGER::lr(nlev)
    INTEGER::lrh
    INTEGER::maxlev
    INTEGER::a1       ! pivot key
    INTEGER::a2       ! pivot key
    INTEGER::at       ! pivot key

    INTEGER, INTENT(IN OUT) :: a(2,*)
    INTEGER, INTENT(IN)     :: n
    !     ...
    maxlev=0
    lev=0
    l=1
    r=n
10  IF(r-l == 1) THEN     ! sort two elements L and R
        IF(a(1,l) > a(1,r).OR.( a(1,l) == a(1,r).AND.a(2,l) > a(2,r))) THEN
            at=a(1,l)       ! exchange L <-> R
            a(1,l)=a(1,r)
            a(1,r)=at
            at=a(2,l)
            a(2,l)=a(2,r)
            a(2,r)=at
        END IF
        r=l
    END IF
    IF(r == l) THEN
        IF(lev <= 0) THEN
            WRITE(*,*) 'SORT2K (quicksort): maxlevel used/available =', maxlev,'/64'
            RETURN
        END IF
        lev=lev-2
        l=lr(lev+1)
        r=lr(lev+2)
    ELSE
        !        LRH=(L+R)/2
        lrh=(l/2)+(r/2)          ! avoid bit overflow
        IF(MOD(l,2) == 1.AND.MOD(r,2) == 1) lrh=lrh+1
        a1=a(1,lrh)      ! middle
        a2=a(2,lrh)
        i=l-1            ! find limits [J,I] with [L,R]
        j=r+1
20      i=i+1
        IF(a(1,i) < a1) GO TO 20
        IF(a(1,i) == a1.AND.a(2,i) < a2) GO TO 20
30      j=j-1
        IF(a(1,j) > a1) GO TO 30
        IF(a(1,j) == a1.AND.a(2,j) > a2) GO TO 30
        IF(i <= j) THEN
            at=a(1,i)     ! exchange I <-> J
            a(1,i)=a(1,j)
            a(1,j)=at
            at=a(2,i)
            a(2,i)=a(2,j)
            a(2,j)=at
            GO TO 20
        END IF
        IF(lev+2 > nlev) STOP 'SORT2K (quicksort): stack overflow'
        IF(r-i < j-l) THEN
            lr(lev+1)=l
            lr(lev+2)=j
            l=i
        ELSE
            lr(lev+1)=i
            lr(lev+2)=r
            r=j
        END IF
        lev=lev+2
        maxlev=MAX(maxlev,lev)
    END IF
    GO TO 10
END SUBROUTINE sort2k

!> Chi2/ndf cuts.
!!
!! Return limit in Chi^2/ndf for N sigmas (N=1, 2 or 3).
!!
!! \param[in] n  number of sigmas
!! \param[in] nd ndf
!! \return Chi2/ndf cut value

REAL FUNCTION chindl(n,nd)
    IMPLICIT NONE
    INTEGER :: m
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: nd
    !
    REAL:: sn(3)
    REAL::table(30,3)
    !     REAL PN(3)
    !     DATA PN/0.31731,0.0455002785,2.69985E-3/         ! probabilities
    DATA sn/0.47523,1.690140,2.782170/
    DATA table/ 1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630,  &
        1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321,  &
        1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136,  &
        1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040,  &
        4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736,  &
        1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658,  &
        1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346,  &
        1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742,  &
        9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468,  &
        2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635,  &
        2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901,  &
        1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681/
    SAVE sn,table
    !     ...
    IF(nd < 1) THEN
        chindl=0.0
    ELSE
        m=MAX(1,MIN(n,3))         ! 1, 2 or 3 sigmas
        IF(nd <= 30) THEN
            chindl=table(nd,m)     ! from table
        ELSE                      ! approximation for ND > 30
            chindl=(sn(m)+SQRT(FLOAT(nd+nd-1)))**2/FLOAT(nd+nd)
        END IF
    END IF
END FUNCTION chindl

!> LLT decomposition.
!!
!! Decomposition: C = L L^T.
!!
!! Variable-band matrix row-Doolittle decomposition of pos. def. matrix.
!! A variable-band NxN symmetric matrix, is stored row by row in the
!! array C(.). For each row all coefficients from the first
!! non-zero element in the row to the diagonal is stored.
!! The pointer array INDIA(N) contains the indices in C(.) of the
!! diagonal elements. INDIA(1) is always 1, and INDIA(N) is equal
!! to the total number of coefficients stored, called the profile.
!! The form of a variable-band matrix is preserved in the L D L^T
!! decomposition. No fill-in is created ahead in any row or ahead of the
!! first entry in any column, but existing zero-values will become
!! non-zero. The decomposition is done "in-place".
!!
!!  - NRKD = 0   no component removed
!!
!!  - NRKD < 0   1 component removed, negative index
!!
!!  - NRKD > 1   number of
!!
!! The matrix C is assumed to be positive definite, e.g. from the
!! normal equations of least squares. The (positive) diagonal elements
!! are reduced during decomposition. If a diagonal element is reduced
!! by about a word length (see line "test for linear dependence"),
!! then the pivot is assumed as zero and the entire row/column is
!! reset to zero, removing the corresponding element from the solution.
!!
!! \param [in]      n      size of matrix
!! \param [in,out]  c      variable-band matrix, replaced by L
!! \param [in]      india  pointer array
!! \param [out]     nrkd   removed components

SUBROUTINE lltdec(n,c,india,nrkd)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: kj
    INTEGER :: mj
    INTEGER :: mk
    DOUBLE PRECISION::diag

    INTEGER, INTENT(IN)              :: n
    DOUBLE PRECISION, INTENT(IN OUT) :: c(*)
    INTEGER, INTENT(IN)              :: india(n)
    INTEGER, INTENT(OUT)             :: nrkd
    
    !     ...
    nrkd=0
    diag=0.0D0
    IF(c(india(1)) > 0.0) THEN
        c(india(1))=1.0D0/SQRT(c(india(1))) ! square root
    ELSE
        c(india(1))=0.0D0
        nrkd=-1
    END IF

    DO k=2,n
        mk=k-india(k)+india(k-1)+1    ! first index in row K
        DO j=mk,k                     ! loop over row K with index J
            mj=1
            IF(j > 1) mj=j-india(j)+india(j-1)+1   ! first index in row J
            kj=india(k)-k+j              ! index kj
            diag=c(india(j))             ! j-th diagonal element
    
            DO i=MAX(mj,mk),j-1
                !        L_kj = L_kj - L_ki           *D_ii       *L_ji
                c(kj)=c(kj) - c(india(k)-k+i)*c(india(j)-j+i)
            END DO ! I
    
            IF(j /= k) c(kj)=c(kj)*diag
        END DO ! J
  
        IF(diag+c(india(k)) > diag) THEN      ! test for linear dependence
            c(india(k))=1.0D0/SQRT(c(india(k))) ! square root
        ELSE
            DO j=mk,k                  ! reset row K
                c(india(k)-k+j)=0.0D0
            END DO ! J
            IF(nrkd == 0) THEN
                nrkd=-k
            ELSE
                IF(nrkd < 0) nrkd=1
                nrkd=nrkd+1
            END IF
        END IF
  
    END DO ! K
    RETURN
END SUBROUTINE lltdec


!> Forward solution.
!!
!! The matrix equation  A X = B  is solved by forward + backward
!! solution. The matrix is assumed to
!! decomposed before using LLTDEC. The array X(N) contains on entry
!! the right-hand-side B(N); at return it contains the solution.
!!
!! \param [in]      n      size of matrix
!! \param [in,out]  c      decomposed variable-band matrix
!! \param [in]      india  pointer array
!! \param [in,out]  x      r.h.s vector B, replaced by solution vector X

SUBROUTINE lltfwd(n,c,india,x)
    IMPLICIT NONE
    INTEGER :: j
    INTEGER :: k

    INTEGER, INTENT(IN)              :: n
    DOUBLE PRECISION, INTENT(IN)     :: c(*)
    INTEGER, INTENT(IN)              :: india(n)
    DOUBLE PRECISION, INTENT(IN OUT) :: x(n)

    x(1)=x(1)*c(india(1))
    DO k=2,n                       ! forward loop
        DO j=k-india(k)+india(k-1)+1,k-1
            x(k)=x(k)-c(india(k)-k+j)*x(j)  ! X_k := X_k - L_kj * B_j
        END DO ! J
        x(k)=x(k)*c(india(k))
    END DO ! K
    RETURN
END SUBROUTINE lltfwd

!> Backward solution.
!!
!! The matrix equation  A X = B  is solved by forward + backward
!! solution. The matrix is assumed to
!! decomposed before using LLTDEC. The array X(N) contains on entry
!! the right-hand-side B(N); at return it contains the solution.
!!
!! \param [in]      n      size of matrix
!! \param [in,out]  c      decomposed variable-band matrix
!! \param [in]      india  pointer array
!! \param [in,out]  x      r.h.s vector B, replaced by solution vector X

SUBROUTINE lltbwd(n,c,india,x)
    IMPLICIT NONE
    INTEGER :: j
    INTEGER :: k

    INTEGER, INTENT(IN)              :: n
    DOUBLE PRECISION, INTENT(IN)     :: c(*)
    INTEGER, INTENT(IN)              :: india(n)
    DOUBLE PRECISION, INTENT(IN OUT) :: x(n)

    DO k=n,2,-1                    ! backward loop
        x(k)=x(k)*c(india(k))
        DO j=k-india(k)+india(k-1)+1,k-1
            x(j)=x(j)-c(india(k)-k+j)*x(k)  ! X_j := X_j - L_kj * X_k
        END DO ! J
    END DO ! K
    x(1)=x(1)*c(india(1))
END SUBROUTINE lltbwd

!> Decomposition of equilibrium systems.
!!
!!     N x N matrix C:     starting with sym.pos.def. matrix (N)
!!     length  of array C: INDIA(N) + N*M + (M*M+M)/2
!!     Content of array C: band matrix, as described by INDIA(1)...INDIA(N)
!!            followed by: NxM elements of constraint matrix A
!!            followed by: (M*M+M)/2 unused elements
!!                         INDIA(N+1)...INDIA(N+M) defined internally
!!
!! \param [in]      n      size of symmetric matrix
!! \param [in]      m      number of constrains
!! \param [in,out]  c      combined variable-band + constraints matrix, replaced by decomposition
!! \param [in,out]  india  pointer array
!! \param [out]     nrkd   removed components
!! \param [out]     nrkd2  removed components
!!
SUBROUTINE equdec(n,m,c,india,nrkd,nrkd2)
    IMPLICIT NONE
    REAL:: eps
    INTEGER:: i
    INTEGER:: j
    INTEGER:: jk
    INTEGER:: k
    INTEGER:: ntotal

    INTEGER, INTENT(IN)              :: n
    INTEGER, INTENT(IN)              :: m
    DOUBLE PRECISION, INTENT(IN OUT) :: c(*)
    INTEGER, INTENT(IN OUT)          :: india(n+m)
    INTEGER, INTENT(OUT)             :: nrkd
    INTEGER, INTENT(OUT)             :: nrkd2

    
    PARAMETER (eps=2.22044605E-16)
    !     ...
    ntotal=n+n*m+(m*m+m)/2

    CALL lltdec(n,c,india,nrkd)                  ! decomposition G G^T
    DO i=1,m
        CALL lltfwd(n,c,india,c(india(n)+(i-1)*n+1)) ! forward solution K
    END DO

    jk=india(n)+n*m
    DO j=1,m
        DO k=1,j
            jk=jk+1
            c(jk)=0.0D0                                 ! product K K^T
            DO i=1,n
                c(jk)=c(jk)+c(india(n)+(j-1)*n+i)*c(india(n)+(k-1)*n+i)
            END DO
        END DO
    END DO

    india(n+1)=1
    DO i=2,m
        india(n+i)=india(n+i-1)+MIN(i,m)              ! pointer for K K^T
    END DO

    CALL lltdec(m,c(india(n)+n*m+1),india(n+1),nrkd2)  ! decomp. H H^T

    ntotal=n+n*m+(m*m+m)/2

    RETURN
END SUBROUTINE equdec

!> Solution of equilibrium systems (after decomposition).
!!
!!     N x N matrix C:     starting with sym.pos.def. matrix (N)
!!     length  of array C: INDIA(N) + N*M + (M*M+M)/2
!!     Content of array C: band matrix, as described by INDIA(1)...INDIA(N)
!!            followed by: NxM elements of constraint matrix A
!!            followed by: (M*M+M)/2 unused elements
!!                         INDIA(N+1)...INDIA(N+M) defined internally
!!
!! \param [in]      n      size of symmetric matrix
!! \param [in]      m      number of constrains
!! \param [in]      c      decomposed combined variable-band + constraints matrix
!! \param [in]      india  pointer array
!! \param [in,out]  x      r.h.s vector B, replaced by solution vector X
!!
SUBROUTINE equslv(n,m,c,india,x)                   ! solution vector
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: j

    INTEGER, INTENT(IN)              :: n
    INTEGER, INTENT(IN)              :: m
    DOUBLE PRECISION, INTENT(IN)     :: c(*)
    INTEGER, INTENT(IN)              :: india(n+m)
    DOUBLE PRECISION, INTENT(IN OUT) :: x(n+m)

    CALL lltfwd(n,c,india,x)                           ! result is u
    DO i=1,m
        DO j=1,n
            x(n+i)=x(n+i)-x(j)*c(india(n)+(i-1)*n+j)         ! g - K u
        END DO
    END DO
    CALL lltfwd(m,c(india(n)+n*m+1),india(n+1),x(n+1)) ! result is v


    CALL lltbwd(m,c(india(n)+n*m+1),india(n+1),x(n+1)) ! result is -y
    DO i=1,m
        x(n+i)=-x(n+i)                                    ! result is +y
    END DO

    DO i=1,n
        DO j=1,m
            x(i)=x(i)-x(n+j)*c(india(n)+(j-1)*n+i)           ! u - K^T y
        END DO
    END DO
    CALL lltbwd(n,c,india,x)                           ! result is x
END SUBROUTINE equslv

!> Constrained preconditioner, decomposition.
!!
!!     Constrained preconditioner, e.g for GMRES solution:
!!
!!                                                intermediate
!!        (            )  (   )      (   )           (   )
!!        (   C    A^T )  ( x )   =  ( y )           ( u )
!!        (            )  (   )      (   )           (   )
!!        (   A     0  )  ( l )      ( d )           ( v )
!!
!!     input:
!!        C(N) is diagonal matrix and remains unchanged
!!             may be identical to CU(N), then it is changed
!!        A(N,P) is modified
!!        Y(N+P) is rhs vector, unchanged
!!             may be identical to X(N), then it is changed
!!
!!     result:
!!        CU(N) is 1/sqrt of diagonal matrix C(N)
!!        X(N+P) is result vector
!!        S((P*P+P)/2) is Cholesky decomposed symmetric (P,P) matrix
!!
!! \param [in]     p     number of constraints
!! \param [in]     n     size of diagonal matrix
!! \param [in]     c     diagonal matrix (changed if c=cu as actual parameters)
!! \param [out]    cu    1/sqrt(c)
!! \param [in,out] a     constraint matrix (size n*p), modified
!! \param [out]    s     Cholesky decomposed symmetric (P,P) matrix

SUBROUTINE precon(p,n,c,cu,a,s)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ii
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: k
    INTEGER :: kk

    INTEGER, INTENT(IN)              :: p
    INTEGER, INTENT(IN)              :: n
    DOUBLE PRECISION, INTENT(IN)     :: c(n)
    DOUBLE PRECISION, INTENT(OUT)    :: cu(n)
    DOUBLE PRECISION, INTENT(IN OUT) :: a(n,p)
    DOUBLE PRECISION, INTENT(OUT)    :: s((p*p+p)/2)

    DOUBLE PRECISION :: div
    DOUBLE PRECISION :: ratio
    
    DO i=1,(p*p+p)/2
        s(i)=0.0D0
    END DO
    DO i=1,n
        jk=0
        div=c(i)                          ! copy
        IF (div > 0.0D0) THEN
            cu(i)=1.0D0/DSQRT(div)
        ELSE
            cu(i)=0.0D0
        END IF
        DO j=1,p
            a(i,j)=a(i,j)*cu(i)              ! K = A C^{-1/2}
            DO k=1,j
                jk=jk+1
                s(jk)=s(jk)+a(i,j)*a(i,k)       ! S = symmetric matrix K K^T
            END DO
        END DO
    END DO

    ii=0
    DO i=1,p                           ! S -> H D H^T (Cholesky)
        ii=ii+i
        IF(s(ii) /= 0.0D0) s(ii)=1.0D0/s(ii)
        jj=ii
        DO j=i+1,p
            ratio=s(i+jj)*s(ii)
            kk=jj
            DO k=j,p
                s(kk+j)=s(kk+j)-s(kk+i)*ratio
                kk=kk+k
            END DO ! K
            s(i+jj)=ratio
            jj=jj+j
        END DO ! J
    END DO ! I
    RETURN
END SUBROUTINE precon

!> Constrained preconditioner, solution.
!!
!! \param [in]     p     number of constraints
!! \param [in]     n     size of diagonal matrix
!! \param [in]     cu    1/sqrt(c)
!! \param [in]     a     modified constraint matrix (size n*p)
!! \param [in]     s     Cholesky decomposed symmetric (P,P) matrix
!! \param [out]    x     result vector
!! \param [in]     y     rhs vector (changed if x=y as actual parameters)

SUBROUTINE presol(p,n,cu,a,s,x,y) ! solution
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: k
    INTEGER :: kk

    INTEGER, INTENT(IN)              :: p
    INTEGER, INTENT(IN)              :: n
    
    DOUBLE PRECISION, INTENT(IN)     :: cu(n)
    DOUBLE PRECISION, INTENT(IN)     :: a(n,p)
    DOUBLE PRECISION, INTENT(IN)     :: s((p*p+p)/2)
    DOUBLE PRECISION, INTENT(OUT)    :: x(n+p)
    DOUBLE PRECISION, INTENT(IN)     :: y(n+p)

    DOUBLE PRECISION :: dsum

    DO i=1,n+p
        x(i)=y(i)
    END DO
    DO i=1,n
        x(i)=x(i)*cu(i)                   ! u =C^{-1/2} y
        DO j=1,p
            x(n+j)=x(n+j)-a(i,j)*x(i)        ! d - K u
        END DO
    END DO

    jj=0
    DO j=1,p                           ! Cholesky solution for v
        dsum=x(n+j)
        DO k=1,j-1
            dsum=dsum-s(k+jj)*x(n+k)           ! H v = d - K u
        END DO
        x(n+j)=dsum                        ! -> v
        jj=jj+j
    END DO

    DO j=p,1,-1                        ! solution for lambda
        dsum=x(n+j)*s(jj)
        kk=jj
        DO k=j+1,p
            dsum=dsum+s(kk+j)*x(n+k)           ! D H^T lambda = -v
            kk=kk+k
        END DO
        x(n+j)=-dsum                       ! -> lambda
        jj=jj-j
    END DO

    DO i=1,n                           ! u - K^T lambda
        DO j=1,p
            x(i)=x(i)-a(i,j)*x(n+j)
        END DO
    END DO
    DO i=1,n
        x(i)=x(i)*cu(i)                   ! x = C^{-1/2} u
    END DO

END SUBROUTINE presol


!                                                 090817 C. Kleinwort, DESY-FH1
!> Bordered band matrix.
!!
!! Obtain solution of a system of linear equations with symmetric
!! bordered band matrix (V * X = B), on request inverse is calculated.
!! For band part root-free Cholesky decomposition and forward/backward
!! substitution is used.
!!
!! Use decomposition in border and band part for block matrix algebra:
!!
!!     | A  Ct |   | x1 |   | b1 |        , A  is the border part
!!     |       | * |    | = |    |        , Ct is the mixed part
!!     | C  D  |   | x2 |   | b2 |        , D  is the band part
!!
!! Explicit inversion of D is avoided by using solution X of D*X=C (X=D^-1*C,
!! obtained from Cholesky decomposition and forward/backward substitution)
!!
!!     | x1 |   | E*b1 - E*Xt*b2 |        , E^-1 = A-Ct*D^-1*C = A-Ct*X
!!     |    | = |                |
!!     | x2 |   |  x   - X*x1    |        , x is solution of D*x=b2 (x=D^-1*b2)
!!
!! Inverse matrix is:
!!
!!     |  E   -E*Xt          |
!!     |                     |            , only band part of (D^-1 + X*E*Xt)
!!     | -X*E  D^-1 + X*E*Xt |              is calculated for inv=1
!!
!! \param [in,out] v symmetric N-by-N matrix in symmetric storage mode
!!                   (V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, ...),
!!                   replaced by inverse matrix
!! \param [in,out] b N-vector, replaced by solution vector
!! \param [in]     n size of V, B
!! \param [in]     nbdr   border size
!! \param [in]     nbnd   band width
!! \param [in]     inv    =1 calculate band part of inverse (for pulls),
!!                        >1 calculate complete inverse
!! \param [out]    nrank  rank of matrix V
!! \param [out]    vbnd   band part of V
!! \param [out]    vbdr   border part of V
!! \param [out]    aux    solutions for border rows
!! \param [out]    vbk    matrix for border solution
!! \param [out]    vzru   border solution
!! \param [out]    scdiag workspace (D)
!! \param [out]    scflag workspace (I)
!!
SUBROUTINE sqmibb(v,b,n,nbdr,nbnd,inv,nrank,vbnd,vbdr,aux,vbk,vzru,scdiag,scflag)
    ! Double precision scratch arrays:
    !     VBND(N*(NBND+1)) = storage of band   part
    !     VBDR(N* NBDR)    = storage of border part
    !     AUX (N* NBDR)    = intermediate results

    ! cost[dot ops] ~= (N-NBDR)*(NBDR+NBND+1)**2 + NBDR**3/3 (leading term, solution only)

    IMPLICIT NONE
    INTEGER:: i
    INTEGER:: ib
    INTEGER:: ij
    INTEGER:: ioff
    INTEGER:: ip
    INTEGER:: ip1
    INTEGER:: ip2
    INTEGER:: is
    INTEGER:: j
    INTEGER:: j0
    INTEGER:: jb
    INTEGER:: joff
    INTEGER:: mp1
    INTEGER:: nb1
    INTEGER:: nmb
    INTEGER:: npri
    INTEGER:: nrankb

    DOUBLE PRECISION, INTENT(IN OUT)         :: v(*)
    DOUBLE PRECISION, INTENT(OUT)            :: b(n)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(IN)                      :: nbdr
    INTEGER, INTENT(IN)                      :: nbnd
    INTEGER, INTENT(IN)                      :: inv
    INTEGER, INTENT(OUT)                     :: nrank

    DOUBLE PRECISION, INTENT(OUT) :: vbnd(n*(nbnd+1))
    DOUBLE PRECISION, INTENT(OUT) :: vbdr(n*nbdr)
    DOUBLE PRECISION, INTENT(OUT) :: aux(n*nbdr)
    DOUBLE PRECISION, INTENT(OUT) :: vbk((nbdr*nbdr+nbdr)/2)
    DOUBLE PRECISION, INTENT(OUT) :: vzru(nbdr)
    DOUBLE PRECISION, INTENT(OUT) :: scdiag(nbdr)
    INTEGER, INTENT(OUT)          :: scflag(nbdr)

    SAVE npri
    DATA npri / 100 /
    !           ...
    nrank=0
    nb1=nbdr+1
    mp1=nbnd+1
    nmb=n-nbdr
    !     copy band part
    DO i=nb1,n
        ip=(i*(i+1))/2
        is=0
        DO j=i,MIN(n,i+nbnd)
            ip=ip+is
            is=j
            ib=j-i+1
            vbnd(ib+(i-nb1)*mp1)=v(ip)
        END DO
    END DO
    !     copy border part
    IF (nbdr > 0) THEN
        ioff=0
        DO i=1,nbdr
            ip=(i*(i+1))/2
            is=0
            DO j=i,n
                ip=ip+is
                is=j
                vbdr(ioff+j)=v(ip)
            END DO
            ioff=ioff+n
        END DO
    END IF

    CALL dbcdec(vbnd,mp1,nmb,aux)
    ! use? CALL DBFDEC(VBND,MP1,NMB) ! modified decomp., numerically more stable
    !      CALL DBCPRB(VBND,MP1,NMB)
    ip=1
    DO i=1, nmb
        IF (vbnd(ip) <= 0.0D0) THEN
            npri=npri-1
            IF (npri >= 0) THEN
                IF (vbnd(ip) == 0.0D0) THEN
                    PRINT *, ' SQMIBB matrix singular', n, nbdr, nbnd
                ELSE
                    PRINT *, ' SQMIBB matrix not positive definite', n, nbdr, nbnd
                END IF
            END IF
            !           return zeros
            DO ip=1,n
                b(ip)=0.0D0
            END DO
            DO ip=1,(n*n+n)/2
                v(ip)=0.0D0
            END DO
            RETURN
        END IF
        ip=ip+mp1
    END DO
    nrank=nmb

    IF (nbdr == 0) THEN ! special case NBDR=0
  
        CALL dbcslv(vbnd,mp1,nmb,b,b)
        IF (inv > 0) THEN
            IF (inv > 1) THEN
                CALL dbcinv(vbnd,mp1,nmb,v)
            ELSE
                CALL dbcinb(vbnd,mp1,nmb,v)
            END IF
        END IF
  
    ELSE ! general case NBDR>0
  
        ioff=nb1
        DO ib=1,nbdr
            !           solve for aux. vectors
            CALL dbcslv(vbnd,mp1,nmb,vbdr(ioff),aux(ioff))
            !           zT ru
            vzru(ib)=b(ib)
            DO i=0,nmb-1
                vzru(ib)=vzru(ib)-b(nb1+i)*aux(ioff+i)
            END DO
            ioff=ioff+n
        END DO
        !        solve for band part only
        CALL dbcslv(vbnd,mp1,nmb,b(nb1),b(nb1))
        !        Ck - cT z
        ip=0
        ioff=nb1
        DO ib=1,nbdr
            joff=nb1
            DO jb=1,ib
                ip=ip+1
                vbk(ip)=v(ip)
                DO i=0,nmb-1
                    vbk(ip)=vbk(ip)-vbdr(ioff+i)*aux(joff+i)
                END DO
                joff=joff+n
            END DO
            ioff=ioff+n
        END DO
        !        solve border part
        CALL sqminv(vbk,vzru,nbdr,nrankb,scdiag,scflag)
        IF (nrankb == nbdr) THEN
            nrank=nrank+nbdr
        ELSE
            npri=npri-1
            IF (npri >= 0) PRINT *, ' SQMIBB undef border ', n, nbdr, nbnd, nrankb
            DO ib=1,nbdr
                vzru(ib)=0.0D0
            END DO
            DO ip=(nbdr*nbdr+nbdr)/2,1,-1
                vbk(ip)=0.0D0
            END DO
        END IF
        !        smoothed data points
        ioff=nb1
        DO ib=1, nbdr
            b(ib) = vzru(ib)
            DO i=0,nmb-1
                b(nb1+i)=b(nb1+i)-b(ib)*aux(ioff+i)
            END DO
            ioff=ioff+n
        END DO
        !        inverse requested ?
        IF (inv > 0) THEN
            IF (inv > 1) THEN
                CALL dbcinv(vbnd,mp1,nmb,v)
            ELSE
                CALL dbcinb(vbnd,mp1,nmb,v)
            END IF
            !           expand/correct from NMB to N
            ip1=(nmb*nmb+nmb)/2
            ip2=(n*n+n)/2
            DO i=nmb-1,0,-1
                j0=0
                IF (inv == 1) j0=MAX(0,i-nbnd)
                DO j=i,j0,-1
                    v(ip2)=v(ip1)
                    ioff=nb1
                    DO ib=1,nbdr
                        joff=nb1
                        DO jb=1,nbdr
                            ij=MAX(ib,jb)
                            ij=(ij*ij-ij)/2+MIN(ib,jb)
                            v(ip2)=v(ip2)+vbk(ij)*aux(ioff+i)*aux(joff+j)
                            joff=joff+n
                        END DO
                        ioff=ioff+n
                    END DO
                    ip1=ip1-1
                    ip2=ip2-1
                END DO
                ip1=ip1-j0
                ip2=ip2-j0
      
                DO ib=nbdr,1,-1
                    v(ip2)=0.0D0
                    joff=nb1
                    DO jb=1,nbdr
                        ij=MAX(ib,jb)
                        ij=(ij*ij-ij)/2+MIN(ib,jb)
                        v(ip2)=v(ip2)-vbk(ij)*aux(i+joff)
                        joff=joff+n
                    END DO
                    ip2=ip2-1
                END DO
            END DO
    
            DO ip=(nbdr*nbdr+nbdr)/2,1,-1
                v(ip2)=vbk(ip)
                ip2=ip2-1
            END DO
    
        END IF
    END IF

END SUBROUTINE sqmibb
