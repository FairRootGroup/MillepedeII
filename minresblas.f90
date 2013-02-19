
! Code converted using TO_F90 by Alan Miller
! Date: 2012-03-16  Time: 11:07:10

!> \file
!! BLAS routines for MINRES.
!!\verbatim
!!     minresblas.f
!!
!!     This file contains Level 1 BLAS from netlib, Thu May 16 1991
!!     (with declarations of the form dx(1) changed to dx(*)):
!!        daxpy    dcopy    ddot
!!     Also
!!        dnrm2    (from NAG,I think).
!!
!!     Also a few utilities to avoid some of the
!!     loops in MINRES (so the debugger can step past them quickly):
!!        daxpy2   dload2   dscal2
!!
!! 15 Jul 2003: dnrm2  is now the NAG version.
!!\endverbatim

!> Constant times a vector plus a vector.
!!
!! Uses unrolled loops for increments equal to one.
!! (jack dongarra, linpack, 3/11/78)
!!
!! \param[in]  n     size of vectors
!! \param[in]  da    scalar constant
!! \param[in]  dx    input vector
!! \param[in]  incx  increment for dx
!! \param[out] dy    output vector
!! \param[in]  incy  increment for dy

SUBROUTINE daxpy(n,da,dx,incx,dy,incy)

    DOUBLE PRECISION :: dx(*),dy(*),da
    INTEGER :: i,incx,incy,ix,iy,m,mp1,n

    IF(n <= 0)RETURN
    IF (da == 0.0D0) RETURN
    IF(incx == 1.AND.incy == 1)GO TO 20

    !        code for unequal increments or equal increments
    !          not equal to 1

    ix = 1
    iy = 1
    IF(incx < 0)ix = (-n+1)*incx + 1
    IF(incy < 0)iy = (-n+1)*incy + 1
    DO  i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
    END DO
    RETURN

    !        code for both increments equal to 1


    !        clean-up loop

20  m = MOD(n,4)
    IF( m == 0 ) GO TO 40
    DO  i = 1,m
        dy(i) = dy(i) + da*dx(i)
    END DO
    IF( n < 4 ) RETURN
40  mp1 = m + 1
    DO  i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
    END DO

END SUBROUTINE daxpy

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Copies a vector, x, to a vector, y.
!!
!! Uses unrolled loops for increments equal to one.
!! (jack dongarra, linpack, 3/11/78)
!!
!! \param[in]  n     size of vectors
!! \param[in]  dx    input vector
!! \param[in]  incx  increment for dx
!! \param[out] dy    output vector
!! \param[in]  incy  increment for dy

SUBROUTINE  dcopy(n,dx,incx,dy,incy)

    DOUBLE PRECISION :: dx(*),dy(*)
    INTEGER :: i,incx,incy,ix,iy,m,mp1,n

    IF(n <= 0)RETURN
    IF(incx == 1.AND.incy == 1)GO TO 20

    !        code for unequal increments or equal increments
    !          not equal to 1

    ix = 1
    iy = 1
    IF(incx < 0)ix = (-n+1)*incx + 1
    IF(incy < 0)iy = (-n+1)*incy + 1
    DO  i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
    END DO
    RETURN

    !        code for both increments equal to 1


    !        clean-up loop

20  m = MOD(n,7)
    IF( m == 0 ) GO TO 40
    DO  i = 1,m
        dy(i) = dx(i)
    END DO
    IF( n < 7 ) RETURN
40  mp1 = m + 1
    DO  i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
    END DO

END SUBROUTINE dcopy

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Forms the dot product of two vectors.
!!
!! Uses unrolled loops for increments equal to one.
!! (jack dongarra, linpack, 3/11/78)
!!
!! \param[in]  n     size of vectors
!! \param[in]  dx    input vector
!! \param[in]  incx  increment for dx
!! \param[in]  dy    input vector
!! \param[in]  incy  increment for dy
!! \return           dot porduct dx*dy

DOUBLE PRECISION FUNCTION ddot(n,dx,incx,dy,incy)

    DOUBLE PRECISION :: dx(*),dy(*),dtemp
    INTEGER :: i,incx,incy,ix,iy,m,mp1,n

    ddot = 0.0D0
    dtemp = 0.0D0
    IF(n <= 0)RETURN
    IF(incx == 1.AND.incy == 1)GO TO 20

    !        code for unequal increments or equal increments
    !          not equal to 1

    ix = 1
    iy = 1
    IF(incx < 0)ix = (-n+1)*incx + 1
    IF(incy < 0)iy = (-n+1)*incy + 1
    DO  i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
    END DO
    ddot = dtemp
    RETURN

    !        code for both increments equal to 1


    !        clean-up loop

20  m = MOD(n,5)
    IF( m == 0 ) GO TO 40
    DO  i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
    END DO
    IF( n < 5 ) GO TO 60
40  mp1 = m + 1
    DO  i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +  &
            dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
    END DO
60  ddot = dtemp

END FUNCTION ddot

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Euclidean vector norm.
!!
!! dnrm2 returns the Euclidean norm of a vector via the function
!! name, so that dnrm2 := sqrt( x'*x ).
!!
!! 15 Jul 2003: dnrm2  obtained from SNOPT src (probably from NAG).
!! s1flmx replaced by safe large number.
!!
!! \param[in]  n     size of vectors
!! \param[in]   x    input vector
!! \param[in]  incx  increment for x
!! \return           vector norm

DOUBLE PRECISION   FUNCTION dnrm2 ( n, x, incx )

    IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
    INTEGER :: incx, n
    DOUBLE PRECISION :: x(*)

    !!!   double precision   s1flmx
    PARAMETER         (one = 1.0D+0, zero = 0.0D+0 )
    DOUBLE PRECISION :: norm
    INTRINSIC          ABS
    !     ------------------------------------------------------------------
    !     flmax = s1flmx( )
    flmax = 1.0D+50

    IF (     n < 1) THEN
        norm  = zero
  
    ELSE IF (n == 1) THEN
        norm  = ABS( x(1) )
  
    ELSE
        scale = zero
        ssq   = one
  
        DO   ix = 1, 1+(n-1)*incx, incx
    
            IF (x(ix) /= zero) THEN
                absxi = ABS( x(ix) )
      
                IF (scale < absxi) THEN
                    ssq   = one + ssq*(scale/absxi)**2
                    scale = absxi
                ELSE
                    ssq   = ssq +     (absxi/scale)**2
                END IF
            END IF
        END DO

        sqt = SQRT( ssq )
        IF (scale < flmax/sqt) THEN
            norm = scale*sqt
        ELSE
            norm = flmax
        END IF
    END IF

    dnrm2  = norm

END FUNCTION dnrm2

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Set  z = a*x + y.
!!
!! 31 May 1999: First version written for MINRES.
!!
!! \param[in]  n    size of vectors
!! \param[in]  a    scalar constant
!! \param[in]  x    input vector
!! \param[in]  y    input vector
!! \param[out] z    output vector

SUBROUTINE daxpy2( n, a, x, y, z )

    IMPLICIT           NONE
    INTEGER :: n
    DOUBLE PRECISION :: a, x(n), y(n), z(n)
    INTEGER :: i

    DO i = 1, n
        z(i) = a*x(i) + y(i)
    END DO

END SUBROUTINE daxpy2

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Set  x = constant.
!!
!! Load all elements of x with const.
!!
!! \param[in]  n      size of vectors
!! \param[in]  const  scalar constant
!! \param[out] x      output vector

SUBROUTINE dload2( n, const, x )

    IMPLICIT           NONE
    INTEGER :: n
    DOUBLE PRECISION :: const, x(n)

    !     ------------------------------------------------------------------
    !     dload2
    !     ------------------------------------------------------------------

    INTEGER :: i

    DO i = 1, n
        x(i) = const
    END DO

END SUBROUTINE dload2

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Set y = a*x.
!!
!! \param[in]  n    size of vectors
!! \param[in]  a    scalar constant
!! \param[in]  x    input vector
!! \param[out] y    output vector

SUBROUTINE dscal2( n, a, x, y )

    IMPLICIT           NONE
    INTEGER :: n
    DOUBLE PRECISION :: a, x(n), y(n)

    INTEGER :: i

    DO i = 1, n
        y(i) = a*x(i)
    END DO

END SUBROUTINE dscal2
