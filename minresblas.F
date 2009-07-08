
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     minresblas.f
*
*     This file contains Level 1 BLAS from netlib, Thu May 16 1991
*     (with declarations of the form dx(1) changed to dx(*)):
*        daxpy    dcopy    ddot
*     Also
*        dnrm2    (from NAG,I think).
*
*     Also a few utilities to avoid some of the
*     loops in MINRES (so the debugger can step past them quickly):
*        daxpy2   dload2   dscal2   
*
* 15 Jul 2003: dnrm2  is now the NAG version.
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue

      end ! subroutine daxpy

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue

      end ! subroutine dcopy

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp

      end ! function ddot

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision   function dnrm2 ( n, x, incx )

      implicit           double precision (a-h,o-z)
      integer            incx, n
      double precision   x(*)

*     ==================================================================
*     dnrm2 returns the Euclidean norm of a vector via the function
*     name, so that dnrm2 := sqrt( x'*x ).
*
*     15 Jul 2003: dnrm2  obtained from SNOPT src (probably from NAG).
*                  s1flmx replaced by safe large number.
*     ==================================================================
!!!   double precision   s1flmx
      parameter         (one = 1.0d+0, zero = 0.0d+0 )
      double precision   norm
      intrinsic          abs
*     ------------------------------------------------------------------
*     flmax = s1flmx( )
      flmax = 1.0d+50

      if (     n .lt. 1) then
         norm  = zero

      else if (n .eq. 1) then
         norm  = abs( x(1) )

      else
         scale = zero
         ssq   = one

         do 10, ix = 1, 1+(n-1)*incx, incx

            if (x(ix) .ne. zero) then
               absxi = abs( x(ix) )

               if (scale .lt. absxi) then
                  ssq   = one + ssq*(scale/absxi)**2
                  scale = absxi
               else
                  ssq   = ssq +     (absxi/scale)**2
               end if
            end if
   10    continue

         sqt = sqrt( ssq )
         if (scale .lt. flmax/sqt) then
            norm = scale*sqt
         else
            norm = flmax
         end if
      end if

      dnrm2  = norm

      end ! function dnrm2

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine daxpy2( n, a, x, y, z )

      implicit           none
      integer            n
      double precision   a, x(n), y(n), z(n)

*     ------------------------------------------------------------------
*     daxpy2 sets  z = a*x + y.
*     31 May 1999: First version written for MINRES.
*     ------------------------------------------------------------------

      integer            i

      do i = 1, n
         z(i) = a*x(i) + y(i)
      end do

      end ! subroutine daxpy2

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dload2( n, const, x )

      implicit           none
      integer            n
      double precision   const, x(n)

*     ------------------------------------------------------------------
*     dload2 loads all elements of x with const.
*     ------------------------------------------------------------------

      integer            i

      do i = 1, n
         x(i) = const
      end do

      end ! subroutine dload2

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dscal2( n, a, x, y )

      implicit           none
      integer            n
      double precision   a, x(n), y(n)

*     ------------------------------------------------------------------
*     dscal2 sets y = a*x.
*     ------------------------------------------------------------------

      integer            i

      do i = 1, n
         y(i) = a*x(i)
      end do

      end ! subroutine dscal2
