
! Code converted using TO_F90 by Alan Miller
! Date: 2012-03-16  Time: 11:06:43

!> \file
!! MINRES algorithm.

!> Solution of linear equation system.
!! \verbatim
!!     ------------------------------------------------------------------
!!
!!     MINRES  is designed to solve the system of linear equations
!!
!!                Ax = b
!!
!!     or the least-squares problem
!!
!!         min || Ax - b ||_2,
!!
!!     where A is an n by n symmetric matrix and b is a given vector.
!!     The matrix A may be indefinite and/or singular.
!!
!!     1. If A is known to be positive definite, the Conjugate Gradient
!!        Method might be preferred, since it requires the same number
!!        of iterations as MINRES but less work per iteration.
!!
!!     2. If A is indefinite but Ax = b is known to have a solution
!!        (e.g. if A is nonsingular), SYMMLQ might be preferred,
!!        since it requires the same number of iterations as MINRES
!!        but slightly less work per iteration.
!!
!!     The matrix A is intended to be large and sparse.  It is accessed
!!     by means of a subroutine call of the form
!!     SYMMLQ development:
!!
!!                call Aprod ( n, x, y )
!!
!!     which must return the product y = Ax for any given vector x.
!!
!!
!!     More generally, MINRES is designed to solve the system
!!
!!                (A - shift*I) x = b
!!     or
!!         min || (A - shift*I) x - b ||_2,
!!
!!     where  shift  is a specified scalar value.  Again, the matrix
!!     (A - shift*I) may be indefinite and/or singular.
!!     The work per iteration is very slightly less if  shift = 0.
!!
!!     Note: If  shift  is an approximate eigenvalue of  A
!!     and  b  is an approximate eigenvector,  x  might prove to be
!!     a better approximate eigenvector, as in the methods of
!!     inverse iteration and/or Rayleigh-quotient iteration.
!!     However, we're not yet sure on that -- it may be better
!!     to use SYMMLQ.
!!
!!     A further option is that of preconditioning, which may reduce
!!     the number of iterations required.  If M = C C' is a positive
!!     definite matrix that is known to approximate  (A - shift*I)
!!     in some sense, and if systems of the form  My = x  can be
!!     solved efficiently, the parameters precon and Msolve may be
!!     used (see below).  When  precon = .true., MINRES will
!!     implicitly solve the system of equations
!!
!!             P (A - shift*I) P' xbar  =  P b,
!!
!!     i.e.                  Abar xbar  =  bbar
!!     where                         P  =  C**(-1),
!!                                Abar  =  P (A - shift*I) P',
!!                                bbar  =  P b,
!!
!!     and return the solution       x  =  P' xbar.
!!     The associated residual is rbar  =  bbar - Abar xbar
!!                                      =  P (b - (A - shift*I)x)
!!                                      =  P r.
!!
!!     In the discussion below, eps refers to the machine precision.
!!     eps is computed by MINRES.  A typical value is eps = 2.22d-16
!!     for IEEE double-precision arithmetic.
!!
!!     Parameters
!!     ----------
!!
!!     n       input      The dimension of the matrix A.
!!
!!     b(n)    input      The rhs vector b.
!!
!!     r1(n)   workspace
!!     r2(n)   workspace
!!     v(n)    workspace
!!     w(n)    workspace
!!     w1(n)   workspace
!!     w2(n)   workspace
!!
!!     x(n)    output     Returns the computed solution  x.
!!
!!     y(n)    workspace
!!
!!     Aprod   external   A subroutine defining the matrix A.
!!                        For a given vector x, the statement
!!
!!                              call Aprod ( n, x, y )
!!
!!                        must return the product y = Ax
!!                        without altering the vector x.
!!
!!     Msolve  external   An optional subroutine defining a
!!                        preconditioning matrix M, which should
!!                        approximate (A - shift*I) in some sense.
!!                        M must be positive definite.
!!                        For a given vector x, the statement
!!
!!                              call Msolve( n, x, y )
!!
!!                        must solve the linear system My = x
!!                        without altering the vector x.
!!
!!                        In general, M should be chosen so that Abar has
!!                        clustered eigenvalues.  For example,
!!                        if A is positive definite, Abar would ideally
!!                        be close to a multiple of I.
!!                        If A or A - shift*I is indefinite, Abar might
!!                        be close to a multiple of diag( I  -I ).
!!
!!                        NOTE.  The program calling MINRES must declare
!!                        Aprod and Msolve to be external.
!!
!!     checkA  input      If checkA = .true., an extra call of Aprod will
!!                        be used to check if A is symmetric.  Also,
!!                        if precon = .true., an extra call of Msolve
!!                        will be used to check if M is symmetric.
!!
!!     precon  input      If precon = .true., preconditioning will
!!                        be invoked.  Otherwise, subroutine Msolve
!!                        will not be referenced; in this case the
!!                        actual parameter corresponding to Msolve may
!!                        be the same as that corresponding to Aprod.
!!
!!     shift   input      Should be zero if the system Ax = b is to be
!!                        solved.  Otherwise, it could be an
!!                        approximation to an eigenvalue of A, such as
!!                        the Rayleigh quotient b'Ab / (b'b)
!!                        corresponding to the vector b.
!!                        If b is sufficiently like an eigenvector
!!                        corresponding to an eigenvalue near shift,
!!                        then the computed x may have very large
!!                        components.  When normalized, x may be
!!                        closer to an eigenvector than b.
!!
!!     nout    input      A file number.
!!                        If nout .gt. 0, a summary of the iterations
!!                        will be printed on unit nout.
!!
!!     itnlim  input      An upper limit on the number of iterations.
!!
!!     rtol    input      A user-specified tolerance.  MINRES terminates
!!                        if it appears that norm(rbar) is smaller than
!!                              rtol * norm(Abar) * norm(xbar),
!!                        where rbar is the transformed residual vector,
!!                              rbar = bbar - Abar xbar.
!!
!!                        If shift = 0 and precon = .false., MINRES
!!                        terminates if norm(b - A*x) is smaller than
!!                              rtol * norm(A) * norm(x).
!!
!!     istop   output     An integer giving the reason for termination...
!!
!!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!!                        second Lanczos vector is zero.  This means the
!!                        rhs is very special.
!!                        If there is no preconditioner, b is an
!!                        eigenvector of A.
!!                        Otherwise (if precon is true), let My = b.
!!                        If shift is zero, y is a solution of the
!!                        generalized eigenvalue problem Ay = lambda My,
!!                        with lambda = alpha1 from the Lanczos vectors.
!!
!!                        In general, (A - shift*I)x = b
!!                        has the solution         x = (1/alpha1) y
!!                        where My = b.
!!
!!               0        b = 0, so the exact solution is x = 0.
!!                        No iterations were performed.
!!
!!               1        Norm(rbar) appears to be less than
!!                        the value  rtol * norm(Abar) * norm(xbar).
!!                        The solution in  x  should be acceptable.
!!
!!               2        Norm(rbar) appears to be less than
!!                        the value  eps * norm(Abar) * norm(xbar).
!!                        This means that the residual is as small as
!!                        seems reasonable on this machine.
!!
!!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!!                        which should indicate that x has essentially
!!                        converged to an eigenvector of A
!!                        corresponding to the eigenvalue shift.
!!
!!               4        Acond (see below) has exceeded 0.1/eps, so
!!                        the matrix Abar must be very ill-conditioned.
!!                        x may not contain an acceptable solution.
!!
!!               5        The iteration limit was reached before any of
!!                        the previous criteria were satisfied.
!!
!!               6        The matrix defined by Aprod does not appear
!!                        to be symmetric.
!!                        For certain vectors y = Av and r = Ay, the
!!                        products y'y and r'v differ significantly.
!!
!!               7        The matrix defined by Msolve does not appear
!!                        to be symmetric.
!!                        For vectors satisfying My = v and Mr = y, the
!!                        products y'y and r'v differ significantly.
!!
!!               8        An inner product of the form  x' M**(-1) x
!!                        was not positive, so the preconditioning matrix
!!                        M does not appear to be positive definite.
!!
!!                        If istop .ge. 5, the final x may not be an
!!                        acceptable solution.
!!
!!     itn     output     The number of iterations performed.
!!
!!     Anorm   output     An estimate of the norm of the matrix operator
!!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!!
!!     Acond   output     An estimate of the condition of Abar above.
!!                        This will usually be a substantial
!!                        under-estimate of the true condition.
!!
!!     rnorm   output     An estimate of the norm of the final
!!                        transformed residual vector,
!!                           P (b  -  (A - shift*I) x).
!!
!!     ynorm   output     An estimate of the norm of xbar.
!!                        This is sqrt( x'Mx ).  If precon is false,
!!                        ynorm is an estimate of norm(x).
!!     ------------------------------------------------------------------
!!
!!
!!     MINRES is an implementation of the algorithm described in
!!     the following reference:
!!
!!     C. C. Paige and M. A. Saunders (1975),
!!     Solution of sparse indefinite systems of linear equations,
!!     SIAM J. Numer. Anal. 12(4), pp. 617-629.
!!     ------------------------------------------------------------------
!!
!!
!!     MINRES development:
!!            1972: First version, similar to original SYMMLQ.
!!                  Later lost @#%*!
!!        Oct 1995: Tried to reconstruct MINRES from
!!                  1995 version of SYMMLQ.
!!     30 May 1999: Need to make it more like LSQR.
!!                  In middle of major overhaul.
!!     19 Jul 2003: Next attempt to reconstruct MINRES.
!!                  Seems to need two vectors more than SYMMLQ.  (w1, w2)
!!                  Lanczos is now at the top of the loop,
!!                  so the operator Aprod is called in just one place
!!                  (not counting the initial check for symmetry).
!!     22 Jul 2003: Success at last.  Preconditioning also works.
!!                  minres.f added to http://www.stanford.edu/group/SOL/.
!!
!!     FUTURE WORK: A stopping rule is needed for singular systems.
!!                  We need to estimate ||Ar|| as in LSQR.  This will be
!!                  joint work with Sou Cheng Choi, SCCM, Stanford.
!!                  Note that ||Ar|| small => r is a null vector for A.
!!
!!
!!     Michael A. Saunders           na.msaunders@na-net.ornl.gov
!!     Department of MS&E            saunders@stanford.edu
!!     Stanford University
!!     Stanford, CA 94305-4026       (650) 723-1875
!!     ------------------------------------------------------------------
!!
!!
!!     Subroutines and functions
!!
!!     USER       Aprod , Msolve
!!     BLAS1      daxpy , dcopy , ddot  , dnrm2  } These are all in
!!     Utilities  daxpy2, dload2, dscal2         } the file minresblas.f
!!
!!     Functions
!! \endverbatim

SUBROUTINE minres( n, b, r1, r2, v, w, w1, w2, x, y,  &
    aprod, msolve, checka, precon, shift, nout , itnlim, rtol,  &
    istop, itn, anorm, acond, rnorm, ynorm )

    IMPLICIT NONE
    EXTERNAL aprod, msolve
    INTEGER :: n, nout, itnlim, istop, itn
    LOGICAL :: checka, precon
    DOUBLE PRECISION :: shift, rtol, anorm, acond, rnorm, ynorm,  &
        b(n), r1(n), r2(n), v(n), w(n), w1(n), w2(n), x(n), y(n)

    EXTERNAL ddot  , dnrm2
    DOUBLE PRECISION :: ddot  , dnrm2

    !     Local variables

    DOUBLE PRECISION :: alfa  , beta  , beta1 , cs    ,  &
        dbar  , delta , denom , diag  , eps   , epsa  , epsln , epsr  , epsx  ,  &
        agamma, gbar  , gmax  , gmin  , oldb  , oldeps, qrnorm, phi   , phibar,  &
        rhs1  , rhs2  , s     , sn    , t     , tnorm2, ynorm2, z
    INTEGER :: i
    LOGICAL :: debug, prnt

    DOUBLE PRECISION :: zero, one, two, ten
    PARAMETER        ( zero = 0.0D+0,  one =  1.0D+0,  &
        two  = 2.0D+0,  ten = 10.0D+0 )

    CHARACTER (LEN=16) :: enter, EXIT
    CHARACTER (LEN=52) :: msg(-1:8)

    DATA               enter /' Enter MINRES.  '/, EXIT  /' Exit  MINRES.  '/

    DATA               msg  &
        / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',  &
        'beta1 = 0.  The exact solution is  x = 0',  &
        'Requested accuracy achieved, as determined by rtol',  &
        'Reasonable accuracy achieved, given eps',  &
        'x has converged to an eigenvector', 'Acond has exceeded 0.1/eps',  &
        'The iteration limit was reached',  &
        'Aprod  does not define a symmetric matrix',  &
        'Msolve does not define a symmetric matrix',  &
        'Msolve does not define a pos-def preconditioner' /
    !     ------------------------------------------------------------------

    debug = .false.

    !     ------------------------------------------------------------------
    !     Compute eps, the machine precision.  The call to daxpy is
    !     intended to fool compilers that use extra-length registers.
    !     31 May 1999: Hardwire eps so the debugger can step thru easily.
    !     ------------------------------------------------------------------
    eps    = 2.22D-16    ! Set eps = zero here if you want it computed.

    eps  = 0.0D0                             !!!!!!!!!!!!!!!
    gmin = 0.0D0
    gmax = 0.0D0

    IF (eps <= zero) THEN
        eps    = two**(-12)
        DO
            eps    = eps / two
            x(1)   = eps
            y(1)   = one
            CALL daxpy ( 1, one, x, 1, y, 1 )
            IF (y(1) <= one) EXIT
        END DO
        eps    = eps * two
    END IF
    !     WRITE(*,*) 'computed epsilon is ',eps   !!!!!!!!!!!!!!!

    !     ------------------------------------------------------------------
    !     Print heading and initialize.
    !     ------------------------------------------------------------------
    IF (nout > 0) THEN
        WRITE(nout, 1000) enter, n, checka, precon, itnlim, rtol, shift
    END IF
    istop  = 0
    itn    = 0
    anorm  = zero
    acond  = zero
    rnorm  = zero
    ynorm  = zero
    CALL dload2( n, zero, x )

    !     ------------------------------------------------------------------
    !     Set up y and v for the first Lanczos vector v1.
    !     y  =  beta1 P' v1,  where  P = C**(-1).
    !     v is really P' v1.
    !     ------------------------------------------------------------------
    CALL dcopy ( n, b, 1, y , 1 )         ! y  = b
    CALL dcopy ( n, b, 1, r1, 1 )         ! r1 = b
    IF ( precon ) CALL msolve( n, b, y )
    beta1  = ddot  ( n, b, 1, y, 1 )

    IF (beta1 < zero) THEN    ! m must be indefinite.
        istop = 8
        GO TO 900
    END IF

    IF (beta1 == zero) THEN    ! b = 0 exactly.  Stop with x = 0.
        istop = 0
        GO TO 900
    END IF

    beta1  = SQRT( beta1 )       ! Normalize y to get v1 later.

    !     ------------------------------------------------------------------
    !     See if Msolve is symmetric.
    !     ------------------------------------------------------------------
    IF (checka  .AND.  precon) THEN
        CALL msolve( n, y, r2 )
        s      = ddot  ( n, y, 1, y, 1 )
        t      = ddot  ( n,r1, 1,r2, 1 )
        z      = ABS( s - t )
        epsa   = (s + eps) * eps**0.33333D+0
        IF (z > epsa) THEN
            istop = 7
            GO TO 900
        END IF
    END IF

    !     ------------------------------------------------------------------
    !     See if Aprod  is symmetric.
    !     ------------------------------------------------------------------
    IF (checka) THEN
        CALL aprod ( n, y, w )
        CALL aprod ( n, w, r2 )
        s      = ddot  ( n, w, 1, w, 1 )
        t      = ddot  ( n, y, 1,r2, 1 )
        z      = ABS( s - t )
        epsa   = (s + eps) * eps**0.33333D+0
        IF (z > epsa) THEN
            istop = 6
            GO TO 900
        END IF
    END IF

    !     ------------------------------------------------------------------
    !     Initialize other quantities.
    !     ------------------------------------------------------------------
    oldb   = zero
    beta   = beta1
    dbar   = zero
    epsln  = zero
    qrnorm = beta1
    phibar = beta1
    rhs1   = beta1
    rhs2   = zero
    tnorm2 = zero
    ynorm2 = zero
    cs     = - one
    sn     = zero
    CALL dload2( n, zero, w  )        ! w  = 0
    CALL dload2( n, zero, w2 )        ! w2 = 0
    CALL dcopy ( n, r1, 1, r2, 1 )    ! r2 = r1

    IF (debug) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'b    ', b
        WRITE(*,*) 'beta ', beta
        WRITE(*,*) ' '
    END IF

    !     ------------------------------------------------------------------
    !     Main iteration loop.
    !     ------------------------------------------------------------------
    DO
        itn    = itn   +  1               ! k = itn = 1 first time through
        IF (istop /= 0) EXIT

        !-----------------------------------------------------------------
        ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
        ! The general iteration is similar to the case k = 1 with v0 = 0:
        !
        !   p1      = Operator * v1  -  beta1 * v0,
        !   alpha1  = v1'p1,
        !   q2      = p2  -  alpha1 * v1,
        !   beta2^2 = q2'q2,
        !   v2      = (1/beta2) q2.
        !
        ! Again, y = betak P vk,  where  P = C**(-1).
        ! .... more description needed.
        !-----------------------------------------------------------------
        s      = one / beta            ! Normalize previous vector (in y).
        CALL dscal2( n, s, y, v )      ! v = vk if P = I

        CALL aprod ( n, v, y )
        CALL daxpy ( n, (- shift), v, 1, y, 1 )
        IF (itn >= 2) THEN
            CALL daxpy ( n, (- beta/oldb), r1, 1, y, 1 )
        END IF

        alfa   = ddot  ( n, v, 1, y, 1 )     ! alphak

        CALL daxpy ( n, (- alfa/beta), r2, 1, y, 1 )
        CALL dcopy ( n, r2, 1, r1, 1 )
        CALL dcopy ( n,  y, 1, r2, 1 )
        IF ( precon ) CALL msolve( n, r2, y )

        oldb   = beta                        ! oldb = betak
        beta   = ddot  ( n, r2, 1, y, 1 )    ! beta = betak+1^2
        IF (beta < zero) THEN
            istop = 6
            WRITE(*,*) 'unfortunately beta < 0 ',beta,zero
            EXIT
        END IF

        beta   = SQRT( beta )                ! beta = betak+1
        tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

        IF (itn == 1) THEN                 ! Initialize a few things.
            IF (beta/beta1 <= ten*eps) THEN ! beta2 = 0 or ~ 0.
                istop = -1                     ! Terminate later.
            END IF
            !tnorm2 = alfa**2
            gmax   = ABS( alfa )              ! alpha1
            gmin   = gmax                     ! alpha1
        END IF

        ! Apply previous rotation Qk-1 to get
        !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
        !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

        oldeps = epsln
        delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
        gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
        epsln  =               sn * beta ! epsln2 = 0         epslnk+1
        dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

        ! Compute the next plane rotation Qk

        agamma = SQRT( gbar**2 + beta**2 )   ! gammak
        cs     = gbar / agamma               ! ck
        sn     = beta / agamma               ! sk
        phi    = cs * phibar                 ! phik
        phibar = sn * phibar                 ! phibark+1

        IF (debug) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'v    ', v
            WRITE(*,*) 'alfa ', alfa
            WRITE(*,*) 'beta ', beta
            WRITE(*,*) 'gamma', agamma
            WRITE(*,*) 'delta', delta
            WRITE(*,*) 'gbar ', gbar
            WRITE(*,*) 'epsln', epsln
            WRITE(*,*) 'dbar ', dbar
            WRITE(*,*) 'phi  ', phi
            WRITE(*,*) 'phiba', phibar
            WRITE(*,*) ' '
        END IF

        ! Update  x.

        denom = one/agamma

        DO i = 1, n
            w1(i) = w2(i)
            w2(i) = w(i)
            w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
            x(i)  =   x(i) +   phi * w(i)
        END DO

        ! Go round again.

        gmax   = MAX( gmax, agamma)
        gmin   = MIN( gmin, agamma)
        z      = rhs1 / agamma
        ynorm2 = z**2  +  ynorm2
        rhs1   = rhs2  -  delta * z
        rhs2   =       -  epsln * z

        ! Estimate various norms and test for convergence.

        anorm  = SQRT( tnorm2 )
        ynorm  = SQRT( ynorm2 )
        epsa   = anorm * eps
        epsx   = anorm * ynorm * eps
        epsr   = anorm * ynorm * rtol
        diag   = gbar
        IF (diag == zero) diag = epsa

        qrnorm = phibar
        rnorm  = qrnorm

        ! Estimate  cond(A).
        ! In this version we look at the diagonals of  R  in the
        ! factorization of the lower Hessenberg matrix,  Q * H = R,
        ! where H is the tridiagonal matrix from Lanczos with one
        ! extra row, beta(k+1) e_k^T.

        acond  = gmax / gmin

        ! See if any of the stopping criteria are satisfied.
        ! In rare cases, istop is already -1 from above (Abar = const*I).

        IF (istop == 0) THEN
            IF (itn    >= itnlim    ) istop = 5
            IF (acond  >= 0.1D+0/eps) istop = 4
            IF (epsx   >= beta1     ) istop = 3
            IF (qrnorm <= epsx      ) istop = 2
            IF (qrnorm <= epsr      ) istop = 1
        END IF


        ! See if it is time to print something.

        IF (nout > 0) THEN
            prnt   = .false.
            IF (n      <= 40         ) prnt = .true.
            IF (itn    <= 10         ) prnt = .true.
            IF (itn    >= itnlim - 10) prnt = .true.
            IF (MOD(itn,10)  ==     0) prnt = .true.
            IF (qrnorm <=  ten * epsx) prnt = .true.
            IF (qrnorm <=  ten * epsr) prnt = .true.
            IF (acond  >= 1.0D-2/eps ) prnt = .true.
            IF (istop  /=  0         ) prnt = .true.
  
            IF ( prnt ) THEN
                IF (    itn     == 1) WRITE(nout, 1200)
                WRITE(nout, 1300) itn, x(1), qrnorm, anorm, acond
                IF (MOD(itn,10) == 0) WRITE(nout, 1500)
            END IF
        END IF

    END DO
    !     ------------------------------------------------------------------
    !     End of main iteration loop.
    !     ------------------------------------------------------------------

    ! Display final status.

900 IF (nout  > 0) THEN
        WRITE(nout, 2000) EXIT, istop, itn, EXIT, anorm, acond,  &
            EXIT, rnorm, ynorm
        WRITE(nout, 3000) EXIT, msg(istop)
    END IF

    RETURN


1000 FORMAT(// 1P,    a, 5X, 'Solution of symmetric   Ax = b'  &
        / ' n      =', i7, 5X, 'checkA =', l4, 12X, 'precon =', l4  &
        / ' itnlim =', i7, 5X, 'rtol   =', e11.2, 5X, 'shift  =', e23.14)
1200 FORMAT(// 5X, 'itn', 8X, 'x(1)', 10X,  &
        'norm(r)', 3X, 'norm(A)', 3X, 'cond(A)')
1300 FORMAT(1P, i8, e19.10, 3E10.2)
1500 FORMAT(1X)
2000 FORMAT(/ 1P, a, 5X, 'istop =', i3,   14X, 'itn   =', i8  &
        /     a, 5X, 'Anorm =', e12.4, 5X, 'Acond =', e12.4  &
        /     a, 5X, 'rnorm =', e12.4, 5X, 'ynorm =', e12.4)
3000 FORMAT(      a, 5X, a )

END SUBROUTINE minres
