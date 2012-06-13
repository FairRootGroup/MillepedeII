
! Code converted using TO_F90 by Alan Miller
! Date: 2012-03-16  Time: 11:09:33

!> \file
!! Random numbers.
!!
!! Random number generators for Uniform and Normal distribution:
!!
!!     URAN() for U(0,1)
!!     GRAN() for N(0,1)

!> F.Gutbrod random number generator.
!!
!! Return N random numbers U(0,1) in array A(N).
!! Initialization by entry GBRVIN.
!!
!! \param[in]   n  number of requested random number
!! \param[out]  a  array of requested random number

SUBROUTINE gbrshi(n,a)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: ian
    INTEGER :: iboost
    INTEGER :: ic
    INTEGER :: idum
    INTEGER :: irotor
    INTEGER :: iseed
    INTEGER :: it
    INTEGER :: iwarm
    INTEGER :: j
    INTEGER :: jseed
    INTEGER :: jwarm
    INTEGER :: k
    INTEGER :: m
    INTEGER :: mbuff

    INTEGER, INTENT(IN)                      :: n
    REAL, INTENT(OUT)                        :: a(*)
    INTEGER, PARAMETER :: nb=511
    INTEGER, PARAMETER :: ia=16807
    INTEGER, PARAMETER :: im=2147483647
    INTEGER, PARAMETER :: iq=127773
    INTEGER, PARAMETER :: ir=2836
    REAL, PARAMETER :: aeps=1.0E-10
    REAL, PARAMETER :: scalin=4.6566125E-10
    COMMON/ranbuf/mbuff(0:nb),ian,ic,iboost

    INTEGER :: istart

    irotor(m,n)=IEOR(ior(ishft(m,17),ishft(m,-15)),n)
    DATA istart/0/,iwarm/10/,iseed/4711/
    IF(istart /= 0) GO TO 20
    WRITE(*,*) ' Automatic GBRSHI initialization using:'
    !     initialize buffer
10  idum=iseed+9876543          ! prevent damage, if iseed=0
    WRITE(*,*) '           ISEED=',iseed,'   IWARM=',iwarm
    DO j=0,nb+1                 ! fill buffer
        k=idum/iq                  ! minimal standard generator
        idum=ia*(idum-k*iq)-ir*k   !    with Schrages method
        IF(idum < 0) idum=idum+im !
        mbuff(j)=ishft(idum,1)     ! fill in leading bit
    END DO
    ian=IAND(ian,nb)            ! mask angle
    ic=1                        ! set pointer
    iboost=0
    DO j=1,iwarm*nb             ! warm up a few times
        it=mbuff(ian)              ! hit ball angle
        mbuff(ian)=irotor(it,ic)   ! new spin
        ic=it                      ! replace red spin
        ian=IAND(it+iboost,nb)     ! boost and mask angle
        iboost=iboost+1            ! increment boost
    END DO
    IF(istart < 0) RETURN      ! return for RBNVIN
    istart=1                    ! set done-flag
    !     generate array of r.n.
    20   DO i=1,n
        it=mbuff(ian)              ! hit ball angle
        mbuff(ian)=irotor(it,ic)   ! new spin
        ic=it                      ! replace red spin
        ian=IAND(it+iboost,nb)     ! boost and mask angle
        a(i)=FLOAT(ishft(it,-1))*scalin+aeps ! avoid zero output
        iboost=iboost+1            ! increment boost
    END DO
    iboost=IAND(iboost,nb)
    RETURN

    ENTRY gbrvin(jseed,jwarm)   ! initialize, but only once
    IF(istart == 0) THEN
        WRITE(*,*) ' Gbrshi initialization by GBRVIN-call using:'
        iseed=jseed              ! copy seed and
        iwarm=jwarm              ! warm-up parameter
        istart=-1                ! start flag
        GO TO 10
    END IF
END SUBROUTINE gbrshi

!> GBRSHI initialization using TIME().
SUBROUTINE gbrtim
    IMPLICIT NONE
    INTEGER :: jseed
    REAL :: time

    LOGICAL :: done
    DATA    done/.FALSE./
    IF(done) RETURN
    jseed=time()
    WRITE(*,*) ' Gbrshi initialialization using Time()'
    CALL gbrvin(jseed,10)
    done=.TRUE.
END SUBROUTINE gbrtim

!> Random number U(0,1) using RANSHI.
!!
!! \return   random number U(0,1)

REAL FUNCTION uran()     ! U(0,1)
    IMPLICIT NONE
    INTEGER :: indx
    INTEGER :: ndim

    PARAMETER (ndim=100)
    REAL :: buffer(ndim)
    DATA indx/ndim/
    SAVE indx,buffer
    indx=MOD(indx,ndim)+1
    IF(indx == 1) CALL gbrshi(ndim,buffer)
    uran=buffer(indx)
END FUNCTION uran

!> Gauss random number.
!!
!! \return   random number N(0,1)

REAL FUNCTION gran()     ! N(0,1)
    IMPLICIT NONE
    REAL :: al
    REAL :: cs
    INTEGER :: indx
    INTEGER :: kn
    INTEGER :: ndim
    REAL :: radsq
    REAL :: rn1
    REAL :: rn2
    REAL :: sn

    PARAMETER (ndim=100)
    REAL:: buffer(ndim)
    DATA indx/ndim/,kn/1/
    SAVE indx,buffer,kn,cs,al
    !     ...
    IF(kn <= 1) THEN
        !        two U(-1,+1) random numbers
10      indx=MOD(indx,ndim)+2
        IF(indx == 2) CALL gbrshi(ndim,buffer)
        rn1=buffer(indx-1)-1.0+buffer(indx-1)
        rn2=buffer(indx  )-1.0+buffer(indx)
        radsq=rn1*rn1+rn2*rn2
        IF(radsq > 1.0) GO TO 10 ! test point inside circle?
        !        sine and cosine for random phi
        sn=rn1/SQRT(radsq)
        cs=rn2/SQRT(radsq)
        !        transform to gaussians
        al=SQRT(-2.0*ALOG(radsq))
        kn =2
        gran=sn*al
    ELSE
        kn =1
        gran=cs*al
    END IF
END FUNCTION gran
