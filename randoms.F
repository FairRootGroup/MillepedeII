
*     URAN() for U(0,1)
*     GRAN() for N(0,1)


      SUBROUTINE GBRSHI(N,A)
*     F.Gutbrod random number generator
*     return N random numbers U(0,1) in array A(N)
      PARAMETER (NB=511)
      PARAMETER (IA=16807,IM=2147483647,IQ=127773,IR=2836)
      PARAMETER (AEPS=1.0E-10,SCALIN=4.6566125E-10)
      COMMON/RANBUF/MBUFF(0:NB),IAN,IC,IBOOST
      REAL A(*)
      INTEGER ISTART
      IROTOR(M,N)=IEOR(IOR(ISHFT(M,17),ISHFT(M,-15)),N)
      DATA ISTART/0/,IWARM/10/,ISEED/4711/
      IF(ISTART.NE.0) GOTO 20
      WRITE(*,*) ' Automatic GBRSHI initialization using:'
*     initialize buffer
 10   IDUM=ISEED+9876543          ! prevent damage, if iseed=0
      WRITE(*,*) '           ISEED=',ISEED,'   IWARM=',IWARM
      DO J=0,NB+1                 ! fill buffer
       K=IDUM/IQ                  ! minimal standard generator
       IDUM=IA*(IDUM-K*IQ)-IR*K   !    with Schrages method
       IF(IDUM.LT.0) IDUM=IDUM+IM !
       MBUFF(J)=ISHFT(IDUM,1)     ! fill in leading bit
      END DO
      IAN=IAND(IAN,NB)            ! mask angle
      IC=1                        ! set pointer
      IBOOST=0
      DO J=1,IWARM*NB             ! warm up a few times
       IT=MBUFF(IAN)              ! hit ball angle
       MBUFF(IAN)=IROTOR(IT,IC)   ! new spin
       IC=IT                      ! replace red spin
       IAN=IAND(IT+IBOOST,NB)     ! boost and mask angle
       IBOOST=IBOOST+1            ! increment boost
      END DO
      IF(ISTART.LT.0) RETURN      ! return for RBNVIN
      ISTART=1                    ! set done-flag
*     generate array of r.n.
 20   DO I=1,N
       IT=MBUFF(IAN)              ! hit ball angle
       MBUFF(IAN)=IROTOR(IT,IC)   ! new spin
       IC=IT                      ! replace red spin
       IAN=IAND(IT+IBOOST,NB)     ! boost and mask angle
       A(I)=FLOAT(ISHFT(IT,-1))*SCALIN+AEPS ! avoid zero output
       IBOOST=IBOOST+1            ! increment boost
      END DO
      IBOOST=IAND(IBOOST,NB)
      RETURN

      ENTRY GBRVIN(JSEED,JWARM)   ! initialize, but only once
      IF(ISTART.EQ.0) THEN
         WRITE(*,*) ' Gbrshi initialization by GBRVIN-call using:'
         ISEED=JSEED              ! copy seed and
         IWARM=JWARM              !    warm-up parameter
         ISTART=-1                ! start flag
         GOTO 10
      END IF
      END

      SUBROUTINE GBRTIM
*     GBRSHI initialization using TIME()
      LOGICAL DONE
      DATA    DONE/.FALSE./
      IF(DONE) RETURN
      JSEED=TIME()
      WRITE(*,*) ' Gbrshi initialialization using Time()'
      CALL GBRVIN(JSEED,10)
      DONE=.TRUE.
      END


      FUNCTION URAN()     ! U(0,1)
*     random number U(0,1) using RANSHI
      PARAMETER (NDIM=100)
      REAL BUFFER(NDIM)
      DATA INDEX/NDIM/
      SAVE INDEX,BUFFER
      INDEX=MOD(INDEX,NDIM)+1
      IF(INDEX.EQ.1) CALL GBRSHI(NDIM,BUFFER)
      URAN=BUFFER(INDEX)
      END

      FUNCTION GRAN()     ! N(0,1)  
*     Gauss random number
      PARAMETER (NDIM=100)
      REAL BUFFER(NDIM)
      DATA INDEX/NDIM/,KN/1/
      SAVE INDEX,BUFFER,KN,CS,AL
*     ...
      IF(KN.LE.1) THEN
*        two U(-1,+1) random numbers
 10      INDEX=MOD(INDEX,NDIM)+2
         IF(INDEX.EQ.2) CALL GBRSHI(NDIM,BUFFER)
         RN1=BUFFER(INDEX-1)-1.0+BUFFER(INDEX-1)
         RN2=BUFFER(INDEX  )-1.0+BUFFER(INDEX)
         RADSQ=RN1*RN1+RN2*RN2
         IF(RADSQ.GT.1.0) GOTO 10 ! test point inside circle?
*        sine and cosine for random phi
         SN=RN1/SQRT(RADSQ)
         CS=RN2/SQRT(RADSQ)
*        transform to gaussians
         AL=SQRT(-2.0*ALOG(RADSQ))
         KN =2
         GRAN=SN*AL
      ELSE
         KN =1
         GRAN=CS*AL
      END IF
      END
