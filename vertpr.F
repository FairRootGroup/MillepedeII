


      SUBROUTINE PZVERT(N,X)                     ! print vertical
*     ZVERT prints the array X of dimension N  (MAX  120) in 6 lines
      PARAMETER (NN=6)
      REAL X(N)
      CHARACTER*66 PX(10),CH(10)*1
      SAVE
      DATA CH/'0','1','2','3','4','5','6','7','8','9'/
*     ...
      IF(N.LE.0) RETURN
      JM=0
      DO I=1,10
       PX(I)=' '
      END DO

      M=MIN(60,N)
      JL=0
      XM=0.0
      DO J=1,M
       IF(ABS(X(J)).GT.XM) THEN
          XM=ABS(X(J))
          MX=J                             ! index of max
       END IF
       IF(X(J).LT.0.0) PX(1)(6+J:6+J)='-'  ! negative columns
       IF(X(J).NE.0.0) JL=J                ! last non-zero column
      END DO
      IF(XM.EQ.0.0.OR.JL.LE.0) RETURN      ! empty array
      JL=60

      KN=MIN(6,MAX(2,IABS(NN)))
      KE=ALOG10(XM*1.0001)
      IF(XM.LT.1.0) KE=KE-1
   22 FAC=10.0**(KN-1-KE)
      IJ=FAC*XM+0.5
      IF(IJ.GE.10**KN) THEN
         KE=KE+1
         GOTO 22
      END IF
      IA=2+KN

      DO J=1,JL
       IJ=FAC*ABS(X(J))+0.5   ! convert to integer
       IM=0
       IF(IJ.NE.0) THEN
          DO I=1,KN
           IF(IJ.NE.0) THEN
              IN=MOD(IJ,10)   ! last digit
              IJ=IJ/10        ! reduce
              IF(IN.NE.0.AND.IM.EQ.0) IM=IA-I+1
              PX(IA-I)(6+J:6+J)=CH(IN+1)
           END IF
          END DO
       END IF
       JM=MAX(IM,JM)
       END DO

      KL=KE
   50 IF(KE.GE.KN) THEN
         KE=KE-3
         GOTO 50
      END IF
   55 IF(KE.LT.0) THEN
         KE=KE+3
         GOTO 55
      END IF

      IN=KE+2                          ! exponent
      IZ=KL-KE
      PX(IN)(6:6)='.'
      PX(IN)(1:1)='E'
      IF(IZ.LT.0) THEN
         PX(IN)(2:2)='-'
         IZ=-IZ
      END IF
      I1=IZ/10                         ! insert exponent
      I2=MOD(IZ,10)
      PX(IN)(3:3)=CH(I1+1)
      PX(IN)(4:4)=CH(I2+1)
      JM=MIN(2+KN,JM)
      JM=MAX(IN+1,JM)
      DO J=1,JL    ! '0' for small nonzero values
      IF(X(J).NE.0.0.AND.PX(JM-1)(6+J:6+J).EQ.' ') PX(JM-1)(6+J:6+J)='0'
      END DO
      DO I=JM,8
       PX(I)=' '
      END DO

      DO J=1,((JL+9)/10)*10            ! index line below
       IF(PX(JM-1)(6+J:6+J).EQ.' ') PX(JM-1)(6+J:6+J)='_'
       IF(MOD(J,2).NE.1) THEN
          I=MOD(J,10)+1
          PX(JM+1)(6+J:6+J)=CH(I)      ! last digit of even bin numbers
          IF(I.EQ.1) THEN              ! ten'th column
             I=MOD(J/10,10)+1
             PX(JM)(6+J:6+J)=CH(I)

          END IF
       END IF
      END DO

      DO J=1,JL
       IF(X(J).EQ.X(MX)) THEN
          PX(JM)(6+J:6+J)='*'          ! * in max bin
       END IF
      END DO

      JM=JM+1
      IF(NN.LT.0) JM=JM-2              ! no index line
      LC=((JL+9)/10)*10+6
      DO J=1,JM
       WRITE(*,*)  PX(J)(1:LC)       ! print
C       WRITE(*,101)  PX(J)(1:LC)       ! print
      END DO
  100 RETURN
  101 FORMAT(A)
      END


      SUBROUTINE PIVERT(N,LIST)                  ! vert. print integer data
      INTEGER LIST(N)
      REAL X(*),Y(60),SC(7)
      DOUBLE PRECISION SUM
      SAVE 
*     ...
      LL=(N+59)/60 ! compression factor
      M=0
      I=0
 10   NHIST=0
      DO L=1,LL
       IF(I+L.LE.N) NHIST=NHIST+LIST(I+L)
      END DO      
      I=I+LL 
      M=M+1
      Y(M)=NHIST
      IF(I.LT.N) GOTO 10
      CALL PZVERT(M,Y)
      RETURN

      ENTRY PFVERT(N,X)                       ! vert. print fltpt data
            LL=(N+59)/60 ! compression factor
      M=0
      I=0
 20   SUM=0.0
      DO L=1,LL
       IF(I+L.LE.N) SUM=SUM+X(I+L)
      END DO
      I=I+LL
      M=M+1
      Y(M)=SUM
      IF(I.LT.N) GOTO 20
      CALL PZVERT(M,Y)
      RETURN 

      ENTRY PSVERT(XA,XB)                     ! print scale 
*     print scale from XA ... XB
      XC=XB
      DO I=1,7
       SC(I)=(FLOAT(7-I)*XA+FLOAT(I-1)*XC)/6.0
      END DO
      WRITE(*,101) SC
 101  FORMAT(3X,7G10.3)
      END 
 
