
 

      SUBROUTINE GMRES(A,X,B,NX, OMEGA,EPSI,AUX,NAUX,LIMITI,NREPT)
*        GMRES - iterative solution of systems of linear equations 
*     ------------------------------------------------------------------
*
*     A(.)      = matrix of equation A . X = B, to be solved for X 
*     X(.)      = vector of unknowns
*     B(.)      = r.h.s. vector
*     NX        = dimension parameter of vectors A(.) and (.)
*     OMEGA     = facto          0.1 <  OMEGA  < 2
*     EPSI      = epsilon with 1E-16 < epsilon < 1 
*     AUX(NAUX) = auxiliary array
*     NAUX      = dimension of auxiliary array AUX(.)  
*     LIMITI    = max number of iterations
*     NREPT     = actual number of iterations
*
*     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NX,NAUX,NREPT,LIMITI
      DOUBLE PRECISION A(*),X(NX),B(NX),AUX(NAUX),OMEGA,EPSI
*     internal
      INTEGER I,MAXLOP,N,INDR,INDW,INDA,INDV,LOOP,LIMIT
      PARAMETER (MAXLOP=37)    
      DOUBLE PRECISION EPS  
*     ...
      NREPT=0
      N=0
      DO I=1,NX
       IF(A((I*I+I)/2).NE.0.0D0) N=I
      END DO

      N=NX !!!

      DO I=N+1,NX
       X(I)=0.0D0 ! additional elements are zero
      END DO
*               N 
      INDR=1   +N               ! define indices for AUX(.)
      INDW=INDR+N
      INDA=INDW+N
      INDV=INDA+N 
      LOOP=(NAUX+1-INDV)/N      ! max number of loops
      IF(LOOP.LT.10) THEN
         WRITE(*,*) 'GMRES: length of AUX array is',NAUX
         WRITE(*,*) 'GMRES: number of loops',LOOP,' too small - stop'
         STOP
      END IF
      LOOP=MIN(LOOP,N,MAXLOP)  
      EPS=MAX(1.0D-16,MIN(EPSI,1.0D0))
      LIMIT=MAX(1,MIN(LIMITI,2000))

      WRITE(*,*) 'GMRES: LOOP =',LOOP
      WRITE(*,*) 'GMRES: EPS  =',EPS
      WRITE(*,*) 'GMRES: OMEGA=',OMEGA
      WRITE(*,*) 'GMRES: N    =',N   
 
      CALL CGMRES(A,X,B,N,
     +            AUX(1),AUX(INDR),AUX(INDW),AUX(INDA),AUX(INDV),
     +            LOOP,OMEGA,EPS,LIMIT,NREPT)

      WRITE(*,*) 'GMRES: NREPT=',NREPT      
  
      END 


      SUBROUTINE CGMRES(A,X,B,N, AX,R,WJ,AV,V, MLOOP,OMEGAI,EPS,
     +                  LIMIT,NREPT)
*     GMRES - iterative solution of systems of linear equations
      IMPLICIT NONE

c      INTEGER NAUX2                     !!!!
c      PARAMETER (NAUX2=1000)
c      DOUBLE PRECISION AUX2(NAUX2)


      INTEGER N,MLOOP,MXLOOP,NREPT,LIMIT
      PARAMETER (MXLOOP=37)  ! dimension parameter
      DOUBLE PRECISION A(*),X(N),B(N),
     +                 AX(N),R(N),WJ(N),AV(N),V(N,MXLOOP),
     +                 OMEGAI,EPS 
*     internal
      INTEGER LOOP,MALOOP 
*
      DOUBLE PRECISION ALPHA(MXLOOP),BETA,GAMMA(MXLOOP),SUM,CS(2,MXLOOP)
      PARAMETER (MALOOP=(MXLOOP*MXLOOP+MXLOOP)/2)
      DOUBLE PRECISION H(MALOOP)

      DOUBLE PRECISION HIJ,DVXY,SUM2,GAMMA1,EXPN,EXPS,OMEGA 
      INTEGER I,J,K,IJ,JJ,MJ
*     ...
      OMEGA=OMEGAI
      IF(OMEGA.NE.0.0D0) OMEGA=MAX(0.1D0,MIN(OMEGA,2.0D0))
      DO I=1,N
       X(I)=B(I)
      END DO

      IF(OMEGA.NE.0.0D0) THEN 
         WRITE(*,*) 'Input to first MXSOLV:'
         WRITE(*,111) (X(I),I=1,N)
         WRITE(*,*) ' '       

*         CALL KKT(A,X,X,N, AUX2,NAUX2)

         CALL MXSOLV(A,X,N,OMEGA)       ! approx. solution
         WRITE(*,*) 'Result of first MXSOLV::'
         WRITE(*,111) (X(I),I=1,N)
         WRITE(*,*) ' '
      END IF

      LOOP=0                         ! number of restarts
 01   LOOP=LOOP+1
*     calculate residual R, and V_1, Gamma_1
      CALL PDVX(A,X,AX,N)            ! AX := A * X
      DO K=1,N
       R(K)=B(K)-AX(K)               ! r  := b - A * x
      END DO
 
      IF(LOOP.EQ.1.AND.OMEGA.NE.0.0D0) THEN     
       WRITE(*,*) 'Input to first MXSOLV:'
       WRITE(*,111) (R(I),I=1,N)
 111   FORMAT(4G15.5)
       WRITE(*,*) ' '
       END IF
      IF(OMEGA.NE.0.0D0) CALL MXSOLV(A,R,N,OMEGA)       ! preconditioning 
*      IF(OMEGA.NE.0.0D0) CALL KKT(A,R,R,N, AUX2,NAUX2)
  

      IF(LOOP.EQ.1.AND.OMEGA.NE.0.0D0) THEN 
       WRITE(*,*) 'Result of first MXSOLV::'
       WRITE(*,111) (R(I),I=1,N)
       WRITE(*,*) ' '
      END IF

      GAMMA(1)=SQRT(DVXY(R,R,N))     ! gamma = ||r||_2 
      IF(GAMMA(1).EQ.0.0D0) RETURN   ! finished
      DO K=1,N
       V(K,1)=R(K)/GAMMA(1)          ! V_1
      END DO

*     ------------------------------------------------------------------

      DO J=1,MLOOP                   ! up to M vectors
       NREPT=NREPT+1                 ! count    
       CALL PDVX(A,V(1,J),AV,N)      ! AV_j

       IF(OMEGA.NE.0.0D0) THEN      
*        CALL KKT(A,AV,AV,N, AUX2,NAUX2)
        CALL MXSOLV(A,AV,N,OMEGA)     ! preconditioning
       END IF

       DO I=1,J
        IJ=(J*J-J)/2+I
        H(IJ)=DVXY(V(1,I),AV,N)      ! element H_ij = V_i^T  A  V_j
       END DO

       SUM2=0.0D0
       DO K=1,N
        SUM=AV(K)                   ! vector W_j 
        DO I=1,J
         IJ=(J*J-J)/2+I
         SUM=SUM-H(IJ)*V(K,I)
        END DO
        WJ(K)=SUM
        SUM2=SUM2+SUM*SUM 
       END DO 
       H(J+(J*J+J)/2)=SQRT(SUM2)

       DO I=1,J-1                               ! rotation
        IJ=I+(J*J-J)/2
        HIJ=CS(1,I)*H(IJ)+CS(2,I)*H(IJ+1)
        H(IJ+1)=-CS(2,I)*H(IJ)+CS(1,I)*H(IJ+1) 
        H(IJ)=HIJ
       END DO

       JJ=(J*J+J)/2
       BETA=SQRT(H(JJ)**2+H(JJ+J)**2)
       CS(1,J)=H(JJ  )/BETA            ! cosine
       CS(2,J)=H(JJ+J)/BETA            ! sine

       H((J*J+J)/2)=BETA

       GAMMA(J+1)=-CS(2,J)*GAMMA(J)  ! gamma_j+1   
       GAMMA(J  )=+CS(1,J)*GAMMA(J)     ! gamma_j         

       MJ=J ! ... at index J
       GAMMA1=GAMMA(J+1)

       EXPN=LOG10(ABS(GAMMA1)+1.0D-16) ! fuer die Konvergenzerkennung

c       WRITE(*,103) LOOP,J,GAMMA1 
       IF(ABS(GAMMA1).LT.EPS) GOTO 10

       DO I=1,N
        V(I,J+1)=WJ(I)/H((J*J+J)/2+J)      ! V_j+1         
       END DO 
      END DO     ! of J=1,M loop
      MJ=MLOOP   ! end of loop reached  

*     ------------------------------------------------------------------
 10   IF(LOOP.EQ.1) THEN
         EXPS=EXPN
      ELSE
         EXPS=0.5D0*(EXPN+EXPS)
      END IF    
      
      WRITE(*,*) 'Loop:',(LOOP-1)*MLOOP+MJ,EXPN,EXPS 
      DO I=MJ,1,-1
       SUM=0.0D0
       DO K=I+1,MJ
        SUM=SUM+ALPHA(K)*H(I+(K*K-K)/2)
       END DO
       ALPHA(I)=(GAMMA(I)-SUM)/H((I*I+I)/2)  ! ALPHA_i
      END DO     
      
      DO K=1,N 
       SUM=X(K)
       DO I=1,MJ 
        SUM=SUM+ALPHA(I)*V(K,I)              ! X update
       END DO
       X(K)=SUM
      END DO
      IF(ABS(GAMMA1).GE.EPS.AND.LOOP*MLOOP.LT.LIMIT) GOTO 01
 102  FORMAT(I5,I4,2G13.4)      
 103  FORMAT(2I9,13X,G13.4)
      END


      SUBROUTINE MXSOLV(A,X,N,OMEGA) ! solve approximately A * X_out = X_in
      DOUBLE PRECISION A(*),X(N),OMEGA
*     ...
      DO I=N,1,-1                   ! backward
       II=(I*I+I)/2
       SUM=X(I)
       IK=II+I
       DO K=I+1,N
        SUM=SUM-OMEGA*A(IK)*X(K)
        IK=IK+K
       END DO
       X(I)=SUM/A(II)  /1.0001D0
      END DO

      DO I=1,N                      ! scaling
       X(I)=X(I)*A((I*I+I)/2)  
      END DO

      IK=0                          ! forward
      DO I=1,N
       SUM=X(I)
       DO K=1,I-1
        SUM=SUM-OMEGA*A(IK+K)*X(K)
       END DO
       X(I)=SUM/A(IK+I)   /1.0001D0
       IK=IK+I
      END DO
      INIT=1
      END 

      DOUBLE PRECISION FUNCTION DVXY(X,Y,N)
      IMPLICIT NONE  
*     Scalar product of two vectors
*                - - -                 T
*       S = DVXY(X,Y,N)           S = X  * Y (scalar product)
*
      INTEGER N,M,I,J
      DOUBLE PRECISION X(*),Y(*)
      DOUBLE PRECISION SUM

*     ...
      SUM=0.0D0
      M=N/4
      I=0   
      DO J=1,M
       SUM=SUM+X(I+1)*Y(I+1)+X(I+2)*Y(I+2)+X(I+3)*Y(I+3)+X(I+4)*Y(I+4)
       I=I+4      
      END DO      
      DO J=I+1,N  
       SUM=SUM+X(J)*Y(J)
      END DO      
      DVXY=SUM     
      END
  
      SUBROUTINE PDVX(V,X,VX,N)
*     VX = V X
      IMPLICIT NONE
      INTEGER N,I,J,IJ 
      DOUBLE PRECISION V(*),X(N),VX(N)
      DOUBLE PRECISION SUM

*     ...
      DO I=1,N
       SUM=0.0D0
       IJ=(I*I-I)/2
       DO J=1,N
        IF(J.LE.I) THEN
           IJ=IJ+1
        ELSE
           IJ=IJ+J-1
        END IF
        SUM=SUM+V(IJ)*X(J)
       END DO
       VX(I)=SUM
      END DO
      END


