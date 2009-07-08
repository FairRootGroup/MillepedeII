
*
*     Line search
*        MTLINE, CSTEPM  
*    


      SUBROUTINE MTLINE(N,X,F,G,S,STP, WA, DGOUT,INFO)
*     Linesearch Algorithms wih Guaranteed Sufficient Decrease"
*     by JJ More and D Thuente, ACM Transactions on Mathenatical
*     Software 20 (1994), pp 286-307 - modified by V.B.
*
*        ===============================================================
*        The algorithm is designed to find a step which satisfies 
*        the sufficient decrease condition
*              F(X+STEP*S) <= F(X) + FTOL*STP*[gradF(X).S]
*        and the curvature condition
*              ABS[gradF(X+STP*S).S] <= GTOL*ABS[gradF(X).S]
*
*        Typical values for Constants: 
*                               FTOL = 1.0E-4
*                               GTOL = 0.9  for Newton, Quasi-Newton
*                                    = 0.1  for CG
*
*        Conditions: 0 < FTOL < GTOL <  1   for Newton, Quasi-Newton
*                    0 < FTOL < GTOL < 1/2  for CG
*        ===============================================================
*      
*     CALL MTLDEF(...)      ! initialization, see end of routine
*     NC=0                  ! reset
*     define STP =  (in units of S(.), usually 1.0) 
* 01  calculate     (note: first value before first MTLINE call) 
*             F    = function value,
*             G(.) = gradient,
*             S(.) = search vector, at X(.) 
*     NC=NC+1               ! count function calls  
*     CALL MTLINE(N,X,F,G,S,STP, WA, DGOUT,INFO)
*        WA(.) =work array, of dimension N
*        DGOUT =value of directional derivative when the Wolfe 
*               conditions hold; may be used to check descent     
*        INFO  =returned flag (not to be modified) 
*     IF(INFO.LT.0) GOTO 01 ! recalculate function ...  
*     IF(INFO.NE.1)         ! error condition
*     IF(INFO.EQ.1)         ! sufficient decrease and deriv. condition
*                           ! i.e. Wolfe conditions hold 
*     CALL MTINFO(INFO)     ! to print end condition
*     ------------------------------------------------------------------ 

      INTEGER N,MAXF,INFO,NFEV
      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
      SAVE

      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO,DGOUT,STPAC
c     DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,2.618D0,0.0D0/
*     ...

*     -------------------------  initialization  ----------------------
      IF(NFEV.GE.0) GOTO 10
      NFEV =0     ! reset  
      INFOC=1
      INFO =0 
      IF(N.LE.0.OR.STP.LE.ZERO) THEN 
         WRITE(*,*) 'MTLINE call with N=',N,', STP=',STP
         RETURN ! check N and step value
      END IF 
      DGINIT=ZERO
      DO J=1,N
       DGINIT=DGINIT+G(J)*S(J) ! dot product  gradient * search_direction 
      END DO
c      WRITE(*,*) 'Initial DG is',DGINIT,N
      IF(DGINIT.GE.ZERO) THEN
         WRITE(*,*) 'Initial DG is',DGINIT,N
         RETURN ! check initial search direction
      END IF 

      BRACKT= .FALSE.         ! local variables
      STAGE1= .TRUE.

      FINIT = F               ! initial function value
      DGTEST= FTOL*DGINIT     ! initial direct. derivative
      WIDTH = STPMAX - STPMIN ! width 
      WIDTH1= WIDTH/P5        
      DO J = 1, N
       WA(J) = X(J)           ! save initial parameter values
      END DO

      STX=ZERO    ! step is zero             at best step
      FX =FINIT   ! function value
      DGX=DGINIT  ! directional derivative

      STY=ZERO    ! step is zero             at other endpoint
      FY =FINIT   ! function value
      DGY=DGINIT  ! directional derivative

#ifdef graphics
*     graphic
      CALL QFINIT(0.0,0.0)
      CALL QUAT('red full line')
      FUNG=F
      DELD=FTOL*DGINIT
      CALL QD2XY(0.0,FUNG,1.0,FUNG+DELD,'l') ! limit line

      CALL QUAT('green full line') 
      DELD=GTOL*DGINIT
      CALL QD2XY(0.0,FUNG,1.0,FUNG+DELD,'l') ! slope line

      CALL QUAT('blue bars line') 
      FUNG=F 
      DELA=-0.1
      DELB=+0.5
      DEL1=DELA*DGINIT
      DEL2=DELB*DGINIT
      CALL QD2XY(DELA,FUNG+DEL1,DELB,FUNG+DEL2,'l')
      CALL QUAT('black bars line')
      CALL QUAT('black plus marker 2.0')  
#endif 

      GOTO 30

*     ------------------ new function and gradient ---------------------  

 10   IF(INFO.GT.0) GOTO 20
      INFO=0
      NFEV=NFEV+1
c      WRITE(*,*) NFEV,'.th function value',F,' with step STP',STP
      STPAC=STP  ! actual step
      DG=ZERO
      DO J=1,N
       DG=DG+G(J)*S(J) ! dot product  gradient * search_direction 
      END DO
c      WRITE(*,*) 'Actual DG is',DG,STP 

#ifdef graphics
*     graphic
      DELG=0.1
      DELD=DELG*DG
      XELG=STP
      FUNG=F
      CALL QD2XY(XELG-DELG,FUNG-DELD,
     +           XELG+DELG,FUNG+DELD,'lm')
      WRITE(*,*) 'Point', XELG-DELG,FUNG-DELD,XELG+DELG,FUNG+DELD
#endif

*     STP=        ! current step             at current step
*     F=          ! function value
*     DG=         ! directional derivative       
      FTEST1 = FINIT + STP*DGTEST

      IF((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     +      .OR. INFOC .EQ. 0)
     +                            INFO = 6 ! rounding error limitation 
      IF(STP .EQ. STPMAX .AND.
     +       F .LE. FTEST1 .AND. DG .LE. DGTEST)
     +                            INFO = 5 ! step at the upper bound
      IF(STP .EQ. STPMIN .AND.
     +       (F .GT. FTEST1 .OR. DG .GE. DGTEST))
     +                            INFO = 4 ! step at the lower bound
      IF(NFEV .GE. MAXF)          INFO = 3 ! max number of calls reached
      IF(BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)
     +                            INFO = 2 ! relative interval width limit

c      WRITE(*,*) 'Test for INFO=1',F,FTEST1,DG,GTOL*(-DGINIT),NFEV,MINF 

      IF(F.LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT).
     *       AND.NFEV.GE.MINF)    INFO = 1 ! min conditions hold ! 

      IF(INFO.NE.0) THEN ! termination 
         DGOUT=DG
         RETURN
      END IF

*     ------------------------ define STP ------------------------------

*     first stage" seek a step for which the modified function has
*     a non-positive value and non-negative derivative
 20   IF(STAGE1.AND.F.LE.FTEST1.AND.
     +   DG.GE.MIN(FTOL,GTOL)*DGINIT) STAGE1=.FALSE.

c      WRITE(*,*) 'TEST STAGE1 F,FTEST1,DG,FTOL,GTOL,DGINIT,M*DGINIT'
c      WRITE(*,*) STAGE1,F,FTEST1,DG,FTOL,GTOL,DGINIT,
c     +           MIN(FTOL,GTOL)*DGINIT
c      WRITE(*,*) 'STAGE1 FX FTEST1',STAGE1,FX,FTEST1

      IF(STAGE1.AND.F.LE.FX.AND.F.GT.FTEST1) THEN

*        use modified function to predict a step, if no step was
*        obtained for which the modified function has a non-positive
*        function value and non-negative derivative, and if   
*        a lower function value has been obtained, but with insufficient
*        decrease

         FM = F - STP*DGTEST    ! modify function value ...
         FXM = FX - STX*DGTEST
         FYM = FY - STY*DGTEST
         DGM = DG - DGTEST      ! ... and gradient values
         DGXM = DGX - DGTEST
         DGYM = DGY - DGTEST

*        update the interval and compute new step

         CALL CSTEPM(STX,FXM,DGXM,
     +               STY,FYM,DGYM,
     +               STP,FM ,DGM,
     +               BRACKT,STMIN,STMAX,INFOC)

         FX = FXM + STX*DGTEST  ! reset function value ...
         FY = FYM + STY*DGTEST
         DGX = DGXM + DGTEST    ! ... and gradient values
         DGY = DGYM + DGTEST

      ELSE

*        update the interval and compute new step

         CALL CSTEPM(STX,FX,DGX,
     +               STY,FY,DGY,
     +               STP,F ,DG,                       ! update interval
     +               BRACKT,STMIN,STMAX,INFOC)        ! and get new step
      END IF

      IF(BRACKT) THEN ! force sufficient decrease of size of interval
         IF(ABS(STY-STX).GE.P66*WIDTH1) STP=STX+P5*(STY-STX)
         WIDTH1 = WIDTH
         WIDTH = ABS(STY-STX)
      END IF

*     ---------- prepare calculation of next function value ------------

   30 IF(BRACKT) THEN               ! interval of uncertainty defined      
         STMIN=MIN(STX,STY)             ! minimum step
         STMAX=MAX(STX,STY)             ! maximum step
      ELSE                          ! interval of uncertainty undefined 
         STMIN=STX                      ! minimum step
         STMAX=STP+XTRAPF*(STP-STX)     ! maximum step
      END IF

      STP=MAX(STPMIN,MIN(STP,STPMAX)) ! STPMIN < STP < STPMAX

      IF((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     +   .OR. NFEV .GE. MAXF-1 .OR. INFOC .EQ. 0
     +   .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) THEN
*        unusual termination
c         WRITE(*,*) 'NFEV MAXF INFOC',NFEV,MAXF,INFOC
c         WRITE(*,*) 'Unusual termination, new,best STP  is',STP,STX
         STP = STX
c         WRITE(*,*) 'STX STPAC',STX,STPAC
         IF(NFEV.GE. MAXF-1.AND.STX.EQ.STPAC) THEN
            INFO=3
            STP=STPAC
            DGOUT=DG
c            WRITE(*,*) 'Special return !'
            RETURN
         END IF
      END IF

      DO J = 1, N          ! calculate the parameter values
       X(J)=WA(J)+STP*S(J) ! for step STP
      END DO
      INFO=-1 ! return to compute function value
      RETURN

      ENTRY MTINFO(INFO)  ! print information about errors
      IF(INFO.GE.0) WRITE(*,*) ' '
      IF(INFO.EQ.0) THEN
         WRITE(*,101) INFO,'N or STP lt 0 or step no descending'
      ELSE IF(INFO.EQ.1) THEN
         WRITE(*,101) INFO,'Linesearch convergence conditions reached'
      ELSE IF(INFO.EQ.2) THEN
         WRITE(*,101) INFO,'interval of uncertainty at lower limit'
      ELSE IF(INFO.EQ.3) THEN
         WRITE(*,101) INFO,'max nr of line search calls reached'
      ELSE IF(INFO.EQ.4) THEN
         WRITE(*,101) INFO,'step at the lower bound' 
      ELSE IF(INFO.EQ.5) THEN
         WRITE(*,101) INFO,'step at the upper bound'
      ELSE IF(INFO.EQ.6) THEN
         WRITE(*,101) INFO,'rounding error limitation'  
      END IF 
      IF(INFO.GE.0) WRITE(*,*) ' '
 101  FORMAT('     Info =',I3,3X,A)
      RETURN

*     MTLINE search initialization -------------------------------------

      ENTRY MTLDEF(FTOLE, GTOLE, XTOLE,  STMINE, STMAXE, MINFE,MAXFE)
*        example:  1.0E-4 1.0E-1 1.0E-17 1.0E-20 1.0E+20 1     40  
*     initialize "Linesearch Algorithm with Sufficient Decrease" 

      FTOL=FTOLE    ! tolerance C1  for sufficient decrease of function value
      GTOL=GTOLE    ! tolerance C2  for directional derivative condition
      XTOL=XTOLE    ! tolerance for width of intercal
      STPMIN=STMINE ! lower bound for the step
      STPMAX=STMAXE ! upper bound for the step
      MINF=MINFE    ! min number of linesearch function calls
      MAXF=MAXFE    ! max number of linesearch function calls (>1)
      IF(FTOL.LT.ZERO.OR.GTOL.LT.ZERO.OR.XTOL.LT.ZERO.OR.
     +   STMIN.LT.ZERO.OR.STPMAX.LT.STPMIN.OR.
     +   MINF.GT.MAXF.OR.MAXF.LE.1) THEN
         WRITE(*,*) 'Improper MTLDEF arguments - stop'
         STOP
      END IF 
      NFEV=-1               ! reset counter
      END


      SUBROUTINE CSTEPM(STX,FX,DX,  ! step, fvalue, der at best step
     +                  STY,FY,DY,  !                   at other endpt
     +                  STP,FP,DP,  !                   at current step 
     +                  BRACKT,        ! bracketed true/false
     +                  STPMIN,STPMAX, ! bounds for ste step
     +                  INFO)          ! returned: INFO=0 means error 
*                       INFO=1 ... 4   ! used method of interpolation
*     compute a safeguarded step for a linesearch
*     update an interval of uncertainty for the function minimum
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
*     ...
      INFO=0
      IF((BRACKT.AND.(STP.LE.MIN(STX,STY).OR.      ! check errors
     *    STP.GE.MAX(STX,STY))).OR.
     *    DX*(STP-STX).GE.0.0.OR.STPMAX.LT.STPMIN) RETURN

      SGND=DP*SIGN(1.0D0,DX) 
      IF(FP.GT.FX) THEN
         INFO = 1
*        higher function value: minimum is bracketed.
*        If the cubic step is closer to STX than the quadratic step,
*        the cubic step is taken, else the average of the two steps. 
         BOUND=.TRUE.
         THETA= 3.0*(FX - FP)/(STP - STX) + DX + DP
         S=MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA=S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF(STP.LT.STX) GAMMA = -GAMMA
         P=(GAMMA-DX)+THETA
         Q=((GAMMA-DX)+GAMMA)+DP
         R=P/Q
         STPC=STX+R*(STP-STX)
         STPQ=STX+((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP-STX)
         IF(ABS(STPC-STX).LT.ABS(STPQ-STX)) THEN
            STPF=STPC
         ELSE
            STPF=STPC+(STPQ-STPC)/2
         END IF
         BRACKT = .TRUE.
      ELSE IF(SGND.LT.0.0) THEN
         INFO = 2
*        lower function value and opp. derivative signs: minimum bracketed.
*        If the cubic step is closer to STX than the quadratic step,
*        the cubic step is taken, else the quadratic step.          
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF(ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
         END IF
         BRACKT = .TRUE.
      ELSEIF(ABS(DP).LT.ABS(DX)) THEN
         INFO = 3
*        lower function value and same derivative signs, der. decreases:
*        The cubic step is only used, if the cubic tends to infinity in
*        the driection of the step, or if the minimum of the cuic is beyond
*        STP. Otherwise the cubic step is defined to be either STPMIN or
*        STPMAX. The quadratic step also computed and if the minimum is
*        bracketed, then the step closest to STX is taken, else the step
*        farthest away is taken.
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF(R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
         END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF(BRACKT) THEN
            IF(ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         ELSE
            IF(ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         END IF
      ELSE
         INFO = 4
*        lower function value and same derivative signs, no der. decrease:
*        If the minimum is not bracketed, trhe step is either STPMIN or
*        STPMAX, else the cubic step is taken. 
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
         END IF
      END IF

*     update the interval of uncertainty (does not depend on thenew step) 
      IF(FP.GT.FX) THEN
         STY=STP
         FY =FP
         DY =DP
      ELSE
         IF(SGND.LT.0.0) THEN
            STY=STX
            FY =FX
            DY =DX
         END IF
         STX=STP
         FX =FP
         DX =DP
      END IF

      STPF=MAX(STPMIN,MIN(STPF,STPMAX))      
      STP =STPF                       ! new step
      IF(BRACKT.AND.BOUND) THEN       ! safeguard the step
         IF(STY.GT.STX) THEN
            STP= MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP= MAX(STX+0.66*(STY-STX),STP)
         END IF
      END IF
      END


