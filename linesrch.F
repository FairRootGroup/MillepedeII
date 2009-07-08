

      SUBROUTINE PTLINE(N,X,F,G,S,STEP, INFO) ! - 2 arguments
      IMPLICIT NONE
*     ------------------------------------------------------------------ 
*     "Line search routine with sufficient decrease of slope"
*
*     In many minimization problems the objective function is close to 
*     quadratic, except far from the solution. Close to the minimum the
*     behaviour may be almost quadratic or, due to round-off errors, 
*     it may have a non-smooth behaviour, which often complicates any
*     further progress and the recognition of convergence.  
*     Round-off errors affect the function value, which may be large and
*     small parameter changes result in small relative changes of the 
*     function value. Close to the minimum the gradient becomes small 
*     and the behaviour is not so much affected by Round-off errors.
*
*     ------------------------------------------------------------------
*
*      CALL PTLDEF(0.0,0.0, 0,0) ! init line search 
*      N=...
*      X(.)=...
*      D(.)=...
*      ALPHA=1.0D0
* 10   F(X)=...
*      G(X)=...
*      IF(.) S(X)=..
*      CALL PTLINE(N,X,F,G,D,ALPHA,INFO)
*      IF(INFO.LT.0) GOTO 10
*
*     ------------------------------------------------------------------
*     N       = dimension of problem
*     X(N)    = current iterate
*     F       = associated function value
*     G(N)    = associated gradient
*     S(N)    = search vector
*
*     ALPHA   = step factor (initially = 1.0) 
*
*     INFO    =     information 
*             = -1  repeat function evaluation
*             =  0  input error (e.g. gradient not negative)
*             =  1  convergence reached
*             =  2  convergence assumed, but round-off errors 
*             =  3  too many function calls
*             =  4  ALPHA to small (ALPHA <= TOL)
*
*     ------------------------------------------------------------------
      INTEGER N,INFO,MINFE,MAXFE                    ! arguments
      DOUBLE PRECISION X(N),F,G(N),S(N),STEP        ! arguments 
      REAL GTOLE,STMAX                              ! arguments
*     ------------------------------------------------------------------

C      DOUBLE PRECISION AUX(1000)

      INTEGER I,IM                                  ! internal
      DOUBLE PRECISION ALPHA                        ! internal
      DOUBLE PRECISION DGINIT,DG,FSAVED,TOT,FP1,FP2 ! internal
*     DGINIT = initial slope
*     DG     = current slope
*     FSAVE  = initla function value
*     TOT    = total step
*     -----------------------------------------------------------------
      INTEGER MSFD,NSFD,IDGL,IDGR,IDGM,MINF,MAXF,LSINFO 
      PARAMETER (MSFD=20)
      DOUBLE PRECISION SFD,STMX,GTOL
      COMMON/CLNSEA/SFD(4,MSFD),STMX,GTOL,NSFD,IDGL,IDGR,IDGM,MINF,MAXF
     +             ,LSINFO
*     -----------------------------------------------------------------
      INTEGER I1,I2
C      DOUBLE PRECISION CENTER,REVAL,COMP
*     -----------------------------------------------------------------
*     -----------------------------------------------------------------
      INTEGER LUNLOG
      COMMON/CLULOG/LUNLOG
*     -----------------------------------------------------------------
      SAVE
* 
*     initialization ---------------------------------------------------
*

      INFO=0             ! reset INFO flag
      DG=0.0D0
      DO I=1,N           !
       DG=DG-G(I)*S(I)   ! DG = scalar product: grad x search
      END DO!

      IF(NSFD.EQ.0) THEN    ! initial call
         DGINIT=DG          ! DG = initial directional gradient
         IF(DGINIT.GE.0.0D0) GOTO 100 ! error: step not decreasing
         STEP=1.0D0         ! initial step factor is one 
         ALPHA=STEP         ! get initial step factor
         TOT=0.0D0          ! reset total step
         IDGL=1             ! index of smallest negative slope 
         IDGR=0             ! index of smallest positive slope 
         FSAVED=F           ! initial Function value
         NSFD=1             ! starting point of iteration
         SFD(1,1)=0.0       ! abscissa
         SFD(2,1)=0.0       ! referencefunction value 
         SFD(3,1)=DGINIT    ! slope
         SFD(4,1)=0.0       ! predicted zero 
         IM=1               ! optimum
      ELSE                  ! subsequent call
         NSFD=NSFD+1
         SFD(1,NSFD)=TOT          ! abscissa
         SFD(2,NSFD)=F-FSAVED     ! function value difference to reference
         SFD(3,NSFD)=DG           ! slope
         SFD(4,NSFD)=0.0          ! predicted zero (see below)  
         IF(DG.LT.SFD(3,IM)) THEN
            IM=NSFD
C            DO I=1,N
C             AUX(I)=G(I)
C            END DO 
         END IF

*        define interval indices IDGL and IDGR 
         IF(DG.LE.0.0D0) THEN
            IF(DG.GE.SFD(3,IDGL)) IDGL=NSFD 
         END IF
         IF(DG.GE.0.0D0) THEN     ! limit to the right
            IF(IDGR.EQ.0) IDGR=NSFD 
            IF(DG.LE.SFD(3,IDGR)) IDGR=NSFD
         END IF 

         IF(IDGR.EQ.0) THEN
            I1=NSFD-1
            I2=NSFD
         ELSE
            I1=IDGL
            I2=IDGR
         END IF
         FP1=SFD(3,I1)
         FP2=SFD(3,I2)                       ! interpolation
         SFD(4,NSFD)=(SFD(1,I1)*FP2-SFD(1,I2)*FP1)/(FP2-FP1)
c         CENTER=0.5D0*(SFD(1,I1)+SFD(1,I2))
c         REVAL=ABS(2.0D0*MIN(ABS(FP1),ABS(FP2))/ABS(FP1-FP2)-1.0D0) 
c         IF(ABS(FP1).GT.ABS(FP2)) THEN
c            FP2=FP2+(0.01*FP2-0.01*FP1)*REVAL**2  
c         ELSE
c            FP1=FP1+0.05D0*(FP1-FP2)*REVAL**2
c         END IF 
c         COMP=(SFD(1,I1)*FP2-SFD(1,I2)*FP1)/(FP2-FP1)
c         WRITE(LUNLOG,*) 'compare',SFD(4,NSFD),COMP,' at',REVAL
c         IF(ABS(DG)/ABS(DGINIT).GT.0.1D0.AND.FP1*FP2.LT.0.0D0) THEN
c            WRITE(LUNLOG,*) 'changed'
c            SFD(4,NSFD)=COMP
c         END IF

c         IF(IDGR.NE.0) THEN
c            CENTER=0.5D0*(SFD(1,I1)+SFD(1,I2))
c            SFD(4,NSFD)=0.25D0*(CENTER+3.0D0*SFD(4,NSFD)) 
c         END IF

*        convergence tests
         IF(NSFD.GE.MINF.AND.ABS(DG).LE.ABS(DGINIT)*GTOL) THEN
*           normal convergence return with INFO=1 ----------------------
            ALPHA=TOT+ALPHA          ! total ALPHA is returned
            STEP =ALPHA
            IDGM=IDGL
            IF(IDGR.NE.0) THEN
               IF(SFD(3,IDGR)+SFD(3,IDGL).LT.0.0D0) IDGM=IDGR
            END IF
#ifdef graphics 
            CALL PLSFD
            CALL PLDERL
#endif
            GOTO 101
         END IF
         IF(NSFD.GE.MAXF) GOTO 102 ! max number of function calls
         ALPHA=MIN(SFD(4,NSFD),STMX)-TOT     ! new step from previous 
         IF(ABS(ALPHA).LT.1.0D-3.AND.SFD(4,NSFD).GT.STMX) GOTO 103   
         IF(ABS(ALPHA).LT.1.0D-3) GOTO 104
      END IF

*     prepare next function call ---------------------------------------

      DO I=1,N
       X(I)=X(I)+ALPHA*S(I)    ! step by ALPHA -> new X
      END DO
      TOT=TOT+ALPHA            ! 
      STEP=TOT
      INFO=-1                  ! recalculate function and gradient
      LSINFO=INFO 
      RETURN

*     error exits ------------------------------------------------------    
 104  INFO=INFO+1              ! 4: step small
 103  INFO=INFO+1              ! 3: maximum reached
 102  INFO=INFO+1              ! 2: too many function calls
 101  INFO=INFO+1              ! 1: normal convergence
      LSINFO=INFO
      IM=1
      DO I=1,NSFD 
       IF(ABS(SFD(3,I)).LT.ABS(SFD(3,IM))) IM=I
      END DO 
      ALPHA=SFD(1,IM)-SFD(1,NSFD)
      IF(IM.EQ.NSFD) RETURN    ! already at minimum
      DO I=1,N
       X(I)=X(I)+ALPHA*S(I)    ! step by ALPHA to slope minimum
C      G(I)=AUX(I)             ! copy gradient at IM  
      END DO
      F=SFD(2,IM)+FSAVED       ! F at minimum
      STEP=SFD(1,IM)           ! total step at convergence 
      IF(IM.NE.1) RETURN       ! improvement
      INFO=5                   ! no improvement
 100  STEP=0.0D0               ! 0: initial slope not negative
      LSINFO=INFO  
      RETURN

*     ================================================================== 

      ENTRY PTLDEF(GTOLE,STMAX,MINFE,MAXFE)
*
*     initialize "Linesearch Algorithm with Sufficient Slope Decrease"
*
*                        --- range ----       default
*     slope ratio        1.0E-4 ... 0.9         0.9 
*     min. F-calls       1 ... 2                 1 
*     max. F-calls       2 ... 10                5 
*
      GTOL=MAX(1.0E-4,MIN(GTOLE,0.9E0))  ! slope ration
      IF(GTOLE.EQ.0.0) GTOL=0.9D0   ! default slope ratio
      STMX=STMAX                    ! maximum total step
      IF(STMX.EQ.0.0D0) STMX=10.0D0 ! default limit
      MINF=MAX(1,MIN(MINFE,MSFD-2)) ! minimum number of evaluations
      MAXF=MAX(2,MIN(MAXFE,MSFD-1)) ! maximum number of evaluations
      IF(MAXFE.EQ.0) MAXF=5         ! default max number of values
      NSFD=0                        ! reset
      END


      SUBROUTINE PTLOPT(NF,M,SLOPES,STEPS)
*     return NF     = number of function values
*            M      = index of function value with smallest slope
*            SLOPES(1) = initial slope
*            SLOPES(2) = current slope         
*            SLOPES(3) = smallest slope
*            STEPS(1)  = initial position
*            STEPS(2)  = current step
*            STEPS(3)  = smallest step

*     -----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NF,M,I
      REAL SLOPES(3),STEPS(3)
*     -----------------------------------------------------------------
      INTEGER MSFD,NSFD,IDGL,IDGR,IDGM,MINF,MAXF,LSINFO
      PARAMETER (MSFD=20)
      DOUBLE PRECISION SFD,STMX,GTOL
      COMMON/CLNSEA/SFD(4,MSFD),STMX,GTOL,NSFD,IDGL,IDGR,IDGM,MINF,MAXF
     +             ,LSINFO
*     -----------------------------------------------------------------
*     ...
      NF=NSFD
      IF(NSFD.EQ.0) THEN  ! no values
         M=0
         DO I=1,3
          SLOPES(I)=0.0
          STEPS(I) =0.0
         END DO 
      ELSE                ! values exist
         M=1
         DO I=1,NSFD
          IF(ABS(SFD(3,I)).LT.ABS(SFD(3,M))) M=I
         END DO
         SLOPES(1)=SFD(3,1)       
         SLOPES(2)=SFD(3,NSFD)
         SLOPES(3)=SFD(3,M)
         STEPS(1) =SFD(1,1)
         STEPS(2) =SFD(1,NSFD)
         STEPS(3) =SFD(1,M)
      END IF
      END

      SUBROUTINE PTLPRT
*     print line search data
      IMPLICIT NONE
*     -----------------------------------------------------------------
      INTEGER MSFD,NSFD,IDGL,IDGR,IDGM,MINF,MAXF,LSINFO
      PARAMETER (MSFD=20)
      DOUBLE PRECISION SFD,STMX,GTOL
      COMMON/CLNSEA/SFD(4,MSFD),STMX,GTOL,NSFD,IDGL,IDGR,IDGM,MINF,MAXF
     +             ,LSINFO
*     -----------------------------------------------------------------
      INTEGER LUN,LUNLOG
      COMMON/CLULOG/LUNLOG
*     -----------------------------------------------------------------
      INTEGER I,J,IM
      REAL RATIO
      CHARACTER*2 TLR
*     ...
      LUN=LUNLOG
      IF(LUN.EQ.0) LUN=6
      IF(NSFD.LE.0) RETURN
      WRITE(LUN,*) ' '
      WRITE(LUN,*) 'PTLINE: line-search method based on slopes',
     +             ' with sufficient slope-decrease'
      WRITE(LUN,*) 'PTLDEF: slope ratio limit=',GTOL
      WRITE(LUN,*) 'PTLDEF: maximum step =',STMX
      WRITE(LUN,*) 'PTLDEF:',MINF,' <= nr of calls <=',MAXF
      WRITE(LUN,101)
      IM=1
      DO I=1,NSFD
       IF(ABS(SFD(3,I)).LT.ABS(SFD(3,IM))) IM=I
      END DO
      DO I=1,NSFD
       TLR='  '
       IF(I.EQ.IM)   TLR='**'
       IF(I.EQ.IDGL) TLR(1:1)='L'
       IF(I.EQ.IDGR) TLR(2:2)='R'
       IF(I.EQ.1) THEN
          WRITE(LUN,102) I-1, SFD(1,I),TLR,(SFD(J,I),J=2,4)
       ELSE
          RATIO=ABS(SFD(3,I)/SFD(3,1))
          WRITE(LUN,103) I-1, SFD(1,I),TLR,(SFD(J,I),J=2,4),RATIO
       END IF 

      END DO
      IF(LSINFO.EQ.0) WRITE(LUN,*)
     +    'PTLINE: INFO=0  input error (e.g. gradient not negative)'
      IF(LSINFO.EQ.1) WRITE(LUN,*)
     +    'PTLINE: INFO=1  convergence reached' 
      IF(LSINFO.EQ.2) WRITE(LUN,*)
     +    'PTLINE: INFO=2  too many function calls'
      IF(LSINFO.EQ.3) WRITE(LUN,*)
     +    'PTLINE: INFO=3  maximum step reached'
      IF(LSINFO.EQ.4) WRITE(LUN,*)
     +    'PTLINE: INFO=4  step too small (< 0.001)' 
      WRITE(LUN,*) ' '  

#ifdef graphics
            CALL PLSFD
            CALL PLDERL
#endif

 101  FORMAT('  i      x              F(x)         F''(x)',
     +       '          minimum     F''(x)')
 102  FORMAT(I3,F12.6,1X,A2,G15.6,G14.6,F12.6,'     ratio')
 103  FORMAT(I3,F12.6,1X,A2,G15.6,G14.6,F12.6,F10.3)

c 101  FORMAT('  i          x              F(x)         F''(x)',
c     +       '              mininum F''(x)')
c 102  FORMAT(I3,F16.10,1X,A2,G15.6,G14.6,F16.10,' ratio')
c 103  FORMAT(I3,F16.10,1X,A2,G15.6,G14.6,F16.10,F6.2)
      END

#ifdef graphics
      SUBROUTINE PLSFD
*     plot line search data
      IMPLICIT NONE
*     -----------------------------------------------------------------
      INTEGER MSFD,NSFD,IDGL,IDGR,IDGM,MINF,MAXF,LSINFO
      PARAMETER (MSFD=20)
      DOUBLE PRECISION SFD,STMX,GTOL
      COMMON/CLNSEA/SFD(4,MSFD),STMX,GTOL,NSFD,IDGL,IDGR,IDGM,MINF,MAXF
     +             ,LSINFO
*     -----------------------------------------------------------------
      REAL XG,FUNG,XMAX,DELD,DELG,XELG,ARG1,ARG2
      INTEGER I,J 
*     -----------------------------------------------------------------
*     ...
      CALL PFINIT
      XMAX=0.0D0
      DO I=2,NSFD
       XMAX=MAX(XMAX,SFD(1,I))
      END DO

*     parabola

      CALL PLATT('green full line')   ! parabola
      DO I=-10,110,5
       XG=0.01*FLOAT(I)*XMAX
       FUNG=SFD(3,1)*XG*(1.0-0.5*XG)
       CALL PDNXY(1,XG,FUNG,'f') ! parabola
      END DO

      CALL PLATT('cyan full line')   ! parabola
      DO J=2,NSFD
       DO I=-10,110,5
        XG=0.01*FLOAT(I)*XMAX
        IF(SFD(4,J).NE.0.0D0) THEN
           FUNG=SFD(3,1)*XG-0.5*SFD(3,1)/SFD(4,J)*XG**2
           CALL PDNXY(1,XG,FUNG,'f') ! parabola
        END IF   
       END DO
       CALL PLATT('cyan full line')  
      END DO


*     next attributes

      CALL PLATT('black bars line 2.0')

*     function value and slopes 

      DELG=0.1
      DO I=1,NSFD
       XELG=SFD(1,I)
       FUNG=SFD(2,I)
       DELD=DELG*SFD(3,I)
       CALL PD2XY(XELG-DELG,FUNG-DELD,
     +            XELG+DELG,FUNG+DELD,'l')
      END DO 



      CALL PLATT('red fcircle symbol')
      DO I=1,NSFD
       XELG=SFD(1,I)
       FUNG=SFD(2,I)
       CALL PDSYMB(XELG,FUNG,0.0,0.0)
      END DO
      CALL PLSHOW
      RETURN

*     ------------------------------------------------------------------      
      ENTRY PLDERL
*     plot derivative dependence 
      CALL PFINIT
      CALL PTITLE(0,'slope dependence')
c      CALL PLATT('black full line')
      CALL PLATT('black plus marker 4.0') 
      XMAX=0.0D0
      DO I=1,NSFD
       XMAX=MAX(XMAX,SFD(1,I))
       ARG1=SFD(1,I)
       ARG2=SFD(3,I)
       CALL PDNXY(1,ARG1,ARG2,'m')
      END DO
      CALL PLATT('red full line')
      CALL PLSHOW
      END
#endif


