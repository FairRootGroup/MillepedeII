
      PROGRAM MPTWO                              ! Millepede II main program

*     Millepede II main program

      CHARACTER*24 CHDATE,CHOST 
      REAL RST,RSTP
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
* MPTEST      
      PARAMETER (NPLAN=100)                ! test data ...
      COMMON/TDETEC/DEL(NPLAN),DVD(NPLAN)  ! ...of misalignment
* MPTST2      
      PARAMETER (NLYR=10)                  ! #layers
      PARAMETER (NMX=10)                    ! #module in X
      PARAMETER (NMY=5)                    ! #module in Y 
      PARAMETER (NTOT=NLYR*NMX*NMY)      
      COMMON/SDEVEC/SDEVX(NTOT),SDEVY(NTOT)
*         
      DOUBLE PRECISION SUMS(9)
!$    INTEGER OMP_GET_NUM_PROCS,OMP_GET_MAX_THREADS

      SAVE
*     ...
      CALL ETIME(TA,RSTP)
      CALL FDATE(CHDATE)

*     millepede.log file
      LUNLOG=8
      LVLLOG=1
      CALL MVOPEN(LUNLOG,'millepede.log')
      CALL GETENV('HOSTNAME',CHOST)
      IF (CHOST(1:1).EQ.' ') CALL GETENV('HOST',CHOST)
      WRITE(*,*) '($Rev$)'
!$    WRITE(*,*) 'using OpenMP (TM)'
      WRITE(*,*) ' '
      WRITE(*,*) '  <  Millepede II-P starting ... ',CHDATE 
      WRITE(*,*) '                                 ',CHOST
      WRITE(*,*) ' '
      CALL MEGNIT             ! initialize dynal memory management
      WRITE(*,*) ' '
            
      WRITE(8,*) ' '
      WRITE(8,*) 'Log-file Millepede II-P                        ',
     +            CHDATE
      WRITE(8,*) '                                               ',
     +            CHOST
C     CALL MEGDEB  ! debug switched on

*     read command line and text files

      CALL FILETC   ! command line and steering file analysis
      CALL FILETX   ! read text files
      LVLLOG=MPRINT ! export print level
!$    WRITE(*,*)
!$    NPROC=1
!$    MXTHRD=1
!$    NPROC=OMP_GET_NUM_PROCS()         ! number of processors available
!$    CALL OMP_SET_NUM_THREADS(MTHRD)   ! set max number of threads to MTHRD
!$    MXTHRD=OMP_GET_MAX_THREADS()      ! get max number of threads back
!$    WRITE(*,*) 'Number of processors available:   ', NPROC
!$    WRITE(*,*) 'Maximum number of OpenMP threads: ', MXTHRD
!$    WRITE(*,*) 'Number of threads for processing: ', MTHRD
!$    IF (MXREC.GT.0) MTHRDR=1          ! to get allways the same MXREC records
!$    WRITE(*,*) 'Number of threads for reading:    ', MTHRDR
!$POMP INST INIT                        ! start profiling with ompP

*     histogram file
      LUN=7
      CALL MVOPEN(LUN,'millepede.his')
      CALL HMPLUN(LUN) ! unit for histograms
      CALL GMPLUN(LUN) ! unit for xy data  

*     debugging
      IF(NRECPR.NE.0.OR.NRECP2.NE.0) THEN
         CALL MVOPEN(1,'mpdebug.txt')
      END IF 

      CALL ETIME(TA,RSTEXT)
      TIMES(0)=RSTEXT-RSTP ! time for text processing
      IF(INTRAC.EQ.1) THEN 
         CALL WAITST(ISTOP,'command line and text files read')
         IF(ISTOP.NE.0) STOP  
      END IF
*
c      METSOL=4

*     preparation of data sub-arrays

      CALL LOOP1
      CALL ETIME(TA,RLOOP1)
      TIMES(1)=RLOOP1-RSTEXT ! time for LOOP1 
      IF(INTRAC.EQ.1) THEN 
         CALL WAITST(ISTOP,'Loop1 finished')
         IF(ISTOP.NE.0) STOP
      END IF


      CALL LOOP2
      IF(CHICUT.NE.0.0) THEN
         WRITE(8,*) 'Chi square cut equiv 3 st.dev applied ...'
         WRITE(8,*) ' in  first iteration with factor',CHICUT
         WRITE(8,*) ' in second iteration with factor',CHIREM
         WRITE(8,*) ' (reduced by sqrt in next iterations)'
      END IF
      IF(LHUBER.NE.0) THEN
         WRITE(8,*) 'Down-weighting of outliers in',
     +               LHUBER,' iterations'
         WRITE(8,*) 'Cut on downweight fraction',DWCUT
      END IF


      CALL ETIME(TA,RLOOP2)
      TIMES(2)=RLOOP2-RLOOP1 ! time for LOOP2 
      IF(INTRAC.EQ.1) THEN  
         CALL WAITST(ISTOP,'Loop2 finished')
         IF(ISTOP.NE.0) STOP 
      END IF 

*     use different solution methods

      CALL MSTART('Iteration')   ! Solution module starting

c      CALL XEXECS               ! all methods

      CALL XLOOPN                ! all methods

c      IF(NREJEC(1)+NREJEC(2)+NREJEC(3).NE.0) THEN
c         WRITE(*,*) 'Data rejected in last loop:   ',NREJEC(1),
c     +   ' (Ndf=0)   ',NREJEC(2),' (huge)   ',NREJEC(3),' (large)'
c       END IF



*     ------------------------------------------------------------------

      IF(NLOOPN.GT.2.AND.NHISTP.NE.0) THEN       ! last iteration
         CALL HMPRNT(3)  ! scaled residual of single measurement (with global deriv.)
         CALL HMPRNT(12) ! scaled residual of single measurement (no global deriv.)
         CALL HMPRNT(4)  ! chi^2/Ndf
      END IF
      IF(NLOOPN.GT.2) THEN
         CALL HMPWRT(3)
         CALL HMPWRT(12)
         CALL HMPWRT(4)
         CALL GMPWRT(4) ! location, dispersion (res.) as a function of record nr
         IF (NLOOPN.LE.LFITNP) THEN
            CALL HMPWRT(13)
            CALL HMPWRT(14)
            CALL GMPWRT(5)
         ENDIF
      END IF
      IF(NHISTP.NE.0) THEN
         CALL GMPRNT(1)       
         CALL GMPRNT(2)
      END IF
      CALL GMPWRT(1)             ! output of xy data 
      CALL GMPWRT(2)             ! output of xy data

      IF(ICTEST.EQ.1) THEN
         WRITE(*,*) ' '
         WRITE(*,*) 'Misalignment test wire chamber'
         WRITE(*,*) ' '

         CALL HMPDEF( 9,-0.0015,+0.0015,'True - fitted displacement')
         CALL HMPDEF(10,-0.0015,+0.0015,'True - fitted Vdrift')
         DO I=1,4
          SUMS(I)=0.0D0
         END DO
         DO I=1,100
          DIFF=SNGL(-DEL(I)-DQ(IND0/2+I))
          SUMS(1)=SUMS(1)+DIFF
          SUMS(2)=SUMS(2)+DIFF*DIFF
          DIFF=SNGL(-DVD(I)-DQ(IND0/2+100+I))
          SUMS(3)=SUMS(3)+DIFF
          SUMS(4)=SUMS(4)+DIFF*DIFF
         END DO
         SUMS(1)=0.01D0*SUMS(1)
         SUMS(2)=SQRT(0.01D0*SUMS(2)) 
         SUMS(3)=0.01D0*SUMS(3)
         SUMS(4)=SQRT(0.01D0*SUMS(4))
         WRITE(*,143) 'Parameters   1 - 100: mean =',SUMS(1),
     +               'rms =',SUMS(2)                       
         WRITE(*,143) 'Parameters 101 - 200: mean =',SUMS(3),
     +               'rms =',SUMS(4)                       
 143     FORMAT(6X,A28,F9.6,3X,A5,F9.6)
         WRITE(*,*) ' '
         WRITE(*,*) ' ' 
         WRITE(*,*) '    I '
         WRITE(*,*) '   --- '             
         DO I=1,100
          WRITE(*,102) I,-DEL(I),DQ(IND0/2+I),-DEL(I)-DQ(IND0/2+I),
     +                -DVD(I),DQ(IND0/2+100+I),-DVD(I)-DQ(IND0/2+100+I)
          DIFF=SNGL(-DEL(I)-DQ(IND0/2+I))
          CALL HMPENT( 9,DIFF)
          DIFF=SNGL(-DVD(I)-DQ(IND0/2+100+I))
          CALL HMPENT(10,DIFF)
         END DO 
         IF(NHISTP.NE.0) THEN
            CALL HMPRNT( 9)
            CALL HMPRNT(10)
         END IF  
         CALL HMPWRT( 9)
         CALL HMPWRT(10)
      END IF
      IF(ICTEST.GT.1) THEN
         WRITE(*,*) ' '
         WRITE(*,*) 'Misalignment test Si tracker'
         WRITE(*,*) ' '

         CALL HMPDEF( 9,-0.0025,+0.0025,'True - fitted displacement X')
         CALL HMPDEF(10,-0.025,+0.025,'True - fitted displacement Y')
         DO I=1,9
          SUMS(I)=0.0D0
         END DO
        NMXY=NMX*NMY
        IX=0
        IY=NTOT
        DO I=1,NLYR
         DO K=1,NMXY
          IX=IX+1
          DIFF=SNGL(-SDEVX((I-1)*NMXY+K)-DQ(IND0/2+IX))
          SUMS(1)=SUMS(1)+1.0D0
          SUMS(2)=SUMS(2)+DIFF
          SUMS(3)=SUMS(3)+DIFF*DIFF
          IXV=JVGBI(IX)
          IF (IXV.GT.0.AND.METSOL.EQ.1.OR.METSOL.EQ.2) THEN
           II=(IXV*IXV+IXV)/2
           GMATI=SNGL(DQ(IGMAT/2+II))
           ERR=SQRT(ABS(GMATI))
           DIFF=DIFF/ERR
           SUMS(7)=SUMS(7)+1.0D0
           SUMS(8)=SUMS(8)+DIFF
           SUMS(9)=SUMS(9)+DIFF*DIFF
          END IF
         END DO
         IF (MOD(I,3).EQ.1) THEN
          DO K=1,NMXY
           IY=IY+1
           DIFF=-SNGL(SDEVY((I-1)*NMXY+K)-DQ(IND0/2+IY))
           SUMS(4)=SUMS(4)+1.0D0
           SUMS(5)=SUMS(5)+DIFF
           SUMS(6)=SUMS(6)+DIFF*DIFF
           IXV=JVGBI(IY)
           IF (IXV.GT.0.AND.METSOL.EQ.1.OR.METSOL.EQ.2) THEN
            II=(IXV*IXV+IXV)/2
            GMATI=SNGL(DQ(IGMAT/2+II))
            ERR=SQRT(ABS(GMATI))
            DIFF=DIFF/ERR
            SUMS(7)=SUMS(7)+1.0D0
            SUMS(8)=SUMS(8)+DIFF
            SUMS(9)=SUMS(9)+DIFF*DIFF
           END IF
          END DO       
         END IF
         END DO
         SUMS(2)=SUMS(2)/SUMS(1)
         SUMS(3)=SQRT(SUMS(3)/SUMS(1)) 
         SUMS(5)=SUMS(5)/SUMS(4)
         SUMS(6)=SQRT(SUMS(6)/SUMS(4))
         WRITE(*,143) 'Parameters   1 - 500: mean =',SUMS(2),
     +               'rms =',SUMS(3)                       
         WRITE(*,143) 'Parameters 501 - 700: mean =',SUMS(5),
     +               'rms =',SUMS(6)                       
         IF (SUMS(7).GT.0.5D0) THEN
          SUMS(8)=SUMS(8)/SUMS(7)
          SUMS(9)=SQRT(SUMS(9)/SUMS(7))
          WRITE(*,143) 'Parameter pulls, all: mean =',SUMS(8),
     +               'rms =',SUMS(9)               
        END IF
        WRITE(*,*) ' '
         WRITE(*,*) ' ' 
         WRITE(*,*) '    I '
         WRITE(*,*) '   --- '             
        IX=0
        IY=NTOT
         DO I=1,NLYR
         DO K=1,NMXY
          IX=IX+1
           DIFF=SNGL(-SDEVX((I-1)*NMXY+K)-DQ(IND0/2+IX))
           CALL HMPENT( 9,DIFF)                  
           WRITE(*,102) IX,-SDEVX((I-1)*NMXY+K),DQ(IND0/2+IX),-DIFF
         END DO
        END DO 
         DO I=1,NLYR
         IF (MOD(I,3).EQ.1) THEN
          DO K=1,NMXY
           IY=IY+1
            DIFF=SNGL(-SDEVY((I-1)*NMXY+K)-DQ(IND0/2+IY))
            CALL HMPENT(10,DIFF)                 
            WRITE(*,102) IY,-SDEVY((I-1)*NMXY+K),DQ(IND0/2+IY),-DIFF
          END DO       
         END IF
         END DO         
         IF(NHISTP.NE.0) THEN
            CALL HMPRNT( 9)
            CALL HMPRNT(10)
         END IF  
         CALL HMPWRT( 9)
         CALL HMPWRT(10)
      END IF
           
      IF(NREC1+NREC2.GT.0) THEN
         WRITE(8,*) ' '
         IF(NREC1.GT.0) THEN
            WRITE(8,*) 'Record',NREC1,' has largest residual:',VALUE1
         END IF
         IF(NREC2.GT.0) THEN
            WRITE(8,*) 'Record',NREC2,' has largest Chi^2/Ndf:',VALUE2
         END IF
      END IF
      IF(NREC3.LT.MAXI4) THEN
         WRITE(8,*) 'Record',NREC3,
     +               ' is first with error (rank deficit/NaN)'
      ENDIF
      WRITE(8,*) ' '
      WRITE(8,*) 'In total 2 +',NLOOPN,' loops through the data files'
      IF (MNRSIT.GT.0) THEN
         WRITE(8,*) ' '
         WRITE(8,*) 'In total    ',MNRSIT,' internal MINRES iterations'
      ENDIF

      WRITE(8,104) GBM,UBM,UBM/GBM*100.0
 104  FORMAT(/' Memory space: total   ',F13.6,' GB'/
     +        '                used   ',F13.6,' GB = ',F9.2,' %')


      WRITE(8,103) TIMES(0),TIMES(1),TIMES(2),TIMES(4),TIMES(7),   
     +             TIMES(5),TIMES(8),TIMES(3),TIMES(6)
      CALL ETIME(TA,RST) 
      DELTAT=RST-RSTP
      NTSEC=NINT(DELTAT)
      CALL SECHMS(DELTAT,NHOUR,MINUT,SECND)
      NSECND=NINT(SECND)  ! round
      WRITE(8,*) 'Total time =',NTSEC,' seconds =',NHOUR,' h',MINUT,
     +           ' m',NSECND,' seconds'
      CALL FDATE(CHDATE)
      WRITE(8,*) 'end                                            ',
     +            CHDATE


*     Rejects ----------------------------------------------------------

      IF(NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3).NE.0) THEN
         WRITE(8,*) ' '
         WRITE(8,*) 'Data rejected in last iteration:   '
         WRITE(8,*) '   ',
     +   NREJEC(0), ' (rank deficit/NaN) ',NREJEC(1),' (Ndf=0)   ',
     +   NREJEC(2), ' (huge)   ',NREJEC(3),' (large)'     
         WRITE(8,*) ' '
      END IF
      CALL EXPLFC(8)

      WRITE(*,*) ' '
      WRITE(*,*) '  <  Millepede II-P ending   ... ', CHDATE ! with exit code',ITEXIT,' >'
      WRITE(*,*) ' '
C      CALL EXIT(ITEXIT)



 102  FORMAT(2X,I4,2X,3F10.5,2X,3F10.5)
 103  FORMAT(' Times [in sec] for     text processing',F12.3/
     +       '                                  LOOP1',F12.3/
     +       '                                  LOOP2',F12.3/
     +       '   func. value                         ',F12.3,' *',F4.0/
     +       '   func. value, global matrix, solution',F12.3,' *',F4.0/
     +       '                           new solution',F12.3,' *',F4.0/)
      END 

      SUBROUTINE SOLGLO(IVGBI) ! single error from minres
#include "largeint.inc"
#include "dynal.inc"
      DOUBLE PRECISION SHIFT,RTOL,ANORM,ACOND,RNORM,YNORM,GMATI,DIAG
      EXTERNAL AVPROD, MVSOLV
      SAVE
#include "mpinds.inc"
      DATA IPH/0/
*     ...
      IF(IPH.EQ.0) THEN
         IPH=1
         WRITE(*,101)
      END IF
      ITGBI=JTGBI(IVGBI)
      ITGBL=JTGBL(ITGBI)

      CALL MEGZER('f') ! reset rhs vector IGVEC 
      DQ(IGVEC/2+IVGBI)=1.0D0

c      NOUT  =6
      NOUT  =0
      ITNLIM=200
      SHIFT =0.0D0
      RTOL  = MRESTL ! from steering

      NVAR=NAGB
      CALL MINRES(NVAR,DQ(IGVEC/2+1),
     +   DQ(IDUX1/2+1),       DQ(IDUX1/2+1+NAGB),  DQ(IDUX1/2+1+2*NAGB),
     +   DQ(IDUX1/2+1+3*NAGB),DQ(IDUX1/2+1+4*NAGB),DQ(IDUX1/2+1+5*NAGB),
     +                 DQ(ISOLV/2+1),DQ(IDUX7/2+1),
     +                 AVPROD, MVSOLV, .TRUE. , .FALSE. , SHIFT,
     +                 NOUT , ITNLIM, RTOL,
     +                 ISTOP, ITN, ANORM, ACOND, RNORM, YNORM)

c      subroutine MINRES( n, b, r1, r2, v, w, w1, w2, x, y,
c     $                   AVPROD, Msolve, .TRUE. , .FALSE. , SHIFT,
c     $                   NOUT , ITNLIM, rtol,
c     $                   ISTOP, ITN, Anorm, Acond, rnorm, ynorm )

      PAR=SNGL(DQ(IND0/2+ITGBI))
      DPA=PAR-QM(IND2+ITGBI)
      GMATI=DQ(ISOLV/2+IVGBI)
      ERR=SQRT(ABS(SNGL(GMATI)))
      IF(GMATI.LT.0.0D0) ERR=-ERR
      IF(MQ(IGMAT-2).EQ.1) THEN ! normal matrix               ! ???
         JK=IJSYM(IVGBI,IVGBI)
      ELSE IF(MQ(IGMAT-2).EQ.2) THEN ! sparse matrix
         JK=IJADD(IVGBI,IVGBI)
      END IF
      IF (JK.GT.0) THEN
         DIAG=DQ(IGMAT/2+JK)
      ELSE
         DIAG=DBLE(QM(IDOT1-JK))
      ENDIF
      GCOR2=SNGL(1.0D0-1.0D0/(GMATI*DIAG)) ! global correlation (squared)
      WRITE(*,102) ITGBL,PAR,QM(IND1+ITGBI),DPA,ERR,GCOR2,ITN
 101  FORMAT(1X,'    label     parameter    presigma      differ',
     +                   '       Error gcor^2   iit'/
     +       1X,'---------',2X,5('-----------'),2X,'----')
 102  FORMAT(I10,2X,4G12.4,F7.4,I6,I4)
      END

      SUBROUTINE ADDCST
*     add constraint information to matrix and vector
#include "largeint.inc"
#include "dynal.inc"
      DOUBLE PRECISION RHS,DRHS(4)
      INTEGER             IDRH (4)
      SAVE
#include "mpinds.inc"
*     ...
      NOP=0
      HUGE=1.0 !  E6
      IF(IND8.EQ.0) RETURN  ! no constraints
      CLIMIT=1.0E-5         ! limit for printout 
c      WRITE(*,102) ' Constraints: only equation values >',
c     +               CLIMIT,' are printed'
      IRHS=0 ! number of values in DRHS(.), to be printed

      I=0
      DO ICGB=1,NCGB
c       WRITE(*,*) 'Constraint ',ICGB 
       RHS=QM(IND8+I+2) ! right hand side
c       RHSI=RHS
       SGM=QM(IND8+I+4) ! sigma parameter
c       WRITE(*,*) 'ADDCST   ',ICGB,RHS,SGM
       I=I+4
 10    LABEL=MQ(IND8+I+1)
       FACTR=QM(IND8+I+2)
       ITGBI=INONE(LABEL) ! -> ITGBI= index of parameter label
Cno       ITGBI=MQ(IND8+I+1) ! already translated 
       IVGB =JVGBI(ITGBI) ! -> variable-parameter index
c       IF(IVGB.LE.0) WRITE(*,*) 'Error in ADDCST1: LABEL ITGBI IVGB',
c     +                LABEL,ITGBI,IVGB

       IF(ICALCM.EQ.1.AND.IVGB.GT.0) THEN   
          CALL MUPDAT(NVGB+ICGB,IVGB,DBLE(HUGE*FACTR)) ! add to matrix
       END IF 

       RHS=RHS-FACTR*DQ(IND0/2+ITGBI)     ! reduce residuum
       I=I+2
       IF(I.LT.LND8) THEN
          IF(MQ(IND8+I+1).NE.0) GOTO 10
       END IF
c       WRITE(*,*) 'RHS is reduced to ',RHS,' from',RHSI 
       IF(ABS(RHS).GT.CLIMIT) THEN
          IRHS=IRHS+1
          IDRH(IRHS)=ICGB 
          DRHS(IRHS)=RHS
          NOP=1
          IF(IRHS.EQ.4) THEN
             WRITE(*,101) (IDRH(L),DRHS(L),L=1,IRHS)
             IRHS=0
          END IF 
       END IF 
       IF(RHS.EQ.0.0) RHS=1.0E-9  
       DQ(IGVEC/2+NVGB+ICGB)=RHS*HUGE      
      END DO

      IF(IRHS.NE.0) THEN
         WRITE(*,101) (IDRH(L),DRHS(L),L=1,IRHS)
      END IF
      IF(NOP.EQ.0) RETURN
      WRITE(*,102) ' Constraints: only equation values >',
     +               CLIMIT,' are printed'
 101  FORMAT('            ',4(I4,G11.3))
 102  FORMAT(A,G11.2,A)
      END 


      SUBROUTINE FEASMA ! matrix for feasible solution
#include "largeint.inc"
#include "dynal.inc"
      DOUBLE PRECISION RHS,SUM1,SUM2,SUM3
      SAVE
#include "mpinds.inc"
*     ...

      CALL MEGRM('f')
      CALL MEGARR('f D',2*NAGB,'rhs vector') ! double precision vector

      IF(IND8.EQ.0) RETURN  ! no constraints

      CALL MEGARR('k S',  NAGB,'auxiliary array')  ! single aux 1
      CALL MEGARR('l D',2*NAGB,'auxiliary array')  ! double aux 1

*     constraint matrix A and product A A^T (A is stored as transposed) 

      CALL MEGARR('y D',2*NCGB*NVGB,'transposed matrix of constraints')
      CALL MEGZER('f') ! reset rhs vector IGVEC

      I=0
      DO ICGB=1,NCGB
       RHS=QM(IND8+I+2) ! right hand side of constraint
       SGM=QM(IND8+I+4) ! sigma parameter
       I=I+4
 10    LABEL=MQ(IND8+I+1)
       FACTR=QM(IND8+I+2)
       ITGBI=INONE(LABEL) ! -> ITGBI= index of parameter label
       IVGB =JVGBI(ITGBI) ! -> variable-parameter index
c        IF(IVGB.LE.0) WRITE(*,*) 'Error in FEASMA1: LABEL ITGBI IVGB',
c     +                LABEL,ITGBI,IVGB

       IF(IVGB.GT.0) DQ(INDY/2+ICGB+(IVGB-1)*NCGB)=FACTR ! matrix element
       RHS=RHS-FACTR*DQ(IND0/2+ITGBI)     ! reduce residuum
       I=I+2
       IF(I.LT.LND8) THEN
          IF(MQ(IND8+I+1).NE.0) GOTO 10
       END IF
       DQ(IGVEC/2+NVGB+ICGB)=RHS        ! constraint discrepancy      
      END DO

      CALL MEGARR('x D',(NCGB*NCGB+NCGB),
     +            'product matrix of constraints') 

      DO L=1,NVGB
       IJ=0
       DO I=1,NCGB
        DO J=1,I 
         IJ=IJ+1  
         DQ(INDX/2+IJ)=DQ(INDX/2+IJ)+
     +                 DQ(INDY/2+(L-1)*NCGB+I)*DQ(INDY/2+(L-1)*NCGB+J)
        END DO
       END DO 
      END DO

*     inversion of product matrix of constraints      
      CALL SQMINV(DQ(INDX/2+1),
     +            DQ(IGVEC/2+NVGB+1),NCGB,NRANK,
     +            DQ(IDUX1/2+1),MQ(IAUX1+1))

      NMISS1=NCGB-NRANK
      
      WRITE(*,*) ' '
      WRITE(*,*) 'Rank of product matrix of constraints is',NRANK,
     +           ' for',NCGB,' constraint equations'
      WRITE(8,*) 'Rank of product matrix of constraints is',NRANK,
     +           ' for',NCGB,' constraint equations'
      IF(NRANK.LT.NCGB) THEN
         WRITE(*,*) 'Warning: insufficient constraint equations!' 
         WRITE(8,*) 'Warning: insufficient constraint equations!' 
      END IF

      CALL MEGRM('y')
      CALL MEGARR('y D',2*NCGB,'constraint correction factor')

      CALL MEGRM('k')
      CALL MEGRM('l')

      RETURN


      ENTRY FEASIB(CONCUT,IACT) ! make parameters feasible
*                        ! acts on subarray "0"
      IACT=0
      IF(IND8.EQ.0) RETURN  ! no constraints

      DO ITER=1,2
       CALL MEGZER('r') ! reset rhs vector IDUX7

*      calculate right constraint equation discrepancies vrom "0" ------

       I=0
       DO ICGB=1,NCGB
        RHS=QM(IND8+I+2) ! right hand side of constraint
        SGM=QM(IND8+I+4) ! sigma parameter
        I=I+4
 20     LABEL=MQ(IND8+I+1)
        FACTR=QM(IND8+I+2)
        ITGBI=INONE(LABEL) ! -> ITGBI= index of parameter label
        RHS=RHS-FACTR*DQ(IND0/2+ITGBI)     ! reduce residuum
        I=I+2
        IF(I.LT.LND8) THEN
           IF(MQ(IND8+I+1).NE.0) GOTO 20
        END IF
        DQ(IDUX7/2+NVGB+ICGB)=RHS        ! constraint discrepancy
       END DO

*      constraint equation discrepancies -------------------------------

       SUM1=0.0D0
       SUM2=0.0D0
       SUM3=0.0D0
       DO ICGB=1,NCGB
        SUM1=SUM1+DQ(IDUX7/2+NVGB+ICGB)**2 
        SUM2=SUM2+ABS(DQ(IDUX7/2+NVGB+ICGB))
        SUM3=MAX(SUM3,ABS(DQ(IDUX7/2+NVGB+ICGB)))
       END DO
       SUM1=SQRT(SUM1/DFLOAT(NCGB))
       SUM2=SUM2/DFLOAT(NCGB)

       IF(ITER.EQ.1.AND.SUM1.LT.CONCUT) RETURN  ! do nothing if correction small

       IF(ITER.EQ.1.AND.NCGB.LE.12) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Constraint equation discrepancies:'      
          WRITE(*,101) (ICGB,DQ(IDUX7/2+NVGB+ICGB),ICGB=1,NCGB)
 101      FORMAT(4X,4(I5,G12.4))
          WRITE(*,103) CONCUT
 103      FORMAT(10X,' Cut on rms value is',G8.1)
       END IF

       IF(IACT.EQ.0) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Improve constraints'
       END IF
       IACT=1               

       WRITE(*,102) ITER,SUM1,SUM2,SUM3
 102   FORMAT(I6,'   rms',G12.4,'  avrg_abs',G12.4,'  max_abs',G12.4)

*      multiply inverse matrix and constraint vector
       CALL DBSVX(DQ(INDX/2+1),DQ(IDUX7/2+NVGB+1),DQ(INDY/2+1),NCGB)
 
       CALL MEGZER('s') ! use "s" as aux vector 
      
       I=0
       DO ICGB=1,NCGB
        RHS=QM(IND8+I+2) ! right hand side of constraint
        SGM=QM(IND8+I+4) ! sigma parameter
        I=I+4
 30     LABEL=MQ(IND8+I+1)
        FACTR=QM(IND8+I+2)
        ITGBI=INONE(LABEL) ! -> ITGBI= index of parameter label
        IVGB =JVGBI(ITGBI) ! -> variable-parameter index
c        IF(IVGB.LE.0) WRITE(*,*) 'Error in FEASMA2: LABEL ITGBI IVGB',
c     +                LABEL,ITGBI,IVGB
        IF(IVGB.GT.0) THEN   
           DQ(ISOLV/2+IVGB)=DQ(ISOLV/2+IVGB)+DQ(INDY/2+ICGB)*FACTR
        END IF
        I=I+2
        IF(I.LT.LND8) THEN
           IF(MQ(IND8+I+1).NE.0) GOTO 30
        END IF
       END DO

       DO I=1,NVGB ! add corrections
        ITGBI=JTGBI(I)
        DQ(IND0/2+ITGBI)=DQ(IND0/2+ITGBI)+DQ(ISOLV/2+I)
       END DO 

      END DO ! iteration 1 and 2

      END 


      SUBROUTINE PEREAD(MORE)
      LOGICAL PRINT,FLOOP
      DOUBLE PRECISION DS0,DS1,DS2,DW
!$    INTEGER OMP_GET_THREAD_NUM
      SAVE 
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
*
      INDER(I)=MQ(IDOT3+NOFF+I)
      GLDER(I)=QM(IDOT4+NOFF+I)

      DATA PRINT/.TRUE./
      DATA FLOOP/.TRUE./
      DATA NPRI / 0 /, MPRI / 1000 /
*     ...
      IF(IFILE.EQ.0) THEN  ! start/restart
         NREC=0
         NTOT=0
         NRPR=1
         DO K=3,6
           MQ(IDOT3-K )=0 ! reset statistics
         ENDDO
         MQ(IDOT3-5)=LDOT3
         NTHR=MQ(IDOT3)
         DO K=0,NTHR*5-1
           MQ(IDOT2+LDOT2-K)=0  ! reset header (located at end)
         ENDDO
*
         NACT=0   ! active threads (have something still to read)
         DO K=0,NTHR-1
           IF (IFILE.LT.NFILB) THEN
             IFILE=IFILE+1
             MQ(IDOT2+LDOT2     -K)=IFILE 
             MQ(IDOT2+LDOT2-NTHR-K)=NACT
             NACT=NACT+1 
           ENDIF
         ENDDO
         MQ(IDOT4)=NACT ! active threads
      END IF 
      NACT=MQ(IDOT4)
      MQ(IDOT2-1)=MQ(IDOT2-2)/NACT
      MQ(IDOT3-1)=MQ(IDOT3-2)/NACT
      MORE=-1
      DO K=0,NTHR-1
         IACT=MQ(IDOT2+LDOT2-NTHR-K)
         MQ(IDOT2+LDOT2-NTHR*3-K)=0                 ! reset counter
         MQ(IDOT2+LDOT2-NTHR*4-K)=IACT*MQ(IDOT3-1)  ! reset offset
      ENDDO
      MQ(IDOT3-3)=MQ(IDOT3-3)+1 ! new block
*
!$OMP  PARALLEL
!$OMP& DEFAULT(PRIVATE)
!$OMP& SHARED(MQ,QM,IDOT2,IDOT3,IDOT4,LDOT2,NTHR,NFILF,FLOOP,
!$OMP&       IFD,KFD,IFILE,NFILB,WFD)  ! 'global' variables 
!$OMP& NUM_THREADS(NTHR)
*
      ITHR=0
!$    ITHR=OMP_GET_THREAD_NUM()         ! thread number 
      JFILE=MQ(IDOT2+LDOT2       -ITHR)  ! file index
      IACT =MQ(IDOT2+LDOT2-NTHR  -ITHR)  ! active thread number
      JREC =MQ(IDOT2+LDOT2-NTHR*2-ITHR)  ! records read
      IOFFP=IACT*MQ(IDOT2-1)
      IF (JFILE.LE.0) GOTO 50
      KFILE=KFD(2,JFILE)
*
 01   CONTINUE  ! read next record
         NBUF=MQ(IDOT2+LDOT2-NTHR*3-ITHR)+1
          NOFF=MQ(IDOT2+LDOT2-NTHR*4-ITHR)+2     ! 2 header words per record
         IF(KFILE.LE.NFILF) THEN ! Fortran file
            LUN=KFILE+10
            READ(LUN,END=10) N,(QM(IDOT4+NOFF+I),I=1,N/2),
     &                         (MQ(IDOT3+NOFF+I),I=1,N/2)
            NR=N/2 
         ELSE         ! C file
            LUN=KFILE-NFILF
            NR=NDIMB
#ifdef READ_C_FILES
            CALL READC(QM(IDOT4+NOFF+1),MQ(IDOT3+NOFF+1),NR,LUN,IERRC)
            N=NR+NR
#else
            IERRC=0
#endif
            IF(IERRC.LE.0) GOTO 11   ! end-of-files or error
         END IF 
         GOTO 02 ! read oK

 10      REWIND LUN                    ! Fortran eof
 11      MQ(IDOT2+LDOT2-ITHR)=-JFILE   ! flag eof
         IF (KFD(1,JFILE).EQ.0) THEN
            print *, 'PEREAD: file ', KFILE, 
     &        'read the first time, found',JREC,' records'
            KFD(1,JFILE)=-JREC
         ENDIF
*        take next file
!$OMP CRITICAL
         IF (IFILE.LT.NFILB) THEN
            IFILE=IFILE+1
            JREC=0
            MQ(IDOT2+LDOT2       -ITHR)=IFILE 
            MQ(IDOT2+LDOT2-NTHR*2-ITHR)=JREC
         ENDIF
!$OMP END CRITICAL
         JFILE=MQ(IDOT2+LDOT2-ITHR)
         IF (JFILE.GT.0) THEN
           KFILE=KFD(2,JFILE)
           GOTO 07 ! continue with new file
         ENDIF
         GOTO 50
*
 02      CONTINUE
         JREC=JREC+1
         MQ(IDOT2+LDOT2-NTHR*2-ITHR)=JREC
         MQ(IDOT2+LDOT2-NTHR*3-ITHR)=NBUF
         MQ(IDOT2+LDOT2-NTHR*4-ITHR)=NOFF+NR
*
         MQ(IDOT2+IOFFP+NBUF)=NOFF            ! pointer to start of buffer
         MQ(IDOT3+NOFF  )    =NOFF+NR         ! pointer to  end  of buffer
         MQ(IDOT3+NOFF-1)    =IFD(KFILE)+JREC ! global record number (available with LOOP2)
         MQ(IDOT4+NOFF  )    =KFILE           ! file number
         QM(IDOT4+NOFF-1)    =WFD(KFILE)      ! weight
*
         IF(FLOOP.AND.ITHR.EQ.0) THEN
            CALL HMPLNT(1,N)
            IF(INDER(1).NE.0) CALL HMPENT(8,FLOAT(INDER(1)))
         END IF 
         IF(N.GT.2*NDIMB) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'Histogram of record length'
            CALL HMPRNT(1) 
            WRITE(*,*) 'PEREAD record length ',N,', buffer is',2*NDIMB 
            STOP 'recompile with larger buffer size (NDIMB)'
         END IF
*
 07      CONTINUE
         IF ((NOFF+NR+2+NDIMB.LT.MQ(IDOT3-1)*(IACT+1)).AND.
     &       (NBUF.LT.MQ(IDOT2-1))) GOTO 01 ! buffer not yet full
*
 50   CONTINUE
!$OMP END PARALLEL
*     compress pointers
      NRD=MQ(IDOT2+LDOT2-3*NTHR) ! buffers from 1 .thread
      DO K=1,NTHR-1
         IACT =MQ(IDOT2+LDOT2-NTHR-K)
         IOFFP=IACT*MQ(IDOT2-1)
         NBUF=MQ(IDOT2+LDOT2-3*NTHR-K)
         DO L=1,NBUF
            MQ(IDOT2+NRD+L)=MQ(IDOT2+IOFFP+L)
         ENDDO
         NRD=NRD+NBUF
      ENDDO
*
      MORE=0
      DO K=0,NTHR-1
         JFILE=MQ(IDOT2+LDOT2-K)
         IF (JFILE.GT.0) THEN ! no eof yet
            MQ(IDOT2+LDOT2-NTHR-K)=MORE
            MORE=MORE+1
         ELSE
* no more files, thread retires
            MQ(IDOT2+LDOT2       -K)=0
            MQ(IDOT2+LDOT2-NTHR  -K)=-1
            MQ(IDOT2+LDOT2-NTHR*2-K)=0
         ENDIF
      ENDDO
*     record limit ?
      IF (MXREC.GT.0.AND.(NTOT+NRD).GE.MXREC) THEN
         NRD=MXREC-NTOT
         MORE=-1
         DO K=0,NTHR-1
            JFILE=MQ(IDOT2+LDOT2-K)
            IF (JFILE.GT.0) THEN   ! rewind files
               NRC=MQ(IDOT2+LDOT2-NTHR*2-K)
               IF (KFD(1,JFILE).EQ.0) KFD(1,JFILE)=-NRC
               IF (KFILE.LE.NFILF) THEN
                  LUN=KFILE+10
                  REWIND LUN
               ELSE
                  LUN=KFILE-NFILF
#ifdef READ_C_FILES
                  CALL RESETC(LUN)
#endif
               ENDIF
            ENDIF
         ENDDO
      ENDIF 
*
      NTOT=NTOT+NRD
      NREC=NTOT
      MQ(IDOT2)=NRD
      MQ(IDOT4)=MORE  ! number of active threads
*
      MQ(IDOT3-4)=MQ(IDOT3-4)+NRD
      MQ(IDOT3-5)=MIN(MQ(IDOT3-5),NRD)
      MQ(IDOT3-6)=MAX(MQ(IDOT3-6),NRD)
*
 55   CONTINUE
      IF (NLOOPN.EQ.0.AND.NTOT.GE.NRPR) THEN
         WRITE(*,*) ' Record ',NRPR
         IF (NRPR.LT.100000) THEN
            NRPR=NRPR*10
         ELSE
            NRPR=NRPR+100000
         ENDIF
         GOTO 55
      ENDIF
      IF (NCACHE.GT.0.AND.NLOOPN.LE.1.AND.
     +    NPRI.LT.MPRI.AND.MPRINT.GT.1) THEN
         NPRI=NPRI+1
         IF (NPRI.EQ.1) WRITE(*,100)
         WRITE(*,101) NREC, NRD, MORE ,IFILE
  100    FORMAT(/' PeRead        records         active      file'
     +          /'            total     block   threads    number')
  101    FORMAT(' PeRead',4I10)
      ENDIF
*
      IF (MORE.LE.0) THEN
         IFILE=0
         IF (FLOOP) THEN
*           check for file weights
            DS0=0.0D0
            DS1=0.0D0
            DS2=0.0D0
            DO K=1,NFILB
               DW=DFLOAT(-KFD(1,K))
               IF (WFD(K).NE.1.0) NFILW=NFILW+1
               DS0=DS0+DW
               DS1=DS1+DW*DBLE(WFD(K))
               DS2=DS2+DW*DBLE(WFD(K)**2)
            ENDDO
            IF (NFILW.GT.0.AND.DS0.GT.0.0D0) THEN
               DS1=DS1/DS0
               DS2=DS2/DS0-DS1*DS1
               DO LUN=6,LUNLOG,2
                  WRITE(LUN,1177) NFILW,SNGL(DS1),SNGL(DS2)
 1177      FORMAT(/' !!!!!',I4,' weighted binary files',
     &            /' !!!!! mean, variance of weights =',2G12.4)
               ENDDO
            ENDIF
*           integrate record numbers
            DO K=2,NFILB 
               IFD(K)=IFD(K-1)-KFD(1,K-1)
            ENDDO
*           sort
            IF (NTHR.GT.1) CALL SORT2K(KFD,NFILB)
         ENDIF
         PRINT=.FALSE.
         FLOOP=.FALSE.
         IF (NCACHE.GT.0.AND.NLOOPN.LE.1.AND.MPRINT.GT.0)
     &     WRITE(*,1179) (MQ(IDOT3-K),K=3,6)
 1179      FORMAT(/' Read  cache usage (#blocks, #records, ',
     &            'min,max records/block'/17X,4I10)
      ENDIF
      END

      SUBROUTINE PEPREP(MODE)
*     prepare records
*     MODE <=0: build index table (INONE) for global variables
*     MODE  >0: use index table, can be parallelized
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
*
      ISFRST(IBUF)=MQ(IDOT2+IBUF)+1
      ISLAST(IBUF)=MQ(IDOT3+MQ(IDOT2+IBUF))
      INDER(I)=MQ(IDOT3+NOFF+I)
      GLDER(I)=QM(IDOT4+NOFF+I)
*
      IF (MODE.GT.0) THEN
        ICHUNK=MIN((MQ(IDOT3)+MTHRD-1)/MTHRD/8+1,256)
!$OMP  PARALLEL DO                       ! parallelize record loop
!$OMP&  DEFAULT(PRIVATE) 
!$OMP&  SHARED(MQ,QM,DQ,IDOT2,IDOT3,IDOT4,ICHUNK)
!$OMP&  SCHEDULE(DYNAMIC,ICHUNK)
        DO IBUF=1,MQ(IDOT2)           ! buffer for current record
          IPROC=0
!$        IPROC=OMP_GET_THREAD_NUM()         ! thread number
          IST=ISFRST(IBUF)
          NST=ISLAST(IBUF)
 10       CALL ISJAJB(NST,IST,JA,JB,JSP)
          IF(JB.NE.0) THEN
            DO J=1,IST-JB
              MQ(IDOT3+JB+J)=INONE( MQ(IDOT3+JB+J) ) ! translate to index
            END DO
            GOTO 10
          ENDIF
        ENDDO
!$OMP  END PARALLEL DO
      ENDIF

!$POMP INST BEGIN(peprep)
      IF (MODE.LE.0) THEN
        DO IBUF=1,MQ(IDOT2)           ! buffer for current record
          IST=ISFRST(IBUF)
          NST=ISLAST(IBUF)
 20       CALL ISJAJB(NST,IST,JA,JB,JSP)
          IF(JB.NE.0) THEN
            DO J=1,IST-JB
              MQ(IDOT3+JB+J)=INONE( MQ(IDOT3+JB+J) ) ! translate to index
            END DO
            GOTO 20
          ENDIF
        ENDDO
      ENDIF
!$POMP INST END(peprep)
*
      END

      SUBROUTINE ISJAJB(NST,IS,JA,JB,JSP)
*
*     Get indices JA, JB, IS for next measurement within record
*     NST = index of last  word of record 
*     IS  = index of first word of record at the first call of ISJAJB
*     ISJAJB returns:  JA  = index of measured value       (=0 at end)
*                          = pointer to local derivatives
*                      JB  = index of standard deviation   (=0 at end)
*                          = pointer to global derivatives
*                      IS  = index of last global derivative
*                      JSP = index to special data
*     Local derivatives are:
*        INDER(JA+J),GLDER(JA+J),J=1,JB-JA-1)  i.e. JB-JA-1 derivatives
*     Global derivatives are:
*        (INDER(JB+J),GLDER(JB+J),J=1,IS-JB)  i.e. IST-JB derivatives
*
*     End_of_data is indicated by returned values JA=0 and JB=0
*     Special data are ignored. At end_of_data the info to the
*     special data is returned: IS = pointer to special data;
*        number of words is NSP=-GLDER(IS)
*
#include "largeint.inc"
#include "dynal.inc"
      SAVE
*     ...
      INDER(I)=MQ(IDOT3+I)
      GLDER(I)=QM(IDOT4+I)
*
      JSP=0
 01   JA=0
      JB=0
      IF(IS.GE.NST) RETURN
 10   IS=IS+1
      IF(INDER(IS).NE.0) GOTO 10
      JA=IS
 20   IS=IS+1
      IF(INDER(IS).NE.0) GOTO 20
      JB=IS
 30   IF(INDER(IS+1).NE.0.AND.IS.LT.NST) THEN
         IS=IS+1
         GOTO 30
      END IF
      IF(JA+1.EQ.JB.AND.GLDER(JB).LT.0.0) THEN
*        special data
         JSP=JB          ! pointer to special data
         IS=IS+IFIX(-GLDER(JB)+0.5) ! skip NSP words
         GOTO 01
      END IF
      END



************************************************************************
*                        LOOPN ...

      SUBROUTINE LOOPN                 ! loop with fits and sums 
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      DOUBLE PRECISION ADDER,FUNREF,DCHI2S,SNDF
      INTEGER INIRIS
      COMMON/CIRIS/INIRIS
      SAVE
*     ...
      ISFRST(IBUF)=MQ(IDOT2+IBUF)+1
      ISLAST(IBUF)=MQ(IDOT3+MQ(IDOT2+IBUF))
      INDER(I)=MQ(IDOT3+I)
      GLDER(I)=QM(IDOT4+I)
*     ----- book and reset ---------------------------------------------   
      IF(NLOOPN.EQ.0) THEN      ! first call
         LASTIT=-1
         IITERA=0
      END IF

      INIRIS=0 !!!!!!!!!!

      NLOOPN=NLOOPN+1           ! increase loop counter
      NDFSUM=0
      SUMNDF=0.0D0

      IF(NLOOPN.EQ.1) THEN      ! book histograms for 1. iteration
         CALL GMPDEF(1,4,'Function value in iterations') 
         IF (METSOL.EQ.3) THEN ! extend to GMRES, i.e. 4?
            CALL GMPDEF(2,3,'Number of MINRES steps vs iteration nr')
         END IF
         CALL HMPDEF( 5,0.0,0.0,'Number of degrees of freedom')
         CALL HMPDEF(11,0.0,0.0,'Number of local parameters')
         CALL HMPDEF(23,0.0,0.0,
     +   'SQRT of diagonal elements without presigma')
         CALL HMPDEF(24,0.0,0.0,
     +   'Log10 of off-diagonal elements')
         CALL HMPDEF(25,0.0,0.0,
     +   'Relative individual pre-sigma')
         CALL HMPDEF(26,0.0,0.0,
     +   'Relative global pre-sigma')
      END IF


      CALL HMPDEF(3,0.0,0.0,          ! book  
     +   'Normalized residuals of single (global) measurement')
      CALL HMPDEF(12,0.0,0.0,          ! book  
     +   'Normalized residuals of single (local) measurement')
      CALL HMPDEF(13,0.0,0.0,          ! book  
     +   'Pulls of single (global) measurement')
      CALL HMPDEF(14,0.0,0.0,          ! book  
     +   'Pulls of single (local) measurement')
      CALL HMPDEF(4,0.0,0.0,'Chi^2/Ndf after local fit')
      CALL GMPDEF(4,5,'location, dispersion (res.) vs record nr') 
      CALL GMPDEF(5,5,'location, dispersion (pull) vs record nr') 

C      WRITE(*,*) 'LOOPN ', NLOOPN, ' executing ICALCM=', ICALCM 

*     reset

      CALL MEGZER('f') ! reset rhs vector              IGVEC
      IF(ICALCM.EQ.1) THEN
         CALL MEGZER('j') ! reset matrix (double)      IGMAT
         CALL MEGZER('.') ! reset matrix (single)
         CALL MEGZER('v') ! reset variable-band matrix INDV
      END IF 
 
      IF(NLOOPN.EQ.2) CALL HMPDEF(6,0.0,0.0,'Down-weight fraction')

      NEWITE=.FALSE.
      IF(ITERAT.NE.LASTIT) THEN ! new iteration
         NEWITE=.TRUE.
         FUNREF=FVALUE
         IF(NLOOPN.GT.1) THEN
            NREJ=NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3)
c            CALL MEND
            IF(ITERAT.EQ.1) THEN
               CHICUT=CHIREM
            ELSE IF(ITERAT.GE.1) THEN 
               CHICUT=SQRT(CHICUT)
               IF(CHICUT.NE.0.0.AND.CHICUT.LT.1.5) CHICUT=1.0
               IF(CHICUT.NE.0.0.AND.NREJ.EQ.0)     CHICUT=1.0
            END IF
         END IF
c         WRITE(*,111) ! header line 
      END IF 

      DO I=0,3
       NREJEC(I)=0   ! reset reject counter
      END DO
      DO K=3,6
       MQ(IND5-K)=0  ! cache usage
       MQ(INDH-K)=0
      ENDDO
*     ----- read next data ----------------------------------------------
 10   CALL PEREAD(NR)  ! read records
      CALL PEPREP(1)   ! prepare records
      NDFS  =0
      SNDF  =0.0D0
      DCHI2S=0.0D0
      CALL LOOPBF(NREJEC,
     &            NAGBN*NALCN,         DQ(IGLMA/2+1),      ! 'g'
     &            NAEQN,               QM(INDW+1),         ! 'w'
     &            NDFS,SNDF,DCHI2S) 
      NDFSUM=NDFSUM+NDFS
      SUMNDF=SUMNDF+SNDF
      CALL ADDSUM(DCHI2S)
      IF (NR.GT.0) GOTO 10 ! -> read next block of events 
*     sum up RHS (over threads)
      IOFFB=0
      DO IPR=2,MTHRD
        IOFFB=IOFFB+MQ(INDF-1)
        DO K=1,MQ(INDF-1)
          DQ(INDF/2+K)=DQ(INDF/2+K)+DQ(INDF/2+IOFFB+K)
        ENDDO
      ENDDO
*
      IF (ICALCM.EQ.1) THEN
C        print *, ' cache/w ',(MQ(INDH-K),K=3,6),(MQ(IND5-K),K=3,6)
        NPARL=MQ(IND5-3)
        NCRIT=MQ(IND5-4)
        USED=FLOAT(MQ(INDH-5))/FLOAT(MQ(INDH-3))*0.1
        USEI=FLOAT(MQ(IND5-5))/FLOAT(MQ(IND5-3))*0.1
        PEAKD=FLOAT(MQ(INDH-6))*0.1
        PEAKI=FLOAT(MQ(IND5-6))*0.1
        WRITE(*,111) NPARL,NCRIT,USEI,USED,PEAKI,PEAKD
  111   FORMAT(' Write cache usage (#flush,#overrun,<levels>,',
     &         'peak(levels))'/2I7,',',4(F6.1,'%'))
      ENDIF
*     ----- after end-of-data add contributions from pre-sigma ---------

      IF(NLOOPN.EQ.1) THEN
*        plot diagonal elements
         ELMT=0.0
         DO I=1,NVGB        ! diagonal elements
          II=0
          IF(MATSTO.EQ.1) II=(I*I+I)/2
          IF(MATSTO.EQ.2) II=I
          IF(MATSTO.EQ.3) II=I
          IF(II.NE.0) THEN
             ELMT=SNGL(DQ(IGMAT/2+II))
             IF(ELMT.GT.0.0) CALL HMPENT(23,1.0/SQRT(ELMT))
          END IF
         END DO 
      END IF  



*     add pre-sigma contributions to matrix diagonal     

c      WRITE(*,*) 'Adding to diagonal ICALCM IND6',ICALCM,IND6

      IF(ICALCM.EQ.1) THEN
         DO IVGB=1,NVGB                        ! add evtl. pre-sigma
c          WRITE(*,*) 'Index ',IVGB,IVGB,QM(IND6+IVGB)
          IF(QM(IND6+IVGB).NE.0.0) THEN 
             IF(IVGB.GT.0) CALL MUPDAT(IVGB,IVGB,DBLE(QM(IND6+IVGB)))
          END IF
         END DO
      END IF

      CALL HMPWRT(23)
      CALL HMPWRT(24) 
      CALL HMPWRT(25)
      CALL HMPWRT(26) 


*     add regularization term to F and to rhs --------------------------
 
c      WRITE(*,*) 'NREGUL ',NREGUL,NLOOPN

      IF(NREGUL.NE.0) THEN ! add regularization term to F and to rhs
         DO IVGB=1,NVGB
          ITGBI=JTGBI(IVGB) ! global parameter index
          DQ(IGVEC/2+IVGB)=DQ(IGVEC/2+IVGB)
     +                    -DQ(IND0/2+ITGBI)*QM(IND6+IVGB)
          ADDER=QM(IND6+IVGB)*DQ(IND0/2+ITGBI)**2
          CALL ADDSUM(ADDER)
         END DO
      END IF


*     ----- add contributions from "measurement" -----------------------

      IF(IND9.EQ.0) GOTO 60

      I=0
 40   RHS=QM(IND9+I+2) ! right hand side
      SGM=QM(IND9+I+4) ! sigma parameter
      I=I+4
      WEIGHT=0.0
      IF(SGM.GT.0.0) WEIGHT=1.0/SGM**2

      DSUM=-RHS 

*     loop over label/factor pairs
      IA=I
 50   LABEL=MQ(IND9+I+1)
      FACTR=QM(IND9+I+2)
      ITGBI=INONE(LABEL)
      IF(ITGBI.NE.0) THEN   
         DSUM=DSUM+FACTR*SNGL(DQ(IND0/2+ITGBI))     ! residuum
      END IF
      I=I+2
      IF(I.LT.LND9) THEN
         IF(MQ(IND9+I+1).NE.0) GOTO 50
      END IF
      IB=I-2     
        
      ADDER=WEIGHT*DSUM**2
      CALL ADDSUM(ADDER)
      FVALUE=FVALUE+WEIGHT*DSUM**2 ! accumulate chi-square

      DO J=IA,IB,2
       LABELJ=MQ(IND9+J+1)
       FACTRJ=QM(IND9+J+2)
       ITGBIJ=INONE(LABELJ)
*      add to vector
       IVGBIJ=0
       IF(ITGBIJ.NE.0) IVGBIJ=JVGBI(ITGBIJ) ! variable-parameter index
       IF(IVGBIJ.GT.0) THEN
           DQ(IGVEC/2+IVGBIJ)=DQ(IGVEC/2+IVGBIJ)
     +                      -WEIGHT*DSUM*FACTRJ ! vector
       END IF 

       IF(ICALCM.EQ.1.AND.IVGBIJ.GT.0) THEN
          DO K=IA,J,2
           LABELK=MQ(IND9+K+1)
           FACTRK=QM(IND9+K+2)
           ITGBIK=INONE(LABELK)
*          add to matrix
           IVGBIK=0
           IF(ITGBIK.NE.0) IVGBIK=JVGBI(ITGBIK) ! variable-parameter index
           IF(IVGBIJ.GT.0.AND.IVGBIK.GT.0) THEN    ! 
              CALL MUPDAT(IVGBIJ,IVGBIK,DBLE(WEIGHT*FACTRJ*FACTRK)) 
           END IF
          END DO
       END IF 
      END DO
      IF(I.LT.LND9) GOTO 40        

*     ----- printout ---------------------------------------------------

 60   CONTINUE


      CALL GETSUM(FVALUE)   ! get accurate sum (Chi^2)
      FLINES=0.5D0*FVALUE   ! Likelihood function value
      RLOOP=ITERAT+0.01*NLOOPN
      ACTFUN=SNGL(FUNREF-FVALUE)
      IF(NLOOPN.EQ.1) ACTFUN=0.0
      NGRAS=NINT(ANGRAS)
      RATAE=0.0                              !!!
      IF(DELFUN.NE.0.0) THEN
         RATAE=MIN(99.9,ACTFUN/DELFUN)  !!!
         RATAE=MAX(-99.9,RATAE)
      END IF 
      NSECND=NINT(SECND)

*     rejects ...

      NREJ  =NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3)
      IF(NLOOPN.EQ.1) THEN
         IF(NREJ.NE.0) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'Data rejected in initial loop:'
            WRITE(*,*) '   ',
     +      NREJEC(0), ' (rank deficit/NaN) ',NREJEC(1),' (Ndf=0)   ',
     +      NREJEC(2), ' (huge)   ',NREJEC(3),' (large)'           
         END IF 
      END IF
c      IF(NREJEC(1)+NREJEC(2)+NREJEC(3).NE.0) THEN
c         WRITE(LUNLOG,*) 'Data rejected in initial loop:',NREJEC(1),
c     +   ' (Ndf=0)   ',NREJEC(2),' (huge)   ',NREJEC(3),' (large)'
c      END IF
     

      IF(NEWITE.AND.ITERAT.EQ.2) THEN 
         IF(NRECPR.NE.0.OR.NRECP2.NE.0) NRECER=NREC3
         IF(NRECPR.LT.0) THEN
            NRECPR=NREC1
         END IF
         IF(NRECP2.LT.0) THEN
            NRECP2=NREC2
         END IF
      END IF

      IF(NLOOPN.LE.2) THEN
         IF(NHISTP.NE.0) THEN 
c            CALL HMPRNT(3)  ! scaled residual of single measurement
c            CALL HMPRNT(12) ! scaled residual of single measurement
c            CALL HMPRNT(4)  ! chi^2/Ndf 
         END IF 
         CALL HMPWRT(3)
         CALL HMPWRT(12)
         CALL HMPWRT(4)
         CALL GMPWRT(4) ! location, dispersion (res.) as a function of record nr
         IF (NLOOPN.LE.LFITNP) THEN
            CALL HMPWRT(13)
            CALL HMPWRT(14)
            CALL GMPWRT(5) ! location, dispersion (pull) as a function of record nr
         ENDIF
      END IF
c      IF(NLOOPN.EQ.2.AND.NHISTP.NE.0) CALL HMPRNT(6)
      IF(NLOOPN.EQ.2) CALL HMPWRT(6) 
      IF(NLOOPN.LE.1) THEN
c         IF(NHISTP.NE.0) CALL HMPRNT(5)  ! number of degrees of freedom  
c         IF(NHISTP.NE.0) CALL HMPRNT(11) ! Nlocal  
         CALL HMPWRT(5)
         CALL HMPWRT(11)
      END IF

*     local fit: band matrix structure !?
      IF (NLOOPN.EQ.1.AND.NBNDR.GT.0) THEN
        DO LUN=6,8,2
           WRITE(LUN,*) ' '
           WRITE(LUN,*) 
     +         ' === local fits have bordered band matrix structure ==='
            WRITE(LUN,101) ' NBNDR',NBNDR,'number of records'
            WRITE(LUN,101) ' NBDRX',NBDRX,'max border size'
            WRITE(LUN,101) ' NBNDX',NBNDX,'max band width'
         ENDDO
      ENDIF
*
      LASTIT=ITERAT

 101  FORMAT(1X,A6,' =',I10,' = ',A)
c 101  FORMAT(' LOOPN',I6,' Function value',F22.8,10X,I6,' records') 
c 102  FORMAT('   incl. constraint penalty',F22.8) 
c 103  FORMAT(I13,3X,A,G12.4)  
      END 


      SUBROUTINE PLOOPA(LUNP)               ! title for iteration
#include "largeint.inc"
#include "dynal.inc"                   
#include "mpinds.inc"
      REAL SLOPES(3),STEPS(3)
      CHARACTER*4 CCALCM(4)
      DATA CCALCM / ' end','   S', ' F  ',' FMS' /
      SAVE
*     .. 
      WRITE(LUNP,*)   ' '
      WRITE(LUNP,101) ! header line
      WRITE(LUNP,102) ! header line
      RETURN

      ENTRY PLOOPB(LUNP)                    ! iteration line 
      NREJ=NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3)    ! rejects
      IF(NREJ.GT.9999999) NREJ=9999999
      CALL ETIME(TA,RSTB)
      DELTIM=RSTB-RSTART
      CALL SECHMS(DELTIM,NHOUR,MINUT,SECND) ! time
      NSECND=NINT(SECND)
      IF(ITERAT.EQ.0) THEN
         WRITE(LUNP,103) ITERAT,NLOOPN,FVALUE,
     +                   CHICUT,NREJ,NHOUR,MINUT,NSECND,CCALCM(LCALCM)
      ELSE
         CALL PTLOPT(NFA,MA,SLOPES,STEPS)  ! slopes steps
         RATAE=ABS(SLOPES(2)/SLOPES(1)) 
         STEPL=STEPS(2)
         WRITE(LUNP,104) ITERAT,NLOOPN,FVALUE,DELFUN,RATAE,ANGRAS,
     +                   IITERA,ISTOPA,LSINFO,STEPL,
     +                   CHICUT,NREJ,NHOUR,MINUT,NSECND,CCALCM(LCALCM)
      END IF 
      RETURN

      ENTRY PLOOPC(LUNP)                    ! sub-iteration line
      NREJ=NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3)    ! rejects
      IF(NREJ.GT.9999999) NREJ=9999999
      CALL ETIME(TA,RSTB)
      DELTIM=RSTB-RSTART
      CALL SECHMS(DELTIM,NHOUR,MINUT,SECND) ! time
      NSECND=NINT(SECND)
      CALL PTLOPT(NFA,MA,SLOPES,STEPS)  ! slopes steps
      RATAE=ABS(SLOPES(2)/SLOPES(1))
      STEPL=STEPS(2)
      WRITE(LUNP,105) NLOOPN,FVALUE, RATAE,LSINFO,
     +                STEPL,NREJ,NHOUR,MINUT,NSECND,CCALCM(LCALCM)
      RETURN

      ENTRY PLOOPD(LUNP)                     ! solution line (ICALCM=-1)
      CALL ETIME(TA,RSTB)
      DELTIM=RSTB-RSTART
      CALL SECHMS(DELTIM,NHOUR,MINUT,SECND) ! time
      NSECND=IFIX(SECND+0.5)

      WRITE(LUNP,106) NHOUR,MINUT,NSECND,CCALCM(LCALCM)
     
 101  FORMAT(' it fc','   fcn_value dfcn_exp  slpr costh  iit st',
     +       ' ls  step cutf',1X,'rejects hmmsec FMS')
 102  FORMAT(' -- --',' ----------- --------  ---- -----  --- --',
     +       ' --  ---- ----',1X,'------- ------ ---')

 103  FORMAT(I3,I3,E12.5,38X,F5.1,
     +      1X,I7,  I2,I2,I3,A4)
 104  FORMAT(I3,I3,E12.5,1X,E8.2,F6.3,F6.3,I5,2I3,F6.2,F5.1,
     +      1X,I7,  I2,I2,I3,A4)
 105  FORMAT(3X,I3,E12.5,9X,     F6.3,14X,I3,F6.2,6X,
     +         I7,  I2,I2,I3,A4)
 106  FORMAT(69X,I2,I2,I3,A4)
      END 



      SUBROUTINE EXPLFC(LUNIT)         ! print table headings
       WRITE(LUNIT,*) ' '
       WRITE(LUNIT,102)
     + 'Explanation of iteration table'
       WRITE(LUNIT,102)
     + '=============================='
       WRITE(LUNIT,101) 'it',
     + 'iteration number. Global parameters are improved for it > 0.'
       WRITE(LUNIT,102)
     + 'First function evaluation is called iteraton 0.'
       WRITE(LUNIT,101) 'fc',
     + 'number of function evaluations.'
       WRITE(LUNIT,101) 'fcn_value',
     + 'value of 2 x Likelihood function (LF).'
       WRITE(LUNIT,102)
     + 'The final value is the chi^2 value of the fit and should'
       WRITE(LUNIT,102)
     + 'be about equal to the NDF (see below).'
       WRITE(LUNIT,101) 'dfcn_exp',
     + 'expected reduction of the value of the Likelihood function (LF)'
       WRITE(LUNIT,101) 'slpr',
     + 'ratio of the actual slope to inital slope.'
        WRITE(LUNIT,101) 'costh',
     + 'cosine of angle between search direction and -gradient'
       WRITE(LUNIT,101) 'iit',
     + 'number of internal iterations in GMRES/MINRES algorithmus'
       WRITE(LUNIT,101) 'st',
     + 'stop code of GMRES/MINRES algorithmus'
       WRITE(LUNIT,102)
     + '< 0:   rhs is very special, with beta2 = 0'
       WRITE(LUNIT,102)
     + '= 0:   rhs b = 0, i.e. the exact solution is  x = 0'
       WRITE(LUNIT,102)
     + '= 1    requested accuracy achieved, as determined by rtol'
       WRITE(LUNIT,102)
     + '= 2    reasonable accuracy achieved, given eps'
       WRITE(LUNIT,102)
     + '= 3    x has converged to an eigenvector'
       WRITE(LUNIT,102)
     + '= 4    matrix ill-conditioned (Acond has exceeded 0.1/eps)'
       WRITE(LUNIT,102)
     + '= 5    the iteration limit was reached'
       WRITE(LUNIT,102)
     + '= 6    Matrix x vector does not define a symmetric matrix'
       WRITE(LUNIT,102)
     + '= 7    Preconditioner does not define a symmetric matrix'
       WRITE(LUNIT,101) 'ls',
     + 'line search info'
       WRITE(LUNIT,102)
     + '< 0    recalculate function'
       WRITE(LUNIT,102)
     + '= 0:   N or STP lt 0 or step not descending'
       WRITE(LUNIT,102)
     + '= 1:   Linesearch convergence conditions reached'
       WRITE(LUNIT,102)
     + '= 2:   interval of uncertainty at lower limit'
       WRITE(LUNIT,102)
     + '= 3:   max nr of line search calls reached'
       WRITE(LUNIT,102)
     + '= 4:   step at the lower bound'
       WRITE(LUNIT,102)
     + '= 5:   step at the upper bound'
       WRITE(LUNIT,102)
     + '= 6:   rounding error limitation'
       WRITE(LUNIT,101) 'step',
     + 'the factor for the Newton step during the line search. Usually'
       WRITE(LUNIT,102)
     + 'a value of 1 gives a sufficient reduction of the LF. Oherwise'
       WRITE(LUNIT,102)
     + 'other step values are tried.'
       WRITE(LUNIT,101) 'cutf',
     + 'cut factor. Local fits are rejected, if their chi^2 value'
       WRITE(LUNIT,102)
     + 'is larger than the 3-sigma chi^2 value times the cut factor.'
       WRITE(LUNIT,102)
     + 'A cut factor of 1 is used finally, but initially a larger'
       WRITE(LUNIT,102)
     + 'factor may be used. A value of 0.0 means no cut.'
       WRITE(LUNIT,101) 'rejects',
     + 'total number of rejected local fits.'
       WRITE(LUNIT,101) 'hmmsec',
     + 'the time in hours (h), minutes (mm) and seconds.'
       WRITE(LUNIT,*) ' '

 101  FORMAT(A9,' =  ',A)
 102  FORMAT(13X,A)
      END

      SUBROUTINE MUPDAT(I,J,ADD)       ! add ADD to matrix elm(I,J)
      DOUBLE PRECISION ADD ! double precision !
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      INTEGER(KIND=LARGE) IJADD,IA,JA,IJ
*     ...
      IF(I.LE.0.OR.J.LE.0) RETURN
      IA=MAX(I,J)          ! larger
      JA=MIN(I,J)          ! smaller
      IJ=0
      IF(MATSTO.EQ.1) THEN                  ! full symmetric matrix
         IJ=JA+(IA*IA-IA)/2                 ! ISYM index
         DQ(IGMAT/2+IJ)=DQ(IGMAT/2+IJ)+ADD
      ELSE IF(MATSTO.EQ.2) THEN             ! sparse symmetric matrix
         IJ=IJADD(I,J) !  inline code requires same time
         IF (IJ.EQ.0) RETURN                ! pair is suppressed
         IF (IJ.GT.0) THEN
            DQ(IGMAT/2+IJ)=DQ(IGMAT/2+IJ)+ADD
         ELSE
            QM(IDOT1-IJ)=QM(IDOT1-IJ)+SNGL(ADD)
         ENDIF
      END IF
      IF(METSOL.GE.3) THEN
         IF(MBANDW.GT.0) THEN     ! for Cholesky decomposition
            IF(IA.LE.NVGB) THEN   ! variable global parameter
               IJ=MQ(INDU+IA)-IA+JA
               IF(IA.GT.1.AND.IJ.LE.MQ(INDU+IA-1)) IJ=0
               IF(IJ.NE.0) DQ(INDV/2+IJ)=DQ(INDV/2+IJ)+ADD
               IF(IJ.LT.0.OR.IJ.GT.LNDV/2) STOP 'MUPDAT: bad index'
            ELSE                  ! Lagrange multiplier
               IJ=MQ(INDU+NVGB)+(IA-NVGB-1)*NVGB+JA
               IF(IJ.NE.0) DQ(INDV/2+IJ)=DQ(INDV/2+IJ)+ADD
            END IF             
         ELSE IF(MBANDW.EQ.0) THEN      ! default preconditioner 
            IF(IA.LE.NVGB) THEN   ! variable global parameter
               IF(JA.EQ.IA) DQ(INDV/2+IA)=DQ(INDV/2+IA)+ADD ! diag
            ELSE                  ! Lagrange multiplier
               IJ=NVGB+(IA-NVGB-1)*NVGB+JA
               IF(IJ.NE.0) DQ(INDV/2+IJ)=DQ(INDV/2+IJ)+ADD
            END IF
         END IF
      END IF
      END 

*
*      loop over events in buffer, fits and sums
*
      SUBROUTINE LOOPBF(NREJ,NLG,VLG,NR,RV,NDFS,SNDF,DCHI2S)
*     contains  
*     SUBROUTINE FITLOC(NALC,NALG) ! number of local, global parameters
*
*     090817 C. Kleinwort, DESY-FH1
*     for local fits with bordered band matrices use band matrix algebra 
*     (implements Broken Line Track fit:
*      A new fast track-fit algorithm based on broken lines, V. Blobel, 
*      Nuclear Instruments and Methods A, 566 (October 2006), pp. 14-17)
*
#include "largeint.inc"
#include "dynal.inc"
*     
*     arrays for local fit
#include "localfit.inc"     
      PARAMETER (MLMSYM=(MLOCAL*MLOCAL+MLOCAL)/2)
      DOUBLE PRECISION BLVEC(MLOCAL),CLMAT(MLMSYM),SCDIAG(MLOCAL)
      DOUBLE PRECISION DCHI2,DVAR,DCHI2S,SNDF,VLG(NLG),DW1,DW2,SUMM
      DIMENSION RV(NR),NREJ(0:3)
      INTEGER SCFLAG(MLOCAL),IBANDH(MLOCAL)
*      scratch arrays for local fit
      DOUBLE PRECISION VBND(MLOCAL*(MXBND+1)),VBDR(MLOCAL*MXBDR),
     +   AUX(MLOCAL*MXBDR),VBK((MXBDR*MXBDR+MXBDR)/2),VZRU(MXBDR)
!$    INTEGER OMP_GET_THREAD_NUM
*       
      LOGICAL LPRNT,LHIST
      CHARACTER*3 CHAST
      DATA CHUBER/1.345/  ! constant for Huber down-weighting
      DATA CAUCHY/2.3849/ ! constant for Cauchy down-weighting 
      SAVE CHUBER,CAUCHY
#include "mpinds.inc"
*     ... 
      ISFRST(IBUF)=MQ(IDOT2+IBUF)+1
      ISLAST(IBUF)=MQ(IDOT3+MQ(IDOT2+IBUF))
      INDER(I)=MQ(IDOT3+I)
      GLDER(I)=QM(IDOT4+I)
*
      ICHUNK=MIN((MQ(IDOT3)+MTHRD-1)/MTHRD/8+1,256)
      DO K=1,MTHRD*9         ! reset header, 3 words per thread:
        MQ(IND5+K)=0         !    number of entries, offset to data, indices
      ENDDO
C      IF (ICALCM.EQ.1) print *, ' ICHUNK ', MQ(IDOT3),MTHRD,ICHUNK
      NPRDBG=0
      IPRDBG=-1
      
*
!$OMP  PARALLEL DO                       ! parallelize record loop
!$OMP&  DEFAULT(PRIVATE) 
!$OMP&  SHARED(MQ,QM,DQ,IND0,INDA,INDC,INDE,INDF,IDOT2,IDOT3,IDOT4,NLG,
!$OMP&     NAGBN,ICALCM,ICHUNK,NLOOPN,NRECER,NPRDBG,IPRDBG,
!$OMP&     NEWITE,CHICUT,LHUBER,CHUBER,ITERAT,NRECPR,IND5,INDH,MTHRD,
!$OMP&     DWCUT,CHHUGE,NRECP2,CAUCHY,LFITNP,LFITBB)  ! 'global' variables 
!$OMP&  REDUCTION(+:NDFS,SNDF,DCHI2S,NREJ,NBNDR) ! private copy of NDFS,.. for each thread, combined at end, init with 0.
!$OMP&  REDUCTION(MAX:NBNDX,NBDRX)
!$OMP&  REDUCTION(MIN:NREC3)
!$OMP&  SCHEDULE(DYNAMIC,ICHUNK)
      DO IBUF=1,MQ(IDOT2)           ! buffer for current record
       NRC=MQ(IDOT3+ISFRST(IBUF)-2)   ! record 
       WRC=QM(IDOT4+ISFRST(IBUF)-2)   ! weight 
       DW1=DBLE(WRC)
       DW2=DSQRT(DW1)
*
       IPROC=0
!$       IPROC=OMP_GET_THREAD_NUM()         ! thread number
       IOFFB=MQ(INDF-1)*IPROC                                    ! offset 'f'.
       IOFFC=MQ(INDC-1)*IPROC                                    ! offset 'c'.
       IOFFE=MQ(INDE-1)*IPROC                                    ! offset 'e'
       IOFFD=MQ(INDH-1)*IPROC+MQ(IND5+MTHRD  +IPROC+1)           ! offset data
       IOFFI=MQ(IND5-1)*IPROC+MQ(IND5+MTHRD*2+IPROC+1)+MTHRD*9+2 ! offset indices
*     ----- reset ------------------------------------------------------
       LPRNT=.FALSE.
       LHIST=(IPROC.EQ.0)
       REC=NRC            ! floating point value
       IF(NLOOPN.EQ.1.AND.MOD(NRC,100000).EQ.0) THEN 
        WRITE(*,*) 'Record',NRC,' ... still reading'
       END IF
*
*      printout/debug only for one thread at a time
*

*      flag for record printout -----------------------------------------

        LPRNT=.FALSE.
        IF(NEWITE.AND.(ITERAT.EQ.1.OR.ITERAT.EQ.3)) THEN
         IF(NRC.EQ.NRECPR) LPRNT=.TRUE.
         IF(NRC.EQ.NRECP2) LPRNT=.TRUE.
         IF(NRC.EQ.NRECER) LPRNT=.TRUE.
        END IF
        IF (LPRNT)THEN
!$OMP ATOMIC
         NPRDBG=NPRDBG+1               ! number of threads with debug
         IF (NPRDBG.EQ.1) IPRDBG=IPROC ! first thread with debug
         IF (IPROC.NE.IPRDBG) LPRNT=.FALSE.
*         print *, ' LPRNT ', NRC, NPRDBG, IPRDBG, IPROC, LPRNT
        ENDIF
        IF(LPRNT) THEN
         WRITE(1,*) ' '
         WRITE(1,*) '------------------ Loop',NLOOPN,
     +              ': Printout for record',NRC,IPROC
         WRITE(1,*) ' '
        END IF  

*     ----- print data -------------------------------------------------

        IF(LPRNT) THEN
         IMEAS=0              ! local derivatives
         IST=ISFRST(IBUF)
         NST=ISLAST(IBUF)
1025     CALL ISJAJB(NST,IST,JA,JB,JSP)
         IF(JA.NE.0) THEN
            IF(IMEAS.EQ.0) WRITE(1,1121)
            IMEAS=IMEAS+1
            WRITE(1,1122) IMEAS,GLDER(JA),GLDER(JB), 
     +      (INDER(JA+J),GLDER(JA+J),J=1,JB-JA-1)
            GOTO 1025
         END IF
1121     FORMAT(/'Measured value and local derivatives'/
     +   '  i measured std_dev  index...derivative ...')
1122     FORMAT(I3,2G12.4,3(I3,G12.4)/(27X,3(I3,G12.4)))

         IMEAS=0              ! global derivatives
         IST=ISFRST(IBUF)
         NST=ISLAST(IBUF)
1026     CALL ISJAJB(NST,IST,JA,JB,JSP)
         IF(JA.NE.0) THEN
           IF(IMEAS.EQ.0) WRITE(1,1123)
           IMEAS=IMEAS+1
           IF (JB.LT.IST) THEN
            IF(IST-JB.GT.2) THEN
            WRITE(1,1124) IMEAS,(JTGBL(INDER(JB+J)),INDER(JB+J),
     +      JVGBI(INDER(JB+J)),GLDER(JB+J),J=1,IST-JB)
            ELSE
            WRITE(1,1125) IMEAS,(JTGBL(INDER(JB+J)),INDER(JB+J),
     +      JVGBI(INDER(JB+J)),GLDER(JB+J),J=1,IST-JB)
            END IF
           ENDIF
           GOTO 1026
         END IF
1123     FORMAT(/'Global derivatives'/
     +   '  i  label gindex vindex derivative ...')   
1124     FORMAT(I3,2(I9,I7,I7,G12.4)/(3X,2(I9,I7,I7,G12.4)))
1125     FORMAT(I3,2(I9,I7,I7,G12.4))
        END IF

*      ----- first loop -------------------------------------------------
*       ------ prepare local fit ------ 
*      count local and global derivates 
*      subtract actual alignment parameters from the measured data 

       IF(LPRNT) THEN
          WRITE(1,*) ' '
          WRITE(1,*) 
     &      'Data corrections using values of global parameters'
          WRITE(1,*) 
     &      '=================================================='
          WRITE(1,101)
       END IF
       NALG=0                         ! count number of global derivatives
       NALC=0                         ! count number of local derivatives
       IMEAS=0
       IST=ISFRST(IBUF)
       NST=ISLAST(IBUF)
 01    CALL ISJAJB(NST,IST,JA,JB,JSP)
       IF(JA.NE.0) THEN       
          RMEAS=GLDER(JA)               ! data
*         subtract global ... from measured value
          DO J=1,IST-JB                 ! global parameter loop
           ITGBI=INDER(JB+J)            ! global parameter label
           RMEAS=RMEAS-GLDER(JB+J)*SNGL(DQ(IND0/2+ITGBI)) ! subtract   !!! reversed
           IF (ICALCM.EQ.1) THEN
            IJ=JVGBI(ITGBI)             ! index of variable global parameter
            IF(IJ.GT.0) THEN
              IJN=MQ(INDE+IOFFE+IJ)          ! get index of index
              IF(IJN.EQ.0) THEN              ! not yet included
                NALG=NALG+1                  ! count
                MQ(INDC+IOFFC+NALG)=IJ       ! store global index
                MQ(INDE+IOFFE+IJ  )=NALG     ! store back index 
              END IF 
             END IF
           ENDIF
          END DO
          IF(LPRNT) THEN
             IMEAS=IMEAS+1
             IF (JB.LT.IST) 
     +         WRITE(1,102) IMEAS,GLDER(JA),RMEAS,GLDER(JB) 
          END IF
          QM(IDOT4+JA)=RMEAS               ! global contribution subtracted
          DO J=1,JB-JA-1              ! local parameter loop
           IJ=INDER(JA+J)
           NALC=MAX(NALC,IJ)          ! number of local parameters
          END DO
          GOTO 01 
       END IF 
 101   FORMAT(' index measvalue   corrvalue          sigma') 
 102   FORMAT(I6,2X,2G12.4,' +-',G12.4)  

       IF(NALC.LE.0) GOTO 90
       IF(NALC.GT.MLOCAL) THEN
          WRITE(*,*) NALC, ' local parameters (=>skip), max. is', 
     &       MLOCAL, NRC
          GOTO 90
       END IF
       IF (ICALCM.EQ.1) THEN
         DO K=1,NALG*NALC
           VLG(K)=0.0D0 ! reset global-local matrix
         ENDDO
         NGG=(NALG*NALG+NALG)/2 
         MQ(IND5+IOFFI-1)=NRC       ! index header:
         MQ(IND5+IOFFI  )=NALG      ! event number, number of global par
         CALL SORT1K(MQ(INDC+IOFFC+1),NALG) ! sort global par.
         DO K=1,NALG
           IEXT=MQ(INDC+IOFFC+K)
           MQ(IND5+IOFFI+K)=IEXT    ! global par indices
           MQ(INDE+IOFFE+IEXT)=K    ! update back index
         ENDDO
         DO K=1,NGG
           DQ(INDH/2+IOFFD+K)=0.0D0 ! reset global-global matrix 
         ENDDO
         MQ(IND5        +IPROC+1)=MQ(IND5        +IPROC+1)+1
         MQ(IND5+MTHRD  +IPROC+1)=MQ(IND5+MTHRD  +IPROC+1)+NGG
         MQ(IND5+MTHRD*2+IPROC+1)=MQ(IND5+MTHRD*2+IPROC+1)+NALG+2
       ENDIF
*      ----- iteration start and check ---------------------------------       

       NTER=1                         ! first loop without down-weighting
       IF(NLOOPN.NE.1.AND.LHUBER.NE.0) NTER=LHUBER

*      check matrix for bordered band structure (MBDR+MBND+1 <= NALC)     
       MBND=-1
       MBDR=NALC
       DO I=1, NALC
          IBANDH(I)=0
       ENDDO      
      
       ITER=0
 05    ITER=ITER+1                  ! outlier suppresssion iteration
       RESMAX=0.0
       IF(LPRNT) THEN
          WRITE(1,*) ' '
          WRITE(1,*) 'Outlier-suppression iteration',ITER,' of',NTER
          WRITE(1,*) '==========================================' 
          WRITE(1,*) ' '
          IMEAS=0 
       END IF

*      ----- second loop ------------------------------------------------
*      accumulate normal equations for local fit and determine solution

       DO I=1,NALC
        BLVEC(I)=0.0D0                  ! reset vector
       END DO
       DO I=1,(NALC*NALC+NALC)/2 ! GF: FIXME - not really, local parameter number...
        CLMAT(I)=0.0D0                  ! reset matrix  
       END DO
       NEQ=0
       NDOWN=0
       NWEIG=0
       IST=ISFRST(IBUF)
       NST=ISLAST(IBUF)
 10    CALL ISJAJB(NST,IST,JA,JB,JSP)
       IF(JA.NE.0) THEN       
          RMEAS=GLDER(JA)               ! data
          RERR =GLDER(JB)               ! ... and the error
          WGHT =1.0/RERR**2             ! weight from error
          NEQ=NEQ+1                     ! count equation
          NWEIG=NWEIG+1
          RESID=RMEAS-RV(NEQ)    ! subtract previous fit
          IF(NLOOPN.NE.1.AND.ITER.NE.1.AND.LHUBER.NE.0) THEN
             IF(ITER.LE.3) THEN
                IF(ABS(RESID).GT.CHUBER*RERR) THEN     ! down-weighting
                   WGHT=WGHT*CHUBER*RERR/ABS(RESID)
                   NDOWN=NDOWN+1
                END IF
             ELSE       ! Cauchy 
                WGHT=WGHT/(1.0+(RESID/RERR/CAUCHY)**2)
             END IF  
          END IF
          CHAST='   '
          IF(ABS(RESID).GT.CHUBER*RERR) CHAST='*  '
          IF(ABS(RESID).GT.3.0*RERR) CHAST='** '
          IF(ABS(RESID).GT.6.0*RERR) CHAST='***'
 
          IF(LPRNT.AND.ITER.NE.1.AND.NTER.NE.1) THEN
             IF(IMEAS.EQ.0) WRITE(1,*) 'Second loop: accumulate' 
             IF(IMEAS.EQ.0) WRITE(1,103) 
             IMEAS=IMEAS+1
             DOWN=1.0/SQRT(WGHT)
             R1=RESID/RERR
             R2=RESID/DOWN
             WRITE(1,104) IMEAS,RMEAS,RESID,RERR,R1,CHAST,R2
          END IF
 103      FORMAT(' index corrvalue    residuum          sigma', 
     +       '     nresid     cnresid')
 104      FORMAT(I6,2X,2G12.4,' +-',G12.4,F7.2,1X,A3,F8.2)

          DO J=1,JB-JA-1 ! normal equations, local parameter loop
           IJ=INDER(JA+J)          ! local parameter index J
           BLVEC(IJ)=BLVEC(IJ)+DBLE(WGHT)*DBLE(RMEAS)*DBLE(GLDER(JA+J))
           DO K=1,J
            IK=INDER(JA+K)         ! local parameter index K
            JK=IJSYM(IJ,IK)        ! index in symmetric matrix
            CLMAT(JK)=CLMAT(JK)    ! force double precision
     +        +DBLE(WGHT)*DBLE(GLDER(JA+J))*DBLE(GLDER(JA+K))
*           check for band matrix substructure
            IF (ITER.EQ.1) THEN
              ID=IABS(IJ-IK)+1
              IM=MIN(IJ,IK)
              IBANDH(ID)=MAX(IBANDH(ID),IM)
            ENDIF
           END DO
          END DO
          GOTO 10
       END IF  
*      for non trivial fits check for bordered band matrix structure     
       IF (ITER.EQ.1.AND.NALC.GT.5.AND.LFITBB.GT.0) THEN
         KX=-1
         KBDR=0
         KBDRX=0
         ICMN=NALC**3 ! cost (*6) should improve by at least factor 2
         DO K=MIN(NALC,MXBND+2),2,-1
            KBND=K-2
            KBDR=MAX(KBDR,IBANDH(K))
            ICOST=6*(NALC-KBDR)*(KBND+KBDR+1)**2+2*KBDR**3
             IF (ICOST.LT.ICMN.AND.KBDR.LE.MXBDR) THEN
               ICMN=ICOST
               KX=K
               KBDRX=KBDR
            ENDIF
         ENDDO
          IF (KX.GT.0) THEN
             MBND=KX-2
             MBDR=KBDRX
          ENDIF          
       ENDIF
*
       IF (MBND.GE.0) THEN
*      fast solution for border banded matrix (inverse for ICALCM>0)
          IF (NLOOPN.EQ.1) THEN
             NBNDR=NBNDR+1
             NBDRX=MAX(NBDRX,MBDR)
             NBNDX=MAX(NBNDX,MBND)
          ENDIF
*        
          INV=0
          IF (NLOOPN.LE.LFITNP.AND.ITER.EQ.1) INV=1 ! band part of inverse (for pulls)
          IF (ICALCM.EQ.1.OR.LPRNT) INV=2     ! complete inverse     
          CALL SQMIBB(CLMAT,BLVEC,NALC,MBDR,MBND,INV,NRANK,
     &                VBND,VBDR,AUX,VBK,VZRU,SCDIAG,SCFLAG)
       ELSE
*      full inversion and solution 
          INV=2
          CALL SQMINV(CLMAT,BLVEC,NALC,NRANK,SCDIAG,SCFLAG)
       ENDIF
*      check for NaNs
       NAN=0
       DO K=1, NALC
          IF ((.NOT.(BLVEC(K).LE.0.0D0)).AND.
     +       (.NOT.(BLVEC(K).GT.0.0D0))) NAN=NAN+1
       ENDDO
            
       IF(LPRNT) THEN 
          WRITE(1,*) ' '
          WRITE(1,*) 'Parameter determination:',NALC,' parameters,',
     +               ' rank=',NRANK
          WRITE(1,*) '-----------------------'
          IF(NDOWN.NE.0) WRITE(1,*) '   ',NDOWN,' data down-weighted'
          WRITE(1,*) ' ' 
       END IF  
      
*      ----- third loop -------------------------------------------------
*      calculate single residuals remaining after local fit and chi^2 
        
       SUMM=0.0D0
       SUWT=0.0
       NEQ=0      
       IMEAS=0
       IST=ISFRST(IBUF)
       NST=ISLAST(IBUF)
 30    CALL ISJAJB(NST,IST,JA,JB,JSP)
       IF(JA.NE.0) THEN
          RMEAS=GLDER(JA)               ! data (global contrib. subtracted)
          RERR =GLDER(JB)               ! ... and the error
          WGHT =1.0/RERR**2             ! weight from error
          NEQ=NEQ+1                     ! count equation
          RMLOC=0.0                     ! local fit result reset
          DO J=1,JB-JA-1                ! local parameter loop
           IJ=INDER(JA+J)
           RMLOC=RMLOC+GLDER(JA+J)*SNGL(BLVEC(IJ)) ! local fit result
          END DO
          RMEAS=RMEAS-RMLOC             ! reduced to residual
 
*         calculate pulls? (needs covariance matrix)
          IF(ITER.EQ.1.AND.INV.GT.0.AND.NLOOPN.LE.LFITNP) THEN
           DVAR=0.0D0
           DO J=1,JB-JA-1
            IJ=INDER(JA+J)
            DO K=1,JB-JA-1
             IK=INDER(JA+K)
             JK=IJSYM(IJ,IK)
             DVAR=DVAR+CLMAT(JK)*DBLE(GLDER(JA+J))*DBLE(GLDER(JA+K))
            ENDDO
           ENDDO
*          some variance left to define a pull?          
           IF (0.999999D0/DBLE(WGHT).GT.DVAR) THEN
            PULL=RMEAS/SNGL(DSQRT(1.0D0/DBLE(WGHT)-DVAR))
            IF (LHIST) THEN
             IF (JB.LT.IST) THEN
              CALL HMPENT(13,PULL) ! histogram pull
              CALL GMPMS(5,REC,PULL) 
             ELSE
              CALL HMPENT(14,PULL) ! histogram pull
             ENDIF
            ENDIF
           ENDIF
          ENDIF

          IF(ITER.EQ.1.AND.JB.LT.IST.AND.LHIST) 
     &      CALL GMPMS(4,REC,RMEAS/RERR) ! residual (with global deriv.)

          DCHI2=WGHT*RMEAS*RMEAS
c          DCHIT=DCHI2
          RESID=RMEAS
          IF(NLOOPN.NE.1.AND.ITER.NE.1.AND.LHUBER.NE.0) THEN
             IF(ITER.LE.3) THEN
                IF(ABS(RESID).GT.CHUBER*RERR) THEN     ! down-weighting
                   WGHT=WGHT*CHUBER*RERR/ABS(RESID)
                   DCHI2=2.0*CHUBER*(ABS(RESID)/RERR-0.5*CHUBER)
                END IF
             ELSE
                WGHT=WGHT/(1.0+(RESID/RERR/CAUCHY)**2)
                DCHI2=LOG(1.0+(RESID/RERR/CAUCHY)**2)*CAUCHY**2
             END IF
          END IF
          CHAST='   '
          IF(ABS(RESID).GT.CHUBER*RERR) CHAST='*  '
          IF(ABS(RESID).GT.3.0*RERR) CHAST='** '
          IF(ABS(RESID).GT.6.0*RERR) CHAST='***'

          RV(NEQ)=RMLOC            ! local fit result
          DOWN=1.0/SQRT(WGHT)
          R1=RESID/RERR
          R2=RESID/DOWN
c          SUWT=SUWT+DCHI2/DCHIT
          SUWT=SUWT+RERR/DOWN
          IF(LPRNT) THEN
             IF(IMEAS.EQ.0) WRITE(1,*) 'Third loop: single residuals'
             IF(IMEAS.EQ.0) WRITE(1,105)
             IMEAS=IMEAS+1
             IF(RESID.LT.0.0) R1=-R1
             IF(RESID.LT.0.0) R2=-R2
             WRITE(1,106) IMEAS,GLDER(JA),RMEAS,RERR,R1,CHAST,R2
          END IF 
 105      FORMAT(' index corrvalue    residuum          sigma', 
     +        '     nresid     cnresid')
 106      FORMAT(I6,2X,2G12.4,' +-',G12.4,F7.2,1X,A3,F8.2)

          IF(ITER.EQ.NTER) THEN
             QM(IDOT4+JA)=RMEAS            ! store remaining residual
             RESMAX=MAX(RESMAX,ABS(RMEAS)/RERR) 
          END IF 

          IF(ITER.EQ.1.AND.LHIST) THEN
             IF (JB.LT.IST) THEN
                CALL HMPENT( 3,RMEAS/RERR) ! histogram norm residual
             ELSE
                CALL HMPENT(12,RMEAS/RERR) ! histogram norm residual
             ENDIF
          ENDIF            
          SUMM=SUMM+DCHI2        ! accumulate chi-square sum
          GOTO 30 
       END IF    
       NDF=NEQ-NRANK 
       RESING=(FLOAT(NWEIG)-SUWT)/FLOAT(NWEIG)
       IF (LHIST) THEN
        IF(ITER.EQ.1) CALL HMPENT( 5,FLOAT(NDF))  ! histogram Ndf      
        IF(ITER.EQ.1) CALL HMPENT(11,FLOAT(NALC)) ! histogram Nlocal      
        IF(NLOOPN.EQ.2.AND.ITER.EQ.NTER) CALL HMPENT(6,RESING)
       ENDIF
       IF(LPRNT) THEN
         WRITE(1,*) ' '
         WRITE(1,*) 'Chi^2=',SUMM,' at',NDF,' degrees of freedom: ',
     +   '3-sigma limit is',CHINDL(3,NDF)*FLOAT(NDF)
         WRITE(1,*) SUWT,' is sum of factors, compared to',NWEIG,
     +             ' Downweight fraction:',RESING 
       END IF
       IF(NRANK.NE.NALC.OR.NAN.GT.0) THEN
          NREJ(0)=NREJ(0)+1         ! count cases
          IF (NREC3.EQ.MAXI4) NREC3=NRC
          IF(LPRNT) THEN
            WRITE(1,*) ' rank deficit/NaN ', NALC, NRANK, NAN
            WRITE(1,*) '   ---> rejected!'
          END IF
          GOTO 90
       END IF 
       IF(NDF.LE.0) THEN
          NREJ(1)=NREJ(1)+1         ! count cases
          IF(LPRNT) THEN
             WRITE(1,*) '   ---> rejected!'
          END IF
          GOTO 90
       END IF 
      
       CHNDF=SNGL(SUMM/DFLOAT(NDF))

       IF(ITER.EQ.1.AND.LHIST) CALL HMPENT(4,CHNDF) ! histogram chi^2/Ndf
       IF(ITER.LT.NTER) GOTO 05     ! outlier iteration loop

       NDFS=NDFS+NDF              ! (local) sum of Ndf
       SNDF=SNDF+DFLOAT(NDF)*DW1  ! (local) weighted sum of Ndf

*      ----- reject eventually ------------------------------------------

       IF(NEWITE.AND.ITERAT.EQ.2) THEN ! find record with largest Chi^2/Ndf
          IF(NRECP2.LT.0.AND.CHNDF.GT.QM(IND5+MTHRD*8+IPROC+1)) THEN
             QM(IND5+MTHRD*8+IPROC+1)=CHNDF
             MQ(IND5+MTHRD*7+IPROC+1)=NRC
          END IF
       END IF

c       IF(CHDFRJ.NE.0.0) THEN ! outlier rejection
c          IF(SUMM/FLOAT(NDF).GT.CHDFRJ) THEN
c             IF(SUMM/FLOAT(NDF).GT.1.0E3) THEN
c                NREJEC(2)=NREJEC(2)+1   ! count cases with huge chi^2 
c                GOTO 90
c             END IF
*            add to FVALUE 
c             DCHI2=CHDFRJ*FLOAT(NDF)    ! total contribution limit
c             CALL ADDSUM(DCHI2)         ! add total contribution
c             NREJEC(3)=NREJEC(3)+1      ! count cases with large chi^2
c             GOTO 90
c          END IF
       CHICHI=CHINDL(3,NDF)*FLOAT(NDF)
* GF       IF(SUMM.GT.50.0*CHICHI) THEN ! huge 
* CHK CHICUT<0: NO cut (1st iteration)
       IF(CHICUT.GE.0) THEN
         IF(SUMM.GT.CHHUGE*CHICHI) THEN ! huge 
             NREJ(2)=NREJ(2)+1    ! count cases with huge chi^2
             IF(LPRNT) THEN
                  WRITE(1,*) '   ---> rejected!'
             END IF
            GOTO 90
         END IF

         IF(CHICUT.GT.0.0) THEN      
            CHLIMT=SNGL(CHICUT)*CHICHI
c          WRITE(*,*) 'chi^2 ',SUMM,CHLIMT,CHICUT,CHINDL(3,NDF),NDF
            IF(SUMM.GT.CHLIMT) THEN
               IF(LPRNT) THEN
                 WRITE(1,*) '   ---> rejected!'   
              END IF
*              add to FVALUE
               DCHI2=CHLIMT               ! total contribution limit
               DCHI2S=DCHI2S+DCHI2*DW1    ! add total contribution
              NREJ(3)=NREJ(3)+1      ! count cases with large chi^2
              GOTO 90
            END IF 
         END IF
       ENDIF

       IF(LHUBER.GT.1.AND.DWCUT.NE.0.0.AND.RESING.GT.DWCUT) THEN
*         add to FVALUE
          DCHI2=SUMM                 ! total contribution
          DCHI2S=DCHI2S+DCHI2*DW1    ! add total contribution
          NREJ(3)=NREJ(3)+1      ! count cases with large chi^2
c          WRITE(*,*) 'Downweight fraction cut ',RESING,DWCUT,SUMM
          IF(LPRNT) THEN
             WRITE(1,*) '   ---> rejected!'
          END IF
          GOTO 90
       END IF 

       IF(NEWITE.AND.ITERAT.EQ.2) THEN ! find record with largest residual
          IF(NRECPR.LT.0.AND.RESMAX.GT.QM(IND5+MTHRD*6+IPROC+1)) THEN
             QM(IND5+MTHRD*6+IPROC+1)=RESMAX
             MQ(IND5+MTHRD*5+IPROC+1)=NRC
          END IF 
       END IF 
 

*      ----- fourth loop ------------------------------------------------
*      update of global matrix and vector according to the "Millepede"
*      principle, from the global/local information

       IST=ISFRST(IBUF)
       NST=ISLAST(IBUF)
   50  CALL ISJAJB(NST,IST,JA,JB,JSP) 
       IF(JA.GT.0) THEN

          RMEAS=GLDER(JA)               ! data residual
          RERR =GLDER(JB)               ! ... and the error
          WGHT =1.0/RERR**2             ! weight from measurement error 
          DCHI2=WGHT*RMEAS*RMEAS        ! least-square contribution

          IF(NLOOPN.NE.1.AND.LHUBER.NE.0) THEN       ! check residual
             RESID=ABS(RMEAS)
             IF(RESID.GT.CHUBER*RERR) THEN
                WGHT=WGHT*CHUBER*RERR/RESID          ! down-weighting                
                DCHI2=2.0*CHUBER*(RESID/RERR-0.5*CHUBER) ! modified contribution
             END IF
          END IF
          DCHI2S=DCHI2S+DCHI2*DW1    ! add to total objective function

*         global-global matrix contribution: add directly to gg-matrix 

          DO J=1,IST-JB
           IVGBJ=JVGBI(INDER(JB+J))     ! variable-parameter index
           IF(IVGBJ.GT.0) THEN
              DQ(INDF/2+IOFFB+IVGBJ)=DQ(INDF/2+IOFFB+IVGBJ)
     &                              +DW1*WGHT*RMEAS*GLDER(JB+J) ! vector  !!! reverse
              IF(ICALCM.EQ.1) THEN
                 IJE=MQ(INDE+IOFFE+IVGBJ)        ! get index of index, non-zero
                 DO K=1,J
                  IVGBK=JVGBI(INDER(JB+K))
                  IF(IVGBK.GT.0) THEN
                    IKE=MQ(INDE+IOFFE+IVGBK)        ! get index of index, non-zero
                    IA=MAX(IJE,IKE)          ! larger
                    IB=MIN(IJE,IKE)          ! smaller
                    IJ=IB+(IA*IA-IA)/2
                    DQ(INDH/2+IOFFD+IJ)=DQ(INDH/2+IOFFD+IJ)
     &                                 -DW1*WGHT*GLDER(JB+J)*GLDER(JB+K)
                  END IF 
                 END DO
              END IF
           END IF
          END DO

*         normal equations - rectangular matrix for global/local pars
*         global-local matrix contribution: accumulate rectangular matrix
          IF (ICALCM.NE.1) GOTO 50
          DO J=1,IST-JB
           IVGBJ=JVGBI(INDER(JB+J))           ! variable-parameter index
           IF(IVGBJ.GT.0) THEN 
              IJE=MQ(INDE+IOFFE+IVGBJ)        ! get index of index, non-zero
              DO K=1,JB-JA-1
               IK=INDER(JA+K)           ! local index
               JK=IK+(IJE-1)*NALC       ! matrix index
               VLG(JK)=VLG(JK)+DW2*WGHT*GLDER(JB+J)*GLDER(JA+K)
              END DO
           END IF
          END DO
          GOTO 50 ! next equation
       END IF 


*      ----- final matrix update ----------------------------------------
*      update global matrices and vectors
       IF(ICALCM.NE.1) GOTO 90 ! matrix update
*      (inverse local matrix) * (rectang. matrix) -> CORM
*                                        T
*      resulting symmetrix matrix  =   G   *   Gamma^{-1}   *   G
*
       CALL DBAVAT(CLMAT,VLG,DQ(INDH/2+IOFFD+1),NALC,-NALG)

*      (rectang. matrix) * (local param vector)   -> CORV
*
*      resulting vector = G * q (q = local parameter)
*
*      CALL DBGAX(DQ(IGLMA/2+1),BLVEC,DQ(ICORV/2+1),NALG,NALC)  ! not done
*
*      the vector update is not done, because after local fit it is zero!
       NFRED=MQ(INDH-1)-MQ(IND5+MTHRD  +IPROC+1)-MQ(INDH-2)
       NFREI=MQ(IND5-1)-MQ(IND5+MTHRD*2+IPROC+1)-MQ(IND5-2)
       IF (NFRED.LT.0.OR.NFREI.LT.0) THEN
        NB=MQ(IND5+IPROC+1)
C        print *, ' flush ', IPROC, NRC, NB, NFRED, NFREI
        JOFFD=MQ(INDH-1)*IPROC           ! offset data
        JOFFI=MQ(IND5-1)*IPROC+MTHRD*9+2 ! offset indices
        USED=FLOAT(MQ(IND5+MTHRD  +IPROC+1))/FLOAT(MQ(INDH-1))
        MQ(IND5+MTHRD*3+IPROC+1)=MQ(IND5+MTHRD*3+IPROC+1)
     &                          +IFIX(1000.0*USED+0.5)
        USED=FLOAT(MQ(IND5+MTHRD*2+IPROC+1))/FLOAT(MQ(IND5-1))
        MQ(IND5+MTHRD*4+IPROC+1)=MQ(IND5+MTHRD*4+IPROC+1)
     &                          +IFIX(1000.0*USED+0.5)
!$OMP CRITICAL
        MQ(INDH-4)=MQ(INDH-4)+1
        MQ(IND5-4)=MQ(IND5-4)+1
*
        DO IB=1,NB
         IJN=0
         DO IN=1,MQ(IND5+JOFFI)
          I=MQ(IND5+JOFFI+IN)
*         DQ(IGVEC/2+I)=DQ(IGVEC/2+I)+DQ(ICORV/2+IN)  ! not done: = zero
          DO JN=1,IN
           IJN=IJN+1
           J=MQ(IND5+JOFFI+JN)
           CALL MUPDAT(I,J,-DQ(INDH/2+JOFFD+IJN))  ! matrix update
          END DO
         END DO
         JOFFD=JOFFD+IJN
         JOFFI=JOFFI+MQ(IND5+JOFFI)+2
        ENDDO
!$OMP END CRITICAL
*       reset counter, pointers
        DO K=0,2
         MQ(IND5+K*MTHRD+IPROC+1)=0
        ENDDO
       ENDIF

 90    CONTINUE

       IF(LPRNT) THEN
          WRITE(1,*) ' '
          WRITE(1,*) '------------------ End of printout for record',NRC
          WRITE(1,*) ' '
       END IF
*
       DO I=1,NALG                 ! reset global index array
        IEXT=MQ(INDC+IOFFC+I) 
        MQ(INDE+IOFFE+IEXT)=0
       END DO
*
      ENDDO
!$OMP END PARALLEL DO
*
      IF (ICALCM.EQ.1) THEN
*     flush remaining matrices
       DO K=1,MTHRD ! update statistics
        MQ(INDH-3)=MQ(INDH-3)+1
        USED=FLOAT(MQ(IND5+MTHRD  +K))/FLOAT(MQ(INDH-1))
        MQ(IND5+MTHRD*3+K)=MQ(IND5+MTHRD*3+K)+IFIX(1000.0*USED+0.5)
        MQ(INDH-5)=MQ(INDH-5)+MQ(IND5+MTHRD*3+K)
        MQ(INDH-6)=MAX(MQ(INDH-6),MQ(IND5+MTHRD*3+K))
        MQ(IND5+MTHRD*3+K)=0
        MQ(IND5-3)=MQ(IND5-3)+1
        USED=FLOAT(MQ(IND5+MTHRD*2+K))/FLOAT(MQ(IND5-1))
        MQ(IND5+MTHRD*4+K)=MQ(IND5+MTHRD*4+K)+IFIX(1000.0*USED+0.5)
        MQ(IND5-5)=MQ(IND5-5)+MQ(IND5+MTHRD*4+K)
        MQ(IND5-6)=MAX(MQ(IND5-6),MQ(IND5+MTHRD*4+K))
        MQ(IND5+MTHRD*4+K)=0
       ENDDO
*
!$OMP  PARALLEL
!$OMP& DEFAULT(PRIVATE)
!$OMP& SHARED(MQ,QM,DQ,IND5,INDH,MTHRD)  ! 'global' variables 
       IPROC=0
!$        IPROC=OMP_GET_THREAD_NUM()         ! thread number
       DO JPROC=0,MTHRD-1
        NB=MQ(IND5+JPROC+1)
C        print *, ' flush end ', JPROC, NRC, NB
        JOFFD=MQ(INDH-1)*JPROC           ! offset data
        JOFFI=MQ(IND5-1)*JPROC+MTHRD*9+2 ! offset indices
        DO IB=1,NB
C         print *, '   buf end ', JPROC,IB,MQ(IND5+JOFFI-1),MQ(IND5+JOFFI)
         IJN=0 
         DO IN=1,MQ(IND5+JOFFI)
          I=MQ(IND5+JOFFI+IN)
!$        IF (MOD(I,MTHRD).EQ.IPROC) THEN
            DO JN=1,IN
             IJN=IJN+1
             J=MQ(IND5+JOFFI+JN)
             CALL MUPDAT(I,J,-DQ(INDH/2+JOFFD+IJN))  ! matrix update
            END DO
!$        ELSE
!$          IJN=IJN+IN
!$        ENDIF
         END DO
         JOFFD=JOFFD+IJN
         JOFFI=JOFFI+MQ(IND5+JOFFI)+2
        ENDDO
       ENDDO
!$OMP END PARALLEL
      ENDIF
*
      IF(NEWITE.AND.ITERAT.EQ.2) THEN ! get worst records (for printrecord -1 -1)
       IF (NRECPR.LT.0) THEN
        DO K=1,MTHRD
         IF (QM(IND5+MTHRD*6+K).GT.VALUE1) THEN
          VALUE1=QM(IND5+MTHRD*6+K)
          NREC1 =MQ(IND5+MTHRD*5+K)
         ENDIF
        ENDDO
       ENDIF
       IF (NRECP2.LT.0) THEN
        DO K=1,MTHRD
         IF (QM(IND5+MTHRD*8+K).GT.VALUE2) THEN
          VALUE2=QM(IND5+MTHRD*8+K)
          NREC2 =MQ(IND5+MTHRD*7+K)
         ENDIF
        ENDDO
       ENDIF
      ENDIF
*
      END  




************************************************************************

      SUBROUTINE PRTGLO    ! print final log file
*
*     label                              I10
*     parameter value                    F12.5
*     presigma                           F12.5
*     difference of parameters values
*     difference at last iteration
*     error (standard deviation)

#include "largeint.inc"
#include "dynal.inc"
      INTEGER LABELE(3)
      REAL    COMPNT(3)
      SAVE
#include "mpinds.inc"
*     ...

      LUP=09
      CALL MVOPEN(LUP,'millepede.res')

      WRITE(*,*) ' '
      WRITE(*,*) '         Result of fit for global parameters'
      WRITE(*,*) '         ==================================='
      WRITE(*,*) ' '

      WRITE(*,101)

      WRITE(LUP,*) 'Parameter   ! first 3 elements per line are',
     +    ' significant (if used as input)'
      IPRLIM=10
      DO ITGBI=1,NTGB  ! all parameter variables
       ITGBL=JTGBL(ITGBI)
       IVGBI=JVGBI(ITGBI)
       PAR=SNGL(DQ(IND0/2+ITGBI))      ! initial value
       IF(IVGBI.GT.0) THEN
          DPA=PAR-QM(IND2+ITGBI)       ! difference
          IF(METSOL.EQ.1.OR.METSOL.EQ.2) THEN
             II=(IVGBI*IVGBI+IVGBI)/2
             GMATI=SNGL(DQ(IGMAT/2+II))
             ERR=SQRT(ABS(GMATI))
             IF(GMATI.LT.0.0D0) ERR=-ERR
          END IF
       END IF
       IF(ITGBI.LE.IPRLIM) THEN 
          IF(IVGBI.LE.0) THEN
             WRITE(*  ,102) ITGBL,PAR,QM(IND1+ITGBI)
          ELSE 
             IF(METSOL.EQ.1.OR.METSOL.EQ.2) THEN
                WRITE(*,102) ITGBL,PAR,QM(IND1+ITGBI),DPA,ERR
             ELSE
                WRITE(*,102) ITGBL,PAR,QM(IND1+ITGBI),DPA
             END IF
          END IF
       ELSE IF(ITGBI.EQ.IPRLIM+1) THEN
          WRITE(*  ,*)
     +    '... (further printout suppressed, but see log file)'
       END IF 

*      file output
       IF(IVGBI.LE.0) THEN
          WRITE(LUP,102) ITGBL,PAR,QM(IND1+ITGBI)
       ELSE
          IF(METSOL.EQ.1.OR.METSOL.EQ.2) THEN
             WRITE(LUP,102) ITGBL,PAR,QM(IND1+ITGBI),DPA,ERR
          ELSE
             WRITE(LUP,102) ITGBL,PAR,QM(IND1+ITGBI),DPA
          END IF
       END IF
      END DO
      REWIND LUP
      CLOSE(UNIT=LUP) 

      IF(METSOL.EQ.2) THEN        ! diagonalisation: write eigenvectors
         CALL MVOPEN(LUP,'millepede.eve')
         DO I=NAGB,1,-1
          IF(DQ(IDUX6/2+I).GT.0.0D0) GOTO 20
         END DO
         I=1
  20     IMIN=I                        ! index of smallest pos. eigenvalue
         IEV=0

         DO ISUB=0,MIN(15,IMIN-1)
          IF(ISUB.LT.10) THEN
             I=IMIN-ISUB
          ELSE
             I=ISUB-9
          END IF 
         
c        DO I=IMIN,MAX(1,IMIN-9),-1    ! backward loop, up to 10 vectors
          WRITE(*,*) 'Eigenvector ',I,' with eigenvalue',DQ(IDUX6/2+I)
          WRITE(LUP,*) 'Eigenvector ',I,' with eigenvalue',DQ(IDUX6/2+I)
          DO J=1,NAGB
           IJ=INDT/2+J+(I-1)*NAGB      ! index with eigenvector array
           IF(J.LE.NVGB) THEN
              ITGBI=JTGBI(J)
              LABEL=JTGBL(ITGBI)
           ELSE
              LABEL=NVGB-J             ! label negative for constraints 
           END IF
           IEV=IEV+1
           LABELE(IEV)=LABEL
           COMPNT(IEV)=SNGL(DQ(IJ))    ! component
           IF(IEV.EQ.3) THEN
              WRITE(LUP,103) (LABELE(IE),COMPNT(IE),IE=1,IEV)
              IEV=0
           END IF 
          END DO
          IF(IEV.NE.0) WRITE(LUP,103) (LABELE(IE),COMPNT(IE),IE=1,IEV)
          IEV=0
          WRITE(LUP,*) ' '
         END DO

      END IF   
      
 101  FORMAT(1X,'    label       parameter      presigma        differ',
     +                   '         error'/
     +       1X,'-----------',4X,4('-------------'))
 102  FORMAT(I10,2X,4G14.5,F8.3)
 103  FORMAT(3(I11,F11.7,2X))
      END



      SUBROUTINE AVPROD(N,X,B)  ! product A(sym) * X => B
#include "largeint.inc"
#include "dynal.inc"
      DOUBLE PRECISION X(N),B(N)
      INTEGER(KIND=LARGE) K,KK,KL,KU,LL,LJ,INDIJ,INDID,I2
      SAVE
#include "mpinds.inc"
*     ...
      IF(MQ(INDJ-2).EQ.2) GOTO 10
      IJ=0     ! full symmetric matrix
      DO I=1,N
       B(I)=DQ(INDJ/2+IJ+I)*X(I)
       DO J=1,I-1
        B(J)=B(J)+DQ(INDJ/2+IJ+J)*X(I)
        B(I)=B(I)+DQ(INDJ/2+IJ+J)*X(J)
       END DO
       IJ=IJ+I
      END DO
      RETURN

 10   IF(LQ(INDZ/LSCALE+2).NE.N+1) 
     +   STOP 'AVPROD: mismatched vector and matrix'
      IENCDB=MQ(INDD-1)
      IENCDM=ISHFT(1,IENCDB)-1
*
!$    DO I=1,N
!$     B(I)=0.0D0             ! reset 'global' B()
!$    END DO
*
      NSPC=MQ(INDZ-1)         ! number of precisions, 1=D, 2=F
!$OMP PARALLEL DO             ! parallelize row loop
!$OMP& PRIVATE(I2,K,KK,LL,KL,KU,INDID,INDIJ,J,JC,JN,LJ) ! variables local to thread
!$OMP& REDUCTION(+:B)         ! private copy of B(N) for each thread, combined at end, init with 0.
!$OMP& SCHEDULE(DYNAMIC,1024) ! slot of 1024 'I' for next idle thread
      DO I=1,N
       B(I)=DQ(INDJ/2+I)*X(I)    ! diagonal elements
*                                ! off-diagonals double precision
       I2=2*I
       KK=LQ(INDZ/LSCALE+I2-1) ! offset in 'd' (column lists) 
       LL=LQ(INDZ/LSCALE+I2)   ! offset in 'j' (matrix)
       KL=0
       KU=LQ(INDZ/LSCALE+I2+1)-1-KK
       INDID=INDD+KK
       INDIJ=INDJ/2+LL
       IF (MQ(INDID).NE.0) THEN  ! no compression
        DO K=KL,KU
         J=MQ(INDID+K)
         B(J)=B(J)+DQ(INDIJ+K)*X(I)
         B(I)=B(I)+DQ(INDIJ+K)*X(J)
        END DO
       ELSE
        LJ=0
        KU=((KU+1)*8)/9-1         ! number of regions (-1)
        INDID=INDID+KU/8+1        ! skip group offsets
        DO KL=0,KU
         JC=MQ(INDID+KL)
         J=ISHFT(JC,-IENCDB)
         JN=IAND(JC, IENCDM)
         DO JJ=1,JN
          B(J)=B(J)+DQ(INDIJ+LJ)*X(I)
          B(I)=B(I)+DQ(INDIJ+LJ)*X(J)
          J=J+1
          LJ=LJ+1
         ENDDO
        ENDDO
       ENDIF
*
       IF (NSPC.GT.1) THEN 
        I2=2*(I+N+1)                  ! off-diagonals single precision
        KK=LQ(INDZ/LSCALE+I2-1) ! offset in 'd' (column lists) 
        LL=LQ(INDZ/LSCALE+I2)   ! offset in '.' (matrix)
        KL=0
        KU=LQ(INDZ/LSCALE+I2+1)-1-KK
        INDID=INDD+KK
        INDIJ=IDOT1+LL
        IF (MQ(INDID).NE.0) THEN  ! no compression
         DO K=KL,KU
          J=MQ(INDID+K)
          B(J)=B(J)+DBLE(QM(INDIJ+K))*X(I)
          B(I)=B(I)+DBLE(QM(INDIJ+K))*X(J)
         END DO
        ELSE
         LJ=0
         KU=((KU+1)*8)/9-1         ! number of regions (-1)
         INDID=INDID+KU/8+1        ! skip group offsets
         DO KL=0,KU
          JC=MQ(INDID+KL)
          J=ISHFT(JC,-IENCDB)
          JN=IAND(JC, IENCDM)
          DO JJ=1,JN
           B(J)=B(J)+DBLE(QM(INDIJ+LJ))*X(I)
           B(I)=B(I)+DBLE(QM(INDIJ+LJ))*X(J)
           J=J+1
           LJ=LJ+1
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDDO
!$OMP END PARALLEL DO
      END

      FUNCTION IJADD(ITEMA,ITEMB)      ! index using "d" and "z"
#include "largeint.inc"
#include "dynal.inc"
      INTEGER(KIND=LARGE) IJADD,K,KK,KL,KU,INDID,ND,LL,K8,ITEM1
*     ...
      IJADD=0
      ND=LQ(INDZ/LSCALE+2)-1          ! dimension of matrix
      ITEM1=MAX(ITEMA,ITEMB)          ! larger index
      ITEM2=MIN(ITEMA,ITEMB)          ! smaller index
      IF(ITEM2.LE.0.OR.ITEM1.GT.ND) RETURN
      IF(ITEM1.EQ.ITEM2) THEN         ! diagonal element
         IJADD=ITEM1
         RETURN  
      ENDIF 
*                                     ! off-diagonal element
      NSPC=MQ(INDZ-1)                 ! number of precisions, 1=D, 2=F
      IENCDB=MQ(INDD-1)               ! encoding info
      IENCDM=ISHFT(1,IENCDB)-1
      ISGN=1
      DO ISPC=1,NSPC
         KK=LQ(INDZ/LSCALE+ITEM1*2-1) ! offset in 'd' (column lists) 
         LL=LQ(INDZ/LSCALE+ITEM1*2)   ! offset in 'j' (matrix)
         KL=0
         KU=LQ(INDZ/LSCALE+ITEM1*2+1)-1-KK
         INDID=INDD+KK
         IF (MQ(INDID).EQ.0) THEN     ! compression ?
*
            KU=((KU+1)*8)/9-1        ! number of regions (-1)
            INDID=INDID+KU/8+1       ! skip group offsets
            KL=0
            IF(KU.LT.KL) THEN
*               WRITE(*,*) ITEM1,ITEM2,' not found',KK,KL,KU
               GOTO 100
            END IF
 10         K=(KL+KU)/2                    ! binary search
            JTEMC=MQ(INDID+K)              ! compressed information
            JTEM =ISHFT(JTEMC,-IENCDB)     ! first column of region
            JTEMN=JTEM+IAND(JTEMC,IENCDM)  ! first column after region
            IF(ITEM2.LT.JTEM.OR.ITEM2.GE.JTEMN) THEN
               IF(ITEM2.LT.JTEM) THEN
                  KU=K-1
               ELSE IF(ITEM2.GE.JTEMN) THEN
                  KL=K+1
               END IF
               IF(KL.LE.KU) GOTO 10
*            WRITE(*,*) ITEM1,ITEM2,' IJADD not found',KK,KL,KU,K
               GOTO 100
            END IF
            K8=K/8                         ! region group (-1)
            LL=LL+MQ(INDD+KK+K8)           ! offset for group of (8) regions 
            DO KL=K8*8,K-1 
               LL=LL+IAND(MQ(INDID+KL),IENCDM) ! add region lengths
            ENDDO
            IJADD=LL+ITEM2-JTEM
* 
         ELSE
*
            IF(KU.LT.KL) THEN
*               WRITE(*,*) ITEM1,ITEM2,' not found',KK,KL,KU
               GOTO 100
            END IF
 20         K=(KL+KU)/2                  ! binary search
            JTEM=MQ(INDID+K)
            JTEMN=JTEM
            IF(ITEM2.NE.JTEM) THEN
               IF(ITEM2.LT.JTEM) THEN
                  KU=K-1
               ELSE IF(ITEM2.GT.JTEM) THEN
                  KL=K+1
               END IF
               IF(KL.LE.KU) GOTO 20
*            WRITE(*,*) ITEM1,ITEM2,' IJADD not found',KK,KL,KU,K
               GOTO 100
            END IF
            IJADD=LL+K
*
         ENDIF
         IJADD=IJADD*ISGN
         RETURN
*
  100    CONTINUE
         ITEM1=ITEM1+ND+1
         ISGN=-ISGN
      ENDDO
      END

      SUBROUTINE SECHMS(DELTAT,NHOUR,MINUT,SECND)
*     DELTAT = time in sec  -> NHOUR,MINUT,SECND 
*     ...
      NSECD=NINT(DELTAT) ! -> integer  
      NHOUR=NSECD/3600
      MINUT=NSECD/60-60*NHOUR
      SECND=DELTAT-60*(MINUT+60*NHOUR) 
      END




*
*     Using indices a,b,d,f,g
*
*     In INONE/UPONE:
*                    a    label & index
*                    b    hash pointer
*
*     In INTWO/UPTWO:                     
*                    f    pairs
*                    g    hash pointer
*
*     In MSPARS:   
*                    d    sparse matrix pointer, column lists (optionally compressed)
*                    z    sparse matrix pointer, row pointers (LARGE integers)



*
*     1- and 2 dim. tables ---------------------------------------------
*
*     Functions INONE(), INSEC(), with sub-arrays "a" and "b"
*     =======================================================
*
*     Functions INONE, INSEC and subroutine UPONE are
*     used to collect items, i.e. labels, and to order and translate them.
*     Two sub-arrays "a" and "b" are used.
*
*     In the first phase items are collected and stored by calling
*        IRES=INONE(ITEM)
*     or by calling
*        IRES=INSEC(ITEM).
*     At the first entry the two sub-arrays "a" and "b" of length 2N
*     are generated with a start length for N=100 entries.
*     In array "a" two words are reserved for each item: (ITEM, count).
*     The function INONE(ITEM) returns the number of the item, the 
*     function INSEC the second value, which is the count at this stage.
*     At each entry the argument is compared with the already stored items,
*     new items are stored. Search
*     for entries is done using hash-indices, stored in sub-array "b".
*     The initial hash-index is
*
*        j = 1 + mod(ITEM, n_prime) + N
*
*     where n_prime is the largest prime number less than N.
*     At each entry the count is increased by one. If N items are stored,
*     the size of the sub-arrays is increased by calling
*        CALL UPONE.
*     The entries in sub-array "a" are sorted, and the hash indices in
*     sub-array "b" are redefined. No new entries are accepted after this
*     call.
*        MQ(INDA-0) = number of words of sub-array
*        MQ(INDA-1) = number of stored items
*        MQ(INDA-2) = 1 (arrays fixed)
*        MQ(INDA-3) = (next dimension value)
*        MQ(INDA-4) = prime number
*        MQ(INDA-5) = number of overflows
*        MQ(INDA-6) = nr of variable parameters
*        MQ(INDA-7) = number of entries  
*
*     The functions return after the UPONE call:  
*     The function INONE(ITEM) returns the number of the item, the
*     function ISEC the second value, which is the count.
*     
*     A loop through all entries is done by the code
*      DO I=1,MQ(INDA-1)
*       MQ(INDA+I+I-1) is ITEM
*       MQ(INDA+I+I  ) is second quantity 
*      END DO
*     In this loop the second quantity (initially the count) can be 
*     redefined by the user (but the first quantity ITEM should never
*     be changed).  
*    
*
*     Function INTWO(,) with sub-arrays "f" and "g"
*     =============================================
*
*     Function INTWO and subroutine UPTWO are used to collect pairs of
*     items, i.e. label pairs, and to order and translate them.
*     Two sub-arrays "f" and "g" are used.  
*
*     In the first phase items are collected and stored by calling
*        IRES=INTWO(ITEMA,ITEMB)
*     Pairs with ITEMA=ITEMB are ignored, pairs (ITEMA,ITEMB) and 
*     (ITEMB,ITEMA) are equivalent. Internal the first item is the larger
*     one of the two items. 
*     
*     After transmisssion of all pairs the final status to the sub-arrays
*     is done by 
*        CALL UPTWO
*     After 
*        MQ(INDF-0) = number of words of sub-array
*        MQ(INDF-1) = number of stored items
*        MQ(INDF-2) = 2 (arrays fixed)
*        MQ(INDF-3) = (next dimension value)
*        MQ(INDF-4) = prime number
*        MQ(INDF-5) = number of overflows
*        MQ(INDF-6) = number 
*        MQ(INDF-7) = number
*  
*     Sparse matrix building:
*         CALL MSPARS(NDIM)
*     Pointer sub-arrays "z" and "d" are created and defined, using the data 
*     from sub-array "f" (from INTWO). All NDIM diagonal elements and
*     all off-diagonal elements from the sub-array "f" are prepared.
*     After the call the corresponding double precision matrix array "j"
*     and with compression the single precision matrix array "." have 
*     to be created by 
*         CALL MEGARR/L('j D',2*(NDIM+NSPARS-NSINGL),'sparse symmetric matrix')
*         MQ(INDJ-2)=2               ! 2=sparse matrix
*         CALL MEGARR/L('. F',NSINGL                ,'sparse symmetric matrix')
*     with NSPARS off-diagonal elements to be used and with compression
*     NSINGL of those in single precision.
*
*     Matrix elements (I,J) are updated by first getting the linear
*     index by
*         IJ=IJADD(I,J)
*     A value zero is returned, if the matrix element is not prepared.
*     Otherwise the element can be updated e.g. by
*        DQ(INDJ/2+IJ)=DQ(INDJ/2+IJ)+PLUS   ! double IJ>0
*        QM(IDOT1 -IJ)=QM(IDOT1 -IJ)+PLUS   ! single IJ<0
*
*     The product of the symmetric matrix with the vector X to give the
*     resulting vector B is done by
*         CALL AVPROD(N,X,B)
*       


      INTEGER FUNCTION INONE(ITEM1)             ! translate 1-D identifier to nrs
*     using subarrays "a" and "b" (remove before}
*
*     length increase by 50 % plus 100 elements
*
      IMPLICIT NONE  
      INTEGER ITEM,ITEM1,ITEM2
#include "largeint.inc"
#include "dynal.inc"
      INTEGER INSEC,J,K,NSTART,IPRIME
      INTEGER(KIND=LARGE) NWORDS,LND
      LOGICAL SEC
*     ...
      SEC=.FALSE.
      ITEM=ITEM1
      GOTO 10
      ENTRY INSEC(ITEM2)
      SEC=.TRUE.
      ITEM=ITEM2
 10   J=0
      IF(ITEM.LE.0) GOTO 40
      IF(INDA.EQ.0) THEN
         NSTART=100                   ! initial number  
         CALL MEGARR('a I',2*NSTART,'INONE: label & index')
         CALL MEGARR('b I',2*NSTART,'INONE: hash pointer')
         MQ(INDA-1)=0                 ! number of stored items
         MQ(INDA-2)=0                 ! =0 during build-up
         MQ(INDA-3)=NSTART            ! next number
         MQ(INDA-4)=IPRIME(NSTART)    ! prime number
         MQ(INDA-5)=0                 ! number of overflows
         MQ(INDA-6)=0                 ! nr of variable parameters 
      END IF

 20   J=1+MOD(ITEM,MQ(INDA-4))+LNDA/2
 30   K=J
      J=MQ(INDB+K)

      IF(J.NE.0) THEN                 ! normal case: find item
         IF(ITEM.NE.MQ(INDA+J+J-1)) GOTO 30
      ELSE                            ! not found
         IF(2*MQ(INDA-1).EQ.LNDA.OR.MQ(INDA-2).NE.0) THEN
            MQ(INDA-5)=MQ(INDA-5)+1 ! overflow
            J=0
            GOTO 40
         END IF
         MQ(INDA-1)=MQ(INDA-1)+1      ! increase number of elements
         MQ(INDA-3)=MQ(INDA-1)           
         J=MQ(INDA-1)
         MQ(INDB+K)=J                 ! hash index
         MQ(INDA+J+J-1)=ITEM          ! add new item  
         IF(2*MQ(INDA-1).EQ.LNDA) THEN
            CALL MEGSPC(NWORDS)       ! nr of free words
            LND=MIN(2*LNDA,LNDA+NWORDS/2)
            IF (.NOT.(LND.GT.LNDA.AND.LND.LE.MAXI4)) THEN
               WRITE(*,*) 'INONE: too many items' 
               STOP 'recompile with larger memory'
            ENDIF
            MQ(INDA-3)=LND/2
            CALL UPONE          ! update with larger dimension
            IF (LVLLOG.GT.1)
     +         WRITE(LUNLOG,*) 'INONE: array increased to',
     +                         LNDA,' words'  
            GOTO 20             ! redefine index 
         END IF   
      END IF
      IF(MQ(INDA-2).EQ.0) THEN
         MQ(INDA+J+J)=MQ(INDA+J+J)+1 ! increase counter
         MQ(INDA-7)=MQ(INDA-7)+1
      END IF  
 40   IF(SEC) THEN
         INSEC=0   
         IF(J.NE.0) INSEC=MQ(INDA+J+J)
      ELSE
         INONE=J
      END IF 
      END


      SUBROUTINE UPONE                  ! update, redefine hash indices
      IMPLICIT NONE 
#include "largeint.inc"
#include "dynal.inc"
      INTEGER I,J,K,MQDIFF,IPRIME  
      SAVE
*     ...
      IF(MQ(INDA-1).EQ.MQ(INDA-3)) THEN
         CALL SORT2K(MQ(INDA+1),MQ(INDA-1)) ! sort items
      END IF 
      CALL MEGARR('a',2*MQ(INDA-3),' ')
      CALL MEGARR('b',2*MQ(INDA-3),' ')
      MQDIFF=MQ(INDA-3)-MQ(INDA-1)
      MQ(INDA-3)=MQ(INDA-1)    
      MQ(INDA-4)=IPRIME(LNDA/2)          ! prime number < LNDA
      CALL MEGZER('b') 

      DO I=1,MQ(INDA-1)                      ! redefine hash 
       J=1+MOD(MQ(INDA+I+I-1),MQ(INDA-4))+LNDA/2
 10    K=J
       J=MQ(INDB+K)
       IF(J.NE.0.AND.J.NE.I) GOTO 10
       MQ(INDB+K)=I 
      END DO
      IF(MQDIFF.NE.0) RETURN
      MQ(INDA-2)=1       ! set flag to inhibit further updates
      IF (LVLLOG.GT.1) THEN
       WRITE(LUNLOG,*) ' '
       WRITE(LUNLOG,*) 'INONE: array reduced to',LNDA,' words' 
       WRITE(LUNLOG,*) 'INONE:',MQ(INDA-1),' items stored.'
      ENDIF
      END 


      INTEGER FUNCTION INTWO(ITEMA,ITEMB)      ! translate 2-D to numbers
      IMPLICIT NONE
      INTEGER ITEMA,ITEMB,IPRIME
#include "largeint.inc"
#include "dynal.inc" 
      INTEGER ITEM1,ITEM2,NSTART,NC1,NC2,NC3,NC4,MIO,J,K
      INTEGER(KIND=LARGE) NWORDS,LND
      SAVE
*     ...
      INTWO=0 
      ITEM1=MAX(ITEMA,ITEMB)
      ITEM2=MIN(ITEMA,ITEMB)
      IF(ITEM1.LE.0.OR.ITEM2.LE.0) RETURN ! no illegal elements

      IF(INDF.EQ.0) THEN
         NSTART=100
         CALL MEGARR('f I',2*NSTART,'INTWO: pairs') 
         CALL MEGARR('g I',2*NSTART,'INTWO: hash pointer') 
         MQ(INDF-1)=0
         MQ(INDF-2)=0
         MQ(INDF-3)=LNDF
         MQ(INDF-4)=IPRIME(LNDF/2)   ! prime number
         MQ(INDF-5)=0                    ! overflows
         MQ(INDF-6)=0                    ! 
         MQ(INDF-7)=0                    ! number of entries

         NC1=0 ! count +1
         NC2=0 ! count +1
         NC3=0 ! count Mio
         NC4=0 ! count Mio
         MIO=1

      END IF
      IF(ITEM1.EQ.ITEM2)           RETURN ! no diagonal elements
      NC1=NC1+1

c 10   J=1+MOD(ITEM2+(MOD(ITEM1,MQ(INDF-4))*ITEM1-ITEM2)/2,MQ(INDF-4))
c     +   +LNDF/2

 10   J=1+LNDF/2
     +   +INT(MOD(DBLE(ITEM1)*DBLE(ITEM2),DBLE(MQ(INDF-4)))+0.5D0)

 20   K=J
      J=MQ(INDG+K)
      IF(J.NE.0) THEN            ! normal case: find item
         NC2=NC2+1
         IF(ITEM1.NE.MQ(INDF+J+J-1)) GOTO 20
         IF(ITEM2.NE.MQ(INDF+J+J  )) GOTO 20
      ELSE                       ! not found
         IF(MQ(INDF-2).NE.0) RETURN ! with result zero
         IF(MQ(INDF-1).EQ.LNDF/2) THEN
            MQ(INDF-5)=MQ(INDF-5)+1      ! overflow
            WRITE(*,*) 'INTWO: ovfl return with J=',
     +                  J,MQ(INDF-1),LNDF/2
            STOP 'recompile with larger memory'
         END IF
         MQ(INDF-1)=MQ(INDF-1)+1         ! increase number of elements
         MQ(INDF-3)=MQ(INDF-1)           ! next 
         J=MQ(INDF-1)
         MQ(INDG+K)=J           ! hash index
         MQ(INDF+J+J-1)=ITEM1   ! add new item
         MQ(INDF+J+J  )=ITEM2
         IF(MQ(INDF-1).EQ.LNDF/2) THEN
            CALL MEGSPC(NWORDS)
            IF(NWORDS.LT.100) THEN
               WRITE(*,*) 'insufficient space for INTWO routine' 
               STOP 'recompile with larger memory'
            END IF
            LND=MIN(2*LNDF,LNDF+NWORDS/2)
            IF (.NOT.(LND.GT.LNDF.AND.LND.LE.MAXI4)) THEN
               WRITE(*,*) 'INTWO: too many item pairs '
               STOP 'recompile with larger memory'
            ENDIF
            MQ(INDF-3)=LND/2 
            CALL UPTWO          ! update with larger dimension
            IF (LVLLOG.GT.1)
     +         WRITE(LUNLOG,*) 'INTWO: array increased to',
     +                          LNDF,' words'
            GOTO 10             ! redefine index
         END IF
      END IF
      MQ(INDF-7)=MQ(INDF-7)+1 
      INTWO=J
      IF(MQ(INDF-2).NE.0) RETURN         ! with result non-zero 
      IF(NC1.GE. 1000000) THEN
         NC3=NC3+1 
         NC1=NC1-1000000
      END IF
      IF(NC2.GE. 1000000) THEN
         NC4=NC4+1 
         NC2=NC2-1000000
      END IF

      IF(NC3.GT.0.AND.MOD(NC3,MIO).EQ.0) THEN
         MIO=2*MIO
         WRITE(*,101) NC3,MQ(INDF-1)
      END IF

 101  FORMAT(I9,' Mio entries  in INTWO',3X,I10,' elements',
     +       F7.2,' tests/entry')
      END


      SUBROUTINE UPTWO                 ! update, redefine hash indices
      IMPLICIT NONE
#include "largeint.inc"
#include "dynal.inc"
      INTEGER MQDIFF,I,J,K,IPRIME
      SAVE 
*     ...
      IF(INDF.EQ.0) RETURN
      IF(MQ(INDF-1).EQ.MQ(INDF-3)) THEN 
         CALL SORT2K(MQ(INDF+1),MQ(INDF-1))  ! sort items
      END IF 
      CALL MEGARR('f',2*MQ(INDF-3),' ')   ! array for pairs
      CALL MEGARR('g',2*MQ(INDF-3),' ')   ! array for hash pointer
      CALL MEGZER('g')                    ! reset hash pointer
      MQDIFF=MQ(INDF-3)-MQ(INDF-1)
      MQ(INDF-3)=MQ(INDF-1)
      MQ(INDF-4)=IPRIME(LNDF/2)       ! prime number < LNDF

      DO I=1,MQ(INDF-1)                   ! redefine hash
c       J=1+MOD(MQ(INDF+I+I)+(MOD(MQ(INDF+I+I-1),MQ(INDF-4))*
c     +         MQ(INDF+I+I-1)-MQ(INDF+I+I))/2,MQ(INDF-4))+LNDF/2 

       J=1+LNDF/2
     +    +INT(MOD(DBLE(MQ(INDF+I+I))*DBLE(MQ(INDF+I+I-1))
     +    ,DBLE(MQ(INDF-4)))+0.5)    

 10    K=J
       J=MQ(INDG+K)
       IF(J.NE.0.AND.J.NE.I) GOTO 10
       MQ(INDG+K)=I
      END DO
      IF(MQDIFF.NE.0) RETURN
      IF (LVLLOG.GT.1) THEN
       WRITE(LUNLOG,*) ' '
       WRITE(LUNLOG,*) 'INTWO: array reduced to',LNDF,' words'
       WRITE(LUNLOG,*) 'INTWO:',MQ(INDF-1),' pairs stored.'
      ENDIF
      END


      INTEGER FUNCTION IPRIME(N)               ! largest prime number < N
      IMPLICIT NONE
      INTEGER N
      INTEGER NPRIME,NSQRT,I  
*     ...
      SAVE
      NPRIME=N                               ! max number
      IF(MOD(NPRIME,2).EQ.0) NPRIME=NPRIME+1 ! ... odd number
  10  NPRIME=NPRIME-2                        ! next lower odd number
      NSQRT=IFIX(SQRT(FLOAT(NPRIME)))
      DO I=3,NSQRT,2                         ! 
       IF(I*(NPRIME/I).EQ.NPRIME) GOTO 10    ! test prime number
      END DO      
      IPRIME=NPRIME
      END 


      SUBROUTINE MSPARS(NAGB)          ! prepare index sparse storage
      IMPLICIT NONE
#include "largeint.inc"
#include "dynal.inc"
      INTEGER NAGB
      INTEGER(KIND=LARGE) K
*     ...
*     New storage scheme for sparse symmetric N*N matrix M
*     (with option 'compress' for lossless compression)
*
*     + subarray "z" contains the row pointers:
*       table with 2 columns and (N+1) rows (2*(N+1) rows for compression)
*       first  column contains for each row offset in column lists "d"
*       second column contains for each row offset in subarray containing
*       the values (row<=N+1: "d", row>N+1: '.')
*       MQ(INDZ-1) contains the number of precisions used 
*         (1= only Double, 2= Double+Single)
*     + subarray "d" contains the column lists:
*       for each row I the list of columns J for the off-diagonal elements
*       Mij is stored
*       - uncompressed: simple list of all columnn indices
*       - compressed  : regions of consecutive columns
*         only first index J1 and number of indices JN is stored
*         encoded in single 4byte word (upper MQ(INDD-1)=int(log2(N))+1
*            bits for J1, lower 31-int(log2(N)) bits for JN,
*            region is split if JN doesn't fit in allocated bits)
*         for each group of 8 regions first the group offsets are stored
*            (sum(JN) of all regions before that group, 
*             allows faster access in IJADD, 
*             first offset is allways 0, tagging compression)
*       - each row is only compressed if the compressed list is
*         shorter than the uncompressed one
*     + subarray "j" contains double precision values:
*       first the N diagonal elements, than frequent off-diagonals
*     + subarray "." contains single precision values:
*       rare off-diagonals (depending on nsingle = 3rd par. of 'pairentries')
*       with nsingle=1 (default from 'compress') only single entry elements
*       are stored as simple floats (lossless as input form binary files is
*       only single precision)
*
*     Compression of column information not used for MSPARS (ISAVES=0),
*     but in NDBITS/SPBITS (ISAVES>0)
* 
      INTEGER L,IROW,JCOL,NELM
      INTEGER(KIND=LARGE) NSPARS,LNROW
      SAVE
*     ...      
      NSPARS=MQ(INDF-1) ! number of off-diagonal elements
      LNROW=NAGB+1
      CALL MEGARL('z L',LNROW*LSCALE*2,'sparse matrix row pointer')
      CALL MEGARL('d I',NSPARS,'sparse matrix column list')
      K    =0
      MQ(INDZ-1)=1   ! no compression, only double precision
      LQ(INDZ/LSCALE+1)=1
      LQ(INDZ/LSCALE+2)=NAGB+1
      L=1
      DO IROW=1,NAGB
 40    IF(L.GT.MQ(INDF-1)) GOTO 50 
       IF(MQ(INDF+L+L-1).NE.IROW) GOTO 50
       JCOL=MQ(INDF+L+L)
       IF(K.EQ.NSPARS) THEN
          STOP 'MSPARS: column list too short'
       END IF
       K=K+1
       MQ(INDD+K)=JCOL
       L=L+1
       GOTO 40
 50    CONTINUE
       LQ(INDZ/LSCALE+IROW*2+1)=K+1
       LQ(INDZ/LSCALE+IROW*2+2)=NAGB+K+1
      END DO
      NELM=K
      END



      SUBROUTINE LOOP1
#include "largeint.inc"
#include "dynal.inc"

      REAL    PLVS(3)    ! vector array: real and ...
      INTEGER LPVS(3)    ! ... integer
      EQUIVALENCE (PLVS(1),LPVS(1))
      CHARACTER*1 CH789(3)
      SAVE
#include "mpinds.inc"
      DATA CH789/'7','8','9'/
*     ...
      WRITE(LUNLOG,*) ' '
      WRITE(LUNLOG,*) 'LOOP1: starting' 
      CALL MSTART('LOOP1') 
*     add labels from parameter, constraints, measurements -------------
      DO K=7,9 ! subarrays 7, 8 and 9
       IF(INDICS(K).NE.0) THEN
          DO ITH=1,MQ(INDICS(K)-1)
           CALL MEGRDV(CH789(K-6),ITH,PLVS,NCOL)
           ITEM=LPVS(1)     ! label
           IDUM=INONE(ITEM)
          END DO
       END IF
      END DO
      IF(INDA.NE.0) THEN
         DO J=1,LNDA/2
          MQ(INDA+J+J)=0   ! reset count 
         END DO
         WRITE(LUNLOG,*) 'LOOP1:',MQ(INDA-1),
     +                   ' labels from txt data stored'
      END IF
      WRITE(LUNLOG,*) 'LOOP1: reading data files' 

*     read all data files and add all labels to global labels table ----
     
      IF(MPRINT.NE.0) THEN
         WRITE(*,*) 'Read all binary data files:'   
      END IF 
      CALL HMPLDF(1,'Number of words/record in binary file')
      CALL HMPDEF(8,0.0,60.0,'not_stored data per record')
*     define read buffer
      NC21=NCACHE/(21*MTHRDR) ! split read cache 1 : 10 : 10 for pointers, ints, floats
      NWRD=NC21+1 
      CALL MEGARR(': I',(NWRD+5)*MTHRDR,'read buffer, pointer')
      MQ(IDOT2-2)=NWRD*MTHRDR     ! number of buffer words (total)
      NWRD=NC21*10+2+NDIMB
      CALL MEGARR('; I',NWRD*MTHRDR,'read buffer, integer')
      MQ(IDOT3  )=MTHRDR          ! number of threads
      MQ(IDOT3-2)=NWRD*MTHRDR     ! number of buffer words 
      CALL MEGARR('! S',NWRD*MTHRDR,'read buffer, float')
      MQ(IDOT4  )=MTHRDR          ! number of (active) threads
      MQ(IDOT4-2)=NWRD*MTHRDR     ! number of buffer words
*
 10   CALL PEREAD(NR)  ! read records
      CALL PEPREP(0)   ! prepare records
      IF(NR.GT.0) GOTO 10 ! no EOD yet
*     release read buffer
      CALL MEGRM('!')
      CALL MEGRM(';')
      CALL MEGRM(':')
*
      IF(NHISTP.NE.0) THEN
         CALL HMPRNT(1)
         CALL HMPRNT(8)
      END IF 
      CALL HMPWRT(1)
      CALL HMPWRT(8) 
      CALL UPONE ! finalize the global label table     
      NTGB = MQ(INDA-1)     ! total number of labels/parameters
      WRITE(LUNLOG,*) 'LOOP1:',MQ(INDA-1),
     +            ' is total number NTGB of labels/parameters'
*     histogram number of entries per label ----------------------------

      CALL HMPLDF(2,'Number of entries per label')
      DO J=1,NTGB
       CALL HMPLNT(2,MQ(INDA+J+J))
      END DO
      IF(NHISTP.NE.0) CALL HMPRNT(2) ! print histogram 
      CALL HMPWRT(2) ! write to his file

      CALL MEGRM('c')                ! remove  

*     three subarrays for all global parameters ------------------------     

      CALL MEGARR('0 D',2*NTGB,'global parameters') 
      CALL MEGARR('1 S',NTGB,'pre-sigmas') ! presigmas 
      CALL MEGARR('2 S',NTGB,'global parameters at start') 
      CALL MEGARR('4 D',2*NTGB,'copy of global parameters') 

      NPRESG=0
      IF(IND7.NE.0) THEN 
         DO I=1,MQ(IND7-1)                  ! subarrays parallel to 'a'
          LABEL=MQ(IND7+(I-1)*MQ(IND7-2)+1)
          PARAM=QM(IND7+(I-1)*MQ(IND7-2)+2)
          PRESG=QM(IND7+(I-1)*MQ(IND7-2)+3)  
          IN=INONE(LABEL)
          IF(IN.NE.0) THEN
             DQ(IND0/2+IN)=PARAM
             IF(PRESG.GT.0.0) NPRESG=NPRESG+1 ! FIXME: check if enough 'entries'?
             QM(IND1+IN)=PRESG              ! insert pre-sigma 0 or > 0
c             IF(PRESG.GT.0.0) QM(IND6+IN)=PRESG ! pre-sigma copy > 0
             QM(IND2+IN)=PARAM
          ELSE
             WRITE(*,*) 'LABEL ...',LABEL,PARAM,PRESG,I,MQ(IND7-1)
             WRITE(*,*) 'Error in LOOP1: label not stored!'
          END IF  
         END DO 
      END IF
      WRITE(LUNLOG,*) 'LOOP1:',NPRESG,' is number of pre-sigmas'
      WRITE(*,*) 'LOOP1:',NPRESG,' is number of pre-sigmas'
      IF(NPRESG.EQ.0) WRITE(*,*) 'Warning: no pre-sigmas defined' 

*     determine flag variable (active) or fixed (inactive) -------------       

      INDAB=0
      DO I=1,MQ(INDA-1) ! = NTGB
       IF(MQ(INDA+I+I).GE.MREQEN.AND.QM(IND1+I).GE.0.0) THEN
          INDAB=INDAB+1
          MQ(INDA+I+I)=INDAB  ! variable, used in matrix (active)
       ELSE 
          MQ(INDA+I+I)=-1     ! fixed, not used in matrix (not active)
       END IF
      END DO   
      MQ(INDA-6)=INDAB ! counted variable
      NVGB=MQ(INDA-6)  ! nr of variable parameters
      WRITE(LUNLOG,*) 'LOOP1:',NVGB,
     +            ' is number NVGB of variable parameters'

*     translation table of length NVGB of total global indices ---------     

      CALL MEGARR('3 I',NVGB,'translation table  var -> total') 
      INDAB=0
      DO I=1,MQ(INDA-1) 
       IF(MQ(INDA+I+I).GT.0) THEN 
          INDAB=INDAB+1
          MQ(IND3+INDAB)=I  
       END IF
      END DO

*     regularization --------------------------------------------------- 
      CALL MEGARR('6 S',NVGB,'pre-sigmas weights') ! presigma weights
      WRITE(*,112) ' Default pre-sigma =',REGPRE,
     +             ' (if no individual pre-sigma defined)'
      WRITE(*,*)   'Pre-sigma factor is',REGULA

      IF(NREGUL.EQ.0) THEN
         WRITE(*,*) 'No regularization will be done'
      ELSE
         WRITE(*,*) 'Regularization will be done, using factor',REGULA
      END IF
 112  FORMAT(A,E9.2,A) 
      IF (NVGB.LE.0) STOP '... no variable global parameters'

      DO IVGBI=1,NVGB         ! IVGBI     = variable parameter index
       ITGBI=JTGBI(IVGBI)     ! ITGBI = global parameter index
       PRESG=QM(IND1+ITGBI)   ! get pre-sigma 
       PREWT=0.0              ! pre-weight
       IF(PRESG.GT.0.0) THEN
         PREWT=1.0/PRESG**2           ! 1/presigma^2
       ELSE IF(PRESG.EQ.0.0.AND.REGPRE.GT.0.0) THEN
         PREWT=1.0/REGPRE**2          ! default 1/presigma^2 
       END IF
       QM(IND6+IVGBI)=REGULA*PREWT    ! weight = factor / presigma^2
      END DO

c      WRITE(*,*) 'GlPa_index  GlPa_label  array1 array6'   
      DO I=1,NTGB
       ITGBL=JTGBL(I)
       IVGBI=JVGBI(I)
       IF(IVGBI.GT.0) THEN
c          WRITE(*,111) I,ITGBL,QM(IND1+I),QM(IND6+IVGBI)
       ELSE
c          WRITE(*,111) I,ITGBL,QM(IND1+I)
       END IF  
      END DO
c 111  FORMAT(I5,I10,F10.5,E12.4)       
      WRITE(*,101) 'NTGB',NTGB,'total number of parameters'
      WRITE(*,101) 'NVGB',NVGB,'number of variable parameters'

*     print overview over important numbers ----------------------------

      NRECAL=NREC
      IF(MPRINT.NE.0) THEN
         WRITE(*,*) ' '
         WRITE(*,101) '  NREC',NREC,'number of records'
         WRITE(*,101) 'MREQEN',MREQEN,'required number of entries'
         IF (MREQPE.GT.1) WRITE(*,101) 
     &      'MREQPE',MREQPE,'required number of pair entries'
         IF (MSNGPE.GE.1) WRITE(*,101) 
     &      'MSNGPE',MSNGPE,'max pair entries single prec. storage'
         WRITE(*,101) 'NTGB',NTGB,'total number of parameters'
         WRITE(*,101) 'NVGB',NVGB,'number of variable parameters'
         IF(MPRINT.GT.1) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'Global parameter labels:'
            MQI=MQ(INDA-1)
            IF(MQI.LE.100) THEN
               WRITE(*,*) (MQ(INDA+I+I-1),I=1,MQI)
            ELSE 
               WRITE(*,*) (MQ(INDA+I+I-1),I=1,30)
               WRITE(*,*) ' ...'
               MQI=((MQI-20)/20)*20+1
               WRITE(*,*) (MQ(INDA+I+I-1),I=MQI,MQ(INDA-1))
            END IF
         END IF 
         WRITE(*,*) ' ' 
         WRITE(*,*) ' '
      END IF
      WRITE(8,*)   ' '
      WRITE(8,101) '  NREC',NREC,'number of records'
      WRITE(8,101) 'MREQEN',MREQEN,'required number of entries'

c      WRITE(*,*) 'Existing sub-arrays:'
c      WRITE(*,*) ' 0   global parameter values'
c      WRITE(*,*) ' 1   pre-sigmas'
c      WRITE(*,*) ' 2   global parameter values at start'
c      WRITE(*,*) ' 3   translation table var -> total'
c      WRITE(*,*) ' 8   translation table var -> total'
c      WRITE(*,*) ' a   INONE: label & index '
c      WRITE(*,*) ' b   INONE: hash table'
c      WRITE(*,*) ' d   label groups'                     no no no no


      WRITE(LUNLOG,*) 'LOOP1: ending'
      WRITE(LUNLOG,*) ' '
      CALL MEND

 101  FORMAT(1X,A6,' =',I10,' = ',A)
      END 


      SUBROUTINE LOOP2
#include "largeint.inc"
#include "dynal.inc"
      SAVE
#include "mpinds.inc"
*
!$    INTEGER OMP_GET_THREAD_NUM
*
      ISFRST(IBUF)=MQ(IDOT2+IBUF)+1
      ISLAST(IBUF)=MQ(IDOT3+MQ(IDOT2+IBUF))
      INDER(I)=MQ(IDOT3+I)
      GLDER(I)=QM(IDOT4+I)
*
      DOUBLE PRECISION DGBM,DUBM,DSTAT(3)
      INTEGER*8 NOFF8,NBITS
      INTEGER(KIND=LARGE) NUSE,NDIMBI,NDIMSA(4),NDGN,MATSIZ(2),LNROW
*     ...
      WRITE(LUNLOG,*) ' '
      WRITE(LUNLOG,*) 'LOOP2: starting'
      CALL MSTART('LOOP2') 

*     two subarrays to get the global parameter indices, used in an event
      CALL MEGARR('c I',NVGB,'global index')
      CALL MEGARR('e I',NVGB,'back index')


*     constraints - determine number of constraints NCGB

      NCGB=0
      IF(IND8.NE.0) THEN
         I=-2
*        find next constraint header and count nr of constraints
 03      IF(I+2.LT.LND8) THEN
            I=I+2
            IF(MQ(IND8+I+1).EQ.0.AND.MQ(IND8+I+3).EQ.(-1)) THEN
               RHS=QM(IND8+I+2)
               SGM=QM(IND8+I+4)
               NCGB=NCGB+1
               I=I+2
            END IF
            GOTO 03
         END IF
      END IF
      WRITE(*,*) 'LOOP2:',NCGB,' constraints'

      NAGB=NVGB+NCGB ! total number of fit parameters
      NOFF8=INT8(NAGB)*INT8(NAGB-1)/2

*     read all data files and add all variable index pairs -------------

      IF(MATSTO.EQ.2) THEN
         IF(ISAVES.EQ.0) THEN  
            IDUM=INTWO(1,1)               ! for sparse storage    
         ELSE
            CALL MEGSPC(NUSE)
            NBITS=INT8(NUSE-11)*32
            MXBITS=NBITS/NOFF8 ! max bit field size fitting into free memory
            IF (MCMPRS.NE.0) ISAVES=MAX(ISAVES,2) ! identify single entries for compression 
            IF (MXBITS.LE.0) THEN
               WRITE(*,*)  'STOP: not enough memory '//
     &            'for bit storage of off-diagonal elements'
               STOP 'recompile with larger memory'
            ELSEIF (ISAVES.GT.MXBITS) THEN
               ISAVES=MXBITS
               MXCNT=2**MXBITS-1
               IF (MSNGPE.GE.MXCNT) MSNGPE=MXCNT-1
               IF (MHISPE.GT.MXCNT) MHISPE=MXCNT
               IF (MREQPE.GT.MXCNT) THEN
                  MREQPE=MXCNT
                  WRITE(*,*) 'Pairentries cut reset to (memory limit)',
     &               MREQPE
               ENDIF
            ENDIF
            CALL CLBITS(NAGB,NDIMBI,NENCDB,ISAVES)   ! get dimension for bit storage, encoding info
            CALL MEGARL('f I',NDIMBI,'INBITS: bit storage')
            CALL INBITS(1,1,1,MQ(INDF+1))
         END IF  
      END IF

*     reading events===reading events===reading events===reading events=
      NRECF =0  ! records with fixed global parameters
      NAEQNG=0  ! count number of equations (with global der.)
      NAEQNF=0  ! count number of equations ( " , fixed)
      NAEQNA=0  ! all
      WRITE(LUNLOG,*) 'LOOP2: start event reading'
*     monitoring for sparse matrix?
      IRECMM=0
      IF (MATSTO.EQ.2.AND.ISAVES.GT.0.AND.MATMON.NE.0) THEN
         NMATMO=0
         IF (MATMON.GT.0) THEN
            NRECMM=MATMON
         ELSE
            NRECMM=1
         ENDIF 
      ENDIF
      DO K=1,3
        DSTAT(K)=0.0D0
      ENDDO
*     define read buffer
      NC21=NCACHE/(21*MTHRDR) ! split read cache 1 : 10 : 10 for pointers, ints, floats
      NWRD=NC21+1 
      CALL MEGARR(': I',(NWRD+5)*MTHRDR,'read buffer, pointer')
      MQ(IDOT2-2)=NWRD*MTHRDR     ! number of buffer words (total)
      NWRD=NC21*10+2+NDIMB
      CALL MEGARR('; I',NWRD*MTHRDR,'read buffer, integer')
      MQ(IDOT3  )=MTHRDR          ! number of threads
      MQ(IDOT3-2)=NWRD*MTHRDR     ! number of buffer words 
      CALL MEGARR('! S',NWRD*MTHRDR,'read buffer, float')
      MQ(IDOT4  )=MTHRDR          ! number of (active) threads
      MQ(IDOT4-2)=NWRD*MTHRDR     ! number of buffer words
*
 10   CALL PEREAD(NR) ! read records
      CALL PEPREP(1)  ! prepare records
      IOFF=0
      DO IBUF=1,MQ(IDOT2)           ! buffer for current record
      NREC=MQ(IDOT3+ISFRST(IBUF)-2)   ! record 
c      IF(NREC.LE.100001) WRITE(*,*) 'Record ',NREC

*     Printout for DEBUG 
      IF(NREC.LE.MDEBUG) THEN
         NDA=0
         KFILE=MQ(IDOT4+ISFRST(IBUF)-1)   ! file
         WREC =QM(IDOT4+ISFRST(IBUF)-2)   ! weight
         WRITE(*,*) ' '
         WRITE(*,*) 'Record number ',NREC,' from file ',KFILE
         IF (WGH.NE.1.0) WRITE(*,*) '       weight ',WREC
         IST=ISFRST(IBUF)
         NST=ISLAST(IBUF)
 12      CALL ISJAJB(NST,IST,JA,JB,JSP)
         IF(JA.NE.0) THEN
            NDA=NDA+1
            IF(NDA.GT.MDEBG2) THEN
               IF(NDA.EQ.MDEBG2+1)  WRITE(*,*) '... and more data'
               GOTO 12
            END IF 
            WRITE(*,*) ' '
            WRITE(*,*) NDA, 
     +       ' Measured value =',GLDER(JA),' +- ',GLDER(JB)
            WRITE(*,*) 'Local derivatives:' 
            WRITE(*,107) (INDER(JA+J),GLDER(JA+J),J=1,JB-JA-1)
 107        FORMAT(6(I3,G12.4))
            IF (JB.LT.IST) THEN
             WRITE(*,*) 'Global derivatives:'
             WRITE(*,108) (JTGBL(INDER(JB+J)),INDER(JB+J),
     +                     JVGBI(INDER(JB+J)),GLDER(JB+J),J=1,IST-JB)
 108         FORMAT(3I11,G12.4)
            ENDIF
            IF(NDA.EQ.1) THEN
               WRITE(*,*)
     +         'total_par_label  __label__   var_par_index   derivative'
            END IF
            GOTO 12
         END IF
         WRITE(*,*) ' '
      END IF

      NAGBN =0                     ! count number of global derivatives
      NALCN =0                     ! count number of local  derivatives
      NAEQN =0                     ! count number of equations
      MAEQNF=NAEQNF
      IST=ISFRST(IBUF)
      NST=ISLAST(IBUF)
      NWRD=NST-IST+1
 15   CALL ISJAJB(NST,IST,JA,JB,JSP)
      IF(JA.EQ.0.AND.JB.EQ.0) GOTO 18
      NAEQN=NAEQN+1
      NAEQNA=NAEQNA+1
      IF(JA.NE.0) THEN
          IF (IST.GT.JB) NAEQNG=NAEQNG+1
          NFIXED=0
          DO J=1,IST-JB
          IJ=INDER(JB+J)           ! index of global parameter
          IJ=JVGBI(IJ)             ! change to variable parameter
          IF(IJ.GT.0) THEN
             IJN=MQ(INDE+IJ)         ! get index of index
             IF(IJN.EQ.0) THEN        ! not yet included
                NAGBN=NAGBN+1         ! count
                MQ(INDC+NAGBN)=IJ    ! store variable index
                MQ(INDE+IJ   )=NAGBN ! store back index 
             END IF
          ELSE
             NFIXED=NFIXED+1   
          END IF 
         END DO
         IF (NFIXED.GT.0) NAEQNF=NAEQNF+1         
      END IF

      IF(JA.NE.0.AND.JB.NE.0) THEN 
         DO J=1,JB-JA-1           ! local parameters
          IJ=INDER(JA+J)
          NALCN=MAX(NALCN,IJ)
         END DO
c         GOTO 15         
      END IF
      GOTO 15

 18   CONTINUE ! end-of-event
      IF (NAEQNF.GT.MAEQNF) NRECF=NRECF+1
      IRECMM=IRECMM+1
*     end-of-event-end-of-event-end-of-event-end-of-event-end-of-event-e

      MQ(INDC-1)=MAX(NAGBN,MQ(INDC-1)) ! maximum number of global parameters
      MQ(INDC-2)=MAX(NALCN,MQ(INDC-2)) ! maximum number of local parameters
      MQ(INDC-3)=MAX(NAEQN,MQ(INDC-3)) ! maximum number of equations

*     sample statistics for caching
      DSTAT(1)=DSTAT(1)+DFLOAT((NWRD+2)*2)               ! record size
      DSTAT(2)=DSTAT(2)+DFLOAT(NAGBN+2)                  ! indices,
      DSTAT(3)=DSTAT(3)+DFLOAT(NAGBN*NAGBN+NAGBN)        ! data for MUPDAT

      CALL SORT1K(MQ(INDC+1),NAGBN) ! sort global par.
      IOFF=IOFF+1
      MQ(IDOT2+IBUF)=IOFF
      MQ(IDOT3+IOFF)=IOFF+NAGBN
      DO I=1,NAGBN                  ! reset global index array
       IEXT=MQ(INDC+I) 
       MQ(INDE+IEXT)=0
       MQ(IDOT3+IOFF+I)=IEXT
      END DO
      IOFF=IOFF+NAGBN
*
      ENDDO
      IOFF=0
*
      IF (MATSTO.EQ.2) THEN
!$OMP  PARALLEL
!$OMP& IF (ISAVES.NE.0)
!$OMP& DEFAULT(PRIVATE)
!$OMP& SHARED(MQ,QM,DQ,INDF,IDOT2,IDOT3,MTHRD,ISAVES)  ! 'global' variables 
        IPROC=0
!$      IPROC=OMP_GET_THREAD_NUM()         ! thread number
        DO IBUF=1,MQ(IDOT2)
         IST=ISFRST(IBUF)
         NST=ISLAST(IBUF)
         DO I=IST,NST                 ! store all combinations     
          IEXT=MQ(IDOT3+I)             ! variable global index
          IF (ISAVES.EQ.0) THEN
           DO L=IST,I
            JEXT=MQ(IDOT3+L)
            IND=INTWO(IEXT,JEXT)       ! add pair to set of pairs
*           inline INTWO version removed - only 3 % faster than function
           END DO
          ELSE
!$         IF (MOD(IEXT,MTHRD).EQ.IPROC) THEN  ! distinct rows per thread
            DO L=IST,I
              JEXT=MQ(IDOT3+L)
              CALL INBITS(IEXT,JEXT,1,MQ(INDF+1)) ! save space
            END DO
!$         ENDIF
          ENDIF
         END DO
        ENDDO
!$OMP END PARALLEL
*        monitoring
         IF (ISAVES.GT.0.AND.MATMON.NE.0.AND.
     &      (IRECMM.GE.NRECMM.OR.IRECMM.EQ.MXREC)) THEN
            IF (NMATMO.EQ.0) THEN
               WRITE(*,*)
               WRITE(*,*) 'Monitoring of sparse matrix construction'
               WRITE(*,*) ' records ........ off-diagonal elements ',
     &                    '....... compression   memory'
               WRITE(*,*) '             non-zero used(double)  used',
     &                    '(float)       [%]       [GB]'
            ENDIF
            NMATMO=NMATMO+1
            JCMPRS=MAX(MCMPRS,MSNGPE)
            CALL CKBITS(NDIMSA,MQ(INDF+1),MREQPE,JCMPRS)
            GBC=4.0E-9*FLOAT(NDIMSA(2)+2*NDIMSA(3)+NDIMSA(4)) ! GB compressed
            GBU=1.2E-8*FLOAT(NDIMSA(3)+NDIMSA(4))             ! GB uncompressed
            CPR=100.0*GBC/GBU
            WRITE(*,1177) IRECMM,NDIMSA(1),NDIMSA(3),NDIMSA(4),CPR,GBC
 1177       FORMAT(I9,3I13,F10.2,F11.6)
   25       CONTINUE
            IF (MATMON.GT.0) THEN
               NRECMM=NRECMM+MATMON
            ELSE
               NRECMM=NRECMM*2
            ENDIF
            IF (IRECMM.GE.NRECMM) GOTO 25
         ENDIF
*
      ENDIF
*
      IF (NR.GT.0) GOTO 10 ! ... to read next block of events
*     release read buffer
      CALL MEGRM('!')
      CALL MEGRM(';')
      CALL MEGRM(':')
* 
      WRITE(LUNLOG,*) 'LOOP2: event reading ended - end of data'
      DO K=1,3
        DSTAT(K)=DSTAT(K)/DFLOAT(NREC)
      ENDDO
*     end=of=data=end=of=data=end=of=data=end=of=data=end=of=data=end=of


*     constraints and index pairs with Lagrange multiplier 


*     constraints - determine number of constraints NCGB and index-pairs  
*        Lagrange multiplier and global parameters


      NCGB=0
      IF(IND8.NE.0) THEN
         I=-2
*        find next constraint header and count nr of constraints
 30      IF(I+2.LT.LND8) THEN
            I=I+2
            IF(MQ(IND8+I+1).EQ.0.AND.MQ(IND8+I+3).EQ.(-1)) THEN
               RHS=QM(IND8+I+2)
               SGM=QM(IND8+I+4)
               NCGB=NCGB+1
                I=I+2
            ELSE IF(MATSTO.EQ.2) THEN ! only for sparse storage mode 
               LABEL=MQ(IND8+I+1)
               ITGBI=INONE(LABEL)
               IJ=JVGBI(ITGBI)         ! change to variable parameter
               IF(IJ.GT.0) THEN
                  IF(ISAVES.EQ.0) THEN
                     IND=INTWO(NVGB+NCGB,IJ) ! add index pair
                  ELSE
                     CALL INBITS(NVGB+NCGB,IJ,MREQPE,MQ(INDF+1)) 
                  END IF
               END IF
            END IF
            GOTO 30
         END IF 
      END IF 

*     measurements - determine index-pairs

      IF(MATSTO.EQ.2.AND.IND9.NE.0) THEN  

         I=0
 40      RHS=QM(IND9+I+2) ! right hand side
         SGM=QM(IND9+I+4) ! sigma parameter
         I=I+4
*        loop over label/factor pairs
         IA=I
 50      LABEL=MQ(IND9+I+1)
         FACTR=QM(IND9+I+2)
         ITGBI=INONE(LABEL)
         I=I+2
         IF(I.LT.LND9) THEN
            IF(MQ(IND9+I+1).NE.0) GOTO 50
         END IF
         IB=I-2     
        
         DO J=IA,IB,2
          LABELJ=MQ(IND9+J+1)
          FACTRJ=QM(IND9+J+2)
          ITGBIJ=INONE(LABELJ)
          MQ(IND9+J+1)=ITGBIJ  !!! label -> total parameter index
*         first index
          IVGBIJ=0
          IF(ITGBIJ.NE.0) IVGBIJ=JVGBI(ITGBIJ) ! variable-parameter index
*         second index
          DO K=IA,J,2
           LABELK=MQ(IND9+K+1)  ! already translated
           FACTRK=QM(IND9+K+2)
           ITGBIK=INONE(LABELK) ! already tramslated
           ITGBIK=MQ(IND9+K+1)  ! take total parameter index

           IVGBIK=0
           IF(ITGBIK.NE.0) IVGBIK=JVGBI(ITGBIK) ! variable-parameter index
           IF(IVGBIJ.GT.0.AND.IVGBIK.GT.0) THEN 
              IF(ISAVES.EQ.0) THEN  
                 IND=INTWO(IVGBIJ,IVGBIK)
              ELSE
                 CALL INBITS(IVGBIJ,IVGBIK,MREQPE,MQ(INDF+1))
              END IF
              IF (MPRINT.GT.1)
     +           WRITE(*,*) 'add index pair ',IVGBIJ,IVGBIK
           END IF
          END DO
         END DO
         IF(I.LT.LND9) GOTO 40
      END IF

*     print numbers ----------------------------------------------------

      IF(MATSTO.EQ.2) THEN
         IF(ISAVES.EQ.0) THEN 
            CALL UPTWO        ! final update
            CALL MEGRM('g')   ! remove hash pointer
            WRITE(LUNLOG,*) 'LOOP2:',MQ(INDF-1),
     +                      ' is total number NTGB of labels pairs'
         END IF
      END IF         
      IF (NAGB.GE.65536) THEN
        NOFF=NOFF8/1000
      ELSE
        NOFF=NOFF8
      ENDIF
      NDGN=0
      NSPC=1 ! number of precision types (double, single) for matrix storage
      IF(MATSTO.EQ.2) THEN
         IF(ISAVES.EQ.0) THEN
            NDGN=MQ(INDF-1)    ! actual number of off-diagonal elements
         ELSE 
            IHIS=0
            IF (MHISPE.GT.0) THEN
               IHIS=15
               CALL HMPDEF(IHIS,0.0,FLOAT(MHISPE),
     &           'NDBITS: #off-diagonal elements')
            ENDIF
*
            JCMPRS=MAX(MCMPRS,MSNGPE)
            IF (JCMPRS.GT.0.AND.ISAVES.GT.1) NSPC=2 ! mixed precision storage
            CALL MEGARR('k I',NAGB*NSPC,'INBITS: row compression') 
            LNROW=NAGB+1 
            CALL MEGARL('z L',LNROW*LSCALE*2*NSPC,
     +                  'sparse matrix row pointer')
            MQ(INDZ-1)=NSPC   ! =1: only double, =2: double + single precision
            CALL NDBITS(NDIMSA,MQ(INDF+1),MQ(INDK+1),LQ(INDZ/LSCALE+1),
     &                  MREQPE,IHIS,JCMPRS)
            CALL MEGARR('f I',NDIMSA(1),'INBITS: bit storage') ! resize to bits
            NDGN=NDIMSA(3)+NDIMSA(4) ! actual number of off-diagonal elements 
*
            IF (MHISPE.GT.0) THEN
               IF (NHISTP.NE.0) CALL HMPRNT(IHIS)
               CALL HMPWRT(IHIS)
            ENDIF
         END IF
      END IF 

      NAGBN=MQ(INDC-1) ! max number of global parameters in one event 
      NALCN=MQ(INDC-2) ! max number of local parameters in one event
      NAEQN=MQ(INDC-3) ! max number of equations in one event
      CALL MEGRM('c')
      CALL MEGRM('e')
*     matrices for event matrices
*       split up cache
      IF (FCACHE(2).EQ.0.0) THEN ! from data (DSTAT) 
        FCACHE(1)=SNGL(DSTAT(1))*FCACHE(1) ! leave some part free for fluctuations
        FCACHE(2)=SNGL(DSTAT(2))
        FCACHE(3)=SNGL(DSTAT(3))
      ENDIF
      FSUM=FCACHE(1)+FCACHE(2)+FCACHE(3)
      DO K=1,3
        FCACHE(K)=FCACHE(K)/FSUM
      ENDDO
      NCACHR=IFIX(FLOAT(NCACHE)*FCACHE(1)+0.5) ! read cache
*     define read buffer
      NC21=NCACHR/(21*MTHRDR) ! split read cache 1 : 10 : 10 for pointers, ints, floats
      NWRD=NC21+1 
      CALL MEGARR(': I',(NWRD+5)*MTHRDR,'read buffer, pointer')
      MQ(IDOT2-2)=NWRD*MTHRDR            ! number of buffer words (per thread)
      NWRD=NC21*10+2+NDIMB
      CALL MEGARR('; I',NWRD*MTHRDR,'read buffer, integer')
      MQ(IDOT3  )=MTHRDR          ! number of threads
      MQ(IDOT3-2)=NWRD*MTHRDR     ! number of buffer words 
      CALL MEGARR('! S',NWRD*MTHRDR,'read buffer, float')
      MQ(IDOT4  )=MTHRDR            ! number of (active) threads
      MQ(IDOT4-2)=NWRD*MTHRDR     ! number of buffer words
*
      NCACHI=IFIX(FLOAT(NCACHE)*FCACHE(2)+0.5) ! index cache
      NCACHD=NCACHE-NCACHR-NCACHI              ! data cache
      NGGD=(NAGBN*NAGBN+NAGBN)/2+NCACHD/(2*MTHRD) ! number of double
      MGGI=9                                      ! header length
      NGGI=2+NAGBN+NCACHI/MTHRD                   ! number of ints
      CALL MEGARR('c I',NAGBN*MTHRD,
     &   'global parameters (dim =max/event)')
      MQ(INDC-1)=NAGBN
      CALL MEGARR('e I',NVGB*MTHRD,'global variable-index array')
      MQ(INDE-1)=NVGB
      CALL MEGARR('g D',2*NAGBN*NALCN,    'local/global matrix') 
      CALL MEGARR('h D',2*NGGD*MTHRD,
     &                 'symmetric update matrices')
      MQ(INDH-1)=NGGD                  ! number of words per thread
      MQ(INDH-2)=(NAGBN*NAGBN+NAGBN)/2 ! min free (double) words
      CALL MEGARR('5 I',(NGGI+MGGI)*MTHRD,
     &                 'symmetric update matrix indices')
      MQ(IND5  )=MGGI
      MQ(IND5-1)=NGGI                  ! number of words per thread
      MQ(IND5-2)=NAGBN+2               ! min free words
      CALL MEGARR('i D',2*NAGBN,          'update vector')

*     print all relevant dimension parameters

      DO LU=6,8,2  ! unit 6 and 8

      WRITE(*,*) ' '
      WRITE(LU,101) 'NTGB',NTGB,'total number of parameters'
      WRITE(LU,102) '(all parameters, appearing in binary files)'
      WRITE(LU,101) 'NVGB',NVGB,'number of variable parameters'
      WRITE(LU,102) '(appearing in fit matrix/vectors)'
      WRITE(LU,101) 'NAGB',NAGB,'number of fit parameters'
      WRITE(LU,102) '(including Lagrange multiplier or reduced)'
      WRITE(LU,101) 'MBANDW',MBANDW,'band width of band matrix'
      WRITE(LU,102) '(if =0, no band matrix)'
      IF (NAGB.GE.65536) THEN
       WRITE(LU,101) 'NOFF/K',NOFF,'max number of off-diagonal elements'
      ELSE
       WRITE(LU,101) 'NOFF',NOFF,'max number of off-diagonal elements'
      ENDIF
      IF(NDGN.NE.0) 
     +WRITE(LU,101) 'NDGN',NDGN,'actual number of off-diagonal elements'
      WRITE(LU,101) 'NCGB',NCGB,'number of constraints' 
      WRITE(LU,101)
     +       'NAGBN',NAGBN,'max number of global parameters in an event'
      WRITE(LU,101)
     +       'NALCN',NALCN,'max number of local parameters in an event' 
      WRITE(LU,101)
     +       'NAEQN',NAEQN,'max number of equations in an event' 
      IF (MPRINT.GT.1) THEN
         WRITE(LU,101) 'NAEQNA',NAEQNA,'number of equations'
         WRITE(LU,101) 'NAEQNG',NAEQNG,
     +   'number of equations with       global derivatives'
         WRITE(LU,101) 'NAEQNF',NAEQNF,             
     +   'number of equations with fixed global derivatives'
         WRITE(LU,101) 'NRECF',NRECF,             
     +   'number of records   with fixed global derivatives'
      ENDIF
      IF (NCACHE.GT.0) THEN
         WRITE(LU,101)
     +       'NCACHE',NCACHE,'number of words for caching' 
         WRITE(LU,111) (FCACHE(K)*100.0,K=1,3)
  111    FORMAT(22X,'cache splitting ',3(F6.1,' %'))
      ENDIF
      WRITE(LU,*) ' '
      IF(MATSTO.EQ.2) THEN
c         FRUSE=FLOAT(MQ(INDF-1))/FLOAT(NOFF)
c         WRITE(LU,123) 100.0*FRUSE
      END IF 
c 123  FORMAT(' Fraction of non-zero off-diagonal elements is',
c     +       F6.1,' %')

*     print information about methods and matrix storage modes 
        
      WRITE(LU,*) ' '
      WRITE(LU,*) 'Solution method and matrix-storage mode:' 
      IF(METSOL.EQ.1) THEN
         WRITE(LU,*) '     METSOL = 1:  matrix inversion'
      ELSE IF(METSOL.EQ.2) THEN
         WRITE(LU,*) '     METSOL = 2:  diagonalization'
      ELSE IF(METSOL.EQ.3) THEN
*GF         WRITE(LU,*) '     METSOL = 3:  MINRES'
         WRITE(LU,*) '     METSOL = 3:  MINRES (rtol', MRESTL,')'
      ELSE IF(METSOL.EQ.4) THEN
         WRITE(LU,*) '     METSOL = 4:  GMRES'
      END IF
      WRITE(LU,*) '                  with',MITERA,' iterations'   
      IF(MATSTO.EQ.1) THEN
         WRITE(LU,*) '     MATSTO = 1:  symmetric matrix, ',
     +              '(n*n+n)/2 elements'
      ELSE IF(MATSTO.EQ.2) THEN
         WRITE(LU,*) '     MATSTO = 2:  sparse matrix'
      END IF 
c      END IF
      IF(DFLIM.NE.0.0) THEN
         WRITE(LU,103) 'Convergence assumed, if expected dF <',DFLIM
      END IF

      END DO ! print loop

*     Wolfe conditions

      IF(0.0.LT.WOLFC1.AND.WOLFC1.LT.WOLFC2.AND.WOLFC2.LT.1.0) GOTO 32
      IF(WOLFC1.EQ.0.0) WOLFC1=1.0E-4      
      IF(WOLFC2.EQ.0.0) WOLFC2=0.9
      IF(0.0.LT.WOLFC1.AND.WOLFC1.LT.WOLFC2.AND.WOLFC2.LT.1.0) GOTO 32
      IF(WOLFC1.LE.0.0) WOLFC1=1.0E-4   
      IF(WOLFC2.GE.1.0) WOLFC2=0.9
      IF(WOLFC1.GT.WOLFC2) THEN ! exchange
         WOLFC3=WOLFC1
         WOLFC1=WOLFC2
         WOLFC2=WOLFC3
      ELSE
         WOLFC1=1.0E-4
         WOLFC2=0.9 
      END IF
      WRITE(*,105) WOLFC1,WOLFC2      
      WRITE(LUN,105) WOLFC1,WOLFC2
 105  FORMAT(' Constants C1, C2 for Wolfe conditions:',G12.4,', ',G12.4)

*     prepare matrix and gradient storage ------------------------------
 32   CONTINUE
      MATSIZ(1)=INT8(NAGB)*INT8(NAGB+1) ! number of words for double precision storage 'j'
      MATSIZ(2)=0                       ! number of words for single precision storage '.'
      IF(MATSTO.EQ.2) THEN     ! sparse matrix
         IF(ISAVES.EQ.0) THEN
            CALL MSPARS(NAGB)  ! pointer for sparse storage -> z, d
            MATSIZ(1)=(LNDD+NAGB)*2
         ELSE
            MATSIZ(1)=(NDIMSA(3)+NAGB)*2
            MATSIZ(2)=NDIMSA(4)
            CALL MEGARL('d I',NDIMSA(2),'sparse matrix column list')
            MQ(INDD-1)=NENCDB ! encoding info (number bits for column counter)
            CALL SPBITS(LQ(INDZ/LSCALE+1),MQ(INDD+1),
     +                  MQ(INDF+1),MQ(INDK+1))
            CALL MEGRM('k')
            CALL MEGRM('f')
         END IF  
      END IF

      CALL FEASMA    ! prepare constraint matrices

      CALL VMPREP(MATSIZ)    ! prepare matrix and gradient storage

      DGBM=4.0D-9*DFLOAT(MEGA)      ! total space of Gbytes 
      CALL MEGSPC(NUSE)            ! unused
      DUBM=4.0D-9*DFLOAT(MEGA-NUSE) ! used space in Gbytes 
      WRITE(*,*) ' '
      WRITE(*,104) DGBM,DUBM,DUBM/DGBM*100.0D0
      GBM=SNGL(DGBM)
      UBM=SNGL(DUBM) 

*     print chi^2 cut tables  

      NDFMAX=NAEQN-1
      WRITE(LUNLOG,*) ' '
      WRITE(LUNLOG,*) '   Cut values of Chi^2/Ndf and Chi2,'
      WRITE(LUNLOG,*) '   corresponding to 2 and 3 standard deviations'
      WRITE(LUNLOG,*) '   Ndf  Chi^2/Ndf(2)  Chi^2(2)   ',
     +                 '  Chi^2/Ndf(3)  Chi^2(3)'
      NDF=0
 48   IF(NDF.GT.NAEQN) GOTO 49
      IF(NDF.LT.10) THEN
         NDF=NDF+1
      ELSE IF(NDF.LT.20) THEN
         NDF=NDF+2
      ELSE IF(NDF.LT.100) THEN
         NDF=NDF+5
      ELSE IF(NDF.LT.200) THEN
         NDF=NDF+10
      ELSE
         GOTO 49 
      END IF
      CHIN2=CHINDL(2,NDF)
      CHIN3=CHINDL(3,NDF)
      WRITE(LUNLOG,106) NDF,CHIN2,CHIN2*FLOAT(NDF),CHIN3,
     +                  CHIN3*FLOAT(NDF)
      GOTO 48    
  49  CONTINUE

      WRITE(LUNLOG,*) 'LOOP2: ending'
      WRITE(LUNLOG,*) ' '
      CALL MEND
 101  FORMAT(1X,A6,' =',I10,' = ',A)
 102  FORMAT(22X,A)
 103  FORMAT(1X,A,G12.4)
 104  FORMAT(' Memory space: total',F13.6,' GB'/
     +       '                used',F13.6,' GB = ',F7.2,' %')
 106  FORMAT(I6,2(3X,F9.3,F12.1,3X))
      END 


      SUBROUTINE VMPREP(MSIZE)

*     prepare storage for vectors and matrices

#include "largeint.inc"
#include "dynal.inc"
      SAVE
#include "mpinds.inc"
      INTEGER(KIND=LARGE) NWORDS,MSIZE(2)
*     ...
*                         Vector/matrix storage
c      IF(MATSTO.EQ.2) THEN      ! sparse matrix  
c         CALL MSPARS(NAGB)  ! pointer for sparse storage -> d
c      END IF 
      CALL MEGRM('f')
      CALL MEGARR('f D',2*NAGB*MTHRD,'rhs vector') ! double precision vector 
      MQ(INDF-1)=NAGB
      CALL MEGARR('w S',NAEQN,'residual vector of one record')

      IF(MATSTO.EQ.1) THEN           ! full symmetric matrix  
         CALL MEGARL('j D',MSIZE(1),'full symmetric matrix')
         MQ(INDJ-2)=1                ! 1=normal symmetric matrix
      ELSE IF(MATSTO.EQ.2  ) THEN      ! sparse matrix 
         CALL MEGARL('j D',MSIZE(1),'sparse symmetric matrix (D)')
         MQ(INDJ-2)=2                ! 2=sparse matrix
         IF (MSIZE(2).GT.0) THEN
            CALL MEGARL('. F',MSIZE(2),'sparse symmetric matrix (F)')
            MQ(INDJ-3)=1
         ENDIF
      END IF      


      IF(METSOL.GE.3) THEN                  ! GMRES/MINRES algorithms
*        array space is:
*           variable-width band matrix or diagonal matrix for parameters 
*           followed by rectangular matrix for constraints
*           followed by symmetric matrix for constraints   
         IF(MBANDW.GT.0) THEN               ! variable-width band matrix
            CALL MEGARR('u I',NAGB,'pointer-array variable-band matrix')
            DO I=1,MIN(MBANDW,NVGB)
             MQ(INDU+I)=(I*I+I)/2           ! increasing number
            END DO
            DO I=MIN(MBANDW,NVGB)+1,NVGB
             MQ(INDU+I)=MQ(INDU+I-1)+MBANDW ! fixed band width
            END DO
            DO I=NVGB+1,NAGB                ! reset
             MQ(INDU+I)=0
            END DO
            NTOTAL=MQ(INDU+NVGB)+NCGB*NVGB+(NCGB*NCGB+NCGB)/2
            CALL MEGARR('v D',2*NTOTAL,'variable-band matrix')
         ELSE                               ! default preconditioner
            NTOTAL=NVGB+NCGB*NVGB+(NCGB*NCGB+NCGB)/2
            CALL MEGARR('v D',2*NTOTAL,'default preconditioner matrix')
         END IF
      END IF   



      CALL MEGARR('s D',2*NAGB,'corrections')      ! double prec corrections

      CALL MEGARR('m D',2*NAGB,'auxiliary array')  ! double aux 2
      CALL MEGARR('k S',  NAGB,'auxiliary array')  ! single aux 1
      CALL MEGARR('r D',2*NAGB,'auxiliary array')  ! double aux 7

      IF(METSOL.EQ.1) THEN
         CALL MEGARR('l D',2*NAGB,'auxiliary array')  ! double aux 1
C         CALL MEGARR('t D',2*NAGB,'auxiliary array')  ! double aux 8
      END IF

      IF(METSOL.EQ.2) THEN
         CALL MEGARR('l D',2*NAGB,'auxiliary array')  ! double aux 1
         CALL MEGARR('n D',2*NAGB,'auxiliary array')  ! double aux 3
         CALL MEGARR('q D',2*NAGB,'auxiliary array')  ! double aux 6
         CALL MEGARR('t D',2*NAGB*NAGB,'matrix U')    ! rotation matrix
      END IF 

      IF(METSOL.GE.3) THEN
         CALL MEGSPC(NWORDS)
         NMATS=NWORDS/(2*NAGB)
         WRITE(*,*) 'Space: double words left ',NWORDS/2,NMATS  
         IF(NMATS.LT.6) THEN
            WRITE(*,*) 'Insufficient space - stop'
            STOP 'recompile with larger memory'
         END IF 
         IF(METSOL.EQ.3) THEN
            CALL MEGARR('l D',6*2*NAGB,'auxiliary array')  ! double aux 1
         ELSE IF(METSOL.EQ.4) THEN
            CALL MEGARR('l D',NMATS*2*NAGB,'auxiliary array') ! double aux 1  
c         CALL MEGARR('n D',2*NAGB,'auxiliary array')  ! double aux 3
c         CALL MEGARR('o D',2*NAGB,'auxiliary array')  ! double aux 4
c         CALL MEGARR('p D',2*NAGB,'auxiliary array')  ! double aux 5
c         CALL MEGARR('q D',2*NAGB,'auxiliary array')  ! double aux 6 
c         CALL MEGARR('t D',2*NAGB,'auxiliary array')  ! double aux 8
         END IF
      END IF 

      END


      SUBROUTINE MINVER 
*     solution by matrix inversion
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      SAVE
*     ...
      LUN=LUNLOG                       ! log file
      IF(LUNLOG.EQ.0) LUNLOG=6

c      WRITE(*,*) 'MINVER ICALCM=',ICALCM
      IF(ICALCM.EQ.1) THEN 
         CALL SQMINV(DQ(IGMAT/2+1),
     +               DQ(ISOLV/2+1),NAGB,NRANK,
     +               DQ(IDUX1/2+1),MQ(IAUX1+1))
         NDEFEC=NAGB-NRANK   ! rank defect
         IF(NDEFEC.NE.0) THEN
            WRITE(*,*)   'The rank defect of the symmetric',NAGB,
     +      '-by-',NAGB,' matrix is ',NDEFEC,' (should be zero).'
            WRITE(LUN,*) 'The rank defect of the symmetric',NAGB,
     +      '-by-',NAGB,' matrix is ',NDEFEC,' (should be zero).'
         ELSE
            WRITE(LUN,*) 'No rank defect of the symmetric matrix'
         END IF 
      ELSE             ! multiply gradient by inverse matrix
         CALL DBSVX(DQ(IGMAT/2+1),DQ(IGVEC/2+1),DQ(ISOLV/2+1),NAGB)
      END IF
      END 


      SUBROUTINE MDIAGS
*     solution by diagonalization
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      SAVE
*     ...


      LUN=LUNLOG                       ! log file
      IF(LUNLOG.EQ.0) LUN=6

      IF(ICALCM.EQ.1) THEN
         NVAR=NAGB
         DO I=1,NVAR      ! used in FEASIB
          DQ(IDUX1/2+I)=DQ(IGMAT/2+(I*I+I)/2) ! save diagonal elements
         END DO

*        IDUX1 = original diagonal elements
*        IDUX8 = matrix of eigen vectors
*        IDUX6 = eigenvalues (ex 2)
*        IDUX3 = work array

*                         eigenvalues   eigenvectors   symm_input  
         CALL DEVROT(NVAR,DQ(IDUX6/2+1),DQ(IDUX8/2+1),DQ(IGMAT/2+1),
     +                    DQ(IDUX3/2+1),MQ(IAUX1+1)) 
*                         work          iwork

*        ---------------------------------------------------------------

*        histogram of positive eigenvalues

         NMAX=INT(1.0+LOG10(SNGL(DQ(IDUX6/2+1)))) ! > log of largest eigenvalue
         DO I=NAGB,1,-1
          IF(DQ(IDUX6/2+I).GT.0.0D0) GOTO 20
         END DO
         I=1
  20     NMIN=INT(LOG10(SNGL(DQ(IDUX6/2+I))))   ! log of smallest pos. eigenvalue 
         IMIN=I                      ! index of smallest pos. eigenvalue
         NTOP=NMIN+6
  21     IF(NTOP.LT.NMAX) THEN
            NTOP=NTOP+3
            GOTO 21
         END IF

         CALL HMPDEF(7,FLOAT(NMIN),FLOAT(NTOP),
     +              'log10 of positive eigenvalues')
         DO IDIA=1,NAGB
          IF(DQ(IDUX6/2+IDIA).GT.0.0D0) THEN ! positive 
             EVALUE=LOG10(SNGL(DQ(IDUX6/2+IDIA)))
             CALL HMPENT(7,EVALUE)
          END IF 
         END DO
         IF(NHISTP.NE.0) CALL HMPRNT(7)
         CALL HMPWRT(7) 

         IAST=MAX(1,IMIN-60)
         CALL GMPDEF(3,2,'low-value end of eigenvalues')
         DO I=IAST,NAGB  
          EVALUE=SNGL(DQ(IDUX6/2+I))
          CALL GMPXY(3,FLOAT(I),EVALUE) 
         END DO
         IF(NHISTP.NE.0) CALL GMPRNT(3)
         CALL GMPWRT(3)
      
         DO I=1,NVAR
          DQ(IDUX3/2+I)=0.0D0
          IF(DQ(IDUX6/2+I).NE.0.0D0) THEN
             DQ(IDUX3/2+I)=MAX(0.0D0,LOG10(ABS(DQ(IDUX6/2+I)))+3.0D0)
             IF(DQ(IDUX6/2+I).LT.0.0D0) DQ(IDUX3/2+I)=-DQ(IDUX3/2+I)
          END IF
         END DO
         WRITE(LUN,*) ' '
         WRITE(LUN,*) 'The first (largest) eigenvalues ...'
         WRITE(LUN,102) (DQ(IDUX6/2+I),I=1,MIN(20,NAGB))
         WRITE(LUN,*) ' '
         WRITE(LUN,*) 'The last eigenvalues ... up to',NVGB
         WRITE(LUN,102) (DQ(IDUX6/2+I),I=MAX(1,NVGB-19),NVGB)
         WRITE(LUN,*) ' '
         WRITE(LUN,*) 'The eigenvalues from',NVGB+1,' to',NAGB
         WRITE(LUN,102) (DQ(IDUX6/2+I),I=NVGB+1,NAGB)         
         WRITE(LUN,*) ' '
         WRITE(LUN,*) 'Log10 + 3 of ',NAGB,' eigenvalues in decreasing',
     +                ' order'
         WRITE(LUN,*) '(for Eigenvalue < 0.001 the value 0.0 is shown)'
         WRITE(LUN,101) (DQ(IDUX3/2+I),I=1,NAGB)
         IF(DQ(IDUX3/2+NVAR).LT.0) WRITE(LUN,*) 'Negative values are ',
     +      'printed for negative eigenvalues'
         CALL DEVSIG(NAGB,DQ(IDUX6/2+1),DQ(IDUX8/2+1),
     +                    DQ(IGVEC/2+1),DQ(IDUX3/2+1))
         WRITE(LUN,*) ' '
         WRITE(LUN,*) NVGB,' significances: insignificant if ',
     +               'compatible with  N(0,1)'
         WRITE(LUN,101) (DQ(IDUX3/2+I),I=1,NVGB)


 101     FORMAT(10F7.1) 
 102     FORMAT(5E14.6)

#ifdef graphics

      DO I=1,NAGB
       QM(IDUX3+I)=DQ(IDUX3/2+I)
      END DO
      CALL PFINIT
      CALL PTITLE(0,'significances')
      CALL PLATT('blue yellow narrow histogram')
      CALL PDIHST(NAGB,1,NAGB,QM(IDUX3+1))
      CALL PLATT('red full line 2.0')
      CALL PD2XY(0.0,1.0,FLOAT(NAGB+1),1.0,'l')
      CALL PLATT('red full line 1.0')
      CALL PD2XY(0.0,3.0,FLOAT(NAGB+1),3.0,'l') 

      CALL PLSHOW


      DO I=1,NAGB 
       QM(IDUX3+I)=DQ(IDUX6/2+I)
      END DO
      CALL PFINIT
      CALL PTITLE(0,'eigenvalues')
      CALL PLATT('blue yellow narrow histogram') 
      CALL PDIHST(NAGB,1,NAGB,QM(IDUX3+1))
      CALL PLSHOW

      CALL PFINIT
      CALL PTITLE(0,'end of eigenvalue spectrum')
      CALL PLATT('blue yellow narrow histogram')
      CALL PDIHST(50,NAGB-49,NAGB,QM(IDUX3+NAGB-49))
      CALL PLSHOW

      CALL PFINIT
      CALL PTITLE(0,'end of eigenvalue spectrum')
      CALL PLATT('blue vertl histogram')
      CALL PDIHST(50,NAGB-49,NAGB,QM(IDUX3+NAGB-49))
      CALL PLSHOW

#endif

      END IF

*     solution ---------------------------------------------------------

*                      eigenvalues   eigenvectors
      CALL DEVSOL(NVAR,DQ(IDUX6/2+1),DQ(IDUX8/2+1),
     +                 DQ(IGVEC/2+1),DQ(ISOLV/2+1),
     +                 DQ(IDUX3/2+1))   
      RETURN

      ENTRY ZDIAGS ! covariance matrix for diagonalization 
*                      eigenvalue    eigenvectors  cov.matrix
      CALL DEVINV(NVAR,DQ(IDUX6/2+1),DQ(IDUX8/2+1),DQ(IGMAT/2+1))  ! inv

      END 

      SUBROUTINE MMINRS  
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      DOUBLE PRECISION SHIFT,RTOL,ANORM,ACOND,RNORM,YNORM
      LOGICAL CHECKA  
      EXTERNAL AVPROD, MVSOLV, MCSOLV
      SAVE
*     ...
      LUN=LUNLOG                       ! log file
      IF(LUNLOG.EQ.0) LUN=6

      NOUT=LUN
      ITNLIM=2000    ! iteration limit
      SHIFT =0.0D0   ! not used
      RTOL = MRESTL ! from steering
      NVAR=NAGB
      CHECKA=.FALSE.

c      WRITE(*,*) 'METSOL MATSTO',METSOL,MATSTO

c      WRITE(*,*) 'Input vector'
c      WRITE(*,*) (DQ(IGVEC/2+I),I=1,NAGB)

      IF(MBANDW.EQ.0) THEN           ! default preconditioner
         IF(ICALCM.EQ.1) THEN
c            WRITE(LUN,*) 'MMINRS: PRECON started'
            CALL PRECON(NCGB,NVGB,DQ(INDV/2+1),DQ(INDV/2+1),
     +                  DQ(INDV/2+1+NVGB),
     +                  DQ(INDV/2+1+NVGB+NCGB*NVGB))
c            WRITE(LUN,*) 'MMINRS: PRECON ended'
         END IF

c         WRITE(*,*) 'NAGB MBANDW ',NAGB,MBANDW
c         WRITE(*,*) IGVEC,IDUX1,ISOLV,IDUX7  

         CALL MINRES(NAGB,DQ(IGVEC/2+1),
     +   DQ(IDUX1/2+1),       DQ(IDUX1/2+1+NAGB),  DQ(IDUX1/2+1+2*NAGB),
     +   DQ(IDUX1/2+1+3*NAGB),DQ(IDUX1/2+1+4*NAGB),DQ(IDUX1/2+1+5*NAGB),        
     +                 DQ(ISOLV/2+1),DQ(IDUX7/2+1),
     +                 AVPROD, MCSOLV, CHECKA ,.TRUE. , SHIFT,
     +                 NOUT , ITNLIM, RTOL,
     +                 ISTOP, ITN, ANORM, ACOND, RNORM, YNORM)
C         WRITE(*,*) 'MINRES ended'   
C    +                 DQ(IDUX1/2+1),DQ(IDUX8/2+1),DQ(IDUX3/2+1),
C    +                 DQ(IDUX4/2+1),DQ(IDUX5/2+1),DQ(IDUX6/2+1),
      ELSE IF(MBANDW.GT.0) THEN                          ! band matrix preconditioner 
         IF(ICALCM.EQ.1) THEN
            WRITE(LUN,*) 'MMINRS: EQUDEC started'
            CALL EQUDEC(NVGB,NCGB,DQ(INDV/2+1),MQ(INDU+1),NRKD,NRKD2)
            WRITE(LUN,*) 'MMINRS: EQUDEC ended'
         END IF 
         CALL MINRES(NAGB,DQ(IGVEC/2+1),
     +   DQ(IDUX1/2+1),       DQ(IDUX1/2+1+NAGB),  DQ(IDUX1/2+1+2*NAGB),
     +   DQ(IDUX1/2+1+3*NAGB),DQ(IDUX1/2+1+4*NAGB),DQ(IDUX1/2+1+5*NAGB),
     +                 DQ(ISOLV/2+1),DQ(IDUX7/2+1),
     +                 AVPROD, MVSOLV, CHECKA ,.TRUE. , SHIFT,
     +                 NOUT,   ITNLIM, RTOL,
     +                 ISTOP,  ITN, ANORM, ACOND, RNORM, YNORM)
      ELSE
c         IF(ICALCM.NE.0) WRITE(LUN,*) 'MMINRS: no preconditioning !!!'
         CALL MINRES(NAGB,DQ(IGVEC/2+1),
     +   DQ(IDUX1/2+1),       DQ(IDUX1/2+1+NAGB),  DQ(IDUX1/2+1+2*NAGB),
     +   DQ(IDUX1/2+1+3*NAGB),DQ(IDUX1/2+1+4*NAGB),DQ(IDUX1/2+1+5*NAGB),
     +                 DQ(ISOLV/2+1),DQ(IDUX7/2+1),
     +                 AVPROD, MVSOLV, CHECKA ,.FALSE., SHIFT,
     +                 NOUT,   ITNLIM, RTOL,
     +                 ISTOP,  ITN, ANORM, ACOND, RNORM, YNORM)
      END IF 
      IITERA=ITN
      ISTOPA=ISTOP
      MNRSIT=MNRSIT+ITN

c      WRITE(*,*) 'Solution vector'
c      WRITE(*,*) (DQ(ISOLV/2+I),I=1,NAGB)

      END 


      SUBROUTINE MCSOLV(N,X,Y)         !  solve M y = x
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      DOUBLE PRECISION X(N),Y(N)
      SAVE
*     ...
      CALL PRESOL(NCGB,NVGB,DQ(INDV/2+1),DQ(INDV/2+1+NVGB),
     +            DQ(INDV/2+1+NVGB+NCGB*NVGB),Y,X)
      END

      SUBROUTINE MVSOLV(N,X,Y)         !  solve M y = x
#include "largeint.inc"
#include "dynal.inc"
#include "mpinds.inc"
      DOUBLE PRECISION X(N),Y(N)
      SAVE
*     ...
      DO I=1,N
       Y(I)=X(I)                    ! copy to output vector
      END DO

      IF(IND8.EQ.0) GOTO 90         ! no constraints
      DO ICGB=1,NCGB                ! copy previous lambda values
c       Y(NVGB+ICGB)=DQ(ISOLV/2+NVGB+ICGB)
      END DO

      I=0
      DO ICGB=1,NCGB
       I=I+4
 10    LABEL=MQ(IND8+I+1)
       FACTR=QM(IND8+I+2)
       ITGBI=INONE(LABEL) ! -> ITGBI= index of parameter label
       IVGB =JVGBI(ITGBI) ! -> variable-parameter index
*                   - A^T * lambda
       IF(IVGB.NE.0) THEN
c          Y(IVGB)=Y(IVGB) - DQ(ISOLV/2+NVGB+ICGB)*FACTR
       END IF
       I=I+2
       IF(I.LT.LND8) THEN
          IF(MQ(IND8+I+1).NE.0) GOTO 10
       END IF
      END DO

 90   CONTINUE
c 90   CALL VABSLV(NVGB,DQ(INDV/2+1),MQ(INDU+1),Y) ! solve for Y
C      CALL VABSLV(NAGB,DQ(INDV/2+1),MQ(INDU+1),Y) ! solve for Y

      CALL EQUSLV(NVGB,NCGB,DQ(INDV/2+1),MQ(INDU+1),Y)
      END



************************************************************************


      SUBROUTINE XLOOPN                ! standard solution algorithm
#include "largeint.inc"
#include "dynal.inc"
      DOUBLE PRECISION STP,DRATIO,DFACIN,DJRAT,DWMEAN
      DOUBLE PRECISION DB,DB1,DB2,DBDOT
      LOGICAL WARNER 
      SAVE
#include "mpinds.inc"
*     ...

*     Printout of algorithm for solution and important parameters ------  

      LUN=LUNLOG                       ! log file 
      IF(LUNLOG.EQ.0) LUNLOG=6  

      DO LUNP=6,LUNLOG,LUNLOG-6
       WRITE(LUNP,*) ' '
       WRITE(LUNP,*) 'Solution algorithm: '
       WRITE(LUNP,121)
     +   '=================================================== '

       IF(METSOL.EQ.1) THEN
          WRITE(LUNP,121) 'solution method:','matrix inversion'
       ELSE IF(METSOL.EQ.2) THEN
          WRITE(LUNP,121) 'solution method:','diagonalization'
       ELSE
          IF(METSOL.EQ.3) THEN
             WRITE(LUNP,121) 'solution method:',
     +                       'minres (Paige/Saunders)'
          ELSE IF(METSOL.EQ.4) THEN
             WRITE(LUNP,121) 'solution method:',
     +       'gmres (generalized minimzation of residuals)'         
          END IF
       END IF
       WRITE(LUNP,123) 'convergence limit at Delta F=',DFLIM
       WRITE(LUNP,122) 'maximum number of iterations=',MITERA
       MATRIT=MIN(MATRIT,MITERA) 
       IF(MATRIT.GT.1) THEN
          WRITE(LUNP,122) 'matrix recalculation up to ',MATRIT,
     +    '. iteration'
       END IF 
       IF(METSOL.GE.3) THEN 
          IF(MATSTO.EQ.1) THEN
             WRITE(LUNP,121) 'matrix storage:','full'
          ELSE IF(MATSTO.EQ.2) THEN
             WRITE(LUNP,121) 'matrix storage:','sparse'
          END IF
          WRITE(LUNP,122) 'pre-con band-width parameter=',MBANDW
          IF(MBANDW.EQ.0) THEN
             WRITE(LUNP,121) 'pre-conditioning:','default'
          ELSE IF(MBANDW.LT.0) THEN
             WRITE(LUNP,121) 'pre-conditioning:','none!' 
          ELSE IF(MBANDW.GT.0) THEN
             WRITE(LUNP,121) 'pre-conditioning=','band-matrix'
          END IF
       END IF
       IF(REGPRE.EQ.0.0D0.AND.NPRESG.EQ.0) THEN
          WRITE(LUNP,121) 'using pre-sigmas:','no'
       ELSE
C FIXME: NPRESG contains parameters that failed the 'entries' cut...
          WRITE(LUNP,124) 'pre-sigmas defined for',
     +       FLOAT(100*NPRESG)/FLOAT(NVGB),' % of variable parameters' 
          WRITE(LUNP,123) 'default pre-sigma=',REGPRE
       END IF 
       IF(NREGUL.EQ.0) THEN
          WRITE(LUNP,121) 'regularization:','no'
       ELSE 
          WRITE(LUNP,121) 'regularization:','yes'
          WRITE(LUNP,123) 'regularization factor=',REGULA
       END IF   

      IF(CHICUT.NE.0.0) THEN
         WRITE(LUNP,121) 'Chi square cut equiv 3 st.dev applied'
         WRITE(LUNP,123) '... in first iteration with factor',CHICUT
         WRITE(LUNP,123) '... in second iteration with factor',CHIREM
         WRITE(LUNP,121) ' (reduced by sqrt in next iterations)'
      END IF

      IF(LHUBER.NE.0) THEN
         WRITE(LUNP,122) 'Down-weighting of outliers in',
     +              LHUBER,' iterations'
         WRITE(LUNP,123) 'Cut on downweight fraction',DWCUT
      END IF


 121  FORMAT(1X,A40,3X,A)
 122  FORMAT(1X,A40,2X,I2,A)
 123  FORMAT(1X,A40,2X,E9.2)
 124  FORMAT(1X,A40,3X,F5.1,A)
      END DO

*     initialization of iterations -------------------------------------

      NTERAT=MITERA             ! max number of iterations
      IITERA=0
      NLOOP =0
      NSOL  =0                  ! counter for solutions  
      INFO  =0
      LSINFO=0
      STP   =0.0D0
      STEPL =SNGL(STP)
      CONCUT=1.0E-12            ! initial constraint accuracy
      CONCU2=1.0E-06            ! constraint accuracy
      ICALCM=1                  ! require matrix calculation
      ITERAT=0                  ! iteration counter
      ITERAT=-1
      LITERA=-2
      NREJ=0                    ! reset number of rejects 
      IF(METSOL.EQ.1) THEN
         WOLFC2=0.5             ! not accurate
         MINF=1
      ELSE IF(METSOL.EQ.2) THEN
         WOLFC2=0.5             ! not acurate
         MINF=2  
      ELSE IF(METSOL.EQ.3) THEN
         WOLFC2=0.1             ! accurate
         MINF=3
      ELSE IF(METSOL.EQ.4) THEN
         WOLFC2=0.1             ! accurate 
         MINF=3
      END IF

*     check initial feasibility of constraint equations ----------------

      WRITE(*,*) ' '
      IF(NOFEAS.EQ.0) THEN        ! make parameter feasible
         WRITE(LUNLOG,*) 'Checking feasibility of parameters:'
         WRITE(*,*) 'Checking feasibility of parameters:'
         CALL FEASIB(CONCUT,IACT) ! check feasibility
         IF(IACT.NE.0) THEN       ! done ...  
            WRITE(*,102) CONCUT
            WRITE(*,*) '   parameters are made feasible'
            WRITE(LUNLOG,102) CONCUT
            WRITE(LUNLOG,*) '   parameters are made feasible'
         ELSE                     ! ... was OK
            WRITE(*,*)
     +      '   parameters are feasible  (i.e. satisfy constraints)'
            WRITE(LUNLOG,*)
     +      '   parameters are feasible  (i.e. satisfy constraints)' 
         END IF 
         CONCUT=CONCU2            ! cut for constraint check
      END IF
      IACT=1                      ! set flag for new data loop  
      NOFEAS=0                    ! set check-feasibility flag

      WRITE(*,*) ' '
      WRITE(*,*)'Reading files and accumulating vectors/matrices ...'
      WRITE(*,*) ' '


*     ICALCM = +1   M, g, F  &  solution 
*     ICALCM =  0   g, F
*     ICALCM = -1   solution
*     ICALCM = -2   end

      CALL ETIME(TA,RSTART) 
      ITERAT=-1
      LITERA= 0
      JCALCM=-1
      IAGAIN= 0 

 09   ICALCM=1

*     Block 1: data loop with vector (and matrix) calculation ---------- 

 10   CONTINUE
      IF(ITERAT.GE.0) THEN
         LCALCM=JCALCM+3 ! mode (1..4) of last loop
         IF(JCALCM+1.NE.0) THEN
            IF(ITERAT.EQ.0) THEN
               CALL PLOOPA(6) ! header
               CALL PLOOPB(6)
               CALL PLOOPA(LUNLOG) ! iteration line
               CALL PLOOPB(LUNLOG)
               ITERAT=1
               CALL GMPXYD(1,FLOAT(NLOOPN),REAL(FVALUE),0.5,0.) ! fcn-value graph (no Delta)
            ELSE
               IF(ITERAT.NE.LITERA) THEN 
                  CALL PLOOPB(6)
C                  CALL PLOOPA(LUNLOG)
                  CALL PLOOPB(LUNLOG)
                  LITERA=ITERAT
                  CALL GMPXYD(1,FLOAT(NLOOPN),REAL(FVALUE),0.5,DELFUN) ! fcn-value (with expected)
                  IF(METSOL.EQ.3) THEN ! extend to 4, i.e. GMRES?
                     CALL GMPXY(2,FLOAT(ITERAT),FLOAT(IITERA)) ! MINRES iterations
                  END IF
               ELSE
                  CALL PLOOPC(6) ! sub-iteration line
                  CALL PLOOPC(LUNLOG)
                  CALL GMPXYD(1,FLOAT(NLOOPN),REAL(FVALUE),0.5,0.) ! fcn-value graph (no Delta)
               END IF 
            END IF       
         ELSE
            CALL PLOOPD(6) ! solution line
            CALL PLOOPD(LUNLOG)      
         END IF
         CALL ETIME(TA,RSTART)
* CHK 
         IF (IABS(JCALCM).LE.1) THEN     
            IDX=JCALCM+4
            TIMES(IDX  )=(TIMES(IDX  )*TIMES(IDX+3)+DELTIM)
     +                  /(TIMES(IDX+3)+1.0)
            TIMES(IDX+3)= TIMES(IDX+3)+1.0
         ENDIF
      ENDIF    
      JCALCM=ICALCM   
*  
      IF(ICALCM.LT.0) GOTO 20          ! ICALCM = +1 & 0 
      CALL LOOPN                       ! data loop
      CALL ADDCST                      ! constraints
      LREJ=NREJ
      NREJ=NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3) ! total number of rejects
      IF(3*NREJ.GT.NRECAL) THEN
         WRITE(*,*) ' '
         WRITE(*,*) 'Data rejected in previous loop:   '
         WRITE(*,*) '   ',
     +   NREJEC(0), ' (rank deficit/NaN) ',NREJEC(1),' (Ndf=0)   ',
     +   NREJEC(2), ' (huge)   ',NREJEC(3),' (large)'    
         WRITE(*,*) 'Too many rejects (>33.3%) - stop'
         STOP 
      END IF 

*     Block 2: new iteration with calculation of solution --------------

 20   IF(ABS(ICALCM).NE.1) GOTO 30     ! ICALCM = +1 & -1
      DO I=1,NAGB
       DQ(ISOLV/2+I)=DQ(IGVEC/2+I)     ! copy rhs
      END DO
      DO I=1,NVGB
       ITGBI=JTGBI(I)
       DQ(IDUX2/2+I)=DQ(IND0/2+ITGBI)  ! copy X for line search
      END DO
      ITERAT=ITERAT+1                  ! increase iteration count
      IF(METSOL.EQ.1) THEN
         CALL MINVER                   ! inversion
      ELSE IF(METSOL.EQ.2) THEN
         CALL MDIAGS                   ! diagonalization
      ELSE IF(METSOL.EQ.3) THEN
         CALL MMINRS                   ! MINRES
      ELSE IF(METSOL.EQ.4) THEN
         WRITE(*,*) '... reserved for GMRES (not yet!)'
         CALL MMINRS                   ! GMRES not yet
      END IF
      

*     check feasibility and evtl. make step vector feasible

      DO I=1,NVGB
       ITGBI=JTGBI(I)
       DQ(IND4/2+ITGBI)=DQ(IND0/2+ITGBI)               ! save
       DQ(IND0/2+ITGBI)=DQ(IND0/2+ITGBI)+DQ(ISOLV/2+I) ! update
      END DO
      CALL FEASIB(CONCUT,IACT)  ! improve constraints
      CONCUT=CONCU2             ! new cut for constraint check
      DO I=1,NVGB
       ITGBI=JTGBI(I)
       DQ(ISOLV/2+I)=DQ(IND0/2+ITGBI)-DQ(IND4/2+ITGBI) ! feasible stp
       DQ(IND0/2+ITGBI)=DQ(IND4/2+ITGBI)               ! restore
      END DO

*     initialize line search based on slopes and prepare next

      CALL PTLDEF(WOLFC2, 10.0, MINF,10)
      IF(METSOL.EQ.1) THEN
         WOLFC2=0.5            ! not accurate
         MINF=3
      ELSE IF(METSOL.EQ.2) THEN
         WOLFC2=0.5            ! not acurate
         MINF=3
      ELSE IF(METSOL.EQ.3) THEN
         WOLFC2=0.1            ! accurate
         MINF=4
      ELSE IF(METSOL.EQ.4) THEN
         WOLFC2=0.1            ! accurate
         MINF=4
      END IF
      DB=DBDOT(NVGB,DQ(ISOLV/2+1),DQ(IGVEC/2+1))
      DB1=DBDOT(NVGB,DQ(ISOLV/2+1),DQ(ISOLV/2+1))
      DB2=DBDOT(NVGB,DQ(IGVEC/2+1),DQ(IGVEC/2+1))
      DELFUN=SNGL(DB)
      ANGRAS=SNGL(DB/DSQRT(DB1*DB2))

      IF(DB.LE.0.0D0) THEN 
         WRITE(*,*) 'Function not decreasing:',DB
         IF(DB.LE.-1.0D-3) THEN ! 100311, VB/CK: allow some margin for numerics
            IAGAIN=IAGAIN+1
            IF (IAGAIN.LE.1) THEN
               WRITE(*,*) '... again matrix calculation' 
               GOTO 09
            ELSE
               WRITE(*,*) '... aborting iterations' 
               GOTO 90
            ENDIF
         ELSE
            WRITE(*,*) '... stopping iterations' 
            IAGAIN=0
            GOTO 90
         ENDIF
      ELSE
         IAGAIN=0 
      END IF
      ICALCM=0                         ! switch 

*     Block 3: line searching ------------------------------------------

 30   CONTINUE
      IF(ICALCM+2.EQ.0) GOTO 91 
      CALL PTLINE(NAGB,DQ(IDUX2/2+1),    ! current parameter values
     +            FLINES,           ! chi^2 dunction value
     +            DQ(IGVEC/2+1),    ! gradient
     +            DQ(ISOLV/2+1),    ! step vector
     +            STP,              ! returned step factor
     +            INFO)             ! returned information
C      WRITE(*,*) 'PTLINE returns INFO, STP=',INFO, STP
      LSINFO=INFO

      STEPL=SNGL(STP)
      NAN=0
      DO I=1,NVGB
       ITGBI=JTGBI(I)
       IF ((.NOT.(DQ(IDUX2/2+I).LE.0.0D0)).AND.
     +     (.NOT.(DQ(IDUX2/2+I).GT.0.0D0))) NAN=NAN+1
       DQ(IND0/2+ITGBI)=DQ(IDUX2/2+I) ! current parameter values
      END DO
      
      IF (NAN.GT.0) THEN
         WRITE(*,*) 'Result vector containes ', NAN,' NaNs - stop'
         STOP
      ENDIF
      
*     subito exit, if required -----------------------------------------

      IF(ISUBIT.NE.0) THEN ! subito
         WRITE(*,*) 'Subito!     Exit after first step.'
         GOTO 90
      END IF
      
      IF(INFO.EQ.0) THEN
         WRITE(*,*) 'INFO=0 should not happen (line search input err)'
         IF (IAGAIN.LE.0) GOTO 09
      END IF
      IF(INFO.LT.0) GOTO 10

*     Block 4: line search convergence ---------------------------------

      CALL PTLPRT
      CALL FEASIB(CONCUT,IACT)               ! check constraints
      IF(IACT.NE.0.OR.CHICUT.GT.1.0) THEN
         ICALCM=-1
         IF(ITERAT.LT.MATRIT) ICALCM=+1
         GOTO 10 ! iterate
      END IF 
      IF(DELFUN.LE.DFLIM)  GOTO 90           ! convergence
      IF(ITERAT.GE.MITERA) GOTO 90           ! ending 
      ICALCM=-1
      IF(ITERAT.LT.MATRIT) ICALCM=+1 
      GOTO 10                                ! next iteration

*     Block 5: iteration ending ----------------------------------------

 90   ICALCM=-2
      GOTO 10
 91   IF(NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3).NE.0) THEN
         WRITE(*,*) ' '
         WRITE(*,*) 'Data rejected in last loop:   '
         WRITE(*,*) '   ',
     +   NREJEC(0), ' (rank deficit/NaN) ',NREJEC(1),' (Ndf=0)   ',
     +   NREJEC(2), ' (huge)   ',NREJEC(3),' (large)'
      END IF

      DWMEAN=SUMNDF/DFLOAT(NDFSUM)
      DRATIO=FVALUE/DWMEAN/DFLOAT(NDFSUM-NAGB)
      CATIO=SNGL(DRATIO)
      IF(LHUBER.NE.0) CATIO=1.07*CATIO  ! correction   
      MRATI=NINT(100.0*CATIO)

      DO LUNP=6,LUNLOG,LUNLOG-6
       WRITE(LUNP,*) ' '
       IF (NFILW.LE.0) THEN
          WRITE(LUNP,*) 'Sum(Chi^2)/Sum(Ndf) =',FVALUE
          WRITE(LUNP,*) '                    / (',NDFSUM,'-',NAGB,')'
          WRITE(LUNP,*) '                    =',DRATIO
       ELSE
          WRITE(LUNP,*) 'Sum(W*Chi^2)/Sum(Ndf)/<W> =',FVALUE
          WRITE(LUNP,*) '                          / (',NDFSUM,'-',
     +                                                  NAGB,')'
          WRITE(LUNP,*) '                          /',DWMEAN
          WRITE(LUNP,*) '                          =',DRATIO
       ENDIF
       WRITE(LUNP,*) ' '
       IF(LHUBER.NE.0) WRITE(LUNP,*)
     + '            with correction for down-weighting   ',CATIO     
      END DO
      NREJ=NREJEC(0)+NREJEC(1)+NREJEC(2)+NREJEC(3) ! total number of rejects

*     ... the end with exit code ???????????????????????????????????????
C
C      WRITE(*,199)              ! write exit code
C     + '-----------------------------------------------------------'
C      IF(ITEXIT.EQ.0) WRITE(*,199)
C     +  'Exit code = 0:  Convergence reached'
C      IF(ITEXIT.EQ.1) WRITE(*,199)
C     +  'Exit code = 1:  No improvement in last iteration'
C      IF(ITEXIT.EQ.2) WRITE(*,199)
C     +  'Exit code = 2:  Maximum number of iterations reached'
C      IF(ITEXIT.EQ.3) WRITE(*,199)
C     +  'Exit code = 3:  Failure'
C      WRITE(*,199)
C     + '-----------------------------------------------------------'
C      WRITE(*,199) ' '


      NRATI=NINT(FLOAT(10000*NREJ)/FLOAT(NRECAL))
      DJRAT=0.01D0*DFLOAT(NRATI)
      NFACI=NINT(100.0*SQRT(CATIO))

      DRATIO=0.01D0*DFLOAT(MRATI)
      DFACIN=0.01D0*DFLOAT(NFACI)

      WARNER=.FALSE.  
      IF(MRATI.LT.90.OR.MRATI.GT.110) WARNER=.TRUE.
      IF(NRATI.GT. 1) WARNER=.TRUE.
      IF(NMISS1.NE.0) WARNER=.TRUE.
      IF(IAGAIN.NE.0) WARNER=.TRUE.

      IF(.NOT.WARNER) GOTO 95 
      WRITE(*,199) ' '
      WRITE(*,199) ' '
      WRITE(*,199) 
     + 'WarningWarningWarningWarningWarningWarningWarningWarningWar'
      WRITE(*,199)
     + 'arningWarningWarningWarningWarningWarningWarningWarningWarn'
      WRITE(*,199)
     + 'rningWarningWarningWarningWarningWarningWarningWarningWarni'
      WRITE(*,199)
     + 'ningWarningWarningWarningWarningWarningWarningWarningWarnin'
      WRITE(*,199)
     + 'ingWarningWarningWarningWarningWarningWarningWarningWarning'
      WRITE(*,199)     
     + 'ngWarningWarningWarningWarningWarningWarningWarningWarningW'
      WRITE(*,199)
     + 'gWarningWarningWarningWarningWarningWarningWarningWarningWa'

      IF(MRATI.LT.90.OR.MRATI.GT.110) THEN
      WRITE(*,199) ' '
      WRITE(*,*) '        Chi^2/Ndf = ',DRATIO,  
     +   '  (should be close to 1)'
      WRITE(*,*) '        => multiply all input standard ',
     +   'deviations by factor',DFACIN
      END IF 

      IF(NRATI.GT. 1) THEN
      WRITE(*,199) ' ' 
      WRITE(*,*) '        Fraction of rejects =',DJRAT,' %',
     +   '  (should be far below 1 %)'  
      WRITE(*,*) '        => please provide correct mille data' 
      END IF 

      IF(IAGAIN.NE.0) THEN
      WRITE(*,199) ' ' 
      WRITE(*,*) '        Matrix not positiv definite '//
     +    '(function not decreasing)'  
      WRITE(*,*) '        => please provide correct mille data' 
      END IF 

      IF(NMISS1.NE.0) THEN 
      WRITE(*,199) ' '
      WRITE(*,*) '        Rank defect =',NMISS1,
     +   '  for constraint equations, should be 0'
      WRITE(*,*) '        => please correct constraint definition' 
      END IF 

      WRITE(*,199) ' '
      WRITE(*,199) 
     + 'WarningWarningWarningWarningWarningWarningWarningWarningWar'
      WRITE(*,199)
     + 'arningWarningWarningWarningWarningWarningWarningWarningWarn'
      WRITE(*,199)
     + 'rningWarningWarningWarningWarningWarningWarningWarningWarni'
      WRITE(*,199)
     + 'ningWarningWarningWarningWarningWarningWarningWarningWarnin'
      WRITE(*,199)
     + 'ingWarningWarningWarningWarningWarningWarningWarningWarning'
      WRITE(*,199)     
     + 'ngWarningWarningWarningWarningWarningWarningWarningWarningW'
      WRITE(*,199)
     + 'gWarningWarningWarningWarningWarningWarningWarningWarningWa'
      WRITE(*,199) ' '
*

*
 95   CONTINUE

      CALL MEND                ! modul ending

*     ------------------------------------------------------------------

      IF(METSOL.EQ.1) THEN

      ELSE IF(METSOL.EQ.2) THEN
         CALL ZDIAGS
      ELSE IF(METSOL.EQ.3) THEN
*        errors and correlations from MINRES
        DO 98 K=1,MNRSEL
           LABELG=LBMNRS(K)
           IF(LABELG.EQ.0) GOTO 98
           ITGBI=INONE(LABELG)
           IVGBI=0
           IF(ITGBI.NE.0) IVGBI=JVGBI(ITGBI)
           IF(IVGBI.LT.0) IVGBI=0
           IF(IVGBI.EQ.0) GOTO 98
*          determine error and global correlation for parameter IVGBI
           CALL SOLGLO(IVGBI)
 98     CONTINUE
*
      ELSE IF(METSOL.EQ.4) THEN

      END IF

      CALL PRTGLO              ! print result

 102  FORMAT(' Call FEASIB with cut=',G10.3)
c 103  FORMAT(1X,A,G12.4)
 199  FORMAT(7X,A)
      END 


      SUBROUTINE FILETC

*     Fetch and interprete command line options
*     if steering file specified, check file existence (calling NUFILE) 

*     Initialize dynal memory management

*     If no steering file specified, check default steering file  

*     Read steering file, print some lines, detect names of text and
*     binary files, check file existence, store all file names

*     open all binary files


      CHARACTER*1024 TEXT,FILNAM 
#include "largeint.inc"
#include "dynal.inc"
      CHARACTER*14 BITE(3)
      CHARACTER*32 KEYSTX
      DOUBLE PRECISION DNUM(100)
      SAVE
#include "mpinds.inc"
      DATA BITE/'C_binary','text  ','Fortran_binary'/
*     ...
      CALL MSTART('FILETC')

      INTRAC=0  ! default is non-interactive
      ICTEST=0  ! no test 
      METSOL=0  ! reset solution method
      MATSTO=2  ! reset matrix storage mode
      MPRINT=1  ! reset print flag (0:minimal, 1:normal, >1: more)
      MDEBUG=0  ! reset debug flag (number of records to print)
      MDEBG2=10 ! number of measurements for record debug printout
      MREQEN=10 ! required number of entries 
      MITERA=1  ! number of iterations
      NLOOPN=0  ! number of data reading, fitting loops
      MBANDW=0  ! width of band matrix
      LUNKNO=0  ! flag for unkown keywords
      LHUBER=0  ! Huber down-weighting flag 
      CHICUT=0  ! cut in terms of 3-sigma cut, first iteration
      CHIREM=0  ! cut in terms of 3-sigma cut, other iterations
      CHHUGE=50.0 ! cut in terms of 3-sigma for unreasonable data, all iterations
      NRECPR=0  ! record number with printout
      NRECP2=0  ! record number with printout
      NREC1 =0  ! record number with largest residual
      NREC2 =0  ! record number with largest chi^2.Ndf
      VALUE1=0.0! largest residual
      VALUE2=0.0! largest chi^2.Ndf
      DWCUT=0.0 ! down-weigth cut
      ISUBIT=0  ! subito flag
      WOLFC1=0.0! C_1 of strong Wolfe condition
      WOLFC2=0.0! C_2 of strong Wolfe condition 
      MRESTL=1.0D-06 ! tolerance criterion for MINRES
      NOFEAS=0  ! default for making par feasible
      NHISTP=0  ! no histogram printout
      DELFUN=0.0! expected function change
      ACTFUN=0.0! actual function change
      ANGRAS=0.0! angle between gradient and search direction 
      ITERAT=0  ! iterations in solution
      NUF   =1  ! C binary is default 
      NREGUL=0  ! no regularization
      REGULA=1.0! regularization parameter
      REGPRE=0.0! default presigma
      MATRIT=0  ! matrix calculation up to iteration MATRIT
      ICALCM=0  ! if nonzero: calculate matrix
      ISAVES=1  ! save space if non-zero
      NBNDR =0  ! number of records with bordered band matrix for local fit
      NBDRX =0  ! max border size for local fit
      NBNDX =0  ! max band width for local fit
      NRECER=0  ! record with error (rank deficit or Not-a-Number) for printout
      NREC3 =MAXI4  ! (1.) record number with error
      MREQPE=1  ! min number of pair entries
      MHISPE=0  ! upper bound for pair entry histogrammimg
      MSNGPE=0  ! upper bound for pair entry single precision storage
      MCMPRS=0  ! compression for sparsity (column indices)
      MTHRD =1  ! no multithreading
!$    MTHRD =1  ! number of (OpenMP) threads
      MXREC =0  ! max number of records
      MATMON=0  ! record interval for monitoring of (sparse) matrix
      LFITNP=MAXI4 ! local fit: number of iteration to calculate pulls
      LFITBB=1  ! local fit: check for bordered band matrix (if >0)
      MNRSEL=0  ! number of MINRES error labels in LBMNRS (calc err, corr with S   OLGLO)
      NCACHE= MIN(MEGA/20,100000000) ! buffer size for caching (5% or 100M)
      FCACHE(1)= 0.8 ! read cache, average fill level
      FCACHE(2)= 0.0 ! write cache
      FCACHE(3)= 0.0 ! dynamic size
      MTHRDR=1  ! number of threads for reading binary files
      MNRSIT=0  ! total number of MINRES internal iterations

      DO I=1,8
       TIMES(I)=0.0
      END DO

*     read command line options ----------------------------------------

      FILNAM=' '   ! print command line options and find steering file
      DO I=1,IARGC()
       IF(I.EQ.1) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Command line options: '
          WRITE(*,*) '--------------------- '
       END IF
       CALL GETARG(I,TEXT)         ! get I.th text from command line
       CALL RLTEXT(TEXT,IA,IB,NAB) ! return indices for non-blank area
       WRITE(*,101) I,TEXT(1:NAB)  ! echo print
       IF(TEXT(IA:IA).NE.'-') THEN 
          NU=NUFILE(TEXT(IA:IB))   ! inquire on file existence
          IF(NU.EQ.2) THEN         ! existing text file
             IF(FILNAM.NE.' ') THEN
                WRITE(*,*) 'Second text file in command line - stop'
                STOP
             ELSE
                FILNAM=TEXT
             END IF
          ELSE
             WRITE(*,*) 'Open error for file:',TEXT(IA:IB),' - stop'
             STOP
          END IF
       ELSE
          IF(INDEX(TEXT(IA:IB),'i').NE.0) INTRAC=1 ! interactive
          IF(INDEX(TEXT(IA:IB),'j').NE.0) INTRAC=2 ! interactive j
          IF(INDEX(TEXT(IA:IB),'b').NE.0) THEN
              MDEBUG=3 ! debug flag
              INTRAC=1 ! ... and interactive
              WRITE(*,*) 'Debugging requested'
          END IF
         IT=INDEX(TEXT(IA:IB),'t')
          IF(IT.NE.0) THEN
            ICTEST=1  ! internal test files
            IEQ=INDEX(TEXT(IA+IT:IB),'=')+IT
            IF (IT.NE.IEQ) THEN
               IF (INDEX(TEXT(IA+IEQ:IB),'SL0' ).NE.0) ICTEST=2
               IF (INDEX(TEXT(IA+IEQ:IB),'SLE' ).NE.0) ICTEST=3
               IF (INDEX(TEXT(IA+IEQ:IB),'BP'  ).NE.0) ICTEST=4
               IF (INDEX(TEXT(IA+IEQ:IB),'BRLF').NE.0) ICTEST=5
               IF (INDEX(TEXT(IA+IEQ:IB),'BRLC').NE.0) ICTEST=6
            ENDIF
         ENDIF   
          IF(INDEX(TEXT(IA:IB),'s').NE.0) ISUBIT=1  ! like "subito"
       END IF
       IF(I.EQ.IARGC()) WRITE(*,*) '--------------------- '
      END DO
      

*     create test files for option -t ----------------------------------

      IF(ICTEST.GE.1) THEN
         WRITE(*,*) ' ' 
        IF (ICTEST.EQ.1) THEN
            CALL MPTEST ! 'wire chamber'
        ELSE
           CALL MPTST2(ICTEST-2) ! 'silicon tracker'
        ENDIF     
         IF(FILNAM.EQ.' ') FILNAM='mp2str.txt'
         WRITE(*,*) ' '  
      END IF   

*     check default steering file with file-name "steerfile" -----------

      IF(FILNAM.EQ.' ') THEN  ! check default steering file  
         TEXT='steerfile'
         CALL RLTEXT(TEXT,IA,IB,NAB) ! return indices for non-blank area
         NU=NUFILE(TEXT(IA:IB))      ! inquire on file existence and type
         IF(NU.GT.0) THEN
            FILNAM=TEXT
         ELSE
            STOP 'in FILETC: no steering file.                      .'
         END IF
      END IF


*     open, read steering file:
*                               end
*                               fortranfiles
*                               cfiles
  

      CALL RLTEXT(FILNAM,IA,IB,NFNAM) ! return indices for non-blank area
      WRITE(*,*) ' '
      WRITE(*,*) 'Listing of steering file: ',FILNAM(1:NFNAM)
      WRITE(*,*) '-------------------------'
      OPEN(10,FILE=FILNAM(1:NFNAM),IOSTAT=IOS)
      IF(IOS.NE.0) THEN
         WRITE(*,*) 'Open error for steering file - stop'
         STOP
      END IF 
      IFILE =0 
      NFILES=0

      NLINE=0
 10   READ(10,102,END=20) TEXT        ! read steering file
      CALL RLTEXT(TEXT,IA,IB,NAB)     ! return indices for non-blank area
      NLINE=NLINE+1 
      IF(NLINE.LE.50) THEN   ! print up to 50 lines
         WRITE(*,101) NLINE,TEXT(1:NAB)
         IF(NLINE.EQ.50)  WRITE(*,*) '     ...'
      END IF

      CALL RLTEXT(TEXT,IA,IB,NAB)        ! test content   'end'    
      IF(IB.EQ.IA+2) THEN
         MAT=MATINT(TEXT(IA:IB),'end',NPAT,NTEXT)
         IF(MAT.EQ.3) THEN
            TEXT=' '
            CALL INTEXT(TEXT,NLINE)
            WRITE(*,*) '    end-statement after',NLINE,' text lines'
            GOTO 20 
         END IF
      END IF

      KEYSTX='fortranfiles'
cGF      MAT=MATINT(TEXT(IA:IB),KEYSTX,NPAT,NTEXT)
cGF      IF(MAT.GE.NTEXT-NTEXT/10) THEN
      IF (TEXT(IA:IB).EQ.KEYSTX) THEN ! exact comparison by GF
         NUF=3
c         WRITE(*,*) 'Fortran files'
         GOTO 10   
      END IF

      KEYSTX='Cfiles'
cGF      MAT=MATINT(TEXT(IA:IB),KEYSTX,NPAT,NTEXT)
cGF      IF(MAT.GE.NTEXT-NTEXT/10) THEN 
      IF (TEXT(IA:IB).EQ.KEYSTX) THEN ! exact comparison by GF
         NUF=1
c         WRITE(*,*) 'Cfiles'
         GOTO 10
      END IF
  
* file names
*     check for file options (' -- ')
      IE=IB
      IOPT=INDEX(TEXT(IA:IB),' -- ')
      IF (IOPT.GT.0) IE=IOPT-1
*
      IF(NAB.EQ.0) GOTO 10
      NU=NUFILE(TEXT(IA:IE))          ! inquire on file existence 
      IF(NU.GT.0) THEN                ! existing file
         IF(NFILES.GE.MFILES) THEN
            STOP 'in FILETC: too many files.                           '
         END IF
         NFILES=NFILES+1              ! count number of files
         IF(NU.EQ.1) NU=NUF           !
         MFD(NFILES)=NU               ! cbinary =1, text =2, fbinary=3
         LFD(NFILES)=IE-IA+1          ! number of chars in file name 
         TFD(NFILES)=TEXT(IA:IE)      ! file name
         NFD(NFILES)=NLINE            ! line number 
         OFD(NFILES)=1.0              ! option for file
         IF (IOPT.GT.0) THEN
            CALL RATEXT(TEXT(IOPT+4:NAB),NUMS,DNUM) ! translate text to DP numbers
            IF (NUMS.GT.0) OFD(NFILES)=SNGL(DNUM(1))
         ENDIF
      ELSE   
c         WRITE(*,*) 'Open error for file ',TEXT(IA:IB)
c         STOP 
      END IF
      GOTO 10
 20   REWIND 10
      WRITE(*,*) '-------------------------'
      WRITE(*,*) ' '

*     print table of files ---------------------------------------------  

      IF (MPRINT.GT.1) THEN
         WRITE(*,*) 'Table of files:'
         WRITE(*,*) '---------------'
      ENDIF
      WRITE(8,*) ' '
      WRITE(8,*) 'Text and data files:' 
      DO I=1,NFILES
         IF (MPRINT.GT.1)
     +   WRITE(*,103) I,BITE(MFD(I)),TFD(I)(1:LFD(I))
         WRITE(8,103) I,BITE(MFD(I)),TFD(I)(1:LFD(I))
      END DO
      IF (MPRINT.GT.1) THEN
         WRITE(*,*) '---------------'
         WRITE(*,*) ' '
      ENDIF

*     open the binary (data) files on unit 11, 12, ...
 
      IOSUM=0
      NFILF=0 
      NFILB=0
      NFILW=0
      DO I=1,NFILES                                 ! Fortran files
       IF(MFD(I).EQ.3) THEN
          WRITE(*,*) 'Opening Fortran file ',10+NFILF+1,
     +               ' ',TFD(I)(1:LFD(I))
          OPEN(10+NFILF+1,FILE=TFD(I)(1:LFD(I)),IOSTAT=IOS,
     +         FORM='UNFORMATTED')
          IF(IOS.NE.0) THEN
             WRITE(*,*) 'Open error for file ',TFD(I)(1:LFD(I)) 
             IOSUM=IOSUM+1
          ELSE 
             NFILF=NFILF+1
             NFILB=NFILB+1
             WFD(NFILB)=OFD(I)
          END IF 
       END IF 
      END DO

*     open the binary C files

      NFILC=-1
      DO I=1,NFILES                                 ! Cfiles
       IF(MFD(I).EQ.1) THEN
#ifdef READ_C_FILES
          IF(NFILC.LT.0) CALL INITC ! uncommented by GF
          IF(NFILC.LT.0) NFILC=0
          WRITE(*,*) 'Opening C file ',NFILC+1, ': ',TFD(I)(1:LFD(I)) ! by GF
          CALL OPENC(TFD(I)(1:LFD(I)),IOS)
          IF(IOS.NE.0) THEN
*GF             WRITE(*,*) 'Open error for file ',TFD(I)(1:LFD(I))
             WRITE(*,*) 'Open error ',IOS,' for file ',TFD(I)(1:LFD(I))
             IOSUM=IOSUM+1 ! typo fixed by GF
          ELSE 
             NFILC=NFILC+1
             NFILB=NFILB+1
             WFD(NFILB)=OFD(I)
          END IF 
#else
          WRITE(*,*) 'Opening of C-files not supported.'
* GF add
          IOSUM=IOSUM+1
* GF add end
#endif
       END IF
      END DO
      DO K=1,NFILB
        KFD(1,K)=0   ! reset record counters
        KFD(2,K)=K   ! set file number
        IFD(K)=0
      ENDDO

      IF(IOSUM.NE.0) THEN
         STOP 'FILETC: open error                                      '
      END IF
      IF(NFILB.EQ.0) THEN
         STOP 'FILETC: no binary files                                 '
      END IF
      WRITE(*,*) NFILB,' binary files opened' ! corrected by GF
      RETURN



      ENTRY FILETX ! ---------------------------------------------------

*     Reset flags and read steering and all other text files
*     Print some lines from each file 
*     Store parameter values and constraints
*     Check flags METSOL (method of solution) and
*     MATSTO (matrix storage mode)
*     set default values for flags, which are undefined

*     Parameter values, stored in subarray '7':
*       1  label
*       2  (initial) parameter value
*       3  pre-sigma
*       4  label
*       5  (initial) parameter value
*       6  pre-sigma  
*       7  label
*     ...  ...
*     (number of words is multiple of 3)
*  

*     Constraint and Wconstrained data, stored in subarray '8':
*       1  0                 ! constraint header of four words:
*       2  right-hand-side   ! 0 and -1 ...
*       3  -1; -2            ! ... indicate ...
*       4  sigma             ! ... header
*       5  label
*       6  factor
*       7  label
*       8  factor
*       9  ...
*     ...  ...
*     (number of words is multiple of 2, at least 6)  

*     Measured data, stored in subarray '9':
*       1  0                 ! constraint header of four words:
*       2  right-hand-side   ! 0 and -1 ...
*       3  -1                ! ... indicate ...
*       4  sigma             ! ... header
*       5  label
*       6  factor
*       7  label
*       8  factor
*       9  ...
*     ...  ...
*     (number of words is multiple of 2, at least 6)

c      CALL MSTART('FILETX')



      WRITE(*,*) ' ' 
      WRITE(*,*) 'Processing text files ...'
      WRITE(*,*) ' '

*      CALL MEGARR('7 S',0,'Parameter values from text file')
*      CALL MEGARR('8 S',0,'Constraints from text file')
*      CALL MEGARR('9 S',0,'Measurements from text file') 

      IOSUM=0
      DO I=0,NFILES
       IF(I.EQ.0) THEN 
          WRITE(*,*) 'File ',FILNAM(1:NFNAM)
          NLINMX=100
       ELSE 
          NLINMX=10 
          IF(MFD(I).NE.2) GOTO 80 ! exclude binary files
          WRITE(*,*) 'File ',TFD(I)(1:LFD(I))
          IF (MPRINT.GT.1) WRITE(*,*) ' '
          OPEN(10,FILE=TFD(I)(1:LFD(I)),IOSTAT=IOS,FORM='FORMATTED')
          IF(IOS.NE.0) THEN
              WRITE(*,*) 'Open error for file ',TFD(I)(1:LFD(I))
              IOSUM=IOSUM+1
              GOTO 80
          END IF
       END IF
 
       NLINE=0
       NFILN=1
*      read text file     
 30    READ(10,102,END=40) TEXT 
       NLINE=NLINE+1
       IF(NLINE.LE.NLINMX.AND.MPRINT.GT.1) THEN ! print first 10 lines of every text file
          CALL RLTEXT(TEXT,IA,IB,NAB)
          NAB=MAX(1,NAB)
          WRITE(*,101) NLINE,TEXT(1:NAB)
          IF(NLINE.EQ.NLINMX) WRITE(*,*) '    ...'
       END IF 

       CALL RLTEXT(TEXT,IA,IB,NAB)        ! test content   'end'    
       IF(IB.EQ.IA+2) THEN
          MAT=MATINT(TEXT(IA:IB),'end',NPAT,NTEXT)
          IF(MAT.EQ.3) THEN
             TEXT=' '
             CALL INTEXT(TEXT,NLINE)
             WRITE(*,*) '    end-statement after',NLINE,' text lines'
             GOTO 50 
          END IF
       END IF

       IF(I.EQ.0) THEN ! first text file - exclude lines with file names
          IF(NFILN.LE.NFILES.AND.NLINE.EQ.NFD(NFILN)) THEN
             NFILN=NFILN+1
             TEXT=' '
C             WRITE(*,*) 'line is excluded ',TEXT(1:10)
          END IF
       END IF
c       WRITE(*,*) TEXT(1:40),'  < interprete text'
       CALL INTEXT(TEXT,NLINE)   ! interprete text
       GOTO 30      
 40    TEXT=' '
       CALL INTEXT(TEXT,NLINE)
       WRITE(*,*) '    end-of-file after',NLINE,' text lines'
 50    WRITE(*,*) ' ' 
       REWIND 10
       CLOSE(UNIT=10)
 80   CONTINUE
      END DO

      IF(IOSUM.NE.0) THEN
         STOP 'FILETX: open error(s) in text files                     '
      END IF
      CALL MEGLEN('7') ! adjust length of sub-array 7
      CALL MEGLEN('8') ! adjust length of sub-array 8
      CALL MEGLEN('9') ! adjust length of sub-array 9 
      IF(LND7.EQ.0) CALL MEGRM('7')
      IF(LND8.EQ.0) CALL MEGRM('8')
      IF(LND9.EQ.0) CALL MEGRM('9')

      WRITE(*,*) '... end of text file processing.' 
      WRITE(*,*) ' '

      IF(LUNKNO.NE.0) THEN
         WRITE(*,*) ' '
         WRITE(*,*) LUNKNO,' unknown keywords in steering files, ',
     +              'or file non-existing,' 
         WRITE(*,*) '   see above!'
         WRITE(*,*) '------------>    stop'
         WRITE(*,*) ' '
         STOP  
      END IF

*     check methods

      IF(METSOL.EQ.0) THEN        ! if undefined
         IF(MATSTO.EQ.0) THEN     ! if undefined
c            METSOL=1 ! default is matrix inversion 
c            MATSTO=1 ! default is symmetric matrix
         ELSE IF(MATSTO.EQ.1) THEN ! if symmetric  
            METSOL=3 ! MINRES
         ELSE IF(MATSTO.EQ.2) THEN ! if sparse
            METSOL=3 ! MINRES
         END IF
      ELSE IF(METSOL.EQ.1) THEN   ! if inversion
         MATSTO=1    !            
      ELSE IF(METSOL.EQ.2) THEN   ! if diagonalization
         MATSTO=1
      ELSE IF(METSOL.EQ.3) THEN   ! if MINRES
*        MATSTO=2 or 1 
      ELSE IF(METSOL.EQ.4) THEN   ! if GMRES
*        MATSTO=2 or 1 
      ELSE
         WRITE(*,*) 'MINRES forced with sparse matrix!'
         WRITE(*,*) ' '
         WRITE(*,*) 'MINRES forced with sparse matrix!'
         WRITE(*,*) ' '
         WRITE(*,*) 'MINRES forced with sparse matrix!'
         METSOL=3 ! forced
         MATSTO=2 ! forced
      END IF
      IF(MATSTO.GT.2) THEN
         WRITE(*,*) 'MINRES forced with sparse matrix!'
         WRITE(*,*) ' '
         WRITE(*,*) 'MINRES forced with sparse matrix!'
         WRITE(*,*) ' '
         WRITE(*,*) 'MINRES forced with sparse matrix!'
         METSOL=3 ! forced
         MATSTO=2 ! forced
      END IF

*     print information about methods and matrix storage modes 
        
      WRITE(*,*) ' '
      WRITE(*,*) 'Solution method and matrix-storage mode:' 
      IF(METSOL.EQ.1) THEN
         WRITE(*,*) '     METSOL = 1:  matrix inversion'
      ELSE IF(METSOL.EQ.2) THEN
         WRITE(*,*) '     METSOL = 2:  diagonalization'
      ELSE IF(METSOL.EQ.3) THEN
         WRITE(*,*) '     METSOL = 3:  MINRES'
      ELSE IF(METSOL.EQ.4) THEN
         WRITE(*,*) '     METSOL = 4:  GMRES'

      END IF

      WRITE(*,*) '                  with',MITERA,' iterations'   

      IF(MATSTO.EQ.1) THEN
         WRITE(*,*) '     MATSTO = 1:  symmetric matrix, ',
     +              '(n*n+n)/2 elements'
      ELSE IF(MATSTO.EQ.2) THEN
         WRITE(*,*) '     MATSTO = 2:  sparse matrix'
      END IF 
      IF(MBANDW.NE.0) THEN
      WRITE(*,*) '                  and band matrix, width',MBANDW
      END IF



      IF(CHICUT.NE.0.0) THEN 
         WRITE(*,*) 'Chi square cut equiv 3 st.dev applied ...'
         WRITE(*,*) ' in  first iteration with factor',CHICUT
         WRITE(*,*) ' in second iteration with factor',CHIREM     
         WRITE(*,*) ' (reduced by sqrt in next iterations)'     
      END IF

      IF(LHUBER.NE.0) THEN
         WRITE(*,*) '   Down-weighting of outliers in',
     +              LHUBER,' iterations'
         WRITE(*,*) '   Cut on downweight fraction',DWCUT
      END IF

      CALL MEND 

 101  FORMAT(I3,2X,A)
 102  FORMAT(A)
 103  FORMAT(I3,2X,A14,3X,A)
      END



      FUNCTION NUFILE(FILNAM) ! inquire on file

*     = 1 for existing binary file, =2 for existing text file, else =0
*     < 0 open error 
*     Text file names are recognized by the filename extension, which
*     should contain xt or tx.

      CHARACTER*(*) FILNAM
      LOGICAL EX
      SAVE
*     ...
      NUFILE=0
      IF(FILNAM(1:5).EQ.'rfio:') NUPRAE=1
      IF(FILNAM(1:5).EQ.'dcap:') NUPRAE=2
      IF(FILNAM(1:5).EQ.'root:') NUPRAE=3
      IF(NUPRAE.EQ.0) THEN
         INQUIRE(FILE=FILNAM,IOSTAT=IOS,EXIST=EX)
         IF(IOS.NE.0) NUFILE=-ABS(IOS)
         IF(IOS.NE.0) RETURN
      ELSE IF(NUPRAE.EQ.1) THEN  ! rfio:
         LL=LEN(FILNAM)
         FILNAM=FILNAM(6:LL)
         EX=.TRUE.
         NUFILE=1
         RETURN 
      ELSE
         EX=.TRUE.   ! assume file existence
      END IF
      IF(EX) THEN
         NUFILE=1                                   ! binary
         LL=LEN(FILNAM) 
         L1=MAX(1,LL-3)
         NE=MATINT('xt',FILNAM(L1:LL),NPAT,NTEXT)
         IF(NE.EQ.2) NUFILE=2                       ! text
         IF(NE.LT.2) THEN
            NE=MATINT('tx',FILNAM(L1:LL),NPAT,NTEXT)
            IF(NE.EQ.2) NUFILE=2                    ! text
         END IF
      END IF      
      END 


      SUBROUTINE INTEXT(TEXT,NLINE) ! interprete text
      CHARACTER*(*) TEXT
      PARAMETER (NKEYS=9,NMETH=4)
      CHARACTER*16 METHXT(NMETH)
      CHARACTER*16 KEYLST(NKEYS)
      CHARACTER*32 KEYWRD,KEYSTX 
      DOUBLE PRECISION DNUM(100)
      PARAMETER (NDIM=100)
      COMMON/CNUMS/LAST,KEYA,KEYB,ICD(2,NDIM)
      REAL    PLVS(3)    ! vector array: real and ... 
      INTEGER LPVS(3)    ! ... integer
      EQUIVALENCE (PLVS(1),LPVS(1))

      DATA KEYLST/'unknown','parameter','constraint','wconstraint',
     +            'measurement',
     +            'method',
     +            'mestimate',
     +            'atleast','option'/

*add  number of iterations
*     wconstraint like constraint
*     measured   r sigma label factor ...(more)    

      SAVE
#include "largeint.inc"
#include "mpinds.inc"
      DATA METHXT/'diagonalization','inversion','fullMINRES',
     +            'sparseMINRES'/
      DATA LKEY/-1/                 ! last keyword

*     ...
      NKEY=-1                       ! new keyword 
      CALL RLTEXT(TEXT,IA,IB,NAB)   ! return indices for non-blank area
      IF(NAB.EQ.0) GOTO 10          
      CALL RATEXT(TEXT(1:NAB),NUMS,DNUM) ! translate text to DP numbers

      IF(NUMS.NE.0) NKEY=0
      IF(KEYB.NE.0) THEN
         KEYWRD=TEXT(KEYA:KEYB) !          text is TEXT(KEYA)...TEXT(KEYB)
c         WRITE(*,*) 'Keyword is ',KEYWRD

*        compare keywords

         DO NKEY=2,NKEYS-1      ! loop over all pede keywords
          KEYSTX=KEYLST(NKEY)   ! copy NKEY.th pede keyword
          MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
          IF(MAT.GE.NTEXT-NTEXT/5) GOTO 10
         END DO   

*        more comparisons 
          
         KEYSTX='print'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
c         WRITE(*,*) KEYSTX,MAT,NTEXT
c         IF(MAT.GE.NTEXT) THEN
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
C         IF(MAT.GE.(NTEXT+NTEXT+3)/3) THEN
            MPRINT=1
            IF(NUMS.GT.0) MPRINT=IDNINT(DNUM(1))
            RETURN 
         END IF

         KEYSTX='debug'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN          
c         IF(MAT.GE.(NTEXT+NTEXT+3)/3) THEN
            MDEBUG=3
* GF            IF(NUMS.GT.0) MPRINT=DNUM(1)
            IF(NUMS.GT.0) MDEBUG=IDNINT(DNUM(1))
            IF(NUMS.GT.1) MDEBG2=IDNINT(DNUM(2))
           RETURN
         END IF

         KEYSTX='entries'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN          
c         IF(MAT.GE.(NTEXT+NTEXT+3)/3) THEN
            MREQEN=0
            IF(NUMS.GT.0) MREQEN=IDNINT(DNUM(1))
            IF (MREQEN.LE.0) MREQEN=1
            RETURN
         END IF

         KEYSTX='printrecord'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            NRECPR=IDNINT(DNUM(1))
            NRECP2=IDNINT(DNUM(2))
            RETURN
         END IF

         KEYSTX='maxrecord'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            IF (NUMS.GT.0.AND.DNUM(1).GT.0.) MXREC=IDNINT(DNUM(1))
            RETURN
         END IF

         KEYSTX='cache'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            IF (NUMS.GT.0.AND.DNUM(1).GE.0.) NCACHE=IDNINT(DNUM(1)) ! cache size, <0 keeps default
            IF (NUMS.EQ.2.AND.DNUM(2).GT.0..AND.DNUM(2).LE.1.0)     ! read cache fill level
     &         FCACHE(1)=SNGL(DNUM(2))
            IF (NUMS.GE.4) THEN                                     ! explicit cache splitting
               DO K=1,3
                 FCACHE(K)=SNGL(DNUM(K+1))
               ENDDO
            ENDIF
            RETURN
         END IF

         KEYSTX='chisqcut'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            IF(NUMS.EQ.0) THEN ! always 3-sigma cut
               CHICUT=1.0 
               CHIREM=1.0
            ELSE
               CHICUT=DNUM(1)
               IF(CHICUT.LT.1.0) CHICUT=-1.0
               IF(NUMS.EQ.1) THEN 
                  CHIREM=1.0  ! 3-sigma cut, if not specified
               ELSE  
                  CHIREM=SNGL(DNUM(2))
                  IF(CHIREM.LT.1.0) CHIREM=1.0
                  IF(CHICUT.GE.1.0) CHIREM=MIN(CHIREM,SNGL(CHICUT))
               END IF
            END IF
            RETURN
         END IF

* GF added:
         KEYSTX='hugecut'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN          
            IF(NUMS.GT.0) CHHUGE=SNGL(DNUM(1))
            IF(CHHUGE.LT.1.0) CHHUGE=1.0 ! at least (!!) 3-sigma
            RETURN
         END IF
* GF added end

         KEYSTX='localfit'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN          
            IF(NUMS.GT.0) LFITNP=IDNINT(DNUM(1))
            IF(NUMS.GT.1) LFITBB=IDNINT(DNUM(2))
            RETURN
         END IF

         KEYSTX='regularization'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            NREGUL=1
            REGULA=SNGL(DNUM(1))
            IF(NUMS.GE.2) REGPRE=SNGL(DNUM(2))
            RETURN
         END IF

         KEYSTX='regularisation'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            NREGUL=1
            REGULA=SNGL(DNUM(1))
            IF(NUMS.GE.2) REGPRE=SNGL(DNUM(2))
            RETURN
         END IF

         KEYSTX='presigma'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            REGPRE=SNGL(DNUM(1))
            RETURN
         END IF

         KEYSTX='matiter'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            MATRIT=IDNINT(DNUM(1))
            RETURN
         END IF

         KEYSTX='matmoni'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            MATMON=-1
            IF (NUMS.GT.0.AND.DNUM(1).GT.0.) MATMON=IDNINT(DNUM(1))
            RETURN
         END IF

         KEYSTX='bandwidth'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
c         IF(MAT.GE.(NTEXT+NTEXT+3)/3) THEN
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            IF(NUMS.GT.0)   MBANDW=IDNINT(DNUM(1))
            IF(MBANDW.LT.0) MBANDW=-1
            RETURN
         END IF

C         KEYSTX='outlierrejection'
C         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
c         WRITE(*,*) KEYSTX,MAT,(NTEXT+NTEXT)/3
c         IF(MAT.GE.(NTEXT+NTEXT+NTEXT-2)/3) THEN
C         IF(MAT.GE.(NPAT-NPAT/5)) THEN
C            CHDFRJ=DNUM(1)
C            IF(CHDFRJ.LT.3.0) CHDFRJ=100.0
C            RETURN
C         END IF

C         KEYSTX='outliersuppression'
C         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
c         WRITE(*,*) KEYSTX,MAT,(NTEXT+NTEXT)/3
c         IF(MAT.GE.(NTEXT+NTEXT+NTEXT-2)/3) THEN
C         IF(MAT.GE.(NPAT-NPAT/5)) THEN
C            LHUBER=DNUM(1)
C            IF(LHUBER.LE.2) LHUBER=2 ! at least 2 Huber iterations
C            RETURN
C         END IF

         KEYSTX='outlierdownweighting'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
c         WRITE(*,*) KEYSTX,MAT,(NTEXT+NTEXT)/3
c         IF(MAT.GE.(NTEXT+NTEXT+NTEXT-2)/3) THEN
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            LHUBER=IDNINT(DNUM(1))
            IF(LHUBER.GT.0.AND.LHUBER.LE.2) LHUBER=2 ! at least 2 Huber iterations (if any)
            RETURN
         END IF

         KEYSTX='dwfractioncut'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            DWCUT=SNGL(DNUM(1))
            IF(DWCUT.GT.0.5) DWCUT=0.5
            RETURN
         END IF

         KEYSTX='subito'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            ISUBIT=1
            RETURN
         END IF

         KEYSTX='threads'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            IOMP=0
!$          IOMP=1
!$          IF (IOMP.GT.0) THEN
!$             IF (NUMS.GE.1.AND.DNUM(1).GT.0.) MTHRD =IDNINT(DNUM(1))
!$             IF (NUMS.GE.2) THEN
!$                MTHRDR=MTHRD
!$                IF (DNUM(2).GT.0.) MTHRDR=IDNINT(DNUM(2))
!$             ENDIF
!$          ELSE
               WRITE(*,*) 'WARNING: multithreading not available'
!$          ENDIF
            RETURN
         END IF

         KEYSTX='compress'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            MCMPRS=1
            RETURN
         END IF

         KEYSTX='errlabels'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5).AND.MNRSEL.LT.100) THEN
            NL=MIN(NUMS,100-MNRSEL)
            DO K=1,NL
               LBMNRS(MNRSEL+K)=IDNINT(DNUM(K))
            ENDDO
            MNRSEL=MNRSEL+NL
            RETURN
         END IF

         KEYSTX='pairentries'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
c     This option could be implemented to get rid of parameter pairs
c     that have very few entries - to save matrix memory size. 
            IF (NUMS.GT.0.AND.DNUM(1).GT.0.0) THEN
               MREQPE=IDNINT(DNUM(1))
               IF (NUMS.GE.2.AND.DNUM(2).GE.DNUM(1)) 
     +            MHISPE=IDNINT(DNUM(2))
               IF (NUMS.GE.3.AND.DNUM(3).GT.0.0) 
     +            MSNGPE=IDNINT(DNUM(3))
               ICOUNT=MAX(MSNGPE+1,MHISPE)
               ICOUNT=MAX(MREQPE,ICOUNT)
               ISAVES=1 ! number of bits needed to count up to ICOUNT
               MXCNT=2
               DO K=1,30
                  IF (ICOUNT.GE.MXCNT) THEN
                     ISAVES=ISAVES+1
                     MXCNT=MXCNT*2
                  ENDIF
               ENDDO
            ENDIF
            RETURN
         END IF

         KEYSTX='wolfe'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            WOLFC1=SNGL(DNUM(1))
            WOLFC2=SNGL(DNUM(2))
            RETURN
         END IF

* GF added:
* convergence tolerance for minres:
         KEYSTX='mrestol'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            IF(NUMS.GT.0) THEN
               IF (DNUM(1).LT.1.0D-10.OR.DNUM(1).GT.1.0D-04) THEN
                  WRITE(*,*) 'ERROR: need 1.0D-10 <= MRESTL ',
     +                 '<= 1.0D-04, but get ', DNUM(1)
               ELSE
                  MRESTL=IDNINT(DNUM(1))
               ENDIF
            ENDIF
            RETURN
         END IF
* GF added end

         KEYSTX='nofeasiblestart'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            NOFEAS=1 ! do not make parameters feasible at start
            RETURN
         END IF 

         KEYSTX='histprint'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.NTEXT-NTEXT/10) THEN
c         IF(MAT.GE.(NPAT-NPAT/5)) THEN
            NHISTP=1 ! print histograms
            RETURN
         END IF

         KEYSTX='fortranfiles'
         MAT=MATINT(TEXT(IA:IB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.NTEXT-NTEXT/10) RETURN

         KEYSTX='Cfiles'
         MAT=MATINT(TEXT(IA:IB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.NTEXT-NTEXT/10) RETURN

         KEYSTX=KEYLST(1)
         NKEY=1                 ! unknown keyword
         IF(NUMS.NE.0) NKEY=0   

         WRITE(*,*) ' '
         WRITE(*,*) '**************************************************'
         WRITE(*,*) ' '
         WRITE(*,*) 'Unknown keyword(s): ',TEXT(1:MIN(NAB,50))
         WRITE(*,*) ' '
         WRITE(*,*) '**************************************************'
         WRITE(*,*) ' '
         LUNKNO=LUNKNO+1

      END IF 
 10   CONTINUE
*     result: NKEY = -1    blank
*             NKEY =  0    numerical data, no text keyword or unknown 
*             NKEY >  0    keyword NKEY from list, keyword = KEYSTX 


*     content/lastcontent
*     -------------------
*     blank            -1
*     data              0
*     keyword
*      unknown          1
*      parameter        2 
*      constraint       3
*      method           4 
*      atleast          5 
*      option           6
*      cut              7

      IF(NKEY.GT.0) THEN     ! new keyword 
         LKEY=NKEY
         IF(LKEY.EQ.2) THEN                                 ! parameter
            IF(NUMS.EQ.3) THEN
               LPVS(1)=IDNINT(DNUM(1)) ! label
               PLVS(2)=SNGL(DNUM(2))   ! start value
               PLVS(3)=SNGL(DNUM(3))   ! pre-sigma
               IF(LPVS(1).NE.0) THEN
                  CALL MEGVEC('7',PLVS,3,0)      ! always 3 words
               ELSE
                  WRITE(*,*) 'Line',NLINE,' error, label=',LPVS(1)
               END IF
            ELSE IF(NUMS.NE.0) THEN
               KKEY=1   ! switch to "unknown"  ?
               WRITE(*,*) 'Wrong text in line',NLINE
               WRITE(*,*) 'Status: new parameter'
               WRITE(*,*) '> ',TEXT(1:NAB)
            END IF 
         ELSE IF(LKEY.EQ.3.OR.LKEY.EQ.4) THEN               ! constraint
c            WRITE(*,*) 'Keyword is constraint!',NUMS,' numerical data'
            IF(NUMS.GE.1.AND.NUMS.LE.2) THEN ! start constraint
c               WRITE(*,*) 'NUMs is ',NUMS
               LPVS(1)=0
               PLVS(2)=SNGL(DNUM(1))         ! r = r.h.s. value
               CALL MEGVEC('8',PLVS,2,0)
c               WRITE(*,*) 'LPVS PLVS ',LPVS,PLVS
               LPVS(1)=-1                    ! constraint  
               IF(LKEY.EQ.4) LPVS(1)=-2      ! wconstraint (weighted)
               PLVS(2)=0.0
c               WRITE(*,*) 'LPVS PLVS ',LPVS,PLVS
               IF(NUMS.EQ.2) PLVS(2)=SNGL(DNUM(2)) ! sigma 
               CALL MEGVEC('8',PLVS,2,0)
c               WRITE(*,*) 'LPVS PLVS ',LPVS,PLVS
            ELSE
               KKEY=1   ! switch to "unknown"
               WRITE(*,*) 'Wrong text in line',NLINE
               WRITE(*,*) 'Status: new keyword (w)constraint'
               WRITE(*,*) '> ',TEXT(1:NAB)
            END IF
         ELSE IF(LKEY.EQ.5) THEN                           ! measurement
            IF(NUMS.EQ.2) THEN ! start measurement
               LPVS(1)=0
               PLVS(2)=SNGL(DNUM(1))         ! r = r.h.s. value
               CALL MEGVEC('9',PLVS,2,0)
               LPVS(1)=-1                    ! constraint
               PLVS(2)=SNGL(DNUM(2))         ! sigma
               CALL MEGVEC('9',PLVS,2,0)
            ELSE
               KKEY=1   ! switch to "unknown"
               WRITE(*,*) 'Wrong text in line',NLINE
               WRITE(*,*) 'Status: new keyword measurement'
               WRITE(*,*) '> ',TEXT(1:NAB)
            END IF

         ELSE IF(LKEY.EQ.6) THEN                                ! method
            MITER=MITERA
            IF(NUMS .GE.1) MITER=IDNINT(DNUM(1))
            IF(MITER.GE.1) MITERA=MITER
            DFLIM=SNGL(DNUM(2))
            LKEY=0
            DO I=1,NMETH  
             KEYSTX=METHXT(I)
             MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
             IF(MAT.GE.NTEXT-NTEXT/5) THEN
               IF(I.EQ.1) THEN            ! diagonalization
                  METSOL=2
                  MATSTO=1
               ELSE IF(I.EQ.2) THEN       ! inversion
                  METSOL=1
                  MATSTO=1 
               ELSE IF(I.EQ.3) THEN       ! fullMINRES
                  METSOL=3
                  MATSTO=1
c                  METSOL=4 !!!
               ELSE IF(I.EQ.4) THEN       ! sparseMINRES
c                  METSOL=4
                  METSOL=3
                  MATSTO=2

               END IF
             END IF 
            END DO
         END IF         
      ELSE IF(NKEY.EQ.0) THEN  ! data for continuation
         IF(LKEY.EQ.2) THEN                                  ! parameter
            IF(NUMS.GE.3) THEN   ! store data from this line
               LPVS(1)=IDNINT(DNUM(1)) ! label
               PLVS(2)=SNGL(DNUM(2))   ! start value
               PLVS(3)=SNGL(DNUM(3))   ! pre-sigma
               IF(LPVS(1).NE.0) THEN
                  CALL MEGVEC('7',PLVS,3,0)      ! always 3 words
               ELSE
                  WRITE(*,*) 'Line',NLINE,' error, label=',LPVS(1)
               END IF
            ELSE IF(NUMS.GT.1.AND.NUMS.LT.3) THEN
               KKEY=1   ! switch to "unknown"  ?
               WRITE(*,*) 'Wrong text in line',NLINE
               WRITE(*,*) 'Status continuation parameter'
               WRITE(*,*) '> ',TEXT(1:NAB)
            END IF

         ELSE IF(LKEY.EQ.3.OR.LKEY.EQ.4) THEN            ! (w)constraint
            IER=0
            DO I=1,NUMS,2 
             LABEL=IDNINT(DNUM(I))
             IF(LABEL.LE.0) IER=1
            END DO
            IF(MOD(NUMS,2).NE.0) IER=1 ! reject odd number
            IF(IER.EQ.0) THEN
               DO I=1,NUMS,2
                LPVS(1)=IDNINT(DNUM(I)) ! label
                PLVS(2)=SNGL(DNUM(I+1)) ! factor
                CALL MEGVEC('8',PLVS,2,0)
               END DO
            ELSE
               KKEY=0
               WRITE(*,*) 'Wrong text in line',NLINE
               WRITE(*,*) 'Status continuation (w)constraint'
               WRITE(*,*) '> ',TEXT(1:NAB)
            END IF

         ELSE IF(LKEY.EQ.5) THEN                           ! measurement
c            WRITE(*,*) 'continuation                          < ',NUMS
            IER=0
            DO I=1,NUMS,2 
             LABEL=IDNINT(DNUM(I))
             IF(LABEL.LE.0) IER=1
            END DO
            IF(MOD(NUMS,2).NE.0) IER=1 ! reject odd number
c            WRITE(*,*) 'IER NUMS ',IER,NUMS
            IF(IER.EQ.0) THEN
               DO I=1,NUMS,2
                LPVS(1)=IDNINT(DNUM(I)) ! label
                PLVS(2)=SNGL(DNUM(I+1)) ! factor
                CALL MEGVEC('9',PLVS,2,0)
               END DO
            ELSE
               KKEY=0
               WRITE(*,*) 'Wrong text in line',NLINE
               WRITE(*,*) 'Status continuation measurement'            
               WRITE(*,*) '> ',TEXT(1:NAB)
            END IF

         END IF
      END IF 
      END

*     ------------------------------------------------------------------ 

      SUBROUTINE MSTART(TEXT)
      CHARACTER*(*) TEXT
      CHARACTER*16  TEXTC
      CHARACTER*74  TEXTL
      SAVE
*     ...
      DO I=1,74
       TEXTL(I:I)='_'
      END DO 
      L=LEN(TEXT)
      KA=(74-L)/2
      KB=KA+L-1
      TEXTL(KA:KB)=TEXT(1:L)
      WRITE(*,*) ' '
      WRITE(*,*) TEXTL
      WRITE(*,*) ' ' 
      TEXTC=TEXT(1:L)//'-end'
*     
      DO I=1,74
       TEXTL(I:I)='_'
      END DO 
      L=L+4
      KA=(74-L)/2
      KB=KA+L-1
      TEXTL(KA:KB)=TEXTC(1:L)
      RETURN

      ENTRY MEND
      WRITE(*,*) ' '
      WRITE(*,*) TEXTL
      CALL PETIME 
      WRITE(*,*) ' '
      END 

      SUBROUTINE MVOPEN(LUN,FILNAM)
*     open file. Evtl. move existing file to ~ version.
      CHARACTER*(*) FILNAM
      CHARACTER*33  NAFILE,NBFILE
      LOGICAL EX
      SAVE
*     ...
      L=LEN(FILNAM)
      IF(L.GT.32) STOP 'File name too long                    '
      NAFILE=FILNAM 
      NAFILE(L+1:L+1)='~'

      INQUIRE(FILE=NAFILE(1:L),EXIST=EX)
      IF(EX) THEN
         INQUIRE(FILE=NAFILE(1:L+1),EXIST=EX)
         IF(EX) THEN
            CALL SYSTEM('rm '//NAFILE)
         END IF
         NBFILE=NAFILE
         NAFILE(L+1:L+1)=' '
         CALL SYSTEM('mv '//NAFILE//NBFILE)
      END IF
      OPEN(UNIT=LUN,FILE=FILNAM)
      END


*     Sparse matrix structure using bits (no counting) or bit fields (some counting)
*
*     CALL CLBITS(N,NDIMB,NBFW)   ! reset for N x N symmetric matrix
*                                 ! N      = matrix dimension
*                                 ! NDIMB  = dimension for bit array 
*                                 !         (NBITS= bit field array: NDIMB*NBFW)
*                                 ! NENCDB = number of bits for encoding column counter 
*                                 ! NBFW   = bit field width
*
*     do 
*       CALL INBITS(I,J,INC,NBITS)! add index pair (i,j)
*                                 ! I, J         = indices
*                                 ! INC          = increment (usually 1)
*     end do                      ! NBITS(NDIMB) = array for bits (bit fields) 
*   
*
*     CALL NDBITS(NDIMS,NBITS,NCMPRS,MNPAIR,IHST,ICMPRS)  ! return dimension parameter for sparse
*                                 ! NDIMS(1)     = (reduced) dimension for NBITS
*                                 ! NDIMS(2)     = dimension for NSPARC
*                                 ! NDIMS(3)     = number of off diagonal elements in double prec.
*                                 ! NDIMS(4)     = number of off diagonal elements in single prec.
*                                 ! NBITS(NDIMB) = array for bits (bit fields)
*                                 ! NCMPRS(N)    = compression mode per row
*                                 ! MNPAIR       = min. entries for pairs
*                                 ! IHST         = >0: histogram number
*                                 ! ICOMPRS      = <>0: compress row information (column indices)
*     
*     CALL SPBITS(NSPARR,NSPARC,NBITS)   ! build sparse structure  
*                                 ! NSPARR(N+1)   = sparse structure  I: row pointers
*                                 ! NSPARC(NDIMS) = sparse structure II: column lists
*                                 ! NBITS(NDIMB)  = array for bits (bit fields)
*     array NBITS(NDIMB) no longer used after SPBITS(NSPARR,NSPARC,NBITS)

*     ------------------------------------------------------------------ 
*     fill bit fields
      SUBROUTINE INBITS(IM,JM,INC,NBITS)        ! include element (I,J)
#include "largeint.inc"
      INTEGER(KIND=LARGE) NDIMB
      COMMON / MPBITS / NDIMB,N,IBFW,MXCNT,NENCDM,NENCDB,NTHRD
      INTEGER NBITS(*)
      INTEGER(KIND=LARGE) L,LL
      INTEGER*8 NOFFI
      LOGICAL BTEST
*
      IF(IM.EQ.JM) RETURN  ! diagonal
      J=MIN(IM,JM)
      I=MAX(IM,JM)
      IF(J.LE.0) RETURN    ! out low
      IF(I.GT.N) RETURN    ! out high

      NOFFI=INT8(I-1)*INT8(I-2)*INT8(IBFW)/2 ! for J=1
      NOFFJ=(J-1)*IBFW
      L=NOFFI/32+I ! row offset
     & +NOFFJ/32   ! column offset
*     add I instead of 1 to keep bit maps of different rows in different words (openMP !)
      M=MOD(NOFFJ,32)
      IF (IBFW.LE.1) THEN
         NBITS(L)=IBSET(NBITS(L),M)
      ELSE
*        get counter from bit field
         LL=L
         MM=M
         ICOUNT=0
         DO IB=0,IBFW-1
            IF (BTEST(NBITS(LL),MM)) ICOUNT=IBSET(ICOUNT,IB)
            MM=MM+1
            IF (MM.GE.32) THEN
               LL=LL+1
               MM=MM-32
            ENDIF
         ENDDO
*        increment
         JCOUNT=ICOUNT
         ICOUNT=MIN(ICOUNT+INC,MXCNT)
*        store counter into bit field
         IF (ICOUNT.NE.JCOUNT) THEN
            LL=L
            MM=M
            DO IB=0,IBFW-1
               IF (BTEST(ICOUNT,IB)) THEN
                  NBITS(LL)=IBSET(NBITS(LL),MM)
               ELSE
                  NBITS(LL)=IBCLR(NBITS(LL),MM)
               ENDIF
               MM=MM+1
               IF (MM.GE.32) THEN
                 LL=LL+1
                 MM=MM-32
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END

*     ------------------------------------------------------------------ 
*     calculate bit field size
      SUBROUTINE CLBITS(IN,IDIMB,IENCDB,JBFW)           ! reset for I,J =< IN
#include "largeint.inc"
      INTEGER(KIND=LARGE) NDIMB,IDIMB
      COMMON / MPBITS / NDIMB,N,IBFW,MXCNT,NENCDM,NENCDB,NTHRD
      INTEGER*8 NOFFD
!$    INTEGER OMP_GET_MAX_THREADS

      N=IN
      IBFW=JBFW
      MXCNT=2**IBFW-1
      NOFFD=INT8(N)*INT8(N-1)*INT8(IBFW)/2
      NDIMB=NOFFD/32+N
      IDIMB=NDIMB
      MB=INT(4.0E-6*FLOAT(NDIMB))
      WRITE(*,*) ' '
      WRITE(*,*) 'CLBITS: symmetric matrix of dimension',N
      WRITE(*,*) 'CLBITS: off-diagonal elements',NOFFD
      WRITE(*,*) 'CLBITS: dimension of bit-array',NDIMB , '(',MB,'MB)'
      IF (IBFW.GT.1) WRITE(*,*) 'CLBITS: bit field width',IBFW
*     encoding for compression
      NBCOL=16    ! 16 bits for column number, 16 bits for column counter
      DO I=16,30
         IF (BTEST(N,I)) NBCOL=I+1
      ENDDO
      NENCDB=32-NBCOL
      IENCDB=NENCDB
      NENCDM=ISHFT(1,NENCDB)-1
      NTHRD=1
!$    NTHRD=OMP_GET_MAX_THREADS()
c      print *, ' encoding ', NENCDB, NENCDM
      RETURN
      END

*     ------------------------------------------------------------------
*     analyze bit fields
      SUBROUTINE NDBITS(NDIMS,NBITS,NCMPRS,NSPARR,MNPAIR,IHST,JCMPRS)
#include "largeint.inc"
      INTEGER(KIND=LARGE) NDIMB
      COMMON / MPBITS / NDIMB,N,IBFW,MXCNT,NENCDM,NENCDB,NTHRD
      INTEGER NBITS(*),NCMPRS(*)
      INTEGER NWCP(0:1),IRGN(2),INR(2)
      INTEGER(KIND=LARGE) NDIMS(4),LL,LB,NIN,NTOT,NSPARR(*)
      INTEGER*8 NOFFI
      LOGICAL BTEST
!$    INTEGER OMP_GET_THREAD_NUM

      NDIMS(1)=NDIMB
      NDIMS(2)=0
      NDIMS(3)=0
      NDIMS(4)=0
      NTOT=0
      ICMPRS=JCMPRS
*     reduce bit field counters to bits
      KBFW=1
      LL=0
      LB=0
      ICHUNK=MIN(N/NTHRD/8+1,256)
*
      IF (IBFW.GT.1.OR.ICMPRS.NE.0) THEN
         IF (IBFW.GT.1.AND.ICMPRS.GT.0) KBFW=2 ! need to tag single precision entries
         JBFW=KBFW ! temp. bit field width
         IF (NTHRD.GT.1) JBFW=IBFW ! don't mix rows in NBITS
*
!$OMP  PARALLEL DO                 ! parallelize row loop
!$OMP& PRIVATE(I,NOFFD,LL,MM,LB,MB,INR,IRGN,LAST,LRGN,
!$OMP&         J,ICOUNT,NEXT,IB,ICP,NWCP,JB,JP,IR)
!$OMP& REDUCTION(+:NDIMS,NTOT)    ! private copy of NDIMS,NTOT for each thread, combined at end, init with 0.
!$OMP& SCHEDULE(DYNAMIC,ICHUNK)    ! some comment
         DO I=1,N
            NOFFI=INT8(I-1)*INT8(I-2)*INT8(IBFW)/2
            LL=NOFFI/32+I 
            MM=0
            NOFFI=INT8(I-1)*INT8(I-2)*INT8(JBFW)/2
            LB=NOFFI/32+I 
            MB=0
            INR(1)=0
            INR(2)=0
            IRGN(1)=0
            IRGN(2)=0
            LAST=0
            LRGN=0
            IPROC=0
!$          IPROC=OMP_GET_THREAD_NUM()         ! thread number
            DO J=1,I-1 
               ICOUNT=0
               NEXT=0
               DO IB=0,IBFW-1
                  IF (BTEST(NBITS(LL),MM)) ICOUNT=IBSET(ICOUNT,IB)
                  MM=MM+1
                  IF (MM.GE.32) THEN
                     LL=LL+1
                     MM=MM-32
                  ENDIF
               ENDDO
               DO JB=0,KBFW-1
                  NBITS(LB)=IBCLR(NBITS(LB),MB+JB)
               ENDDO
*
               IF (ICOUNT.GT.0) THEN
                  NTOT=NTOT+1
                  IF (IPROC.EQ.0.AND.IHST.GT.0) 
     &               CALL HMPENT(IHST,FLOAT(ICOUNT))
               ENDIF
*              keep pair ?
               IF (ICOUNT.GE.MNPAIR) THEN
                  NEXT=1 ! double
                  IF (ICOUNT.LE.ICMPRS.AND.ICMPRS.GT.0) NEXT=2 ! single 
                  INR(NEXT)=INR(NEXT)+1
                  NBITS(LB)=IBSET(NBITS(LB),MB+NEXT-1)
                  IF (NEXT.NE.LAST.OR.LRGN.GE.NENCDM) THEN
                     IRGN(NEXT)=IRGN(NEXT)+1
                     LRGN=0
                  ENDIF
                  LRGN=LRGN+1
               ENDIF
               MB=MB+KBFW
               IF (MB.GE.32) THEN
                  LB=LB+1
                  MB=MB-32
               ENDIF
               LAST=NEXT
            ENDDO
*
            DO JP=1,KBFW
               ICP=0
               NWCP(0)=INR(JP)                       ! list of column indices (default)
               IF (INR(JP).GT.0) THEN
                  NWCP(1)=IRGN(JP)+(IRGN(JP)+7)/8    ! list of regions of consecutive columns
*                 compress row ?
                  IF (NWCP(1).LT.NWCP(0).AND.ICMPRS.NE.0) THEN
                     ICP=1
                     NCMPRS(I+N*(JP-1))=IRGN(JP)
                  ENDIF
                  NDIMS(2)   =NDIMS(2)   +NWCP(ICP)
                  NDIMS(JP+2)=NDIMS(JP+2)+NWCP(0)
               ENDIF
               IR=(I+(N+1)*(JP-1))*2
               NSPARR(IR+1)=NWCP(ICP)
               NSPARR(IR+2)=NWCP(0)
            ENDDO
*
         ENDDO
!$OMP END PARALLEL DO
*        sum up, fill row offsets
         LB=1
         N1=0
         LL=N+1
         DO JP=1,KBFW
            DO I=1,N
               N1=N1+2
               NSPARR(N1-1)=LB
               NSPARR(N1  )=LL
               LB=LB+NSPARR(N1+1)
               LL=LL+NSPARR(N1+2)
            ENDDO
            N1=N1+2
            NSPARR(N1-1)=LB
            NSPARR(N1  )=LL
            LL=1
         ENDDO
*
         IF (JBFW.NE.KBFW) THEN ! move bit fields 
            DO I=1,N
               NOFFI=INT8(I-1)*INT8(I-2)*INT8(JBFW)/2
               LL=NOFFI/32+I 
               NOFFI=INT8(I-1)*INT8(I-2)*INT8(KBFW)/2
               LB=NOFFI/32+I 
               NJ=((I-1)*KBFW)/32
               DO K=0,NJ
                 NBITS(LB+K)=NBITS(LL+K)
               ENDDO
            ENDDO
         ENDIF
*
         IBFW=KBFW
         NOFFI=INT8(N)*INT8(N-1)*INT8(IBFW)/2
         NDIMB=NOFFI/32+N
         NDIMS(1)=NDIMB
*
      ELSE
*
         NIN=0
         NSPARR(1)=1
         NSPARR(2)=N+1
         N1=0
         DO I=1,N
            NOFFI=INT8(I-1)*INT8(I-2)/2
            LL=NOFFI/32+I 
            NJ=((I-1)*KBFW)/32
            DO K=0,NJ
               DO M=0,31
                 IF(BTEST(NBITS(LL+K),M)) NIN=NIN+1
               ENDDO
            ENDDO
            N1=N1+2
            NSPARR(N1+1)=NSPARR(1)+NIN
            NSPARR(N1+2)=NSPARR(2)+NIN
         ENDDO
         NDIMS(2)=NIN
         NDIMS(3)=NIN
         NTOT=NIN
*
      ENDIF
*
      NIN=NDIMS(3)+NDIMS(4)
      FRACZ=200.0*FLOAT(NTOT)/FLOAT(N)/FLOAT(N-1)
      FRACU=200.0*FLOAT(NIN)/FLOAT(N)/FLOAT(N-1)
      WRITE(*,*) ' '
      WRITE(*,*) 'NDBITS: number of diagonal elements',N
      WRITE(*,*) 'NDBITS: number of used off-diagonal elements',NIN
      WRITE(*,1000) 'fraction of non-zero off-diagonal elements', FRACZ
      WRITE(*,1000) 'fraction of used off-diagonal elements', FRACU
      IF (ICMPRS.NE.0) THEN
         CPR=100.0*FLOAT(NDIMS(2)+2*NDIMS(3)+NDIMS(4))/FLOAT(3*NIN)
      WRITE(*,1000) 'compression ratio for off-diagonal elements', CPR
      ENDIF
 1000 FORMAT(' NDBITS: ',A,F6.2,' %')
      RETURN 
      END

*     ------------------------------------------------------------------

      SUBROUTINE CKBITS(NDIMS,NBITS,MNPAIR,JCMPRS)           ! check sparse matrix
#include "largeint.inc"
      INTEGER(KIND=LARGE) NDIMB
      COMMON / MPBITS / NDIMB,N,IBFW,MXCNT,NENCDM,NENCDB,NTHRD
      INTEGER NBITS(*)
      INTEGER NWCP(0:1),IRGN(2),INR(2)
      INTEGER(KIND=LARGE) NDIMS(4),LL
      INTEGER*8 NOFFI
      LOGICAL BTEST

      DO I=1,4
         NDIMS(I)=0
      ENDDO
      ICMPRS=JCMPRS
      KBFW=1
      IF (IBFW.GT.1.AND.ICMPRS.GT.0) KBFW=2
      LL=0
*
      DO I=1,N
         NOFFI=INT8(I-1)*INT8(I-2)*INT8(IBFW)/2
         LL=NOFFI/32+I 
         MM=0
         INR(1)=0
         INR(2)=0
         IRGN(1)=0
         IRGN(2)=0
         LAST=0
         LRGN=0
         DO J=1,I-1 
            ICOUNT=0
            NEXT=0
            DO IB=0,IBFW-1
               IF (BTEST(NBITS(LL),MM)) ICOUNT=IBSET(ICOUNT,IB)
               MM=MM+1
               IF (MM.GE.32) THEN
                  LL=LL+1
                  MM=MM-32
               ENDIF
            ENDDO
*
            IF (ICOUNT.GT.0) NDIMS(1)=NDIMS(1)+1
*           keep pair ?
            IF (ICOUNT.GE.MNPAIR) THEN
               NEXT=1 ! double
               IF (ICOUNT.LE.ICMPRS.AND.ICMPRS.GT.0) NEXT=2 ! single 
               INR(NEXT)=INR(NEXT)+1
               IF (NEXT.NE.LAST.OR.LRGN.GE.NENCDM) THEN
                  IRGN(NEXT)=IRGN(NEXT)+1
                  LRGN=0
               ENDIF
               LRGN=LRGN+1
            ENDIF
            LAST=NEXT
         ENDDO
*
         IF (ICMPRS.NE.0) THEN
            DO JP=1,KBFW
               IF (INR(JP).GT.0) THEN
                  ICP=0
                  NWCP(0)=INR(JP)                    ! list of column indices (default)
                  NWCP(1)=IRGN(JP)+(IRGN(JP)+7)/8    ! list of regions of consecutive columns
*                 compress row ?
                  IF (NWCP(1).LT.NWCP(0)) ICP=1
                  NDIMS(2)   =NDIMS(2)   +NWCP(ICP)
                  NDIMS(JP+2)=NDIMS(JP+2)+NWCP(0)
               ENDIF
            ENDDO
         ELSE
            NDIMS(2)=NDIMS(2)+INR(1)
            NDIMS(3)=NDIMS(3)+INR(1)
         ENDIF
*
      ENDDO
*
      RETURN
      END 

*     ------------------------------------------------------------------
*     create sparsity information
      SUBROUTINE SPBITS(NSPARR,NSPARC,NBITS,NCMPRS)               ! collect elements
#include "largeint.inc"
      INTEGER(KIND=LARGE) NDIMB
      COMMON / MPBITS / NDIMB,N,IBFW,MXCNT,NENCDM,NENCDB,NTHRD
      INTEGER NBITS(*),NSPARC(*),NCMPRS(*)
      INTEGER(KIND=LARGE) NSPARR(*),KL,L,LL,L1,K8,N1
      INTEGER*8 NOFFI
      LOGICAL BTEST
*
      ICHUNK=MIN(N/NTHRD/8+1,256)
*
      DO JB=0,IBFW-1
*
!$OMP  PARALLEL DO                 ! parallelize row loop
!$OMP& PRIVATE(I,N1,NOFFI,L,M,KL,L1,NRGN,NRGN8,K8,
!$OMP&         LAST,LRGN,LL,J1,JN,J,NEXT)
!$OMP& SCHEDULE(DYNAMIC,ICHUNK)
       DO I=1,N
        N1=2*(I+JB*(N+1))
        NOFFI=INT8(I-1)*INT8(I-2)*INT8(IBFW)/2
        L=NOFFI/32+I 
        M=JB
        KL=NSPARR(N1-1)-1  ! pointer to row in NSPARC
        L1=NSPARR(N1  )    ! pointer to row in sparse matrix
        NRGN=NCMPRS(I+N*JB)! compression  (number of consecutive regions)
        NRGN8=(NRGN+7)/8   ! number of groups (1 offset per group)
        K8=KL
        KL=KL+NRGN8        ! reserve space of offsets
        LAST=0
        LRGN=0
        LL=L1-1
        J1=0
        JN=0
*
        DO J=1,I-1         !  loop for off-diagonal elements  
         NEXT=0
         IF(NBITS(L).NE.0) THEN 
          IF(BTEST(NBITS(L),M)) THEN
           LL=LL+1
           IF (NRGN.LE.0) THEN
            KL=KL+1      
            NSPARC(KL)=J ! column index
           ELSE
            NEXT=1
            IF (LAST.EQ.0.OR.JN.GE.NENCDM) THEN
             IF (MOD(LRGN,8).EQ.0) THEN
              K8=K8+1
              NSPARC(K8)=LL-L1
             ENDIF
             LRGN=LRGN+1
             KL=KL+1
             J1=ISHFT(J,NENCDB)
             JN=0
            ENDIF
            JN=JN+1
            NSPARC(KL)=IOR(J1,JN)
           ENDIF 
          END IF
         END IF
         LAST=NEXT
         M=M+IBFW
         IF (M.GE.32) THEN
          M=M-32
          L=L+1
         ENDIF 
*
        END DO
       END DO
!$OMP END PARALLEL DO
*
      ENDDO
*
      N1=2*(N+1)*IBFW
      WRITE(*,*) ' '
      WRITE(*,*) 'SPBITS: sparse structure constructed ',NSPARR(N1-1),
     &           ' words'
      WRITE(*,*) 'SPBITS: dimension parameter of matrix',NSPARR(2)-1
      IF (IBFW.LE.1) THEN
         WRITE(*,*) 'SPBITS: index of last used location',NSPARR(N1)-1
      ELSE
         WRITE(*,*) 'SPBITS: index of last used double',NSPARR(N1/2)-1
         WRITE(*,*) 'SPBITS: index of last used single',NSPARR(N1)-1
      ENDIF
c      WRITE(*,101) (NSPARC(I),I=1,KL)
c 101  FORMAT(10I6)  
      END

* !!! not up to date !!! see IJADD 
*      FUNCTION IJIND(ITEMA,ITEMB,NSPARR,NSPARC)      ! index for element
*#include "largeint.inc"
*      INTEGER(KIND=LARGE) NSPARR(*), IJIND,K,KK,KL,KU,INDIJ,ND
*      INTEGER NSPARC(*)
*      SAVE
*     ...
*      IJIND=0
*      ND=NSPARR(1)-1                  ! dimension of matrix
*      ITEM1=MAX(ITEMA,ITEMB)          ! larger index
*      ITEM2=MIN(ITEMA,ITEMB)          ! smaller index
*      IF(ITEM2.LE.0.OR.ITEM1.GE.NSPARR(1)) RETURN
*      IF(ITEM1.EQ.ITEM2) THEN         ! diagonal element
*         IJIND=ITEM1
*      ELSE                            ! off-diagonal element
*         KK=NSPARR(ITEM1)-ND
*         KL=0
*         KU=NSPARR(ITEM1+1)-ND-1-KK
*         INDIJ=KK
*         IF(KU.LT.KL) THEN
c            WRITE(*,*) ITEM1,ITEM2,' not found',KK,KL,KU
*            RETURN
*         END IF
* 10      K=(KL+KU)/2                  ! binary search
*         IF(ITEM2.EQ.NSPARC(INDIJ+K)) THEN
*            IJIND=ND+K
*         ELSE
*            IF(ITEM2.LT.NSPARC(INDIJ+K)) THEN
*               KU=K-1
*            ELSE IF(ITEM2.GT.NSPARC(INDIJ+K)) THEN
*               KL=K+1
*            END IF
*            IF(KL.LE.KU) GOTO 10
c            WRITE(*,*) ITEM1,ITEM2,' not found',KK,KL,KU
*            RETURN
*         END IF
*         IJIND=ND+K+KK
*      END IF
*      END

