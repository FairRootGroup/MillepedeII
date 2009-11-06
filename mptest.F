

      SUBROUTINE MPTEST                          ! generate test files
*     MC for simple 100-plane chamber, no B-field, straight tracks
*
*     unit  8: textfile mp2str.txt   = steering file
*     unit  9: textfile mp2con.txt   = constraint file
*     unit 51: binary file mp2test.bin, written using CALL MILLE(.)
*     existing file are removed
*

      PARAMETER (NPLAN=100)
      COMMON/TDETEC/DEL(NPLAN),DVD(NPLAN)  ! misalignment
      COMMON/CDETEC/DETX,DISX,THCK,HEIT,
     +              EFF(NPLAN),SGM(NPLAN),
     +              YNULL,SLOPE,IHITS(NPLAN),YDRFT(NPLAN), 
     +              YNULLM,SLOPEM,
     +              NHITS,XHITS(NPLAN),YHITS(NPLAN),    
     +              YSHFT(NPLAN),DYSHFT(NPLAN),
     +              VDCOR(NPLAN),DVDCOR(NPLAN),SIGMA(NPLAN)

      DOUBLE PRECISION S1,S2,SW,SV,SUM1,SUM2     
      REAL DERLC(2),DERGL(2)
      INTEGER LABEL(2)
      LOGICAL EX1,EX2,EX3
*     ...
CCC      CALL RNTIME
      INQUIRE(FILE='mp2str.txt',IOSTAT=IOS,EXIST=EX1) ! keep, if existing
      INQUIRE(FILE='mp2con.txt',IOSTAT=IOS,EXIST=EX2) ! keep, if existing

      INQUIRE(FILE='mp2tst.bin',IOSTAT=IOS,EXIST=EX3) ! remove, if existing

      WRITE(*,*) ' '
      WRITE(*,*) 'Generating test data for mp II...'
      WRITE(*,*) ' '
*     file management
      IF(EX3) CALL SYSTEM('rm mp2tst.bin')   ! remove old file

      IF(.NOT.EX1) 
     +OPEN(UNIT=7,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     +     FILE='mp2str.txt')
      IF(.NOT.EX2) 
     +OPEN(UNIT=9,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     +     FILE='mp2con.txt')
      OPEN(UNIT=51,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     +     FILE='mp2tst.bin')  

*     define detector geometry

      DETX= 10.0            ! x-value of first plane
      DISX= 10.0            ! distance between planes
      THCK=  2.0            ! thickness of plane
      HEIT=100.0            ! height of detector plane
      DO I=1,NPLAN
       EFF(I)=0.90          ! plane efficiency
       SGM(I)=0.0150        ! measurement sigma
       DEL(I)=0.0           ! true shift is zero
       YSHFT(I)=0.0         ! measured shift is zero
       VDCOR(I)=0.0         ! drift velocity correction
      END DO

      IPL=7                 ! modify one plane (7)
      EFF(IPL)=0.1          ! low efficiency  
      SGM(IPL)=0.0400       ! bad resolution

*     misalign detector planes -----------------------------------------

      DISPL=0.1                        ! displacement 1 mm * N(0,1) 
      DRIFT=0.02                       ! Vdrift deviation 2 %  * N(0,1) 
      DO I=1,NPLAN 
       DEL(I)=DISPL*GRAN()             ! shift               
       DVD(I)=DRIFT*GRAN()             ! rel. drift velocitu deviation
      END DO
      DEL(10)=0.0                      ! no shift 
      DEL(90)=0.0                      ! no shift

*     write text files -------------------------------------------------

      IF(EX1) GOTO 10       
      LUNS=7                           ! steerfile
      WRITE(LUNS,101) '*            Default test steering file'
      WRITE(LUNS,101) 'fortranfiles ! following bin files are fortran'
      WRITE(LUNS,101) 'mp2con.txt   ! constraints text file '
      WRITE(LUNS,101) 'mp2tst.bin   ! binary data file'
      WRITE(LUNS,101) 'Cfiles       ! following bin files are Cfiles'  
C      WRITE(LUNS,101) '*outlierrejection 100.0 ! reject if Chi^2/Ndf >'
C      WRITE(LUNS,101) '*outliersuppression 3   ! 3 local_fit iterations'
      
      WRITE(LUNS,101)
     +  '*hugecut 50.0     !cut factor in iteration 0'
      WRITE(LUNS,101)
     +  '*chisqcut 1.0 1.0 ! cut factor in iterations 1 and 2' 
      WRITE(LUNS,101)
     +  '*entries  10 ! lower limit on number of entries/parameter'
      WRITE(LUNS,101)
     +  '*pairentries 10 ! lower limit on number of parameter pairs',
     +  '                ! (not yet!)'
      WRITE(LUNS,101) 
     +  '*printrecord   1  2      ! debug printout for records'
      WRITE(LUNS,101)
     +  '*printrecord  -1 -1      ! debug printout for bad data records' 
      WRITE(LUNS,101) 
     +  '*outlierdownweighting  2 ! number of internal iterations (> 1)'
      WRITE(LUNS,101)
     +  '*dwfractioncut      0.2  ! 0 < value < 0.5'
      WRITE(LUNS,101)
     +  '*presigma           0.01 ! default value for presigma'
      WRITE(LUNS,101)
     +  '*regularisation 1.0      ! regularisation factor'
      WRITE(LUNS,101)  
     +  '*regularisation 1.0 0.01 ! regularisation factor, pre-sigma' 

      WRITE(LUNS,101) ' '
      WRITE(LUNS,101)
     +     '*bandwidth 0         ! width of precond. band matrix'
      WRITE(LUNS,101)
     +     'method diagonalization 3 0.001 ! diagonalization      '
      WRITE(LUNS,101)
     +     'method fullMINRES       3 0.01 ! minimal residual     '
      WRITE(LUNS,101)
     +     'method sparseMINRES     3 0.01 ! minimal residual     '
      WRITE(LUNS,101)
     +     '*mrestol      1.0D-8          ! epsilon for MINRES'
      WRITE(LUNS,101)
     +     'method inversion       3 0.001 ! Gauss matrix inversion'
      WRITE(LUNS,101) '* last method is applied'
      WRITE(LUNS,101)
     +     '*matiter      3  ! recalculate matrix in iterations'
      WRITE(LUNS,101) ' '
      WRITE(LUNS,101) 'end ! optional for end-of-data'
 10   CONTINUE

      LUNT=9                           ! constraint file
      ONE=1.0                          ! shift constraint
      IF(.NOT.EX2) WRITE(LUNT,*) 'Constraint  0.0'
      DO I=1,NPLAN  
       LABELT=10+I*2
       X=DETX+FLOAT(I-1)*DISX+0.5*THCK
       IF(.NOT.EX2) WRITE(LUNT,103) LABELT,ONE
      END DO

      SW=0.0D0                         ! tilt constraint 
      SV=0.0D0
      S1=0.0D0
      S2=0.0D0
      IF(.NOT.EX2) WRITE(LUNT,*) 'Constraint 0.0'   ! write
      DBAR=0.5*FLOAT(NPLAN-1)*DISX
      XBAR=DETX+0.5*FLOAT(NPLAN-1)*DISX! +0.5*THCK
      DO I=1,NPLAN  
       LABELT=10+I*2
       X=DETX+FLOAT(I-1)*DISX          !+0.5*THCK
       WW=(X-XBAR)/DBAR
       IF(.NOT.EX2) WRITE(LUNT,103) LABELT,WW          ! write
       S1=S1+DEL(I)
       S2=S2+WW*DEL(I)
       SW=SW+WW
       SV=SV+WW*WW 
      END DO
 

      DET=FLOAT(NPLAN)*SV-SW*SW
      EPS=(SV*S1-SW*S2)/DET  
      ETA=(FLOAT(NPLAN)*S2-SW*S1)/DET
      DO I=1,NPLAN
       X=DETX+FLOAT(I-1)*DISX 
       WW=(X-XBAR)/DBAR
       DEL(I)=DEL(I)-EPS-ETA*WW        ! correct displacement ...
      END DO                           ! ... for constraints

      SUM1=0.0
      SUM2=0.0
      DO I=1,NPLAN
       SUM1=SUM1+DEL(I) 
       X=DETX+FLOAT(I-1)*DISX          !+0.5*THCK
       WW=(X-XBAR)/DBAR
       SUM2=SUM2+DEL(I)*WW
      END DO
c      WRITE(*,*) '   Check for constraints ',SUM1,SUM2

*     record loop ------------------------------------------------------

      NCOUNT=10000
      NTHITS=0
      NRECDS=0

      DO ICOUNT=1,NCOUNT
       IP=0
       IF(ICOUNT.EQ.8759) IP=1
c       IF(ICOUNT.EQ.6309) IP=1
c       IF(ICOUNT.EQ.7468) IP=1
       CALL GENLIN(IP)      ! generate hits

       DO I=1,NHITS
        DERLC(1)=1.0
        DERLC(2)=XHITS(I)
        DERGL(1)=1.0
        DERGL(2)=YDRFT(I)
        LABEL(1)=10+IHITS(I)*2
        LABEL(2)=500 + IHITS(I)
        CALL MILLE(2,DERLC,2,DERGL,LABEL,YHITS(I),SIGMA(I))
        NTHITS=NTHITS+1  ! count hits
       END DO
       CALL ENDLE
       NRECDS=NRECDS+1   ! count records
      END DO 

*     ------------------------------------------------------------------
      IF(.NOT.EX1) THEN
         REWIND  (7)
         CLOSE   (7)
      END IF
      IF(.NOT.EX2) THEN 
         REWIND  (9)
         CLOSE   (9)
      END IF
      REWIND (51)
      CLOSE  (51)

c      WRITE(*,*) ' '
c      WRITE(*,*) 'Shifts and drift velocity deviations:'
c      DO I=1,NPLAN
c       WRITE(*,102) I,DEL(I),DVD(I)
c      END DO 


      WRITE(*,*) ' '      
      WRITE(*,*) ' '
      WRITE(*,*) NCOUNT,' tracks generated with ',NTHITS,' hits.'
      WRITE(*,*) NRECDS,' records written.' 
      WRITE(*,*) ' '
 101  FORMAT(A)
c 102  FORMAT(I6,2F10.5)
 103  FORMAT(I8,F10.5)
      END

      SUBROUTINE GENLIN(IP)                        ! generate line and meas
      PARAMETER (NPLAN=100)
      COMMON/TDETEC/DEL(NPLAN),DVD(NPLAN)  ! misalignment
      COMMON/CDETEC/DETX,DISX,THCK,HEIT,
     +              EFF(NPLAN),SGM(NPLAN),
     +              YNULL,SLOPE,IHITS(NPLAN),YDRFT(NPLAN), 
     +              YNULLM,SLOPEM,
     +              NHITS,XHITS(NPLAN),YHITS(NPLAN),    
     +              YSHFT(NPLAN),DYSHFT(NPLAN),
     +              VDCOR(NPLAN),DVDCOR(NPLAN),SIGMA(NPLAN)
*     ...
      YNULL=0.5*HEIT+0.1*HEIT*(URAN()-0.5)   ! uniform vertex
      SLOPE=(URAN()-0.5)*HEIT/(FLOAT(NPLAN-1)*DISX) 
      IF(IP.NE.0) THEN
         WRITE(*,*) ' '
C         WRITE(*,*) 'YNULL=',YNULL,'    SLOPE=',SLOPE
      END IF
      NHITS=0
      DO I=1,NPLAN
       X=DETX+FLOAT(I-1)*DISX  !  +0.5*THCK
       IF(URAN().LT.EFF(I)) THEN
          YLIN        =YNULL+SLOPE*X             ! true y value
          YBIAS       =YLIN-DEL(I)               ! biased value
          NWIRE=1.0+YBIAS/4.0                    ! wire number
          IF(NWIRE.LE.0.OR.NWIRE.GT.25) GOTO 10  ! check wire number
          NHITS=NHITS+1                          ! track hits the plane
          XHITS(NHITS)=X
          IHITS(NHITS)=I
          GR=GRAN()
          YMEAS=SGM(I)*GR
          YDVDS=0.0
          YHITS(NHITS)=YBIAS+YMEAS+YDVDS     ! measured
          YWIRE=FLOAT(NWIRE)*4.0-2.0 
          YDRFT(NHITS)=YBIAS-YWIRE           ! signed drift length
          YDVDS=YDRFT(NHITS)*DVD(I)
          YHITS(NHITS)=YBIAS+YMEAS-YDVDS     ! measured
          SIGMA(NHITS)=SGM(I)
          IF(IP.NE.0) THEN
c             WRITE(*,101) NHITS,I,X,YLIN,YBIAS,YMEAS,
c     +       SGM(I),YHITS(NHITS),GR,DEL(I)
          END IF
       END IF
      END DO 
 10   CONTINUE
c 101  FORMAT(2I3,F5.0,7F8.4)
      END

