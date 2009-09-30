

      SUBROUTINE MPTST2(IMODEL)         ! generate test files
*     MC for simple 10-layer silicon strip tracker, 5 modules per layer (2*10cm)
*                   10 cm spacing, no B-field,
*                   layers 1,4,7,10 have additional +/-5deg stereo modules
*                   intrinsic resolution 20mu, 2% X0 per strip module, 
*                   uniform track offsets/slopes, momentum given  
*     different track models IMODEL
*           0: 'straight line', ignoring multiple scattering
*           1: 'straight line', using diagonal of m.s. error matrix
*           2: 'break points' 
*           3: 'broken lines', fine   
*           4: 'broken lines', coarse (stereo layers combined)   
*
*     unit  8: textfile mp2str.txt   = steering file
*     unit  9: textfile mp2con.txt   = constraint file
*     unit 51: binary file mp2test.bin, written using CALL MILLE(.)
*     existing file are removed
*
      PARAMETER (NLYR=10)                  ! #layers
      PARAMETER (NMLYR=14)                 ! #module layers
      PARAMETER (NMX=10)                    ! #module in X
      PARAMETER (NMY=5)                    ! #module in Y 
      PARAMETER (NTOT=NLYR*NMX*NMY) 
      COMMON/SDEVEC/SDEVX(NTOT),SDEVY(NTOT)   
      COMMON/SDETEC/SARC(NMLYR),SSIG(NMLYR),SPRO(2,NMLYR),
     +      ISLYR(NMLYR),SIZE,THE0, 
     +      NHITS,IHITS(NMLYR),XHITS(NMLYR),YHITS(NMLYR),SIGMA(NMLYR) 

      REAL DERLC(NMLYR*2+3),DERGL(NMLYR*2+3)
      INTEGER LABEL(2)
      LOGICAL EX1,EX2,EX3
*     for broken lines: 1=fine, 2=coarse      
      DIMENSION NBRL(2),LBRL(NMLYR,2),SBRL(NMLYR,2),WBRL(NMLYR,2),
     +          CMBBRL(2)
      DATA CMBBRL / 0.0, 1.0 / ! cut for combining layers 
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

      DETS= 10.0            ! arclength of first plane
      DISS= 10.0            ! distance between planes
      THCK= 0.02            ! thickness of plane (X0)
      OFFS=  0.5            ! offset of stereo modules
      STEREO=0.08727        ! stereo angle
      SIZE= 20.0            ! size of layers
      
      S=DETS
      I=0
      DO LAYER=1,10
       I=I+1
       ISLYR(I)=LAYER       ! layer
       SARC(I)=S            ! arclength
       SSIG(I)=0.002        ! resolution
       SPRO(1,I)=1.0        ! module measures 'X'
       SPRO(2,I)=0.0
       IF (MOD(LAYER,3).EQ.1) THEN
        I=I+1
        ISLYR(I)=LAYER       ! layer
        SARC(I)=S+OFFS      ! arclength stereo module
        SSIG(I)=0.002       ! resolution 
       SPRO(1,I)=SQRT(1.0-STEREO**2)
       SPRO(2,I)=STEREO    ! module measures both 'X' and 'Y'
       STEREO=-STEREO     
       ENDIF
       S=S+DISS
      END DO
      
* define broken lines
      SOLD=-1000.
      NBRL(1)=0
      NBRL(2)=0 
      DO K=1,2
       DO I=1, NMLYR
        IF (ABS(SARC(I)-SOLD).GT.CMBBRL(K)) NBRL(K)=NBRL(K)+1 
       LB=NBRL(K)
       LBRL(I,K)=LB
       SBRL(LB,K)=SBRL(LB,K)+SARC(I)
       WBRL(LB,K)=WBRL(LB,K)+1.0 
        SOLD=SARC(I)
       END DO
       DO I=1,NBRL(K)
        SBRL(I,K)=SBRL(I,K)/WBRL(I,K)
       WBRL(I,K)=SQRT(WBRL(I,K))
       END DO 
      END DO
      IBRL=IMODEL-2 
      
*     misalign detector modules -----------------------------------------

      DISPXM=0.01           ! module displacement in X .05 mm * N(0,1) 
      DISPYM=0.01           ! module displacement in Y .05 mm * N(0,1)
             
      DO I=0,NLYR-1
       DO K=0,NMY-1
         DO L=1,NMX
           SDEVX(((I*NMY+K)*NMX+L))=DISPXM*GRAN()     ! shift in x
           SDEVY(((I*NMY+K)*NMX+L))=DISPYM*GRAN()     ! shift in y
        ENDDO
       ENDDO
      END DO
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
*
* constraints: fix center modules in first/last layer
*
      NCX=(NMX+1)/2
      NMXY=NMX*NMY
      LUNT=9
      ONE=1.0
      DO I=1,NLYR,NLYR-1
       IF(.NOT.EX2) WRITE(LUNT,*) 'Constraint  0.0'
       DO K=0,NMY-1
        LABELT=(I*NMY+K)*NMX+NCX-1
        IF(.NOT.EX2) WRITE(LUNT,103) LABELT,ONE
       SDEVX(((I-1)*NMY+K)*NMX+NCX)=0.0      ! fix center modules at 0.
       END DO 
       IF(.NOT.EX2) WRITE(LUNT,*) 'Constraint  0.0'
       DO K=0,NMY-1
        LABELT=(I*NMY+K)*NMX+NCX+1000-1
        IF(.NOT.EX2) WRITE(LUNT,103) LABELT,ONE
       SDEVY(((I-1)*NMY+K)*NMX+NCX)=0.0      ! fix center modules at 0.
       END DO
      END DO 
      
*     record loop ------------------------------------------------------

      NCOUNT=10000
      NTHITS=0
      NRECDS=0

      DO ICOUNT=1,NCOUNT
*      10..100 GeV      
       P=10.0**(1.+URAN())
       THE0=SQRT(THCK)*0.014/P         
       IP=0
C       IF (ICOUNT.LE.3) IP=1
       CALL GENLN2(IP)      ! generate hits

*        
       DO I=1,NHITS 
* simple straight line  
        LYR=IHITS(I)/NMXY+1
        IM =MOD(IHITS(I),NMXY)
        NALC=4
        DERLC(1)=SPRO(1,LYR)
        DERLC(2)=SPRO(2,LYR)
        DERLC(3)=XHITS(I)*SPRO(1,LYR)
        DERLC(4)=XHITS(I)*SPRO(2,LYR)
        DERGL(1)=SPRO(1,LYR)
        DERGL(2)=SPRO(2,LYR)
        LABEL(1)=IM+NMXY*ISLYR(LYR)
        LABEL(2)=IM+NMXY*ISLYR(LYR)+1000
* add multiple scattering errors (no correlations)
        IF (IMODEL.EQ.1) THEN
        DO J=I,NHITS
         SIGMA(J)=SQRT(SIGMA(J)**2+((XHITS(J)-XHITS(I))*THE0)**2)
        END DO
       END IF  
* add 'break points' for multiple scattering
        IF (IMODEL.EQ.2.AND.I.GT.1) THEN
        DO J=1,I-1
* 2 scattering angles from each layer in front of current        
         NALC=NALC+1
         DERLC(NALC)=(XHITS(I)-XHITS(J))*SPRO(1,LYR)
         NALC=NALC+1
         DERLC(NALC)=(XHITS(I)-XHITS(J))*SPRO(2,LYR)
        END DO 
       END IF  
* add 'broken lines' offsets for multiple scattering
        IF (IMODEL.GE.3) THEN
         NALC=2*NBRL(IBRL)
        DO K=1, NALC
         DERLC(K)=0.0
        END DO
* 2 offsets
         LB=LBRL(LYR,IBRL)    
        DERLC(LB*2-1)=SPRO(1,LYR)
        DERLC(LB*2  )=SPRO(2,LYR)
       END IF  
*       
        CALL MILLE(NALC,DERLC,2,DERGL,LABEL,YHITS(I),SIGMA(I))
        NTHITS=NTHITS+1  ! count hits
       END DO
* additional measurements from MS
       IF (IMODEL.EQ.2) THEN
        DO I=1,(NHITS-1)*2
        NALC=I+4
         DO K=1,NALC
         DERLC(K)=0.0
        END DO
        DERLC(NALC)=1.0     
        CALL MILLE(NALC,DERLC,0,DERGL,LABEL,0.0,THE0)
       END DO
       END IF
*       
       IF (IMODEL.GE.3) THEN
        DO I=2,NBRL(IBRL)-1
        DP=1.0/(SBRL(I,IBRL)-SBRL(I-1,IBRL))
        DN=1.0/(SBRL(I+1,IBRL)-SBRL(I,IBRL))
        NALC=(I+1)*2
        DO L=-1,0
          DO K=1,NALC
          DERLC(K)=0.0
         END DO
         DERLC(2*(I-1)+L)= DP     
         DERLC(2* I   +L)=-DP-DN     
         DERLC(2*(I+1)+L)= DN       
         CALL MILLE(NALC,DERLC,0,DERGL,LABEL,0.0,THE0*WBRL(I,IBRL))
        END DO 
       END DO
       END IF             
*             
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

      SUBROUTINE GENLN2(IP)                        ! generate line and meas
*     
      PARAMETER (NLYR=10)                  ! #layers
      PARAMETER (NMLYR=14)                 ! #module layers
      PARAMETER (NMX=10)                    ! #module in X
      PARAMETER (NMY=5)                    ! #module in Y 
      PARAMETER (NTOT=NLYR*NMX*NMY) 
      COMMON/SDEVEC/SDEVX(NTOT),SDEVY(NTOT)   
      COMMON/SDETEC/SARC(NMLYR),SSIG(NMLYR),SPRO(2,NMLYR),
     +      ISLYR(NMLYR),SIZE,THE0, 
     +      NHITS,IHITS(NMLYR),XHITS(NMLYR),YHITS(NMLYR),SIGMA(NMLYR)
*     track parameters      
      XNULL=SIZE*(URAN()-0.5)   ! uniform vertex
      YNULL=SIZE*(URAN()-0.5)   ! uniform vertex   
      XEXIT=SIZE*(URAN()-0.5)   ! uniform exit point
      YEXIT=SIZE*(URAN()-0.5)   ! uniform exit point
      XSLOP=(XEXIT-XNULL)/SARC(NMLYR)
      YSLOP=(YEXIT-YNULL)/SARC(NMLYR)      
      IF(IP.NE.0) THEN
         WRITE(*,*) ' '
         WRITE(*,*) ' Track ', XNULL, YNULL, XSLOP, YSLOP
      END IF
*  
      NHITS=0
      X=XNULL
      Y=YNULL
      DX=XSLOP
      DY=YSLOP
      SOLD=0.0
*      
      DO 10 I=1,NMLYR
       DS=SARC(I)-SOLD
       SOLD=SARC(I)
*      position with parameters 1. hit       
       XS=XNULL+SARC(I)*XSLOP
       YS=YNULL+SARC(I)*YSLOP    
*      true track position       
       X=X+DX*DS
       Y=Y+DY*DS
*      multiple scattering
       DX=DX+GRAN()*THE0
       DY=DY+GRAN()*THE0
*             
       IMX=IFIX((X+SIZE*0.5)/SIZE*FLOAT(NMX))
       IF (IMX.LT.0.OR.IMX.GE.NMX) GOTO 10
       IMY=IFIX((Y+SIZE*0.5)/SIZE*FLOAT(NMY))
       IF (IMY.LT.0.OR.IMY.GE.NMY) GOTO 10
*       
       IHIT=((I-1)*NMY+IMY)*NMX+IMX
       IOFF=((ISLYR(I)-1)*NMY+IMY)*NMX+IMX+1
       NHITS=NHITS+1
       IHITS(NHITS)=IHIT
       XL=X-SDEVX(IOFF)
       YL=Y-SDEVY(IOFF)
       XHITS(NHITS)=SARC(I)   
       YHITS(NHITS)=(XL-XS)*SPRO(1,I)+(YL-YS)*SPRO(2,I)+GRAN()*SSIG(I)
       SIGMA(NHITS)=SSIG(I)
*       
       IF(IP.NE.0) THEN
        WRITE(*,101) NHITS,I,IHIT,X,Y,XHITS(NHITS),
     +     YHITS(NHITS),SIGMA(NHITS)
       END IF
 10   CONTINUE
 101  FORMAT(3I3,5F8.4)
      END
