

*
*     Program to present Millepede histograms
*   
      CHARACTER*60 FILNAM,KEYSTX
      CHARACTER*80 TEXT

      DOUBLE PRECISION DNUM
      COMMON/NUMDUM/DNUM(100),NUMS

      PARAMETER (NDIM=100)
      COMMON/CNUMS/LAST,KEYA,KEYB,ICD(2,NDIM)
 
*     ...
      LUN=7
      IAR=1
      IF(IARGC().EQ.0) IAR=0

      DO I=IAR,IARGC()
       IF(I.EQ.0) THEN
          FILNAM='millepede.his'
       ELSE 
          CALL GETARG(I,FILNAM)
       END IF
       CALL RLTEXT(FILNAM,IA,IB,NAB)
       LUN=9
       WRITE(*,*) ' ' 
       WRITE(*,*) 'File ',FILNAM(1:NAB) 
       WRITE(*,*) ' '
       OPEN(UNIT=LUN,FILE=FILNAM(1:NAB),IOSTAT=IOS)
       IF(IOS.NE.0) THEN
          WRITE(*,*) 'OPEN error ',IOS,'/File ignored'
          GOTO 90 
       END IF

 10    READ(LUN,101,END=80) TEXT
       CALL RLTEXT(TEXT,IA,IB,NAB)        ! indices for non-blank area
       IF(NAB.EQ.0) GOTO 10               ! skip this line

       CALL RATEXT(TEXT(1:NAB),NUMS,DNUM) ! translate text to DP numbers
       IF(KEYB.EQ.0) GOTO 10              ! skip, no text present
       
       KEYSTX='Histogram'
       MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
       IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
          CALL RATEXT(TEXT(1:NAB),NUMS,DNUM)
          CALL HISTER(LUN)
          GOTO 10
       END IF 
      
       KEYSTX='XY-Data'
       MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
       IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
          CALL RATEXT(TEXT(1:NAB),NUMS,DNUM)
          CALL XYDATA(LUN)
          GOTO 10 
       END IF

       DO J=1,4
        IF(J.EQ.1) KEYSTX='Parameter'
        IF(J.EQ.2) KEYSTX='Constraint'
        IF(J.EQ.3) KEYSTX='WConstraint'
        IF(J.EQ.4) KEYSTX='Measurement'
        MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
        IF(MAT.GE.(NPAT+NPAT+3)/3) THEN 
           REWIND LUN
           CALL MPDATA(LUN)
           GOTO 80
        END IF
       END DO
       GOTO 10

 80    CLOSE(UNIT=LUN) 
 90   END DO
 101  FORMAT(A)
      END


      SUBROUTINE HISTER(LUN)
      CHARACTER*80 TEXT, KEYSTX*24
      DOUBLE PRECISION DNUM
      COMMON/NUMDUM/DNUM(100),NUMS

      PARAMETER (NDIM=100)
      COMMON/CNUMS/LAST,KEYA,KEYB,ICD(2,NDIM)

*     Histogram
      COMMON/HISTEX/IH,IV,IT,NBIN,IBIN,INOUT(3),XA,XB,HIST(120),
     +              XMIN,XMAX,XMEAN,XSIGM 
      CHARACTER*60 COMT
      COMMON/HSTEXT/COMT    

*     ...
      IH=DNUM(1) ! histogram number
      IV=DNUM(2) ! version of histogram
      IT=DNUM(3) ! histogram type
      WRITE(*,*) ' '
      WRITE(*,*) 'Histogram ',IH,'   Version',IV,'    (type',IT,')'
      WRITE(*,*) '============'
      NBIN=0     ! reset
      IBIN=0
      XMIN=0.0
      XMAX=0.0
      XMEAN=0.0
      XSIGM=0.0
      READ(LUN,101,END=20) COMT                  ! read comment
      CALL RLTEXT(COMT,IA,IB,NAB)
      WRITE(*,*) COMT(IA:IB)

 10   READ(LUN,101,END=20) TEXT
      CALL RLTEXT(TEXT,IA,IB,NAB)        ! indices for non-blank area
      IF(NAB.EQ.0) GOTO 10
      CALL RATEXT(TEXT(1:NAB),NUMS,DNUM) ! translate text to DP numbers
       
      IF(KEYB.NE.0) THEN        ! text present

         KEYSTX='bins, limits'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN             ! positive 
            NBIN=DNUM(1)            
            IBIN=0
            BINX=0.0
            BANX=0.0
            XA=DNUM(2)
            XB=DNUM(3)
            WRITE(*,*) NBIN,' bins' 
            WRITE(*,103) 'from  ',XA,' to  ',XB
            GOTO 10
         END IF


         KEYSTX='out-low inside out-high'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN             ! positive 
            INOUT(1)=DNUM(1)
            INOUT(2)=DNUM(2)
            INOUT(3)=DNUM(3)
            GOTO 10
         END IF

         KEYSTX='bincontent'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+PAT+3)/3) THEN             ! positive
            IBIN=0
            GOTO 10
         END IF

         KEYSTX='minmax'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN             ! positive 
            XMIN=DNUM(1)
            XMAX=DNUM(2)
            WRITE(*,103) 'Min =',XMIN,'Max =',XMAX
            GOTO 10
         END IF

         KEYSTX='meansigma'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN             ! positive 
            XMEAN=DNUM(1)
            XSIGM=DNUM(2)
            WRITE(*,103) 'Mean =',XMEAN,'sigma =',XSIGM
            GOTO 10
         END IF

         KEYSTX='end of histogram'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN             ! positive
            IF(BINX.EQ.BANX) RETURN
            CALL PFVERT(IBIN,HIST)
            CALL PSVERT(XA,XB)
            WRITE(*,102) ' '
            WRITE(*,102) 'Outlow:',INOUT(1),'inside:',INOUT(2),
     +           'outhigh:',INOUT(3)
            CALL PHISTE
            RETURN
         END IF
      ELSE 
         DO I=1,NUMS
          IBIN=IBIN+1
          HIST(IBIN)=DNUM(I)
          BINX=MIN(BINX,HIST(IBIN))
          BANX=MAX(BANX,HIST(IBIN))  
         END DO
      END IF
      GOTO 10
 20   CONTINUE 
 101  FORMAT(A)      
 102  FORMAT(A8,I10,6X,A8,I10,6X,A8,I10)
 103  FORMAT(3X,A8,G15.5,6X,A8,G15.5)
      END




      SUBROUTINE XYDATA(LUN)
      CHARACTER*80 TEXT,KEYSTX*24
      INTEGER INOUT(2)
      CHARACTER*2 ML 
      DOUBLE PRECISION DNUM
      COMMON/NUMDUM/DNUM(100),NUMS

      PARAMETER (NDIM=100)
      COMMON/CNUMS/LAST,KEYA,KEYB,ICD(2,NDIM)

      CHARACTER*60 COMT
      COMMON/HSTEXT/COMT

*     ...
      ID=DNUM(1) ! XY-Data number
      IV=DNUM(2) ! version of XY-Data
      IT=DNUM(3) ! data type
      WRITE(*,*) ' '  
      WRITE(*,*) 'XY-Data ',ID,'   Version',IV,'    (type',IT,')'
      WRITE(*,*) '=========='
      NPOINT=0
      ML='m'
      IF(IT.EQ.2) ML='l'
      IF(IT.EQ.3) ML='ml'
      READ(LUN,101,END=90) COMT                  ! read comment

 10   READ(LUN,101,END=90) TEXT
      CALL RLTEXT(TEXT,IA,IB,NAB)        ! indices for non-blank area
      IF(NAB.EQ.0) GOTO 10
      CALL RATEXT(TEXT(1:NAB),NUMS,DNUM) ! translate text to DP numbers
      IF(KEYB.EQ.0) THEN
         GOTO 10 
      ELSE IF(KEYB.NE.0) THEN        ! text present
         KEYSTX='stored  not-stored'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
            INOUT(1)=DNUM(1)
            INOUT(2)=DNUM(2)
            WRITE(*,*) 'Stored points:',INOUT(1),
     +                 '      not stored:',INOUT(2)
            IF(INOUT(1).GT.1) THEN
               CALL RLTEXT(COMT,IA,IB,NAB)
               WRITE(*,*) COMT(IA:IB)
               WRITE(*,*) ' '
#ifdef graphics
               CALL QFINIT(0.0,0.0)
               CALL TXTNUM(FLOAT(ID))
               CALL TXTNUM(FLOAT(IV))
               TEXT='XY-Data @/@: '
               TEXT(14:15+IB-IA)=COMT(IA:IB)
               CALL RLTEXT(TEXT,IA,IB,NAB)
               CALL QTITLE(0,TEXT(1:NAB))
               CALL QUAT('blue full 2.0 line')
               CALL QUAT('red circle 2.0 marker') 
               CALL QUAT('blue fcircle symbol')
#endif 
            END IF
            GOTO 10 
         END IF

         KEYSTX='x-y'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
            WRITE(*,*) '      i       X              Y '
            WRITE(*,*) '     -------------------------------'
 20         READ(LUN,101,END=90) TEXT
            CALL RLTEXT(TEXT,IA,IB,NAB) ! indices for non-blank area
            IF(NAB.EQ.0) GOTO 20
            CALL RATEXT(TEXT(1:NAB),NUMS,DNUM) ! translate text
            IF(KEYB.EQ.0.AND.NUMS.GE.2) THEN
               NPOINT=NPOINT+1
C               WRITE(*,*) (DNUM(L),L=1,NUMS)
               X=DNUM(1) 
               Y=DNUM(2)
               IF(NUMS.EQ.2) THEN
#ifdef graphics
                  CALL QDNXY(1,X,Y,ML)
#endif 
                  WRITE(*,103) NPOINT,X,Y
               ELSE IF(NUMS.EQ.4) THEN
                  DX=DNUM(3)
                  DY=DNUM(4)  
#ifdef graphics
                  CALL QDSYMB(X,Y,DX,DY)
                  WRITE(*,103) NPOINT,X,Y,DX,DY
#endif
               END IF
               GOTO 20
            ELSE IF(KEYB.NE.0) THEN
               KEYSTX='end of xy-data'
               MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT)
               IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
                  WRITE(*,*) ' '
#ifdef graphics
                  CALL QUPDAT
#endif
                  RETURN
               END IF
               GOTO 10
            ELSE
               GOTO 20
            END IF 
         END IF 

         KEYSTX='x-y-dx-dy'
         MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT) ! comparison
         IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
            WRITE(*,*) '      i       X              Y ',
     +                 '             DX             DY '
            WRITE(*,*) '     -----------------------------------------',
     +                 '---------------------'              
 30         READ(LUN,101,END=90) TEXT
            CALL RLTEXT(TEXT,IA,IB,NAB) ! indices for non-blank area
            IF(NAB.EQ.0) GOTO 30
            CALL RATEXT(TEXT(1:NAB),NUMS,DNUM) 
            IF(KEYB.EQ.0.AND.NUMS.GE.2) THEN
               NPOINT=NPOINT+1 
               X=DNUM(1) 
               Y=DNUM(2)
               IF(NUMS.EQ.2) THEN
#ifdef graphics
                  CALL QDNXY(1,X,Y,ML)
#endif
                  WRITE(*,103) NPOINT,X,Y 
               ELSE IF(NUMS.EQ.4) THEN
                  DX=DNUM(3)
                  DY=DNUM(4)  
#ifdef graphics
                  CALL QDSYMB(X,Y,DX,DY)
#endif
                  WRITE(*,103) NPOINT,X,Y,DX,DY
               END IF
               GOTO 30
            ELSE IF(KEYB.NE.0) THEN
               KEYSTX='end of xy-data'
               MAT=MATINT(TEXT(KEYA:KEYB),KEYSTX,NPAT,NTEXT)
               IF(MAT.GE.(NPAT+NPAT+3)/3) THEN
                  WRITE(*,*) ' '
#ifdef graphics
                  CALL QUPDAT
#endif
                  RETURN
               END IF
               GOTO 10
            ELSE
               GOTO 30
            END IF 
         END IF                    
      END IF
      RETURN 
 90   STOP    
 101  FORMAT(A)
 103  FORMAT(I8,G15.5,G15.5,G15.5,G15.5)
      END

      SUBROUTINE MPDATA(LUN)
      CHARACTER*80 TEXT, KEYWRD*32,KEYSTX*24
      DOUBLE PRECISION DNUM
      COMMON/NUMDUM/DNUM(100),NUMS

      PARAMETER (NDIM=100)
      COMMON/CNUMS/LAST,KEYA,KEYB,ICD(2,NDIM)
*     ...
      WRITE(*,*) 'MPDATA'
      END




      SUBROUTINE EXCOMD
      END

      SUBROUTINE PHISTE
      CHARACTER*80 TEXT 
*     Histogram
      COMMON/HISTEX/IH,IV,IT,NBIN,IBIN,INOUT(3),XA,XB,HIST(120),
     +              XMIN,XMAX,XMEAN,XSIGM
      CHARACTER*60 COMT
      COMMON/HSTEXT/COMT
*     ...
      IF(NBIN.EQ.0.OR.NBIN.NE.IBIN) RETURN
*
#ifdef graphics
      CALL QFINIT(0.0,0.0)
      CALL TXTNUM(FLOAT(IH))
      CALL TXTNUM(FLOAT(IV)) 
      TEXT='Hist @/@: '
      CALL RLTEXT(COMT,IA,IB,NAB)
      TEXT(11:12+IB-IA)=COMT(IA:IB)
      CALL QTITLE(0,TEXT(1:12+IB-IA))

      IF(IT.EQ.1.AND.XMIN.NE.XMAX) THEN
         CALL TXTNUM(XMIN)
         CALL TXTNUM(XMAX)
         TEXT='Min = @,  Max = @'
         CALL QTITLE(1,TEXT(1:17))
      ELSE IF(IT.EQ.2) THEN 
         CALL TXTNUM(XMIN)
         CALL TXTNUM(XMAX)
         TEXT='Log10 of entries. Min = @,  Max = @'
         CALL QTITLE(1,TEXT(1:35))
      END IF    

      CALL QUAT('normal blue yellow histogram')
      CALL QDHIST(NBIN,XA,XB,HIST)

      IF(XSIGM.NE.0.0) THEN
         CALL QLGAUS(INOUT(2),NBIN,XA,XB,XMEAN,XSIGM)
      END IF
      CALL QUPDAT
#endif
      END

      SUBROUTINE QLGAUS(IN,NX,XA,XB,XMEAN,XSIGM) ! plot Gauss function
      PARAMETER (FGA=0.39894228)
      CHARACTER*5 TEXT 
*     ...
      BIN=(XB-XA)/FLOAT(NX) 
      X1=XMEAN-3.0*XSIGM
      X2=XMEAN+3.0*XSIGM
      I1=1.0+(X1-XA)/(XB-XA)*FLOAT(NX)
      I2=1.0+(X2-XA)/(XB-XA)*FLOAT(NX)
      I1=MAX( 1,I1)
      I2=MIN(NX,I2)
      IF(I2-I1.LT.9) RETURN
#ifdef graphics
      CALL QUAT('full red 2.0 line')
#endif
      DO I=I1-1,I2
       X=XA+FLOAT(I)*BIN
       F=BIN*FLOAT(IN)*FGA/XSIGM*EXP(-0.5*((X-XMEAN)/XSIGM)**2)
#ifdef graphics
       CALL QDNXY(1,X,F,'f') ! gaussian density function
#endif
      END DO

      PEAK=BIN*FLOAT(IN)*FGA/XSIGM
      XTM=XMEAN+0.1*XSIGM
      XTM=XMEAN-5.0*BIN
      YTM=1.03*PEAK
#ifdef graphics
      CALL TXTNUM(XMEAN)
      CALL TXTNUM(XSIGM)
      TEXT='@   @'
      TEXT(3:3)=CHAR(177)  
      CALL QDTEXT(1,XTM,YTM,TEXT)
#endif
      END


