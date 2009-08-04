

*
*                          HMP... and GMP...
*                    Histogram and XY data in text files 
*
*     Booking:
*    
*     CALL HMPDEF(IH,XA,XB,TEXT)               CALL GMPDEF(IG,ITYP,TEXT) 
*     where                                    where 
*        IH    = 1 ... 10                         IG    = 1 ... 10
*        XA,XB = left, right limit                ITYP  = 1 dots
*        TEXT  = explanation                            = 2 line
*                                                       = 3 dots and line
*                                                       = 4 symbols 
*                                                       = 5 mean/sigma 
*                                                 TEXT  = explanation  
*
*     CALL HMPLUN(LUNW)                        CALL GMPLUN(LUNW)   
*     unit for output                          unit for output 
*
*     CALL HMPENT(IH,X)                        CALL GMPXY(IG,X,Y)
*     entry flt.pt. X                          add (X,Y) pair 
*
*                                              CALL GMPXYD(IG,X,Y,DX,DY) 
*                                              add (X,Y,DX,DY) ITYP=4
*
*        new                                   CALL GMPMS(IG,X,Y)
*                                              mean/sigma from x,y
*                                                
*     Booking log integer histogram:  
*
*     CALL HMPLDF(IH,TEXT)   
*     book and reset log integer histogram
*
*     CALL HMPLNT(IH,IX)                             
*     entry integer IX
*
*     Printing and writing:
*
*     CALL HMPRNT(IH)                          CALL GMPRNT(IG) 
*     print histogram IH or all, if 0          print data Ig or all, if 0     
*
*     CALL HMPWRT(IH)                          CALL GMPWRT(IG)
*     write histogram IH or all to file        write data IG or all to file
*                  
*
*     Storage manager for GMP...
*
*      CALL STMARS                 ! init/reset  storage manager
*
*      CALL STMAPR(JFLC,X,Y)       ! store pair (X,Y)
*
*      CALL STMADP(JFLC,FOUR)      ! store double pair 
*
*      CALL STMACP(JFLC,ARRAY,N)   ! copy (cp) all pairs to array
*
*      CALL STMARM(JFLC)           ! remove (rm) stored paiirs
*
*     _____________________________________________________________________


*     *************************** Histograms ******************************

      SUBROUTINE HMPDEF(IH,XA,XB,TEXT)           ! book, reset histogram
*     book millepede histogram, 120 bins
      PARAMETER (NUMHIS=12)    
      INTEGER INHIST(120,NUMHIS),JNHIST(5,NUMHIS),KHIST(NUMHIS)
      REAL    FNHIST(120,NUMHIS)
      EQUIVALENCE (INHIST(1,1),FNHIST(1,1))  
      INTEGER KVERS(NUMHIS)
      REAL XL(6,NUMHIS)
      DOUBLE PRECISION DL(2,NUMHIS)
      CHARACTER*60 HTEXT(NUMHIS),TEXT*(*)
      SAVE
      DATA KHIST/NUMHIS*0/,LUN/7/
*     ...
      IF(IH.LE.0.OR.IH.GT.NUMHIS) RETURN
c      IF(XA.EQ.XB)                RETURN
      DO I=1,120
       INHIST(I,IH)=0 
      END DO
      DO J=1,5
       JNHIST(J,IH)=0
      END DO
      XL(1,IH)=XA
      XL(2,IH)=XB
      XL(3,IH)=0.0
      IF(XA.NE.XB) XL(3,IH)=120.0/(XB-XA)
      XL(6,IH)=0.5*(XA+XB) ! center
      IF(KHIST(IH).EQ.0) THEN
         KVERS(IH)=0
      ELSE
         KVERS(IH)=KVERS(IH)+1
      END IF
      KHIST(IH)=1    ! flt.pt. (lin) 
      HTEXT(IH)=TEXT
      DL(1,IH)=0.0D0
      DL(2,IH)=0.0D0
      RETURN  

      ENTRY HMPLDF(IH,TEXT)                   ! book, reset log histogram
      IF(IH.LE.0.OR.IH.GT.NUMHIS) RETURN
      DO I=1,120
       INHIST(I,IH)=0
      END DO
      DO J=1,5
       JNHIST(J,IH)=0
      END DO
      IF(KHIST(IH).EQ.0) THEN
         KVERS(IH)=0
      ELSE
         KVERS(IH)=KVERS(IH)+1
      END IF
      KHIST(IH)=2    ! integer log
      HTEXT(IH)=TEXT
      XL(1,IH)=0.0
      XL(2,IH)=6.0
      RETURN

      ENTRY HMPENT(IH,X)                      ! entry flt.pt.
      IF(IH.LE.0.OR.IH.GT.NUMHIS) RETURN
      IF(KHIST(IH).NE.1)          RETURN
      JNHIST(4,IH)=JNHIST(4,IH)+1      ! count
      IF(JNHIST(4,IH).LE.120) THEN   
         FNHIST(JNHIST(4,IH),IH)=X     ! store value
         IF(JNHIST(4,IH).EQ.120) THEN
            CALL HMPMAK(INHIST(1,IH),FNHIST(1,IH),JNHIST(1,IH),
     +      XL(1,IH),DL(1,IH))
         END IF
         RETURN
      END IF 
c      IF(JNHIST(1,IH)+JNHIST(2,IH)+JNHIST(3,IH).EQ.0) THEN 
c         XL(4,IH)=X
c         XL(5,IH)=X
c      END IF
      I=1.0+XL(3,IH)*(X-XL(1,IH))   ! X - Xmin  
      J=2   
      IF(I.LT.  1) J=1
      IF(I.GT.120) J=3
      JNHIST(J,IH)=JNHIST(J,IH)+1
      XL(4,IH)=MIN(XL(4,IH),X)
      XL(5,IH)=MAX(XL(5,IH),X)
      IF(J.NE.2) RETURN
      INHIST(I,IH)=INHIST(I,IH)+1
      DL(1,IH)=DL(1,IH)+ X-XL(6,IH)
      DL(2,IH)=DL(2,IH)+(X-XL(6,IH))**2      
      RETURN   

      ENTRY HMPLNT(IH,IX)                     ! entry integer
      IF(IH.LE.0.OR.IH.GT.NUMHIS) RETURN
      IF(KHIST(IH).NE.2)          RETURN
      IF(IX.LE.0) THEN
         JNHIST(1,IH)=JNHIST(1,IH)+1
      ELSE
         IF(JNHIST(4,IH).EQ.0) JNHIST(4,IH)=IX
         IF(JNHIST(5,IH).EQ.0) JNHIST(5,IH)=IX
         JNHIST(4,IH)=MIN(JNHIST(4,IH),IX)
         JNHIST(5,IH)=MAX(JNHIST(5,IH),IX)
         I=1.0+20.0*LOG10(REAL(IX))
         J=2
         IF(I.LT.  1) J=1
         IF(I.GT.120) J=3
         IF(J.EQ.2)  INHIST(I,IH)=INHIST(I,IH)+1
         JNHIST(J,IH)=JNHIST(J,IH)+1
      END IF 
      RETURN

      ENTRY HMPRNT(IH)                        ! print, content vert
      IF(IH.EQ.0) THEN
         IHA=1
         IHB=NUMHIS 
      ELSE
         IF(IH.LE.0.OR.IH.GT.NUMHIS) RETURN 
         IHA=IH
         IHB=IH
      END IF
      DO IHC=IHA,IHB
       IF(KHIST(IHC).NE.0) THEN
        IF(KHIST(IHC).EQ.1) THEN
           CALL HMPMAK(INHIST(1,IHC),FNHIST(1,IHC),JNHIST(1,IHC),
     +                 XL(1,IHC),DL(1,IHC))  
        END IF   
        NN=JNHIST(1,IHC)+JNHIST(2,IHC)+JNHIST(3,IHC)
        IF(NN.NE.0.OR.KHIST(IHC).EQ.3) THEN
           WRITE(*,111)
 111       FORMAT(' ______',2('______________________________')) 
           IF(KVERS(IHC).EQ.1) THEN
              WRITE(*,*) 'Histogram',IHC,': ',HTEXT(IHC)
           ELSE
              WRITE(*,*) 'Histogram',IHC,'/',KVERS(IHC),': ',HTEXT(IHC)
           END IF  
           IF(KHIST(IHC).EQ.1) THEN
              WRITE(*,*) '   Out_low  inside  out_high = ',
     +                   (JNHIST(J,IHC),J=1,3)
           ELSE IF(KHIST(IHC).EQ.2) THEN
              WRITE(*,*) '   0_or_negative  inside  above_10^6 = ',
     +                   (JNHIST(J,IHC),J=1,3)
           END IF
           IF(KHIST(IHC).EQ.3) THEN
              CALL PFVERT(120,FNHIST(1,IHC))
           END IF
           IF(JNHIST(2,IHC).NE.0) THEN        ! integer content
              CALL PIVERT(120,INHIST(1,IHC))
              IF(KHIST(IHC).EQ.1) THEN
                 CALL PSVERT(XL(1,IHC),XL(2,IHC))  
              ELSE IF(KHIST(IHC).EQ.2) THEN
                 CALL PSVERT(0.0,6.0)
              END IF
           END IF
           IF(KHIST(IHC).EQ.1) THEN
              WRITE(*,*) '   Min and Max are',XL(4,IHC),XL(5,IHC)
              IF(JNHIST(2,IHC).GT.1) THEN
                 XMEAN=XL(6,IHC)+DL(1,IHC)/FLOAT(JNHIST(2,IHC))
                 XCENT=0.5D0*(XL(1,IHC)+XL(2,IHC))
                 XSIGM=(DL(2,IHC)-DL(1,IHC)**2/FLOAT(JNHIST(2,IHC)))
                 XSIGM=SQRT(XSIGM/FLOAT(JNHIST(2,IHC)-1))
                 WRITE(*,*) '   Mean and sigma are',
     +           XMEAN,' +-',XSIGM
              END IF
           ELSE IF(KHIST(IHC).EQ.2) THEN
              WRITE(*,*) '   Plot of log10 of entries. Min and Max are',
     +                   JNHIST(4,IHC),JNHIST(5,IHC)
           END IF 
        END IF
       END IF
      END DO
      RETURN

      ENTRY HMPLUN(LUNW)                      ! unit for output
      LUN=LUNW
      RETURN

      ENTRY HMPWRT(IH)                        ! write histogram text file
      IF(LUN.LE.0) RETURN
      IF(IH.EQ.0) THEN
         IHA=1
         IHB=NUMHIS
      ELSE
         IF(IH.LE.0.OR.IH.GT.NUMHIS) RETURN
         IHA=IH
         IHB=IH
      END IF

      DO IHC=IHA,IHB ! histogram loop        
       IF(KHIST(IHC).NE.0) THEN
          IF(KHIST(IHC).EQ.1) THEN
             CALL HMPMAK(INHIST(1,IHC),FNHIST(1,IHC),JNHIST(1,IHC),
     +                   XL(1,IHC),DL(1,IHC))  
          END IF   
          NBIN=120 
          WRITE(LUN,204) ' '
          WRITE(LUN,201) IHC,KVERS(IHC),KHIST(IHC)
          WRITE(LUN,204) HTEXT(IHC)  
          IF (JNHIST(1,IHC)+JNHIST(2,IHC)+JNHIST(3,IHC).EQ.0
     +         .AND.XL(1,IHC).EQ.XL(2,IHC)) THEN
c     hist is empty and hist range makes no sense
c     - cause: hist with 'variable edges' was never filled
c     - workaround: make lower and upper edge of hist differ in output
             WRITE(LUN,202) NBIN,XL(1,IHC)-0.001,XL(2,IHC)+0.001
          ELSE
             WRITE(LUN,202) NBIN,XL(1,IHC),XL(2,IHC)
          END IF
          WRITE(LUN,203) (JNHIST(J,IHC),J=1,3)
          WRITE(LUN,204) 'bincontent'
          IF(KHIST(IHC).EQ.1.OR.KHIST(IHC).EQ.2) THEN
             CALL KPRINT(LUN,INHIST(1,IHC),NBIN)
          ELSE
             WRITE(LUN,219) (FNHIST(I,IHC),I=1,NBIN)
          END IF

          IF(KHIST(IHC).EQ.1) THEN
             WRITE(LUN,205) XL(4,IHC),XL(5,IHC)
          ELSE IF(KHIST(IHC).EQ.2) THEN 
             WRITE(LUN,205) FLOAT(JNHIST(4,IHC)),FLOAT(JNHIST(5,IHC))
          END IF 
          IF(KHIST(IHC).EQ.1) THEN
             IF(JNHIST(2,IHC).GT.1) THEN
                XMEAN=XL(6,IHC)+DL(1,IHC)/FLOAT(JNHIST(2,IHC))
                XCENT=0.5D0*(XL(1,IHC)+XL(2,IHC))
                XSIGM=(DL(2,IHC)-DL(1,IHC)**2/FLOAT(JNHIST(2,IHC)))
                XSIGM=SQRT(XSIGM/FLOAT(JNHIST(2,IHC)-1))
                WRITE(LUN,206) XMEAN,XSIGM
             END IF
          END IF
          WRITE(LUN,204) 'end of histogram'   
       END IF
      END DO

 201  FORMAT('Histogram ',I4,10X,'version ',I4,10X,'type',I2)
 202  FORMAT(10X,' bins, limits ',I4,2G15.5)
 203  FORMAT(10X,'out-low inside out-high ',3I10)
 204  FORMAT(A)
 205  FORMAT('minmax',2E15.7)
 206  FORMAT('meansigma',2E15.7)

 219  FORMAT(4E15.7)
      END 

      SUBROUTINE HMPMAK(INHIST,FNHIST,JNHIST,XL,DL) ! hist scale from data
      REAL CPHIST(120)
      REAL FNHIST(120),XL(6)
      DOUBLE PRECISION DL(2)
      INTEGER INHIST(120),JNHIST(5)
      SAVE
*     ...
      NN=JNHIST(4)
c      WRITE(*,*) 'HMPMAK: NN,JNHIST(5)',NN,JNHIST(5)
      IF(NN.EQ.0.OR.JNHIST(5).NE.0) RETURN
      JNHIST(5)=1
      DO I=1,NN
c       WRITE(*,*) 'copy ',I,FNHIST(I)
       CPHIST(I)=FNHIST(I)
      END DO 
      CALL HEAPF(CPHIST,NN)
      IF(XL(3).EQ.0.0) THEN
         CALL BINTAB(CPHIST,NN,XA,XB)
         XL(1)=XA
         XL(2)=XB
         XL(3)=0.0
         IF(XA.NE.XB) XL(3)=120.0/(XB-XA)
         XL(6)=0.5*(XA+XB) ! center
      END IF 
      XL(4)=CPHIST( 1)
      XL(5)=CPHIST(NN)
c      WRITE(*,*) 'XL ',XL
      DO I=1,NN
       INHIST(I)=0
      END DO
      DO K=1,NN
       X=CPHIST(K)
       I=1.0+XL(3)*(X-XL(1))   ! X - Xmin
c       WRITE(*,*) 'K,I,X ',K,I,X
       J=2
       IF(I.LT.  1) J=1
       IF(I.GT.120) J=3
       JNHIST(J)=JNHIST(J)+1
       IF(J.EQ.2) THEN
          INHIST(I)=INHIST(I)+1
          DL(1)=DL(1)+ X-XL(6)
          DL(2)=DL(2)+(X-XL(6))**2
       END IF  
      END DO
      END          


 

      SUBROUTINE BINTAB(TAB,N,XA,XB)             ! hist scale from data   
*     Bin limits XA and XB from TAB(N)
      REAL TAB(N)
      REAL BIN(10)
      DATA BIN/1.0,1.5,2.0,3.0,4.0,5.0,8.0,10.0,15.0,20.0/
*     ...

      CALL HEAPF(TAB,N) ! reduced statistic
c      WRITE(*,*) ' '  
c      WRITE(*,*) 'Sorted ',(TAB(I),I=1,N)
      IF(N.LT.100) THEN
         X1=TAB(1)
         X2=TAB(N) 
c         WRITE(*,*) 'reduced statistic X1 X2 ',X1,X2
      ELSE              ! large statistic
         M1=1.0+0.05*FLOAT(N)
         M2=1.0+0.16*FLOAT(N)
         X1=TAB(M1)-4.0*(TAB(M2)-TAB(M1)) 
         IF(X1.LT.0.0.AND.TAB(1).GE.0.0) X1=TAB(1)          
         X2=TAB(N+1-M1)+4.0*(TAB(N+1-M1)-TAB(N+1-M2))
         IF(X2.GT.0.0.AND.TAB(N).LE.0.0) X2=TAB(N)
c         WRITE(*,*) 'large statistic ',X1,X2 
c         WRITE(*,*) 'min und max ',TAB(1),TAB(N)
         IF(X1*TAB(1).LE.0.0) X1=0.0
         IF(X2*TAB(N).LE.0.0) X2=0.0
c         WRITE(*,*) 'large statistic zero ',X1,X2
         IF(X1*X2.LT.0.0.AND.MIN(-X1,X2).GT.0.6*MAX(-X1,X2)) THEN
            XX=MAX(-X1,X2) ! symmetry
            X1=-XX
            X2=+XX  
         ELSE IF(X1*X2.GT.0.0.AND.  ! include zero ?
     +           ABS(MIN(X1,X2)).LT.0.4*ABS(MAX(X1,X2))) THEN
            IF(X1.LT.0.0) THEN
               X2=0.0
            ELSE 
               X1=0.0
            END IF
         END IF
c         WRITE(*,*) 'large statistic ',X1,X2 
      END IF
      IF(X1.EQ.X2) THEN
         X1=X1-1.0
         X2=X2+1.0
      END IF  
      DX=X2-X1
c      WRITE(*,*) 'X1,X2,DX ',X1,X2,DX
      RAT=0.0
      DO J=1,11
       I=J
       IF(J.EQ.11) I=II 
       IEXP=101.0+LOG10(DX)-LOG10(6.0*BIN(I))
       IEXP=IEXP-100 
       DD=BIN(I)*10.0**IEXP

       N1=ABS(X1)/DD
       IF(X1.LT.0.0) N1=-N1
       IF(FLOAT(N1)*DD.GT.X1) N1=N1-1
c       WRITE(*,*) 'Bin ',I,N1,N1*DD,X1

       N2=ABS(X2)/DD
       IF(X2.LT.0.0) N2=-N2
       IF(FLOAT(N2)*DD.LT.X2) N2=N2+1
c       WRITE(*,*) 'Bin ',I,N2,N2*DD,X2
 10    IF(N2-N1.LT.6) THEN
          IF(N1.NE.0) N1=N1-1
          IF(N2-N1.LT.6.AND.N2.NE.0) N2=N2+1
          GOTO 10
       END IF
c       WRITE(*,*) 'corrected N1 N2 ',N1,N2
       XA=SIGN(FLOAT(N1)*DD,X1)
       XB=SIGN(FLOAT(N2)*DD,X2)
c       WRITE(*,*) J,' resulting limits XA XB ',XA,XB
       IF((X2-X1)/(XB-XA).GT.RAT) THEN
          II=I
          RAT=(X2-X1)/(XB-XA)
       END IF
      END DO 
c      WRITE(*,*) J,' resulting limits XA XB ',XA,XB
      END 

      SUBROUTINE KPRINT(LUN,LIST,N)              ! print integer array
*     print integer array LIST(N)
      INTEGER LIST(N),LI(7)
      DATA LI/2,3,4,6,8,9,12/ ! number of characters
*     ...
      IB=0
 10   IA=IB+1
      IF(IA.GT.N) RETURN
      DO K=1,7
       NP=72/LI(K)
       IB=MIN(IA-1+NP,N)
       IF(K.LE.6) THEN
          LP=10**(LI(K)-1)-1 ! maximum positive
          LN=-LP/10       ! minimum negative
          DO I=IA,IB
           IF(LIST(I).GT.LP.OR.LIST(I).LT.LN) GOTO 20
          END DO
       END IF 
       IF(K.EQ.1) THEN
          WRITE(LUN,101) (LIST(I),I=IA,IB)
       ELSE IF(K.EQ.2) THEN
          WRITE(LUN,102) (LIST(I),I=IA,IB)
       ELSE IF(K.EQ.3) THEN
          WRITE(LUN,103) (LIST(I),I=IA,IB)
       ELSE IF(K.EQ.4) THEN
          WRITE(LUN,104) (LIST(I),I=IA,IB)
       ELSE IF(K.EQ.5) THEN
          WRITE(LUN,105) (LIST(I),I=IA,IB)
       ELSE IF(K.EQ.6) THEN
          WRITE(LUN,106) (LIST(I),I=IA,IB)
       ELSE IF(K.EQ.7) THEN
          WRITE(LUN,107) (LIST(I),I=IA,IB)
       END IF
       GOTO 10 
  20   CONTINUE
      END DO
 101  FORMAT(36I2)
 102  FORMAT(24I3)
 103  FORMAT(18I4)
 104  FORMAT(12I6)
 105  FORMAT( 9I8)
 106  FORMAT( 8I9)
 107  FORMAT( 6I12)
      END 

*     ***************************** XY data **************************** 

      SUBROUTINE GMPDEF(IG,ITYP,TEXT)            ! book, reset XY storage
*     ITYP = 1  X,Y     as dots
*          = 2  X,Y     as line 
*          = 3  X,Y     as line and dots 
*          = 4  X,Y, DX,DY symbols
      PARAMETER (NARR=1000)
      REAL ARRAY(2,NARR),ARRAY4(4,NARR/2),ARRAY1(NARR+NARR),FOUR(4)
      EQUIVALENCE (ARRAY(1,1),ARRAY4(1,1),ARRAY1(1)) 
      PARAMETER (NUMGXY=10) 
      INTEGER NSTR(NUMGXY),IGTP(NUMGXY),LVERS(NUMGXY),NST(3,NUMGXY)
      REAL XYPLWS(10,NUMGXY)
      INTEGER JFLC(5,NUMGXY),KFLC(5,NUMGXY)
*     JFLC(1,.) = first used index 
*     JFLC(2,.) = last used index
*     JFLC(3,.) = counter of used places
*     JFLC(4,.) = counter of ignored  
*     JFLC(5,.) = limit for JFLC(3) 
      CHARACTER*60 GTEXT(NUMGXY)
      CHARACTER*(*) TEXT  
      LOGICAL START
      SAVE
      DATA START/.TRUE./,LUN/7/
      DATA NSTR/NUMGXY*0/   ! by GF
*     ...
      IF(START) THEN
         START=.FALSE.
         CALL STMARS    ! initialize storage
         DO I=1,NUMGXY
          DO J=1,5
           JFLC(J,I)=0
           KFLC(J,I)=0
          END DO
         END DO 
      END IF        

      IF(IG.LT.1.OR.IG.GT.NUMGXY) RETURN
      IF(ITYP.LT.1.OR.ITYP.GT. 5) RETURN 
      IF(NSTR(IG).EQ.0) THEN
         LVERS(IG)=0
      ELSE
         LVERS(IG)=LVERS(IG)+1
      END IF
      NSTR(IG)=1 ! by GF
*        remove stored elements
         IF(JFLC(1,IG).NE.0) CALL STMARM(JFLC(1,IG))
         IF(KFLC(1,IG).NE.0) CALL STMARM(KFLC(1,IG))
      IGTP(IG)=ITYP
      GTEXT(IG)=TEXT
      DO J=1,5
       JFLC(J,IG)=0
      END DO
      JFLC(5,IG)=200
      IF(ITYP.EQ.5) THEN
         DO J=1,5
          KFLC(J,IG)=0
         END DO
         JFLC(5,IG)=128 ! maximum of 128 values 
         KFLC(5,IG)=NARR
         NST(1,IG)=0
         NST(2,IG)=0
         NST(3,IG)=1
         DO J=1,10
          XYPLWS(J,IG)=0.0 
         END DO
      END IF
      RETURN

      ENTRY GMPXY(IG,X,Y)                     ! add (X,Y) pair
      IF(IG  .LT.1.OR.IG.GT.NUMGXY) RETURN        ! check argument IG
      IF(IGTP(IG).LT.1.OR.IGTP(IG).GT.3) RETURN   ! check type
      CALL STMAPR(JFLC(1,IG),X,Y)
      RETURN

      ENTRY GMPXYD(IG,X,Y,DX,DY)              ! add (X,Y,DX,DY)
      IF(IG  .LT.1.OR.IG.GT.NUMGXY) RETURN        ! check argument IG
      IF(IGTP(IG).NE.4) RETURN 
      FOUR(1)=X
      FOUR(2)=Y
      FOUR(3)=DX
      FOUR(4)=DY
      CALL STMADP(JFLC(1,IG),FOUR)
      RETURN

      ENTRY GMPMS(IG,X,Y)                     ! mean sigma(X) from Y
*     mean sigma from Y, as a function of X    
C      WRITE(*,*) 'GMPMS ',IG,X,Y

      IF(IG  .LT.1.OR.IG.GT.NUMGXY) RETURN        ! check argument IG
      IF(IGTP(IG).NE.5) RETURN
     
      XYPLWS(10,IG)=X  ! last X  coordinate
      IF(NST(1,IG).EQ.0) THEN
         Y1=Y
         NST(1,IG)=1
         IF(KFLC(3,IG).EQ.0) XYPLWS(9,IG)=X        ! start coordinate
      ELSE 
         NST(1,IG)=0 
         CALL STMAPR(KFLC(1,IG),Y1,Y) ! store pair
         IF(KFLC(3,IG).GE.KFLC(5,IG)) THEN  
            CALL STMACP(KFLC(1,IG),ARRAY,N) ! get data
            CALL STMARM(KFLC(1,IG))         ! remove data
            N=N+N 
            CALL RMESIG(ARRAY,N,XYPLWS(2,IG),XYPLWS(4,IG))
            NST(2,IG)=NST(2,IG)+1
            IF(NST(2,IG).EQ.1) XYPLWS(7,IG)=XYPLWS(9,IG)
            XYPLWS(8,IG)=X                   ! end coordinate
            XYPLWS(5,IG)=XYPLWS(5,IG)+XYPLWS(2,IG)
            XYPLWS(6,IG)=XYPLWS(6,IG)+XYPLWS(4,IG)
            IF(NST(2,IG).EQ.NST(3,IG)) THEN 
               XYPLWS(1,IG)=0.5*(XYPLWS(7,IG)+XYPLWS(8,IG))
               XYPLWS(2,IG)=XYPLWS(5,IG)/FLOAT(NST(3,IG))
               XYPLWS(3,IG)=0.5*(XYPLWS(8,IG)-XYPLWS(7,IG)) 
               XYPLWS(4,IG)=XYPLWS(6,IG)/FLOAT(NST(3,IG))
               XYPLWS(5,IG)=0.0
               XYPLWS(6,IG)=0.0
               NST(2,IG)=0 
               CALL STMADP(JFLC(1,IG),XYPLWS(1,IG)) 
               IF(JFLC(3,IG).GE.JFLC(5,IG)) THEN
                  CALL STMACP(JFLC(1,IG),ARRAY4,N)   ! get data
                  N=N/2
                  CALL STMARM(JFLC(1,IG))            ! remove data
                  DO I=1,N,2                   ! average
                   XYPLWS(7,IG)=ARRAY4(1,I  )-ARRAY4(3,  I)
                   XYPLWS(8,IG)=ARRAY4(1,I+1)+ARRAY4(3,I+1)
                   XYPLWS(1,IG)=0.5*(XYPLWS(7,IG)+XYPLWS(8,IG)) 
                   XYPLWS(2,IG)=0.5*(ARRAY4(2,I)+ARRAY4(2,I+1))
                   XYPLWS(3,IG)=0.5*(XYPLWS(8,IG)-XYPLWS(7,IG))
                   XYPLWS(4,IG)=0.5*(ARRAY4(4,I)+ARRAY4(4,I+1))
                   CALL STMADP(JFLC(1,IG),XYPLWS(1,IG))
                  END DO
                  NST(3,IG)=NST(3,IG)+NST(3,IG)
               END IF
            END IF
         END IF 
      END IF
      RETURN

      ENTRY GMPRNT(IG)                        ! print XY data
      IF(IG.EQ.0) THEN
         IGA=1
         IGB=NUMGXY
      ELSE
         IF(IG.LE.0.OR.IG.GT.NUMGXY) RETURN
         IGA=IG
         IGB=IG
      END IF
      DO IGC=IGA,IGB

       IF(IGTP(IGC).GE.1.AND.IGTP(IGC).LE.3) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Store ',IGC,': ',GTEXT(IGC)
          IF(JFLC(4,IGC).EQ.0) THEN
             WRITE(*,*) '      stored n-tuples: ',JFLC(3,IGC)
          ELSE  
             WRITE(*,*) '   stored n-tuples,  not-stored n-tuples: ',
     +                  JFLC(3,IGC),', ',JFLC(4,IGC)
          END IF

          CALL STMACP(JFLC(1,IGC),ARRAY,NA) ! get all data 

          DO N=1,NA
           WRITE(*,102) N, ARRAY(1,N),ARRAY(2,N) 
          END DO

       ELSE IF(IGTP(IGC).EQ.4) THEN

          WRITE(*,*) ' '
          WRITE(*,*) 'Store ',IGC,': ',GTEXT(IGC)
          IF(JFLC(4,IGC).EQ.0) THEN
             WRITE(*,*) '      stored n-tuples: ',JFLC(3,IGC)
          ELSE
             WRITE(*,*) '   stored n-tuples,  not-stored n-tuples: ',
     +                  JFLC(3,IGC),', ',JFLC(4,IGC)
          END IF

          CALL STMACP(JFLC(1,IGC),ARRAY,NA) ! get all data
          NA=NA/2

          DO N=1,NA
           WRITE(*,102) N,(ARRAY4(J,N),J=1,4)
          END DO

        ELSE IF(IGTP(IGC).EQ.5) THEN

          CALL STMACP(KFLC(1,IGC),ARRAY,N) ! get data
          CALL STMARM(KFLC(1,IGC))         ! remove data
          N=N+N 
          IF(NST(1,IGC).EQ.1) THEN
             N=N+1
             ARRAY1(N)=Y1
             NST(1,IGC)=0 ! reset
          END IF
          IF(N.NE.0) THEN 
             XYPLWS(7,IGC)=XYPLWS( 9,IGC)
             XYPLWS(8,IGC)=XYPLWS(10,IGC)
             CALL RMESIG(ARRAY1,N,XYPLWS(2,IGC),XYPLWS(4,IGC))
             WGHT=FLOAT(N)/FLOAT(NST(3,IGC)*KFLC(5,IGC))
             XYPLWS(5,IGC)=XYPLWS(5,IGC)+XYPLWS(2,IGC)*WGHT
             XYPLWS(6,IGC)=XYPLWS(6,IGC)+XYPLWS(4,IGC)*WGHT
             XYPLWS(2,IGC)=XYPLWS(5,IGC)/(FLOAT(NST(2,IGC))+WGHT)
             XYPLWS(4,IGC)=XYPLWS(6,IGC)/(FLOAT(NST(2,IGC))+WGHT)
             XYPLWS(1,IGC)=0.5*(XYPLWS(7,IGC)+XYPLWS(8,IGC))
             XYPLWS(3,IGC)=0.5*(XYPLWS(8,IGC)-XYPLWS(7,IGC)) 
             CALL STMADP(JFLC(1,IGC),XYPLWS(1,IGC))
          END IF

          WRITE(*,*) ' '
          WRITE(*,*) 'Store ',IGC,': ',GTEXT(IGC)
          IF(JFLC(4,IGC).EQ.0) THEN
             WRITE(*,*) '      stored n-tuples: ',JFLC(3,IGC)
          ELSE
             WRITE(*,*) '   stored n-tuples,  not-stored n-tuples: ',
     +                  JFLC(3,IGC),', ',JFLC(4,IGC)
          END IF

          CALL STMACP(JFLC(1,IGC),ARRAY,NA) ! get all data
          NA=NA/2
          DO N=1,NA
           WRITE(*,102) N,(ARRAY4(J,N),J=1,4)
          END DO
       END IF 
      END DO
      RETURN   

      ENTRY GMPLUN(LUNW)                      ! unit for output
      LUN=LUNW
      RETURN

      ENTRY GMPWRT(IG)                        ! write XY text file
      IF(LUN.LE.0) RETURN
      IF(IG.EQ.0) THEN
         IGA=1
         IGB=NUMGXY
      ELSE
         IF(IG.LE.0.OR.IG.GT.NUMGXY) RETURN
         IGA=IG
         IGB=IG
      END IF
      DO IGC=IGA,IGB
       IF(IGTP(IGC).EQ.5) THEN 

          CALL STMACP(KFLC(1,IGC),ARRAY,N) ! get data
          CALL STMARM(KFLC(1,IGC))         ! remove data
          N=N+N 
          IF(NST(1,IGC).EQ.1) THEN
             N=N+1
             ARRAY1(N)=Y1
             NST(1,IGC)=0 ! reset
          END IF
          IF(N.NE.0) THEN 
             XYPLWS(7,IGC)=XYPLWS( 9,IGC)
             XYPLWS(8,IGC)=XYPLWS(10,IGC)
             CALL RMESIG(ARRAY1,N,XYPLWS(2,IGC),XYPLWS(4,IGC))
             WGHT=FLOAT(N)/FLOAT(NST(3,IGC)*KFLC(5,IGC))
             XYPLWS(5,IGC)=XYPLWS(5,IGC)+XYPLWS(2,IGC)*WGHT
             XYPLWS(6,IGC)=XYPLWS(6,IGC)+XYPLWS(4,IGC)*WGHT
             XYPLWS(2,IGC)=XYPLWS(5,IGC)/(FLOAT(NST(2,IGC))+WGHT)
             XYPLWS(4,IGC)=XYPLWS(6,IGC)/(FLOAT(NST(2,IGC))+WGHT)
             XYPLWS(1,IGC)=0.5*(XYPLWS(7,IGC)+XYPLWS(8,IGC))
             XYPLWS(3,IGC)=0.5*(XYPLWS(8,IGC)-XYPLWS(7,IGC)) 
             CALL STMADP(JFLC(1,IGC),XYPLWS(1,IGC))
          END IF

       END IF 
       IF(JFLC(3,IGC)+JFLC(4,IGC).NE.0) THEN
          WRITE(LUN,204) ' '
          WRITE(LUN,201) IGC,LVERS(IGC),IGTP(IGC)
          WRITE(LUN,204) GTEXT(IGC)
          WRITE(LUN,203) JFLC(3,IGC)+JFLC(4,IGC)
          CALL STMACP(JFLC(1,IGC),ARRAY,NA) ! get all data 
          IF(IGTP(IGC).GE.1.AND.IGTP(IGC).LE.3) THEN
             WRITE(LUN,204) 'x-y'
             DO N=1,NA
              WRITE(LUN,205) ARRAY(1,N),ARRAY(2,N)
             END DO 
          ELSE IF(IGTP(IGC).EQ.4.OR.IGTP(IGC).EQ.5) THEN
             WRITE(LUN,204) 'x-y-dx-dy'
             NA=NA/2
             DO N=1,NA
              WRITE(LUN,205) (ARRAY4(J,N),J=1,4)
             END DO
          END IF
          WRITE(LUN,204) 'end of xy-data' 
       END IF
      END DO
 102  FORMAT(I12,4G15.7)
c 103  FORMAT('       Index    ___X___       ___Y___     '/   
c     +       '       ----- -------------- --------------')
c 104  FORMAT('       Index    ___X___       ___Y___     ',
c     +                   '    ___DX__       ___DY__     '/    
c     +       '       ----- -------------- --------------',
c     +                   ' -------------- --------------')    
 201  FORMAT('XY-Data ',I4,10X,'version ',I4,10X,'type',I2)
 203  FORMAT(10X,'stored  not-stored ',2I10)
 204  FORMAT(A)
 205  FORMAT(3X,4G15.7)
      END 


      SUBROUTINE STMARS                          ! init/reset  storage
      PARAMETER (NDIM=3000)   ! storage dimension
      REAL      TK(2,NDIM)    ! pair storage for data pairs
      INTEGER   NEXT(NDIM)    ! pointer 
      INTEGER IFLC1,IFLC2     ! first and last index of free pairs

      REAL FOUR(4),ARRAY(2,*) ! double_pair, copy array
      INTEGER JFLC(5)         ! user array
*     JFLC(1) = first used index 
*     JFLC(2) = last used index
*     JFLC(3) = counter of used places
*     JFLC(4) = counter of ignored  
*     JFLC(5) = limit for JFLC(3) 
*     ...
      DO I=1,NDIM
       NEXT(I)=I+1   ! pointer to next free location
       TK(1,I)=0.0   ! reset
       TK(2,I)=0.0  
      END DO
      NEXT(NDIM)=0   ! ... and end pointer
      IFLC1=1        ! index first free pair
      IFLC2=NDIM     ! index last free pair
      RETURN

      ENTRY STMAPR(JFLC,X,Y)                  ! store pair (X,Y)
      IFRE=IFLC1                   ! index of free place
      IF(IFRE.EQ.0.OR.JFLC(3).GE.JFLC(5)) THEN ! overflow
         JFLC(4)=JFLC(4)+1
      ELSE
         IFLC1=NEXT(IFRE)          ! pointer to new free location
         IF(JFLC(1).EQ.0) THEN     ! first item
            JFLC(1)=IFRE 
         ELSE
            NEXT(JFLC(2))=IFRE
         END IF
         NEXT(IFRE)=0
         JFLC(2)=IFRE              ! last index
         JFLC(3)=JFLC(3)+1         ! counter
         TK(1,IFRE)=X
         TK(2,IFRE)=Y
      END IF
      RETURN 

      ENTRY STMADP(JFLC,FOUR)                 ! store double pair 
      IFREA=IFLC1                  ! index of 1. free place
      IF(IFREA.EQ.0) THEN ! overflow
         JFLC(4)=JFLC(4)+1
      ELSE
         IFREB=NEXT(IFLC1)         ! index of 2. free place
         IF(IFREB.EQ.0.OR.JFLC(3).GE.2*JFLC(5)) THEN ! overflow
            JFLC(4)=JFLC(4)+1
         ELSE   
            IFLC1=NEXT(IFREB)      ! pointer to new free location
            IF(JFLC(1).EQ.0) THEN  ! first item
               JFLC(1)=IFREA
            ELSE
               NEXT(JFLC(2))=IFREA
            END IF
            NEXT(IFREB)=0
            JFLC(2)=IFREB          ! last index
            JFLC(3)=JFLC(3)+1      ! counter
            TK(1,IFREA)=FOUR(1)
            TK(2,IFREA)=FOUR(2)
            TK(1,IFREB)=FOUR(3)
            TK(2,IFREB)=FOUR(4)
         END IF    
      END IF
      RETURN

      ENTRY STMACP(JFLC,ARRAY,N)              ! copy (cp) all pairs to array
      N=0
      IND=JFLC(1)
 10   IF(IND.EQ.0) RETURN
      N=N+1
      ARRAY(1,N)=TK(1,IND)
      ARRAY(2,N)=TK(2,IND)
      IND=NEXT(IND)
      GOTO 10

      ENTRY STMARM(JFLC)                      ! remove (rm) stored paiirs
      NEXT(IFLC2)=JFLC(1)          ! connect to free space
      IFLC2=JFLC(2)                ! new last free index
      DO J=1,4
       JFLC(J)=0
      END DO
      END 

      SUBROUTINE RMESIG(X,N,XLOC,XSCA)           ! robust mean and sigma
*     robust determination of location and scale parameter, 
*        for Gaussian data: location=mean and scale=standard deviation
*     XLOC = median of X_i            (N values in array X(N))
*     XCSA = median of | X_i - XLOC |, times 1.4826   
      REAL X(N) ! input array, modified 
*     ...
      XLOC=0.0
      XSCA=0.0
      IF(N.LE.0) RETURN
      CALL HEAPF(X,N)  ! sort 
      XLOC=0.5*(X((N+1)/2)+X((N+2)/2))         ! location
      DO I=1,N
       X(I)=ABS(X(I)-XLOC)
      END DO
      CALL HEAPF(X,N) ! sort
      XSCA=1.4826*0.5*(X((N+1)/2)+X((N+2)/2))  ! dispersion
      END



