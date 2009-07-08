


      SUBROUTINE RATEXT(TEXT,NUMS,DNUM)
*
*     Translate TEXT into arrays of double precision numbers
*
      CHARACTER*(*) TEXT, CH*1,KEYWRD*16
      PARAMETER (NDIM=100)
      COMMON/CNUMS/LAST,KEYA,KEYB,ICD(2,NDIM) 
      DOUBLE PRECISION DNUM(*),DIC(NDIM),DUMBER
      INTEGER ICDT(NDIM)
*     ...
      NUMS=0
      LAST=0
      KEYA=0
      KEYB=0 
      IF(TEXT(1:1).EQ.'*') RETURN
      NUM=ICHAR('0')
      LENT=0
      LAST=0
      DO I=1,LEN(TEXT)        ! find comment and end
       IF(LENT.EQ.0.AND.(TEXT(I:I).EQ.'!'.OR.TEXT(I:I).EQ.'%')) LENT=I
       IF(TEXT(I:I).NE.' ') LAST=I
      END DO
      IF(LENT.EQ.0) LENT=LAST+1
      ICD(1,1)=LENT
              
      J=1
c      IF(LENT.GT.1) WRITE(*,*) 'comment-free:',TEXT(1:LENT-1) 
      ICDT(1)=0
      ICL=0
      DO I=1,LENT-1
       CH =TEXT(I:I)
       ICH=ICHAR(CH)
       IC=0
       IF(CH.EQ.'.') IC=1
       IF(CH.EQ.'+') IC=2
       IF(CH.EQ.'-') IC=3
       IF(CH.EQ.'E') IC=4
       IF(CH.EQ.'D') IC=4
       IF(IC.GT.0) THEN
          J=J+1
          ICD(1,J)=I
          ICD(2,J)=I
          ICDT(J)=IC
       ELSE  
          IC=6
          IF(ICH.GE.NUM.AND.ICH.LE.NUM+9) IC=5  ! digit
          IF(IC.NE.ICL) THEN
             J=J+1
             ICD(1,J)=I
             ICDT(J)=IC
          END IF
          ICD(2,J)=I
       END IF 
       ICL=IC ! previous IC
      END DO
      ICDT(J+1)=0

      DO I=1,J                  ! define number number 
       IF(ICDT(I).EQ.5) THEN
          DUMBER=0.0D0
          DO K=ICD(1,I),ICD(2,I)
           DUMBER=10.D0*DUMBER+DFLOAT(ICHAR(TEXT(K:K))-NUM) 
          END DO
          DIC(I)=DUMBER
       END IF 
      END DO
      ICDT(J+1)=0

      DO I=2,J                  ! get dots
       IF(ICDT(I).EQ.1) THEN
          ICODE=0
          IF(ICDT(I-1).EQ.5.AND.ICD(2,I-1)+1.EQ.ICD(1,I)) ICODE=1
          IF(ICDT(I+1).EQ.5.AND.ICD(1,I+1)-1.EQ.ICD(2,I)) ICODE=ICODE+2
          IF(ICODE.EQ.1) THEN            ! 123.
             ICD(2,I-1)=ICD(2,I)
             ICDT(I)=0
          ELSE IF(ICODE.EQ.2) THEN       ! .456 
             DIC(I)=10.0D0**(ICD(1,I+1)-ICD(2,I+1)-1)*DIC(I+1)
             ICDT(I)=5
             ICD(2,I)=ICD(2,I+1)
             ICDT(I+1)=0
          ELSE IF(ICODE.EQ.3) THEN       ! 123.456 
             DIC(I-1)=DIC(I-1)+
     +                10.0D0**(ICD(1,I+1)-ICD(2,I+1)-1)*DIC(I+1)
             ICD(2,I-1)=ICD(2,I+1)
             ICDT(I)=0
             ICDT(I+1)=0
          END IF
       END IF
      END DO

      K=1                         ! remove blanks, compress
      DO I=2,J    
       IF(ICDT(I).EQ.6.AND.TEXT(ICD(1,I):ICD(2,I)).EQ.' ') ICDT(I)=0
       IF(ICDT(I).NE.0) THEN 
          K=K+1
          ICD(1,K)=ICD(1,I)
          ICD(2,K)=ICD(2,I)
          ICDT(K)=ICDT(I)
          DIC(K)=DIC(I)
       END IF
      END DO
      J=K 

      DO I=2,J-1
       IF(ICDT(I).EQ.2.OR.ICDT(I).EQ.3) THEN   !  +- 
          IF(ICDT(I+1).EQ.5) THEN        
             ICD(1,I+1)=ICD(1,I)
             IF(ICDT(I).EQ.3) DIC(I+1)=-DIC(I+1)
             ICDT(I)=0
          END IF 
       END IF
      END DO 

      K=1                         ! compress
      DO I=2,J    
       IF(ICDT(I).EQ.6.AND.TEXT(ICD(1,I):ICD(2,I)).EQ.' ') ICDT(I)=0
       IF(ICDT(I).NE.0) THEN 
          K=K+1
          ICD(1,K)=ICD(1,I)
          ICD(2,K)=ICD(2,I)
          ICDT(K)=ICDT(I)
          DIC(K)=DIC(I)
       END IF
      END DO
      J=K 

      DO I=2,J-1      
       IF(ICDT(I).EQ.4) THEN        ! E or D
          IF(ICDT(I-1).EQ.5.AND.ICDT(I+1).EQ.5) THEN
             ICD(2,I-1)=ICD(2,I+1)
             DIC(I-1)=DIC(I-1)*10.0D0**DIC(I+1) 
             ICDT(I)=0
             ICDT(I+1)=0
          END IF
       END IF
      END DO

      NUMS=0                         ! compress
      DO I=1,J    
       IF(ICDT(I).EQ.5) THEN 
          NUMS=NUMS+1
          ICD(1,NUMS)=ICD(1,I)
          ICD(2,NUMS)=ICD(2,I)
          DNUM(NUMS)=DIC(I)
       END IF
      END DO

      KEYWRD=' '                     ! assemble keyword
      IA=0
      IB=-1
      DO I=1,ICD(1,1)-1 
       IF(IA.EQ.0.AND.TEXT(I:I).NE.' ') IA=I 
       IF(TEXT(I:I).NE.' ') IB=I
      END DO
      IF(IB.GE.0) KEYWRD=TEXT(IA:IB)
      KEYA=IA
      KEYB=MAX(0,IB)

c      WRITE(*,*) ' '
c      WRITE(*,*) ICD(1,1),'<-Keyword is: ',KEYWRD
c      WRITE(*,*) ' '
c      DO I=1,NUMS
c       WRITE(*,101) I,ICD(1,I),ICD(2,I),DNUM(I),
c     +              TEXT(ICD(1,I):ICD(2,I))
c      END DO 
c      WRITE(*,*) 'LAST=',LAST
 101  FORMAT(3I6,G20.10,1X,'[',A,']')
      END 



      SUBROUTINE RLTEXT(TEXT,IA,IB,NAB) ! analyse text and return
      CHARACTER*(*) TEXT
*     IA = index of first non-blank character, or =1
*     IB = index of last non-blank character,  or =0 - comment excluded
*     NAB= index of last non-blank character (=0 for blank text)
*     ...
      IA=0
      IB=0
      NAB=0
      LIM=0
      DO I=1,LEN(TEXT)
       IF(TEXT(I:I).NE.' ') NAB=I
       IF((I.EQ.1.AND.TEXT(1:1).EQ.'*').OR.TEXT(I:I).EQ.'!') THEN
           IF(LIM.EQ.0) LIM=I
       END IF
      END DO
      IF(LIM.EQ.0) THEN
         LIM=NAB
      ELSE
         LIM=LIM-1
      END IF
      DO I=1,LIM
       IF(IA.EQ.0.AND.TEXT(I:I).NE.' ') IA=I
       IF(TEXT(I:I).NE.' ') IB=I
      END DO
      END

      
      FUNCTION MATINT(PAT,TEXT,NPAT,NTEXT)
*     approximate string matching - case insensitive 
*     ------------------------------------------------------------------
*     return number of matches of string PAT in string TEXT, 
*     and number NPAT, NTEXT of characters of string PAT and string TEXT
*     strings are considered from first to last non-blank character
*     Example:
*     MATCH = MATCH(' keYs ','keyWO RD',NPAT,NTEXT)
*     returns MATCH=3, NPAT=4, NTEXT=8       
*     ------------------------------------------------------------------
      CHARACTER*(*) PAT,TEXT
      INTEGER ID(0:100,2) 
*
      LOGICAL START                        ! for case conversion
      CHARACTER*26 CHU,CHL
      INTEGER NJ(0:255)
      DATA  CHU/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA  CHL/'abcdefghijklmnopqrstuvwxyz'/
      DATA  START/.TRUE./
      DATA NJ/256*0/
      SAVE
*     ...
      IF(START) THEN
         START=.FALSE.
         DO J=0,255
          NJ(J)=J
         END DO 
         DO I=1,26
          NJ(ICHAR(CHL(I:I)))=ICHAR(CHU(I:I))
         END DO
      END IF
*     ...
      MATINT=0
      NTEXT=0 
      DO I=1,LEN(TEXT)               ! find indices ITA...ITB
       IF(TEXT(I:I).NE.' ') GOTO 10
      END DO
      GOTO 15
 10   ITA=I
      DO I=ITA,LEN(TEXT)
       IF(TEXT(I:I).NE.' ') ITB=I
      END DO
      NTEXT=ITB-ITA+1               ! number of charcaters in TEXT

 15   NPAT=0
      DO I=1,LEN(PAT)               ! find indices IPA...IPB
       IF(PAT(I:I).NE.' ') GOTO 20
      END DO
      RETURN
 20   IPA=I
      DO I=IPA,LEN(PAT)
       IF(PAT(I:I).NE.' ') IPB=I
      END DO
      NPAT=IPB-IPA+1
      IF(NPAT.GT.100) STOP 'MATINT: string PAT too long!   '

      ID(0,1)=0 
      DO I=0,NPAT
       ID(I,2)=I
      END DO
      JOT=2
      
      DO J=1,NTEXT                                                                  
       JOT=3-JOT
       JT=J+ITA-1
       JC=NJ(ICHAR(TEXT(JT:JT)))
       DO I=1,NPAT  
        IP=I+IPA-1
        IDEQ=ID(I-1,3-JOT) 
        IC=NJ(ICHAR(PAT(IP:IP)))           
        IF(IC.NE.JC) IDEQ=IDEQ+1 
        ID(I,JOT)=MIN(IDEQ,ID(I,3-JOT)+1,ID(I-1,JOT)+1)
       END DO
       MATINT=MAX(MATINT,NPAT-ID(NPAT,JOT)) 
      END DO
      END                                                                        




