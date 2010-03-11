
      SUBROUTINE MEGNIT ! ----------------------------- initialize array
*     internal initialization of memory and defaults
      IMPLICIT NONE  
#include "dynal.inc"
      CHARACTER*40 SUTEXT(0:39)
      CHARACTER*3  FMT
      CHARACTER*(*) TEXT
      CHARACTER*(*) CH, CCH*1
      CHARACTER*1 CHBYI
      CHARACTER*36 CHARRA
      LOGICAL FULLPR,DONE,DEBUG
      INTEGER NWORDS,NDL,LENSML,LENDIF,IBYCH,LLT,LUNDEB,MQI,ICI,MQU9
      INTEGER I,ISTOP,IU,IL,ID,IP,IAI,IA,IBI,IFMT,LCH,ICH,ND,IND,NDC  
      DATA CHARRA/'0123456789abcdefghijklmnopqrstuvwxyz'/
      DATA DONE/.FALSE./
      SAVE 

*     ------------------------------------------------------------------
*     This subroutine allows to subdivide a large array dynamically.
*     The total dimension is defined in a PARAMETER statement. The
*     first call has to be 
*        CALL MEGNIT
*     to initialize the subroutine. Further calls of MEGNIT are ignored.
*     To restart
*        CALL MEGRES
*     which removes all arrays.  
*
*     A short overview over sub-arrays is printed by
*        CALL MEGPRN
*
*     A subarray 'a' is created and the index of the sub-array with ND
*     4-byte words is obtained by 
*        CALL MEGARR('a',ND,Text)
*     With the same call the format for printing can be defined by
*        CALL MEGARR('a I',ND,Text)    integer
*        CALL MEGARR('a S',ND,Text)    real (floating-point)
*        CALL MEGARR('a F',ND,Text)    real (floating-point)
*        CALL MEGARR('a D',ND,Text)    double precision (ND/2 DP words)
*          
*     All words of the sub-array are set to zero.        
*
*     The index INDA is defined and can be used by:
*        MQ(INDA   +1) ... MQ(IND    +ND)    Integer
*        QM(INDA   +1) ... QM(IND    +ND)    Single prec floating pt. 
*        DQ(INDA/2 +1) ... DQ(INDA/2 +ND/2)  Double prec floatinf pt.
*     i.e. for double precision the number ND has to be twice the 
*     number of double precision floating pt. words.
*     The number ND should be even (an odd number is increased by 1).  
*
*     Up to 36 sub-arrays are allowed, with the 26 letters as first
*     argument of the call (case-insensitive) or a digit 0...9 or 
*     characters . :   ;  !.    
*     The corresponding index
*     is IND with the same letter, e.g. the sub-array of the call 
*        CALL MEGARR('h',ND,'Normal matrix')   
*     is INDH and for 
*        CALL MEGARR('7',ND,'gradient vector')
*     the index is IND7. One can also use the index INDICS(I) with
*     i from 0 to 39. 
*     Special charcaters:   .   IDOT1  
*                           :   IDOT2
*                           ;   IDOT3
*                           !   IDOT4   
*
*     character:  0  1  2 ...  8  9  a  b  c  d ...  y  z  . :   ;  !
*     index I:    0  1  2 ...  8  9 10 11 12 13 ... 34 35 36 37 38 39 
*     
*     The text given in the call is printed in MEGPRN.
*     The programs stops with an error printout if to many words are
*     requested. The number of free words can be obtained by the call
*        CALL MEGSPC(NWORDS)
*     at any time to avoid an unexpected stop.
*
*     The length ND of the sub-array can be changed by an identical call
*     with the new length ND. If the length is increased, the additional
*     words are set to zero, other words are unchanged
*
*        CALL MEGRM('a')   removes array 'a'
*     There is an option to remove an existing array by an identical 
*     call with ND < 0. Use the call MEGRM('a') instead.
*     The total array is always reorganised to avoid unused space
*     between used space.    
*
*     Each sub-array has 10 extra words in addition to the ND words, 
*     where 2 words are used by the subroutine (and should never be 
*     modified), and 8 words may be freely used:
*      MQ(INDA -9) is index of index (0...35) (not to be changed)
*     Word MQ(INDA -8) may contain a format code: 1=I, 2=S, 3=D, which
*     is used in the call MEGAPR to select the format of the data part.
*                                     
*         MQ(INDA -7) 
*            ...      can be freely used
*         MQ(INDA -1)
*
*     The format code in word MQ(INDA -8) is set in a call of MEGARR,
*     if one of the (upper case) characters I, S or D is given behind
*     the character, which idenfies the sub-array:  
*     CALL MEGARR('h S',ND,'Normal matrix') means single precision fmt.   
*
*      MQ(INDA   ) is ND                      (not to be changed)
*      MQ(INDA +1) is first data word 
*         ...
*      MQ(INDA +ND) is last data word
*
*     The first 10 words in the common are reserved, the first word
*     contains the total number of used words: 
*        MQU= number of words used
*     MEGA-MQU is the number of unused words.
*
*     A sub-array including all data words is printed by
*        CALL MEGAPR('a)     in format as defin4ed by MEGARR.    
*     A different format can be forced by
*        CALL MEGIPR('a')    as integer words
*        CALL MEGFPR('a')    as single precision fl.pt. words
*        CALL MEGDPR('a')    as double precision fl.pt. words
*     for sub-array a. 
*     By the call 
*        CALL MEGZER('a')
*     all data words of subarray a are zeroed
*
*
*        CALL MEGCP('ab')    subarray a is copied to subarray b
*     Nothing is done if there is no subarray a. An existing subarray b
*     is replaced by subarray a.
*
*        CALL MEGMV('ab')    subarray a is moved to subarray b
*     Nothing is done if there is no subarray a. An existing subarray b
*     is replaced by subarray a. Subarray a is removed. 
*     
*        CALL MEGDEB         switch on/off debugging.
*                            Debugging printout is printed to file
*                            dynal.txt on unit 75.
*
*        CALL MEGDB          print dump of whole (!) array
*       
*
*     ------------------------------------------------------------------
      IF(DONE) RETURN
      ENTRY MEGRES    ! reset all = restart ----------------------------
      DONE=.TRUE. 
      DEBUG=.FALSE.
      FULLPR=.FALSE.
      LUNDEB=6
      DO I=1,10       ! reserved words 
       MQ(I)=0
      END DO
      MQU=0           ! number of used words 
      ISTOP=0
      IU=ICHAR('A')-1 ! upper case
      IL=ICHAR('a')-1 ! lower case   
      ID=ICHAR('0')   ! digit
      DO I=0,39
       INDICS(I)=0 
      END DO
      WRITE(*,*) ' '
      WRITE(*,*) 'Dynamic array allocation with dimension ',
     +            1.0E-6*FLOAT(MEGA),' Mio words'
      CALL PETIME
      RETURN

      ENTRY MEGPRN !------------------------------------- print overview
 02   WRITE(*,*) ' '
      WRITE(*,*) 'Dynamic array [allocated with dimension ',
     +            1.0E-6*FLOAT(MEGA),' Mio words]:'
      WRITE(*,*) ' '
      IP=0
      DO I=0,39
       IF(INDICS(I).NE.0) THEN
          IF(IP.EQ.0) WRITE(*,101)
          IF(IP.EQ.0) WRITE(*,103)
          IP=1 
          FMT=' '
          IF(MQ(INDICS(I)-8).EQ.1) FMT='(I)' 
          IF(MQ(INDICS(I)-8).EQ.2) FMT='(S)' 
          IF(MQ(INDICS(I)-8).EQ.3) FMT='(D)' 
          WRITE(*,102) CHBYI(I),FMT,I,INDICS(I),
     +    MQ(INDICS(I)),SUTEXT(I)
 101      FORMAT('    ch fmt   i      index nr_of_words')
 102      FORMAT(5X,A,1X,A3,1X,I3,I11,I12,2X,A)
 103      FORMAT('    ---------------------------------')
 104      FORMAT('    Total number_of_words ',I11) 
       END IF
      END DO
      IF(IP.EQ.0.AND.ISTOP.EQ.0) RETURN 
      WRITE(*,103)
      WRITE(*,104) MQU
      CALL PETIME
      WRITE(*,*) ' '
      IF(ISTOP.EQ.0) RETURN
      WRITE(*,*) 'STOP follows ...'
      STOP
  01  ISTOP=1     ! error with following stop
      WRITE(*,*) ' '
      WRITE(*,*) 'Error: too many words required for array ',CCH
      WRITE(*,*) 'Request for array ',CCH,' with',ND,' words.'
      STOP 
      GOTO 02
   
      ENTRY MEGSPC(NWORDS) ! --------------- return number of free words
      NWORDS=MEGA-MQU-10              ! number of free words  
      RETURN

      ENTRY MEGARR(CH,NDL,TEXT) ! ---------------------------- sub-array
      IF(DEBUG) THEN
         LLT=LEN(TEXT)
         IF(LLT.LT.50) WRITE(LUNDEB,*) 'MEGARR ',CH,NDL,'   ',TEXT
         IF(LLT.GE.50) WRITE(LUNDEB,*) 'MEGARR ',CH,NDL,'   ',TEXT(1:50)
      END IF
      CCH=CH
      IAI=IBYCH(CH)
      IF(IAI.LT.0) RETURN 
      IFMT=0                            ! format code ?
      LCH=LEN(CH)
      ICH=1
 03   IF(ICH.LT.LCH) THEN
         ICH=ICH+1
         IF(CH(ICH:ICH).EQ.'I') IFMT=1   ! integer
         IF(CH(ICH:ICH).EQ.'S') IFMT=2   ! single precision
         IF(CH(ICH:ICH).EQ.'F') IFMT=2   !      -"-
         IF(CH(ICH:ICH).EQ.'D') IFMT=3   ! double precision
         IF(IFMT.EQ.0) GOTO 03 
      END IF                    
      ND=NDL

      IF(ND.GT.0) ND=2*((ND+1)/2)       ! make ND > 0 even
      IF(INDICS(IAI).EQ.0) THEN 
         IF(NDL.LT.0) GOTO 01           ! probably too large
         SUTEXT(IAI)=' '                !
         IF(ND.LT.0) RETURN   
c         WRITE(*,*) 'Compare MQU+10+ND.GT.MEGA ',MQU,ND,MEGA                 
         MQU9=MQU+10+ND
         IF(.NOT.(MQU9.LE.MEGA)) GOTO 01 ! will work for MQU+10+ND > 2**31-1 (NaN)
         IF(MQU9.LT.0) GOTO 01           ! integer overflow into negative numbers
         IND=MQU+10
         INDICS(IAI)=IND                ! new index
         DO I=-9,ND
          MQ(IND+I)=0                   ! reset all
         END DO
         MQ(IND-8)=IFMT                 ! set format
         MQ(IND-9)=IAI                  ! array number
         MQ(IND)=ND                     ! array length   
         MQU=MQU+10+ND                  ! update "used" counter
         SUTEXT(IAI)=TEXT
         IF(DEBUG) THEN
c            WRITE(*,*) 'MEGARR: new array ',CH(1:1),' with',NDL,' words'
c            WRITE(*,*) '        ',SUTEXT(IAI) 
         END IF 
         RETURN
      ELSE IF(ND.LT.0) THEN             ! remove array
         IND=INDICS(IAI) 
         NDC=-10-MQ(IND) 
         IF(DEBUG) THEN
c            WRITE(*,*) 'MEGARR: array ',CH(1:1),' with',MQ(IND),
c     +                 ' words removed'
c            WRITE(*,*) '       ',SUTEXT(IAI)
         END IF
         IA=IND+MQ(IND) 
         DO I=1+IA,MQU
          MQ(I+NDC)=MQ(I)
         END DO

         INDICS(IAI)=0
         SUTEXT(IAI)=' '
         MQU=MQU+NDC
         IA=IA+NDC
      ELSE
         IND=INDICS(IAI)
         IF(DEBUG) THEN
c            WRITE(*,*) 'MEGARR: array ',CH(1:1),', length change from',
c     +      MQ(IND),' to ',ND,' words'
c            WRITE(*,*) '   old  ',SUTEXT(IAI)
c            WRITE(*,*) '   new  ',TEXT
         END IF
         NDC=ND-MQ(IND)                 ! change in length
         IF(NDC.EQ.0) RETURN
        
c         WRITE(*,*)'Compare MQU+10+NDC.GT.MEGA ',MQU,NDC,MEGA,ND,MQ(IND)
         MQU9=MQU+NDC
         IF(.NOT.(MQU9.LE.MEGA)) GOTO 01 ! will work for MQU+NDC > 2**31-1 (NaN)
         IF(MQU9.LT.0) GOTO 01           ! integer overflow into negative numbers
         IA=IND+MQ(IND)                 ! last word old array
         IF(NDC.LT.0) THEN
            DO I=1+IA,MQU         
             MQ(I+NDC)=MQ(I)
            END DO
            MQ(IND)=ND 
         ELSE
            DO I=MQU,1+IA,-1
             MQ(I+NDC)=MQ(I)
            END DO
            DO I=1+IA,IND+ND            ! reset new elements
             MQ(I)=0
            END DO
            MQ(IND)=ND 
         END IF
         IA=IA+NDC                      ! last word new array                     
         MQU=MQU+NDC 
         IF(TEXT.NE.' ') SUTEXT(IAI)=TEXT
      END IF 
*     reassign indices
 10   IF(IA.GE.MQU) THEN
         RETURN
      END IF       
      IAI=MQ(IA+1)
      INDICS(IAI)=IA+10                 ! new index
      IA=IA+10+MQ(INDICS(IAI))
      GOTO 10

      ENTRY MEGRM(CH)                ! remove sub-array
      IF(DEBUG) THEN 
         WRITE(LUNDEB,*) 'MEGRM  ',CH
      END IF 
      IAI=IBYCH(CH)
      IF(IAI.LT.0) RETURN

      IND=INDICS(IAI) 
      IF(IND.EQ.0) THEN
         SUTEXT(IAI)=' '               
         RETURN
      END IF
      IF(DEBUG) THEN
c            WRITE(*,*) 'MEGRM:  array ',CH(1:1),' with',MQ(IND),
c     +      ' words removed'
c            WRITE(*,*) '        ',SUTEXT(IAI)
      END IF

      NDC=-10-MQ(IND) 
      IA=IND+MQ(IND) 
      DO I=1+IA,MQU
       MQ(I+NDC)=MQ(I)
      END DO
      INDICS(IAI)=0
      SUTEXT(IAI)=' '
      MQU=MQU+NDC
      IA=IA+NDC

 30   IF(IA.GE.MQU) RETURN
      IAI=MQ(IA+1)
      INDICS(IAI)=IA+10                 ! new index
      IA=IA+10+MQ(INDICS(IAI))
      GOTO 30

      ENTRY MEGCP(CH)              ! copy from CHA to CHAB
      IF(DEBUG) THEN
         WRITE(LUNDEB,*) 'MEGCP  ',CH
      END IF
      CCH=CH(2:2)
      IBI=IBYCH(CCH)
      IF(IBI.LT.0) RETURN
C      IBI=IAI
      CCH=CH(1:1)
      IAI=IBYCH(CCH)
      IF(IAI.LT.0) RETURN
      IF(INDICS(IAI).EQ.0) RETURN       ! sub-array missing
      IA=0
      IF(INDICS(IBI).NE.0) THEN         ! remove subarray
         IND=INDICS(IBI)
         NDC=-10-MQ(IND)
         IA=IND+MQ(IND)
         DO I=1+IA,MQU
          MQ(I+NDC)=MQ(I)
         END DO
         INDICS(IBI)=0
         SUTEXT(IBI)=' '
         MQU=MQU+NDC
         IA=IA+NDC
 11      IF(IA.LT.MQU) THEN
            ICI=MQ(IA+1)
            INDICS(ICI)=IA+10           ! new index
            IA=IA+10+MQ(INDICS(ICI))
            GOTO 11
         END IF 
      END IF
      ND=MQ(INDICS(IAI))             ! copy-length
      SUTEXT(IBI)='copy of '//SUTEXT(IAI)                !
      MQU9=MQU+10+ND
      IF(.NOT.(MQU9.LE.MEGA)) GOTO 01 ! will work for MQU+10+ND > 2**31-1 (NaN)
      IF(MQU9.LT.0) GOTO 01           ! integer overflow into negative numbers
      IND=MQU+10
      INDICS(IBI)=IND                ! new index
      DO I=-9,ND
       MQ(IND+I)=MQ(INDICS(IAI)+I)     ! copy all
      END DO
      MQ(IND-9)=IBI
      MQU=MQU+10+ND                  ! update "used" counter
      RETURN 

      ENTRY MEGDMP                 ! dump
      WRITE(LUNDEB,*) ' '
      WRITE(LUNDEB,*) 'MQU MEGA',MQU,MEGA
      WRITE(LUNDEB,*) 'INDICS',INDICS
      WRITE(LUNDEB,105) (MQ(I),I=1,MQU) 
      WRITE(LUNDEB,*) 'end of dump'
      WRITE(LUNDEB,*) ' '
      RETURN   

      ENTRY MEGMV(CH)              ! move from CHA to CHAB
      IF(DEBUG) THEN
         WRITE(LUNDEB,*) 'MEGMV  ',CH
      END IF
      CCH=CH(2:2)
      IBI=IBYCH(CCH)
      IF(IBI.LT.0) RETURN      
      CCH=CH(1:1)
      IAI=IBYCH(CCH)
      IF(IAI.LT.0) RETURN 
      IF(INDICS(IAI).EQ.0) RETURN       ! sub-array missing
      IA=0 
      IF(INDICS(IBI).NE.0) THEN         ! remove subarray  
         IND=INDICS(IBI) 
         NDC=-10-MQ(IND) 
         IA=IND+MQ(IND) 
         DO I=1+IA,MQU
          MQ(I+NDC)=MQ(I)
         END DO
         INDICS(IBI)=0
         SUTEXT(IBI)=' '
         MQU=MQU+NDC
         IA=IA+NDC
 12      IF(IA.LT.MQU) THEN
            ICI=MQ(IA+1)
            INDICS(ICI)=IA+10           ! new index
            IA=IA+10+MQ(INDICS(ICI))
            GOTO 12
         END IF 
      END IF
      MQ(INDICS(IBI)-9)=IBI
      INDICS(IBI)=INDICS(IAI)           ! move subarray
      SUTEXT(IBI)=SUTEXT(IAI)
      INDICS(IAI)=0 
      SUTEXT(IAI)=' '
      RETURN

      ENTRY MEGIPR(CH) ! --------------------------------------- integer
      IFMT=1
      GOTO 20
      ENTRY MEGFPR(CH) ! ------------------------------ single precision
      IFMT=2
      GOTO 20
      ENTRY MEGDPR(CH) ! -------------------------------double precision
      IFMT=3
      GOTO 20
      ENTRY MEGAPR(CH) ! ------------------------------format from array 
      IFMT=0           ! format undefined 
 20   CONTINUE
      CCH=CH
      IAI=IBYCH(CCH)
      IF(IAI.LT.0) RETURN
      IND=INDICS(IAI)
      IF(IND.EQ.0) RETURN
c      FULLPR=.TRUE.
c      IF(IFMT.EQ.0) FULLPR=.FALSE.

c      FULLPR=.FALSE. !!!!!!!!!!!!!!!!

      IF(IFMT.EQ.0) THEN
         IFMT=MQ(IND-8)
         IF(IFMT.LT.1.OR.IFMT.GT.3) IFMT=1
      END IF  
      WRITE(LUNDEB,*) ' '
      IF(IFMT.NE.3) THEN 
         WRITE(LUNDEB,*)
     +              'Sub-array ',CHARRA(IAI+1:IAI+1),' at index',IND,
     +              ' with',MQ(IND),' words:'
      ELSE
         WRITE(LUNDEB,*) 
     +              'Double precision sub-array ',CHARRA(IAI+1:IAI+1),
     +              ' at index',IND/2,' with',MQ(IND)/2,' words:'
      END IF   

      WRITE(LUNDEB,*) 'Control-words ',(MQ(IND+I),I=-8,-1)
      LENSML=200
      LENDIF=60
      IF(IFMT.EQ.1) THEN
         IF(MQ(IND).LE.LENSML.OR.FULLPR) THEN
            WRITE(LUNDEB,105) (MQ(IND+I),I=1,MQ(IND))
         ELSE
            WRITE(LUNDEB,105) (MQ(IND+I),I=1,LENDIF)
            WRITE(LUNDEB,*) ' ...'
            MQI=((MQ(IND)-LENDIF)/20)*20+1
            WRITE(LUNDEB,105) (MQ(IND+I),I=MQI,MQ(IND))
         END IF 
      ELSE IF(IFMT.EQ.2) THEN
         IF(MQ(IND).LE.LENSML.OR.FULLPR) THEN
            WRITE(LUNDEB,106) (QM(IND+I),I=1,MQ(IND))
         ELSE
            WRITE(LUNDEB,106) (QM(IND+I),I=1,LENDIF)
            WRITE(LUNDEB,*) ' ...'
            MQI=((MQ(IND)-LENDIF)/20)*20+1
            WRITE(LUNDEB,106) (QM(IND+I),I=MQI,MQ(IND))
         END IF
      ELSE IF(IFMT.EQ.3) THEN
         IF(MQ(IND).LE.LENSML.OR.FULLPR) THEN   
            WRITE(LUNDEB,106) (DQ(IND/2+I),I=1,MQ(IND)/2)
         ELSE
            WRITE(LUNDEB,106) (DQ(IND/2+I),I=1,LENDIF)
            WRITE(LUNDEB,*) ' ...'
            MQI=((MQ(IND)/2-LENDIF)/20)*20+1
            WRITE(LUNDEB,106) (DQ(IND/2+I),I=MQI,MQ(IND)/2)
         END IF

      END IF
      WRITE(LUNDEB,*) ' '
 105  FORMAT(5I14)
 106  FORMAT(5G14.4)
      RETURN

      ENTRY MEGZER(CH)                        ! zero content
      CCH=CH
      IAI=IBYCH(CCH)
      IF(IAI.LT.0) RETURN
      IND=INDICS(IAI)
      IF(IND.EQ.0) RETURN 
      DO I=1,MQ(IND)
       MQ(IND+I)=0 
      END DO
      RETURN

      ENTRY MEGDEB                            ! debug on/off
      IF(DEBUG) THEN
         WRITE(*,*) 'MEG... debugging switched off'
         LUNDEB=6  
      ELSE
         WRITE(*,*) 'MEG... debugging switched on'
         LUNDEB=75
         CALL MVOPEN(LUNDEB,'dynal.txt')  
      END IF 
      DEBUG=.NOT.DEBUG 
      FULLPR=.NOT.FULLPR

      END 



      SUBROUTINE PETIME                          ! print times
      IMPLICIT NONE
*     print the elapsed and total time
      REAL TA(2),RST,DELTA,RSTP,SECND1,SECND2
      INTEGER NCOUNT,MODE
      INTEGER NHOUR1,MINUT1,NSECD1,NHOUR2,MINUT2,NSECD2

      SAVE
      DATA NCOUNT/0/
*     ...
      MODE=1
 10   NCOUNT=NCOUNT+1
      CALL ETIME(TA,RST)
      IF(NCOUNT.GT.1) THEN
         DELTA=RST
               NSECD1=DELTA ! -> integer
               NHOUR1=NSECD1/3600
               MINUT1=NSECD1/60-60*NHOUR1
               SECND1=DELTA-60*(MINUT1+60*NHOUR1)
         DELTA=RST-RSTP
               NSECD2=DELTA ! -> integer
               NHOUR2=NSECD2/3600
               MINUT2=NSECD2/60-60*NHOUR2
               SECND2=DELTA-60*(MINUT2+60*NHOUR2)
         WRITE(*,101) NHOUR1,MINUT1,SECND1, NHOUR2,MINUT2,SECND2
      END IF 

      RSTP=RST
      RETURN

      ENTRY EATIME
      WRITE(*,*) 'EATIME'
      MODE=2
      GOTO 10
 101  FORMAT(I2,' h',I3,' min',F5.1,' sec total',22X,'elapsed',
     +       I2,' h',I3,' min',F5.1,' sec')
      END




      SUBROUTINE WAITST(ISTOP,TEXT)              ! wait to check printout
      IMPLICIT NONE
*     Input from keyboard:
*     blank             continue
*     q                 quit = stop
*     /                 print short overview over sub-arrays (MEGPRN)
*     a ... z           print array
*     A ... Z           print array
*     0 ... 9           print array
*     $                 MEGDMP 
*     -                 switch off
*      else nonblank    wait for keyboard input
*
*     The program is waiting for input from keyboard
*    
      CHARACTER*1 TDUMMY,CCH, TEXT*(*)
      CHARACTER*8 TBLANK
      LOGICAL START
      INTEGER IWAIT,ISTOP,IAI,IBYCH
      SAVE
      DATA START/.TRUE./
*     ...
      IF(.NOT.START) GOTO 05  
      START=.FALSE.
      IWAIT=1
 01   WRITE(*,*) ' '
      WRITE(*,*) 'Wait and continue:'
      WRITE(*,*) '   hit enter key to continue'
      WRITE(*,*) '   enter / to print the array list'
      WRITE(*,*) '   enter one character: 0 or a .. or A or Z to ',
     +              'print an array'
      WRITE(*,*) '   enter ! to sitch to debugging' 
      WRITE(*,*) '   enter $ to call MEGDMP'  
      WRITE(*,*) '   enter ? to get this text'
      WRITE(*,*) '   enter - to switch wait off' 
      WRITE(*,*) '   enter Q to stop immediately' 
      WRITE(*,*) '   enter q to stop in a controlled way'
      WRITE(*,*) ' '
      WRITE(*,*) ' '

 05   WRITE(*,*) TEXT
      ISTOP=0
 10   IF(IWAIT.EQ.0) RETURN 
      TDUMMY=' '
      WRITE(*,*) '[wait for input]  =>'
      READ(*,101) TDUMMY
 101  FORMAT(A)
      TBLANK=TDUMMY
C      WRITE(*,*) 'code=',ICHAR(TBLANK(1:1))
      IF(TBLANK.EQ.' ') TBLANK=' (blank)' 

      WRITE(*,*) '         [input]  =>',TBLANK    
      IF(TDUMMY.EQ.' ') RETURN
      IF(TDUMMY.EQ.'$') THEN
         CALL MEGDMP ! many lines
         GOTO 10
      END IF
      IF(TDUMMY.EQ.'!') THEN
         CALL MEGDEB
         GOTO 10
      END IF 
      IF(TDUMMY.EQ.'Q') THEN
         STOP
      END IF 
      IF(TDUMMY.EQ.'q') THEN
         ISTOP=1
         RETURN
      END IF
      IF(TDUMMY.EQ.'?') THEN
         GOTO 01
      END IF
      IF(TDUMMY.EQ.'-') THEN
         IWAIT=0
         RETURN
      END IF
      IF(TDUMMY.EQ.'/') THEN
         CALL MEGPRN
         GOTO 05
      END IF 
      CCH=TDUMMY
      IAI=IBYCH(CCH)
      IF(IAI.GE.0) THEN
         CALL MEGAPR(CCH)
      ELSE
         CALL MEGPRN
      END IF
      IF(TDUMMY.NE.' ') GOTO 10
      END


      INTEGER FUNCTION IBYCH(CH)                         ! return index from char
      IMPLICIT NONE
      CHARACTER*(*) CH, CCH*1
      INTEGER IL,IU,ID,ICH,IAI
      SAVE 
*     ...
      CCH=CH 
      IU=ICHAR('A')-1 ! upper case
      IL=ICHAR('a')-1 ! lower case
      ID=ICHAR('0')   ! digit
      ICH=ICHAR(CCH)
      IF(ICH-IU.GT.0.AND.ICH-IU.LE.26)      THEN
         IAI=9+ICH-IU         ! upper case A...Z
      ELSE IF(ICH-IL.GT.0.AND.ICH-IL.LE.26) THEN
         IAI=9+ICH-IL         ! lower case a...z
      ELSE IF(ICH-ID.GE.0.AND.ICH-ID.LE.9)  THEN
         IAI=ICH-ID           ! digit 0...9
      ELSE IF(CCH.EQ.'.') THEN
         IAI=36
      ELSE IF(CCH.EQ.':') THEN
         IAI=37
      ELSE IF(CCH.EQ.';') THEN
         IAI=38
      ELSE IF(CCH.EQ.'!') THEN
         IAI=39
      ELSE
         IAI=-1
      END IF
      IF(IAI.GE.40) IAI=-1
      IBYCH=IAI
      IF(IAI.GE.0) RETURN
      WRITE(*,*) 'Illegal character:',CH,ICHAR(CH),' - stop'
      STOP
      END 

      CHARACTER*1 FUNCTION CHBYI(IAI)            ! return array characater
      IMPLICIT NONE
      CHARACTER*40 CHARS
      INTEGER IAI
      SAVE
      DATA CHARS/'0123456789abcdefghijklmnopqrstuvwxyz.:;!'/ 
*     ...
      IF(IAI.LT.0.OR.IAI.GT.39) THEN
         CHBYI=' '
      ELSE
         CHBYI=CHARS(IAI+1:IAI+1)
      END IF
      END 

  
      INTEGER FUNCTION INTVAL(NSEC)                      ! return time diff
      IMPLICIT NONE  
*
      REAL TA(2),RST,RSTI,RSTP
      LOGICAL START
      INTEGER NSEC,NT 
      SAVE
      DATA START/.TRUE./
*     ...
      INTVAL=0
      IF(START) THEN
         START=.FALSE.
         CALL ETIME(TA,RSTI)
         RSTP=RSTI
      ELSE
         CALL ETIME(TA,RST)
         INTVAL=RSTI-RST
         NT=RSTP-RST
         IF(NT+NSEC.GT.0) RETURN
         RSTP=RST
         INTVAL=RST-RSTI
      END IF
      END


*     vector tables ----------------------------------------------------

      SUBROUTINE MEGVEC(CH,VECT,N,NPLUS)         ! add vector
      IMPLICIT NONE 
#include "dynal.inc"
*
*     ind-2   length of vctor
*     ind-1   number of vectors
*     ind     array length
*
*     array length increase by 20 % plus 100 vectors
*
      CHARACTER*(*) CH, CHC*3
      INTEGER N,NPLUS,IAI,IND,NCOL,NROW,NROWDF,I,IBYCH
      REAL VECT(N)
      SAVE
*     ...
      IAI=IBYCH(CH)
      IF(IAI.LT.0) RETURN ! illegal character
      IND=INDICS(IAI)
      NCOL=N+NPLUS
      IF(IND.EQ.0) THEN
         CHC=CH(1:1)//' S'
         NROWDF=100      ! default length
         IF(N.LE.0.OR.NPLUS.LT.0) RETURN
         CALL MEGARR(CHC,NROWDF*NCOL,'vector storage')
         IND=INDICS(IAI)
      END IF
      IF(MQ(IND-1).EQ.0) MQ(IND-2)=NCOL
      IF(NCOL.NE.MQ(IND-2)) THEN
         RETURN
      END IF
      IF((MQ(IND-1)+1)*MQ(IND-2).GE.MQ(IND)) THEN
         CHC=CH(1:1)//' S'
         NROW=1.2*FLOAT(MQ(IND-1))+100.5
         CALL MEGARR(CHC,NROW*NCOL,'vector storage')
         IND=INDICS(IAI)
      END IF
      DO I=1,N
       QM(IND+NCOL*MQ(IND-1)+I)=VECT(I)
      END DO
      MQ(IND-1)=MQ(IND-1)+1
      END

      SUBROUTINE MEGLEN(CH)                      ! correct length
      IMPLICIT NONE
#include "dynal.inc"
      CHARACTER*(*) CH
      INTEGER IAI,IND,NCOL,NROW,IBYCH
      SAVE
*     ... 
      IAI=IBYCH(CH)
      IF(IAI.LT.0) RETURN ! illegal character
      IND=INDICS(IAI)
      IF(IND.EQ.0) RETURN ! no sub-array
      NCOL=MQ(IND-2)
      NROW=MQ(IND-1)
      CALL MEGARR(CH,NROW*NCOL,' ')
      END

      SUBROUTINE MEGRDV(CH,ITH,VECT,N)           ! read I.th vector
      IMPLICIT NONE  
#include "dynal.inc"
      CHARACTER*(*) CH 
      REAL VECT(*)
      INTEGER I,ITH,N,IAI,IND,NCOL,NROW,IBYCH
      SAVE
*     ...
      N=0
      IAI=IBYCH(CH)
      IF(IAI.LT.0) RETURN ! illegal character
      IND=INDICS(IAI)
      NCOL=MQ(IND-2)
      NROW=MQ(IND-1)
      IF(ITH.LE.0.OR.ITH.GT.NROW) RETURN
      N=NCOL
      DO I=1,N
       VECT(I)=QM(IND+NCOL*(ITH-1)+I)
      END DO
      RETURN

      ENTRY MEGSAV(CH,ITH,VECT,N)                ! save I.th vector 
      IAI=IBYCH(CH)
      IF(IAI.LT.0) RETURN ! illegal character
      IND=INDICS(IAI)
      NCOL=MQ(IND-2)
      NROW=MQ(IND-1)
      IF(ITH.LE.0.OR.ITH.GT.NROW) RETURN
      IF(N.NE.NCOL) RETURN
      DO I=1,N
       QM(IND+NCOL*(ITH-1)+I)=VECT(I)
      END DO
      END



 

