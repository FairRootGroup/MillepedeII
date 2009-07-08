


      SUBROUTINE MILLE(NLC,DERLC,NGL,DERGL,LABEL,RMEAS,SIGMA) ! add data
*         CALL MILLE(...)       ! measured value, derivatives (one set)                         
*         CALL ENDLE            ! complete, write record (many sets)
*     (or CALL KILLE            ! stop record)      
*     The data transmitted by MILLE calls are collected in two arrays,
*     a real array and an integer array, of same length. The collected
*     data are written at the ENDLE call. The content of the arrays:
*         real array              integer array    
*     1   0.0                     error count (this record)  
*     2   RMEAS, measured value   0                            JA
*     3   local derivative        index of local derivative
*     4   local derivative        index of local derivative
*     5    ...
*     6   SIGMA, error (>0)       0                            JB
*         global derivative       label of global derivative
*         global derivative       label of global derivative   IST
*         RMEAS, measured value   0
*         local derivative        index of local derivative
*         local derivative        index of local derivative
*         ...
*         SIGMA, error            0
*         global derivative       label of global derivative
*         global derivative       label of global derivative
*         ...
*     NR  global derivative       label of global derivative 
*
*     i.e. the 0's in the integer array allow to recognize the start
*     of a new set, the measured value and the error. The local and 
*     the global derivatives are inbetween, with a positive value in
*     the integer array, the index of the local derivative or the 
*     label of the global derivative.
*     ------------------------------------------------------------------
*     If more than one output unit is needed: duplicate this subroutine
*     change the entry names to e.g. AMILLE, AENDLE, AKILLE and change
*     the value of LUN and evtl. the dimension parameter in the 
*     parameter statements. 
*     -----------------------------------------------------------------
      PARAMETER     (LUN=51)    ! unit for record output
      PARAMETER     (NDIM=5000) ! dimension parameter of arrays 
      REAL     GLDER(NDIM)      ! real data record array
      INTEGER  INDER(NDIM)      ! integer data record array 
*     -----------------------------------------------------------------
      REAL DERLC(NLC),DERGL(NGL)
      INTEGER         LABEL(NGL)   
      SAVE
      DATA NR/0/                ! initial record length
      DATA ICOUNT/0/
*     ...
      IF(SIGMA.LE.0.0) RETURN   ! error zero - no measurement
      IF(NR.EQ.0) THEN
         NR=1
         GLDER(1)=0.0           
         INDER(1)=0             ! error counter
         ISP=0
      END IF   
      IF(NR+NLC+NGL+2.GT.NDIM) THEN
         ICOUNT=ICOUNT+1
         IF(ICOUNT.LE.10) THEN 
            WRITE(*,*) 'Mille warning: data can not be stored'
            IF(ICOUNT.EQ.10) THEN
               WRITE(*,*) 'Mille warning: no further printout'
            END IF
         END IF 
         INDER(1)=INDER(1)+1    ! count errors         
         RETURN                 ! record dimension too small
      END IF
      NR=NR+1
      GLDER(NR)=RMEAS           ! measured value 
      INDER(NR)=0  
      DO I=1,NLC                ! local derivatives
       IF(DERLC(I).NE.0.0) THEN 
          NR=NR+1
          GLDER(NR)=DERLC(I)    ! derivative of local parameter 
          INDER(NR)=I           ! index of local parameter
       END IF 
      END DO

      NR=NR+1
      GLDER(NR)=SIGMA           ! error of measured value
      INDER(NR)=0  
      DO I=1,NGL                ! global derivatives
       IF(DERGL(I).NE.0.0.AND.LABEL(I).GT.0) THEN 
          NR=NR+1
          GLDER(NR)=DERGL(I)    ! derivative of global parameter 
          INDER(NR)=LABEL(I)    ! index of global parameter
       END IF 
      END DO
      RETURN

      ENTRY MILLSP(NSP,DERGL,LABEL)
*     add NSP special words (floating-point and integer)
      
*     0.0            0    
*     -float(NSP)    0   ! indicates special data
*     following NSP floating and NSP integer data
* 
      IF(NSP.LE.0.OR.ISP.NE.0) RETURN
      ISP=NR
      IF(NR.EQ.0) THEN
         NR=1
         GLDER(1)=0.0
         INDER(1)=0             ! error counter
      END IF
      IF(NR+NSP+2.GT.NDIM) THEN
         INDER(1)=INDER(1)+1    ! count errors
         RETURN                 ! record dimension too small
      END IF
      NR=NR+1                   ! zero pair
      GLDER(NR)=0.0
      INDER(NR)=0
      NR=NR+1                   ! nsp and zero 
      GLDER(NR)=-FLOAT(NSP)
      INDER(NR)=0
      DO I=1,NSP
       NR=NR+1
       GLDER(NR)=DERGL(I)       ! floating-point 
       INDER(NR)=LABEL(I)       ! integer 
      END DO
      RETURN

      ENTRY KILLE                             ! stop record
      NR=0    ! reset
      RETURN

      ENTRY ENDLE                             ! end-of-record
      IF(NR.GT.1) THEN
         WRITE(LUN) NR+NR,(GLDER(I),I=1,NR),(INDER(I),I=1,NR)
      END IF
      NR=0    ! reset 
      END 


        
#ifdef wegdamit
      SUBROUTINE MILLE(NLC,DERLC,NGL,DERGL,LABEL,RMEAS,SIGMA) ! add data
*         CALL MILLE(...)       ! measured value, derivatives (one set)                         
*         CALL ENDLE            ! complete, write record (many sets)
*     (or CALL KILLE            ! stop record)      
*     The data transmitted by MILLE calls are collected in two arrays,
*     a real array and an integer array, of same length. The collected
*     data are written at the ENDLE call. The content of the arrays:
*         real array              integer array    
*     1   RMEAS, measured value   0
*     2   local derivative        index of local derivative
*     3   local derivative        index of local derivative
*         ...
*         SIGMA, error            0 
*         global derivative       label of global derivative
*         global derivative       label of global derivative
*         RMEAS, measured value   0
*         local derivative        index of local derivative
*         local derivative        index of local derivative
*         ...
*         SIGMA, error            0
*         global derivative       label of global derivative
*         global derivative       label of global derivative
*         ...
*     NR  global derivative       label of global derivative 
*
*     i.e. the 0's in the integer array allow to recognize the start
*     of a new set, the measured value and the error. The local and 
*     the global derivatives are inbetween, with a positive value in
*     the integer array, the index of the local derivative or the 
*     label of the global derivative.
*     ------------------------------------------------------------------
*     If more than one output unit is needed: duplicate this subroutine
*     change the entry names to e.g. AMILLE, AENDLE, AKILLE and change
*     the value of LUN and evtl. the dimension parameter in the 
*     parameter statements. 
*     -----------------------------------------------------------------
      PARAMETER     (LUN=51)    ! unit for record output
      PARAMETER     (NDIM=5000) ! dimension parameter of arrays 
      REAL     GLDER(NDIM)      ! real data record array
      INTEGER  INDER(NDIM)      ! integer data record array 
*     -----------------------------------------------------------------
      REAL DERLC(NLC),DERGL(NGL)
      INTEGER         LABEL(NGL)   
      SAVE
      DATA NR/0/                ! initial record length
*     ...
      IF(SIGMA.LE.0.0) RETURN   ! error zero - no measurement
      IF(NR.EQ.0) THEN
         NR=1
         GLDER(1)=0.0           
         INDER(1)=0             ! error counter
      END IF   
      IF(NR+NLC+NGL+2.GT.NDIM) THEN
         INDER(1)=INDER(1)+1    ! count errors         
         RETURN                 ! record dimension too small
      END IF
      NR=NR+1
      GLDER(NR)=RMEAS           ! measured value 
      INDER(NR)=0  
      DO I=1,NLC                ! local derivatives
       IF(DERLC(I).NE.0.0) THEN 
          NR=NR+1
          GLDER(NR)=DERLC(I)    ! derivative of local parameter 
          INDER(NR)=I           ! index of local parameter
       END IF 
      END DO

      NR=NR+1
      GLDER(NR)=SIGMA           ! error of measured value
      INDER(NR)=0  
      DO I=1,NGL                ! global derivatives
       IF(DERGL(I).NE.0.0.AND.LABEL(I).GT.0) THEN 
          NR=NR+1
          GLDER(NR)=DERGL(I)    ! derivative of global parameter 
          INDER(NR)=LABEL(I)    ! index of global parameter
       END IF 
      END DO
      RETURN

      ENTRY KILLE                             ! stop record
      NR=0
      RETURN

      ENTRY ENDLE                             ! end-of-record
      IF(NR.GT.1) THEN
         WRITE(LUN) NR+NR,(GLDER(I),I=1,NR),(INDER(I),I=1,NR)
      END IF
      NR=0 
      END 
#endif

