

*     Numerical subprograms used in mp II: matrix equations, 
*        and matrix products, double precision
*      
*     Solution by inversion
*        SQMINV
*
*     Solution by diagonalization
*        DEVROT, DEVPRT, DEFSOL,DEVINV
*
*     Solution by Cholesky decomposition of symmetric matrix  
*        CHOLDC
*
*     Solution by Cholesky decomposition of variable-band matrix 
*        VABDEC
*
*     Matrix/vector products
*        DBDOT     dot vector product
*        DBAXPY    multiplication and addition
*        DBSVX     symmetric matrix vector
*        DBGAX     general matrix vector
*        DBAVAT    AVAT product          
*        DBMPRV    print parameter and matrix 
*     
*     Chi^2 cut values
*        CHINDL    
*
*     Accurate summation
*        ADDSUM
*
*     Sorting
*        HEAPF    heap sort reals direct
*        SORT2K   sort 2-dim key-array 
*

*     matrix inversion -------------------------------------------------

      SUBROUTINE SQMINV(V,B,N,NRANK,DIAG,NEXT)   ! matrix inversion  
*     obtain solution of a system of linear equations with symmetric
*     matrix and the inverse.
*
*                    - - -
*        CALL SQMINV(V,B,N,NRANK,...,...)      solve  V * X = B
*                    - -   -----
*
*           V = symmetric N-by-N matrix in symmetric storage mode
*               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
*               replaced by inverse matrix
*           B = N-vector, replaced by solution vector
*
*     DIAG(N) =  double precision scratch array
*     NEXT(N) =  integer aux array
*
*     Method of solution is by elimination selecting the  pivot  on  the
*     diagonal each stage. The rank of the matrix is returned in  NRANK.
*     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
*     matrix V and the corresponding elements of  B  are  set  to  zero.
*
      DOUBLE PRECISION V(*),B(N),DIAG(N),VKK,VJK,EPS
      INTEGER        NEXT(N)
      PARAMETER     (EPS=1.0D-10)
      SAVE 
*     ...
      NEXT0=1
      DO I=1,N
       NEXT(I)=I+1                ! set "next" pointer
       DIAG(I)=ABS(V((I*I+I)/2))  ! save abs of diagonal elements
      END DO
      NEXT(N)=-1                  ! end flag

      NRANK=0
      DO I=1,N                    ! start of loop
       K  =0
       VKK=0.0D0

       J=NEXT0
       LAST=0 
 05    IF(J.GT.0) THEN
          JJ=(J*J+J)/2
          IF(ABS(V(JJ)).GT.MAX(ABS(VKK),EPS*DIAG(J))) THEN
             VKK=V(JJ)              
             K=J
             L=LAST
          END IF
          LAST=J
          J=NEXT(LAST)
          GOTO 05   
       END IF 

       IF(K.NE.0) THEN            ! pivot found
          KK=(K*K+K)/2
          IF(L.EQ.0) THEN
             NEXT0=NEXT(K)
          ELSE 
             NEXT(L)=NEXT(K)
          END IF
          NEXT(K)=0               ! index is used, reset
          NRANK=NRANK+1           ! increase rank and ...
          VKK    =1.0/VKK
          V(KK)  =-VKK
          B(K)   =B(K)*VKK
          JK     =KK-K
          JL     =0
          DO J=1,N                ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =V(JK)
              V(JK)=VKK*VJK
              B(J) =B(J)-B(K)*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  V(JL)=V(JL)-V(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE
          DO K=1,N
           IF(NEXT(K).NE.0) THEN
              B(K)=0.0D0       ! clear vector element
              DO J=1,K
               IF(NEXT(J).NE.0) V((K*K-K)/2+J)=0.0D0  ! clear matrix row/col
              END DO
           END IF
          END DO
          GOTO 10
       END IF
      END DO             ! end of loop
 10   DO IJ=1,(N*N+N)/2
       V(IJ)=-V(IJ)      ! finally reverse sign of all matrix elements
      END DO
      END


*     diagonalization --------------------------------------------------

      SUBROUTINE DEVROT(N,DIAG,U,V,WORK,IWORK)   ! diagonalization
*     Determination  of  eigenvalues  and  eigenvectors
*     of symmetric matrix V by  Householder
*     method
*                 -        -    
*     CALL DEVROT(N,DIAG,U,V,WORK)
*                   ---- -
*        N       = dimension of matrix 
*        DIAG(I) = diagonal elements
*        U(I,J)  = transformation matrix
*        V(.)    = symmetric matrix, unchanged        
*        WORK(N) = work array  
*
      PARAMETER (ITMAX=30,TOL=1.0E-15,EPS=1.0E-10)
      DOUBLE PRECISION V(*),DIAG(N),U(N,N),WORK(N)
      DOUBLE PRECISION F,G,H,SH,HH,B,P,R,S,C,WORKD
      INTEGER IWORK(N)
      SAVE
*     ...
*     1. part: symmetric matrix V reduced to tridiagonal from
      IJ=0
      DO I=1,N
       DO J=1,I
        IJ=IJ+1
        U(I,J)=V(IJ)  ! copy half of symmetric matirx
       END DO
      END DO
*
      DO I=N,2,-1
       L=I-2
       F=U(I,I-1)
       G=0.0D0
       IF(L.NE.0) THEN
          DO K=1,L
           IF(ABS(U(I,K)).GT.TOL) G=G+U(I,K)*U(I,K)
          END DO
          H=G+F*F
       END IF
       IF(G.LT.TOL) THEN   ! G too small
          WORK(I)=F           ! skip transformation
          H   =0.0D0
       ELSE
          L=L+1
          SH=SQRT(H)
          IF(F.GE.0.0D0) SH=-SH
          G=SH
          WORK(I)=SH
          H=H-F*G
          U(I,I-1)=F-G
          F=0.0D0
          DO J=1,L
           U(J,I)=U(I,J)/H
           G=0.0D0
*          form element of a u
           DO K=1,J
            IF(ABS(U(J,K)).GT.TOL.AND.ABS(U(I,K)).GT.TOL) THEN
               G=G+U(J,K)*U(I,K)
            END IF 
           END DO
           DO K=J+1,L
            IF(ABS(U(K,J)).GT.TOL.AND.ABS(U(I,K)).GT.TOL) THEN
               G=G+U(K,J)*U(I,K)
            END IF
           END DO
           WORK(J)=G/H
           F=F+G*U(J,I)
          END DO
*         form k
          HH=F/(H+H)
*         form reduced a
          DO J=1,L
           F=U(I,J)
           WORK(J)=WORK(J)-HH*F
           G=WORK(J)
           DO K=1,J
            U(J,K)=U(J,K)-F*WORK(K)-G*U(I,K)
           END DO
          END DO
       END IF
       DIAG(I)=H
      END DO
*
      DIAG(1)=0.0D0
      WORK(1)=0.0D0

*     accumulation of transformation matrices
      DO I=1,N
       IF(DIAG(I).NE.0.0) THEN
          DO J=1,I-1
           G=0.0D0
           DO K=1,I-1
            G=G+U(I,K)*U(K,J)
           END DO
           DO K=1,I-1
            U(K,J)=U(K,J)-G*U(K,I)
           END DO
          END DO
       END IF
       DIAG(I)=U(I,I)
       U(I,I)=1.0D0
       DO J=1,I-1
        U(I,J)=0.0D0
        U(J,I)=0.0D0
       END DO
      END DO
*     2. part: diagonalization of tridiagonal matrix
      DO I=2,N
       WORK(I-1)=WORK(I)
      END DO
      WORK(N)=0.0D0
      B=0.0D0
      F=0.0D0
*
      DO L=1,N
       J=0
       H=EPS*(ABS(DIAG(L))+ABS(WORK(L)))
       IF(B.LT.H) B=H
       DO M=L,N
        IF(ABS(WORK(M)).LE.B) GOTO 10 ! look for small sub-diagonal element
       END DO
       M=L
   10  IF(M.EQ.L) GOTO 30          
*      next iteration
   20  IF(J.EQ.ITMAX) RETURN       ! no convergence
       J=J+1
       G=DIAG(L)
       P=(DIAG(L+1)-G)/(2.0D0*WORK(L))
       R=SQRT(1.0D0+P*P)
       DIAG(L)=WORK(L)
       IF(P.LT.0.0D0) DIAG(L)=DIAG(L)/(P-R)
       IF(P.GE.0.0D0) DIAG(L)=DIAG(L)/(P+R)
       H=G-DIAG(L)
       DO I=L+1,N
        DIAG(I)=DIAG(I)-H
       END DO
       F=F+H
*      QL transformation
       P=DIAG(M)
       C=1.0D0
       S=0.0D0
       DO I=M-1,L,-1          ! reverse loop
        G=C*WORK(I)
        H=C*P
        IF(ABS(P).GE.ABS(WORK(I))) THEN
           C=WORK(I)/P
           R=SQRT(1.0D0+C*C)
           WORK(I+1)=S*P*R
           S=C/R
           C=1.0D0/R
        ELSE
           C=P/WORK(I)
           R=SQRT(1.0D0+C*C)
           WORK(I+1)=S*WORK(I)*R
           S=1.0D0/R
           C=C/R
        END IF
        P=C*DIAG(I)-S*G
        DIAG(I+1)=H+S*(C*G+S*DIAG(I))
*       form vector
        DO K=1,N
         H=U(K,I+1)
         U(K,I+1)=S*U(K,I)+C*H
         U(K,I)=C*U(K,I)-S*H
        END DO
       END DO 
       WORK(L)=S*P
       DIAG(L)=C*P
       IF(ABS(WORK(L)).GT.B) GOTO 20 ! next iteration
   30  DIAG(L)=DIAG(L)+F
      END DO 

      DO I=1,N
       IWORK(I)=I
      END DO

      M=1
 40   M=1+3*M    ! determine initial increment
      IF(M.LE.N) GOTO 40
 50   M=M/3
      DO J=1,N-M ! sort with increment M
       L=J
 60    IF(DIAG(IWORK(L+M)).GT.DIAG(IWORK(L))) THEN ! compare
          LL=IWORK(L+M)      ! exchange the two index values
          IWORK(L+M)=IWORK(L)
          IWORK(L)=LL
          L=L-M
          IF(L.GT.0) GOTO 60
       END IF
      END DO
      IF(M.GT.1) GOTO 50

      DO I=1,N
       IF(IWORK(I).NE.I) THEN
*         move vector from position I to the work area
          WORKD=DIAG(I)
          DO L=1,N
           WORK(L)=U(L,I)
          END DO
          K=I
 70       J=K
          K=IWORK(J)
          IWORK(J)=J
          IF(K.NE.I) THEN
*            move vector from position K to the (free) position J
             DIAG(J)=DIAG(K)
             DO L=1,N
              U(L,J)=U(L,K)
             END DO
             GOTO 70
          END IF
*         move vector from the work area to position J
          DIAG(J)=WORKD
          DO L=1,N
           U(L,J)=WORK(L)
          END DO
       END IF
      END DO
      END


      SUBROUTINE DEVPRT(N,DIAG,U,B,WORK,IWORK,X) ! print eigenvectors 
*     print 
      DOUBLE PRECISION B(N),U(N,N),DIAG(N),X(N),WORK(N),SUM
      INTEGER IWORK(N)
      SAVE
*     ...
      WRITE(*,*) 'Eigenvalues of matrix and conditions: '
      DO I=1,N
       COND=0.0D0
       IF(DIAG(I).NE.0.0) COND=DIAG(1)/DIAG(I)
       WRITE(*,101) I,DIAG(I),COND
      END DO
 101  FORMAT(I3,2G16.8)
      WRITE(*,*) ' '
      WRITE(*,*) 'Eigenvectors of matrix'
      DO I=1,N
       WRITE(*,102) I,(U(J,I),J=1,N)
      END DO
 102  FORMAT(I3,10F6.2/(3X,10F6.2))

*     solution

      WRITE(*,*) 'Solution '
 203  FORMAT('  I Coeff     error      Coeff    error')
      DO I=1,N
       SUM=0.0D0
       DO J=1,N
        SUM=SUM+U(J,I)*B(J)
       END DO
       COEFY=0.0
       ERRFY=0.0
       IF(DIAG(I).GT.0.0) THEN
          COEFY=SUM/DIAG(I)
          ERRFY=1.0/SQRT(DIAG(I)) 
       END IF

       COEFZ=0.0
       ERRFZ=0.0
       IF(DIAG(I).GT.0.0) THEN
          COEFZ=SUM/SQRT(DIAG(I))
          ERRFZ=1.0
       END IF
       WRITE(*,103) I,COEFY,CHAR(177),ERRFY,COEFZ,CHAR(177),ERRFZ
      END DO 
 103  FORMAT(I3,G16.8,1X,A1,1X,G12.4,5X,G16.8,1X,A1,1X,G12.4)

*     select significant terms

      DO I=1,N
       IWORK(I)=I
       SUM=0.0
       DO J=1,N
        SUM=SUM+U(J,I)*B(J)
       END DO
       WORK(I)=0.0
       IF(DIAG(I).GT.0.0) THEN 
          WORK(I)=ABS(SUM)/SQRT(DIAG(I))
       END IF 
      END DO
      WRITE(*,*) 'Significances ',(WORK(IWORK(I)),I=1,N)

*     shell sort

      M=4
 10   IF(3*M+1.LT.N) THEN    ! determine initial increment
         M=3*M+1
         GOTO 10
      END IF
 20   M=M/3                  ! next increment
      DO I=1+M,N             ! ... used in insert sort
       IDXT=IWORK(I)
       J=I
 30    IF(WORK(IWORK(J-M)).LT.WORK(IDXT)) THEN
          IWORK(J)=IWORK(J-M)
          J=J-M
          IF(J.GT.M) GOTO 30
       END IF
       IWORK(J)=IDXT
      END DO
      IF(M.GT.1) GOTO 20
      WRITE(*,*) 'Significances ',(WORK(IWORK(I)),I=1,N)
      WRITE(*,*) 'indices  ',(IWORK(I),I=1,N/2)

      DO I=1,N
       X(I)=0.0
      END DO

      DO IS=1,N/2
       I=IWORK(IS)
       SUM=0.0
       DO J=1,N
        SUM=SUM+U(J,I)*B(J)
       END DO
       YI=SUM/DIAG(I)
       DO J=1,N
        X(J)=X(J)+YI*U(J,I)  
       END DO
      END DO
      WRITE(*,*) 'X ',X
      END


      SUBROUTINE DEVSOL(N,DIAG,U,B,X,WORK)       ! solution  
*     Solution of matrix equation   V * X = B 
*     after diagonalization of V by SMDEV
*
*                 - ---- - - 
*     CALL DEVSOL(N,DIAG,U,B,X,WORK)
*                            - ---- 
*        N       = dimension of matrix 
*        DIAG(.) = diagonal elements      (from SMDEV)
*        U(I,J)  = transformation matrix  (   "-"    )
*        B(N)    = r.h.s. of matrix equation (unchanged)
*        X(N)    = solution vector
*        WORK(N) = work array  
*
      DOUBLE PRECISION B(N),U(N,N),DIAG(N),X(N),WORK(N),S
      SAVE
*     ...
      DO J=1,N
        S=0.0D0
        WORK(J)=0.0D0
        IF(DIAG(J).NE.0.0D0) THEN  
           DO I=1,N
*           j-th eigenvector is U(.,J)           
            S=S+U(I,J)*B(I)  
           END DO 
           WORK(J)=S/DIAG(J) 
        END IF 
      END DO 

      DO J=1,N
        S=0.0D0
        DO JJ=1,N 
         S=S+U(J,JJ)*WORK(JJ) 
        END DO 
        X(J)=S
      END DO 
      END


      SUBROUTINE DEVINV(N,DIAG,U,V)              ! inversion
*     get inverse matrix V from DIAG and U
*                 - ---- -   
*     CALL DEVINV(N,DIAG,U,V)
*                          -
*        N       = dimension of matrix 
*        DIAG(I) = diagonal elements
*        U(I,J)  = transformation matrix
*        V(.)    = smmmetric matrix        
*      
      DOUBLE PRECISION U(N,N),DIAG(N),V(*),SUM
      SAVE
*     ...
      IJ=0
      DO I=1,N
       DO J=1,I                          
        IJ=IJ+1
        SUM=0.0D0
        DO K=1,N
         IF(DIAG(K).NE.0.0D0) THEN
            SUM=SUM+U(I,K)*U(J,K)/DIAG(K)
         END IF 
        END DO
        V(IJ)=SUM
       END DO
      END DO
      END


*     Cholesky decomposition -------------------------------------------


      SUBROUTINE CHOLDC(G,N)                     ! decomposition 
*     
*     Cholesky decomposition of the matrix G:      G =  L  D  L^T
*        G = symmetric matrix, in symmetric storage mode  
*        L = unit triangular matrix (1's on diagonal)
*        D = diagonal matrix (elements store on diagonal of L)
*     The sqrts of the usual Cholesky decomposition are avoided by D. 
*     Matrices L and D are stored in the place of matrix G; after the 
*     decomposition, the solution of matrix equations and the computation
*     of the inverse of the (original) matrix G are done by the
*     entries CHOLSL and CHOLIN.
*
      DOUBLE PRECISION G(*),X(*),V(*),SUM,RATIO
      SAVE     
*     ...
      II=0
      DO I=1,N

c       WRITE(*,*) 'vor I=',I,'  ',(G((L*L+L)/2),L=1,N) 

       II=II+I
       IF(G(II).NE.0.0) G(II)=1.0/G(II)  ! (I,I) div !
       JJ=II
       DO J=I+1,N
        RATIO=G(I+JJ)*G(II)              ! (I,J) (I,I)
        KK=JJ
        DO K=J,N
         G(KK+J)=G(KK+J)-G(KK+I)*RATIO   ! (K,J) (K,I)
         KK=KK+K
        END DO ! K
        G(I+JJ)=RATIO                    ! (I,J)
        JJ=JJ+J
       END DO ! J

c       WRITE(*,*) 'nach I=',I,'  ',(G((L*L+L)/2),L=1,N)

      END DO ! I
      RETURN

      ENTRY CHOLSL(G,X,N)                        ! solution
*
*     The matrix equation  G X = B is solved for X, where the matrix
*     G in the argument is already decomposed by CHOLDC. The vector B
*     is called X in the argument and the content is replaced by the 
*     resulting vector X.  
*
      II=0
      DO I=1,N
       SUM=X(I)
       DO K=1,I-1
        SUM=SUM-G(K+II)*X(K)             ! (K,I)
       END DO
       X(I)=SUM
       II=II+I
      END DO
      DO I=N,1,-1
       SUM=X(I)*G(II)                    ! (I,I)
       KK=II
       DO K=I+1,N
        SUM=SUM-G(KK+I)*X(K)             ! (K,I)
        KK=KK+K
       END DO
       X(I)=SUM
       II=II-I
      END DO
      RETURN

      ENTRY CHOLIN(G,V,N)                        ! inverse matrix
*
*     The inverse of the (original) matrix G is computed and stored
*     in symmetric storage mode in matrix V. Arrays G and V must be 
*     different arrays. 
*
      II=(N*N-N)/2
      DO I=N,1,-1
       SUM=G(II+I)                       ! (I,I)
       DO J=I,1,-1
        DO K=J+1,N
         L=MIN(I,K)
         M=MAX(I,K)
         SUM=SUM-G(J+(K*K-K)/2)*V(L+(M*M-M)/2) ! (J,K) (I,K)
        END DO
        V(II+J)=SUM                      ! (I,J)
        SUM=0.0D0
       END DO
       II=II-I+1
      END DO
      END


*     variable band matrix operations ----------------------------------

      SUBROUTINE VABDEC(N,VAL,ILPTR)             ! A = L D L^T  decomposition
*     Variable_band matrix row Doolittle decomposition
*     A variable-band NxN symmetric matrix, also called skyline, is stored
*     row by row in the array VAL(.). For each row every coefficient
*     between the first non-zero element in the row and the diagonal is
*     stored.  
*     The pointer array ILPTR(N) contains the indices in VAL(.) of the
*     diagonal elements. ILPTR(1) is always 1, and ILPTR(N) is equal
*     to the total number of coefficients stored, called the profile.
*     The form of a variable-band matrix is preserved in the L D L^T
*     decomposition no fill-in is created ahead in any row or ahead of the
*     first entry in any column, but existing zero-values will become
*     non-zero. The decomposition is done "in-place".
*    
      DOUBLE PRECISION VAL(*),X(N),GAMMA,XI,EPS,ONE,TWO,PRD,VALKJ
      INTEGER ILPTR(N)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
*     DATA EPS/0.0D0/  
      DATA EPS/2.22044605E-16/ 
      SAVE
*     ...
      IF(EPS.EQ.0.0D0) THEN
         EPS    = TWO**(-12)
   10    EPS    = EPS/TWO
         PRD=ONE
         CALL DBAXPY(1,ONE,EPS,PRD)
         IF(PRD.GT.ONE) GOTO 10
         EPS=EPS*TWO                      ! EPS is machine presision
         WRITE(*,*) 'Machine precision is ',EPS
      END IF 

      WRITE(*,*) 'Variable band matrix Cholesky decomposition'

      GAMMA=0.0D0
      I=1
      DO J=1,ILPTR(N) ! loop thrugh all matrix elements
       IF(ILPTR(I).EQ.J) THEN ! diagonal element
          IF(VAL(J).LE.0.0D0) GOTO 01   ! exit loop for negative diag 
          GAMMA=MAX(GAMMA,ABS(VAL(J)))  ! max diagonal element
          I=I+1
       END IF 
      END DO
      I=N+1       
 01   IN=I-1      ! IN positive diagonal elements
      WRITE(*,*) '  ',IN,' positive diagonal elements' 
      XI=0.0D0
      I=1
      DO J=1,ILPTR(IN)          ! loop for positive diagonal elements
                                ! through all matrix elements
       IF(ILPTR(I).EQ.J) THEN   ! diagonal element
          I=I+1
       ELSE
          XI=MAX(XI,ABS(VAL(J))) ! Xi = abs(max) off-diagonal element 
       END IF 
      END DO

      DELTA=EPS*MAX(1.0D0,GAMMA+XI)  
      SN=1.0D0
      IF(N.GT.1) SN=1.0D0/SQRT(DFLOAT(N*N-1))
      BETA =SQRT(MAX(EPS,GAMMA,XI*SN))            ! beta 
      WRITE(*,*) '   DELTA and BETA ',DELTA,BETA

*     VAL(1)=VAL(1)             ! d_11 = a_11
      DO K=2,N      
       MK=K-ILPTR(K)+ILPTR(K-1)+1

       THETA=0.0

       DO J=MK,K
        MJ=J-ILPTR(J)+ILPTR(J-1)+1
        KJ=ILPTR(K)-K+J        ! index kj

        DO I=MAX(MJ,MK),J-1
         VAL(KJ)=VAL(KJ)       ! L_kj := L_kj - L_ki D_ii L_ji   
     +          -VAL(ILPTR(K)-K+I)*VAL(ILPTR(I))*VAL(ILPTR(J)-J+I)

c         WRITE(*,101)  'K J I  KJ KI II JI',K,J,I, K,J, I,I, J,I 
 101     FORMAT(A,3I4,3(3X,2I4))
        END DO ! 

        THETA=MAX(THETA,ABS(VAL(KJ)))  ! maximum value of row

c        WRITE(*,*) '                           K J JJ         ',K,J, J,J 
 102    FORMAT(A,2I4,3X,2I4)
        IF(J.NE.K) THEN
           IF(VAL(ILPTR(J)).NE.0.0D0) THEN
              VAL(KJ)=VAL(KJ)/VAL(ILPTR(J))
           ELSE
              VAL(KJ)=0.0D0
           END IF 
        END IF                                ! L_kj := L_kj/D_jj ! D_kk

        IF(J.EQ.K) THEN
           VALKJ=VAL(KJ)
c          WRITE(*,*) 'Indices K and J',K,J,VAL(KJ),KJ,ILPTR(K)
c          WRITE(*,*) VAL(KJ),(THETA/BETA)**2,DELTA,THETA
           IF(K.LE.IN) THEN
              VAL(KJ)=MAX(ABS(VAL(KJ)),(THETA/BETA)**2,DELTA)
              IF(VALKJ.NE.VAL(KJ)) THEN
                 WRITE(*,*) '   Index K=',K
                 WRITE(*,*) '   ',VALKJ,VAL(KJ),
     +                      (THETA/BETA)**2,DELTA,THETA 
              END IF
           END IF 
        END IF
       END DO ! J 

      END DO ! K

      DO K=1,N
        OLD=VAL(ILPTR(K))
        IF(VAL(ILPTR(K)).NE.0.0D0) VAL(ILPTR(K))=1.0D0/VAL(ILPTR(K))
c        WRITE(*,*) 'K old new ',OLD,VAL(ILPTR(K)) 
      END DO

      RETURN

      ENTRY VABMMM(N,VAL,ILPTR)                  ! minimum and maximum 
      KR=1        
      KS=1 
      KP=1
      DO K=1,N
       IF(VAL(ILPTR(K)).GT.VAL(ILPTR(KS))) KS=K
       IF(VAL(ILPTR(K)).LT.VAL(ILPTR(KR))) KR=K
       IF(VAL(ILPTR(K)).GT.0.0.AND.VAL(ILPTR(K)).LT.VAL(ILPTR(KP))) KP=K 
      END DO
      WRITE(*,*) '   Index value ',KS,VAL(ILPTR(KS))
      WRITE(*,*) '   Index value ',KP,VAL(ILPTR(KP))
      WRITE(*,*) '   Index value ',KR,VAL(ILPTR(KR))  

      RETURN

      ENTRY VABSLV(N,VAL,ILPTR,X)                ! solution X
*     The matrix equation  A X = B  is solved. The matrix is assumed to
*     decomposed before using VABDEC. The array X(N) contains on entry
*     the right-hand-side B(N); at return it contains the solution.    
*     ...
      DO K=1,N                  ! forward loop 
       MK=K-ILPTR(K)+ILPTR(K-1)+1        
       DO J=MK,K-1
        X(K)=X(K)-VAL(ILPTR(K)-K+J)*X(J)  ! X_k := X_k - L_kj B_j
       END DO
      END DO ! K

      DO K=1,N                  ! divide by diagonal elements
c       X(K)=X(K)/VAL(ILPTR(K))            ! X_k := X_k/D_kk

       X(K)=X(K)*VAL(ILPTR(K))            ! X_k := X_k*D_kk 

      END DO

      DO K=N,1,-1               ! backward loop
       MK=K-ILPTR(K)+ILPTR(K-1)+1
       DO J=MK,K-1
        X(J)=X(J)-VAL(ILPTR(K)-K+J)*X(K)  ! X_j := X_j - L_kj X_k
       END DO
      END DO ! K
      END 


*     matrix/vector products -------------------------------------------

      DOUBLE PRECISION FUNCTION DBDOT(N,DX,DY)   ! dot product
*     Dot product of two vectors: DX * DY
      DOUBLE PRECISION DX(*),DY(*),DTEMP
*     ...
      DTEMP=0.0D0
      DO I = 1,MOD(N,5)
       DTEMP=DTEMP+DX(I)*DY(I)
      END DO
      DO I =MOD(N,5)+1,N,5
       DTEMP=DTEMP+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)
     +                        +DX(I+3)*DY(I+3)+DX(I+4)*DY(I+4)
      END DO
      DBDOT=DTEMP
      END


      SUBROUTINE DBAXPY(N,DA,DX,DY)              ! multiply, addition
*     Constant times vector added to a vector: DY:=DY+DA*DX
      DOUBLE PRECISION DX(*),DY(*),DA
*     ...
      DO I=1,MOD(N,4)
        DY(I)=DY(I)+DA*DX(I)
      END DO
      DO I=MOD(N,4)+1,N,4
        DY(I  )=DY(I  )+DA*DX(I  )
        DY(I+1)=DY(I+1)+DA*DX(I+1)
        DY(I+2)=DY(I+2)+DA*DX(I+2)
        DY(I+3)=DY(I+3)+DA*DX(I+3)
      END DO
      END


      SUBROUTINE DBSVX(V,A,B,N)                  ! symm. matrix, vector
*     multiply symmetric N-by-N matrix and N-vector
*
*                   - -   - 
*        CALL DBSVX(V,A,B,N)
*                       -
*
*         B   :=    V   *    A
*         N        N*N       N
*
      DOUBLE PRECISION V,A,B
      DIMENSION V(*),A(*),B(*)
      DOUBLE PRECISION SUM
      IJS=1
      DO I=1,N
       SUM=0.0
       IJ=IJS
       DO J=1,N
        SUM=SUM+V(IJ)*A(J)
        IF(J.LT.I) THEN
           IJ=IJ+1
        ELSE
           IJ=IJ+J
        END IF
       END DO
       B(I)=SUM
       IJS=IJS+I
      END DO
      END


      SUBROUTINE DBGAX(A,X,Y,M,N)                ! gen. matrix, vector
*     multiply general M-by-N matrix A and N-vector X 
*
*                    - -   - -
*        CALL  DBGAX(A,X,Y,M,N)          Y   :=   A   *    X
*                        -               M       M*N       N
*
*        where A = general M-by-N matrix (A11 A12 ... A1N  A21 A22 ...) 
*              X = N vector 
*              Y = M vector
*
      DOUBLE PRECISION A(*),X(*),Y(*)
      SAVE
*     ...
      IJ=0
      DO I=1,M
       Y(I)=0.0D0
       DO J=1,N
        IJ=IJ+1
        Y(I)=Y(I)+A(IJ)*X(J)
       END DO  
      END DO 
      END


      SUBROUTINE DBAVAT(V,A,W,N,M)               ! A V AT product    
*     multiply symmetric N-by-N matrix from the left with general M-by-N
*     matrix and from the right with the transposed of the same  general
*     matrix  to  form  symmetric  M-by-M   matrix   (used   for   error
*     propagation).
*
*                    - -   - -                                   T
*        CALL DBAVAT(V,A,W,N,M)         W   =   A   *   V   *   A
*                        -             M*M     M*N     N*N     N*M
*
*        where V = symmetric N-by-N matrix
*              A = general N-by-M matrix
*              W = symmetric M-by-M matrix
*
      DOUBLE PRECISION V,A,W,CIK 
      DIMENSION V(*),A(*),W(*)
      SAVE
*     ...
      DO I=1,(M*M+M)/2
       W(I)=0.0                ! reset output matrix
      END DO
      IL=-N
      IJS=0
      DO I=1,M                 ! do I
       IJS=IJS+I-1             !
       IL=IL+N                 !
       LKL=0                   !
       DO K=1,N                !   do K
        CIK=0.0D0              !  
        LKL=LKL+K-1            ! 
        LK=LKL                 !
        DO L=1,K               !     do L 
         LK=LK+1               !     .
         CIK=CIK+A(IL+L)*V(LK) !     .
        END DO                 !     end do L
        DO L=K+1,N             !     do L 
         LK=LK+L-1             !     . 
         CIK=CIK+A(IL+L)*V(LK) !     .
        END DO                 !     end do L 
        JK=K                   !
        IJ=IJS                 !
        DO J=1,I               !     do J 
         IJ=IJ+1               !     .
         W(IJ)=W(IJ)+CIK*A(JK) !     .
         JK=JK+N               !     . 
        END DO                 !     end do J
       END DO                  !   end do K 
      END DO                   ! end do I 
      END

      SUBROUTINE DBMPRV(LUN,X,V,N)
*     Prints the n-vector X and the symmetric N-by-N  covariance  matrix
*     V, the latter as a correlation matrix.
*
*                    - - -
*        CALL DPRPAC(X,V,N)     X = vector of parameters
*                               V = covariance matrix
*
      DOUBLE PRECISION X,V
      DIMENSION X(*),V(*)
      INTEGER MC(15)
*
      WRITE(LUN,103)
      WRITE(LUN,101)
      II=0
      DO I=1,N
       IJ=II
       II=II+I
       ERR=0.0
       IF(V(II).GT.0.0) ERR=SQRT(V(II))
       L=0
       JJ=0
       DO J=1,I
        JJ=JJ+J
        IJ=IJ+1
        RHO=0.0
        PD=V(II)*V(JJ)
        IF(PD.GT.0.0) RHO=V(IJ)/SQRT(PD)
        L=L+1
        MC(L)=100.0*ABS(RHO)+0.5
        IF(RHO.LT.0.0) MC(L)=-MC(L)
        IF(J.EQ.I.OR.L.EQ.15) THEN
           IF(J.LE.15) THEN
              IF(J.EQ.I) THEN
                 WRITE(LUN,102) I,X(I),ERR,(MC(M),M=1,L-1)
              ELSE
                 WRITE(LUN,102) I,X(I),ERR,(MC(M),M=1,L)
              END IF
           ELSE
              IF(J.EQ.I) THEN
                 WRITE(LUN,103) (MC(M),M=1,L-1)
              ELSE
                 WRITE(LUN,103) (MC(M),M=1,L)
              END IF
              L=0
           END IF
        END IF
       END DO
      END DO
      WRITE(LUN,104)
  100 RETURN
  101 FORMAT(9X,'Param',7X,'error',7X,'correlation coefficients'/)
  102 FORMAT(1X,I5,2G12.4,1X,15I5)
  103 FORMAT(31X,15I5)
  104 FORMAT(33X,'(correlation coefficients in percent)')
      END



*     accurate summation -----------------------------------------------

      
      SUBROUTINE ADDSUM(ADD)                     ! accurate summation
      DOUBLE PRECISION ADD,SUM,DSUM
      DATA NSUM/0/,NEXP/0/,DSUM/0.0D0/,NEXP20/1048576/
      SAVE
*     ...
      NADD=ADD
      NSUM=NSUM+NADD
      DSUM=DSUM+(ADD-DFLOAT(NADD))
      IF(DSUM.GT.16.0D0) THEN
         DSUM=DSUM-16.D0
         NSUM=NSUM+16
      END IF
      IF(NSUM.GT.NEXP20) THEN
         NEXP=NEXP+1  
         NSUM=NSUM-NEXP20
      END IF
      RETURN
       
      ENTRY GETSUM(SUM)                       ! get accurate sum
      SUM=(DSUM+DFLOAT(NSUM))+DFLOAT(NEXP)*DFLOAT(NEXP20)
      DSUM=0.0D0
      NSUM=0
      NEXP=0
      END  

*     sort -------------------------------------------------------------

      SUBROUTINE HEAPF(A,N)                      ! heap sort direct (real)
*     real keys A(*), sorted at return
      INTEGER I,J,L,R,N
      REAL A(*),AT    ! array of keys and pivot key value
      SAVE
*     ...
      IF(N.LE.1) RETURN
      L=N/2+1
      R=N
   10 IF(L.GT.1) THEN
         L=L-1
         AT  =A(L)
      ELSE
         AT  =A(R)
         A(R)=A(1)
         R=R-1
         IF(R.EQ.1) THEN
            A(1)=AT
            RETURN
         END IF
      END IF
      I=L
      J=L+L
   20 IF(J.LE.R) THEN
         IF(J.LT.R) THEN
            IF(A(J).LT.A(J+1)) J=J+1
         END IF
         IF(AT.LT.A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
         ELSE
            J=R+1
         END IF
         GOTO 20
      END IF
      A(I)=AT
      GOTO 10
      END

      SUBROUTINE SORT2K(ITEM,N)                  ! sort 2-dim key-array
*     integer sort using Shells method of 2-dim key-array
      INTEGER ITEM(2,N),ITEMS
      LOGICAL EXCHAN
      SAVE
*     ...
      M=1
 10   M=1+3*M    ! determine initial increment
      IF(M.LE.N) GOTO 10
 20   M=M/3
      DO J=1,N-M ! sort with increment M
       L=J
 30    EXCHAN=.FALSE.
       IF(ITEM(1,L+M).LT.ITEM(1,L)) THEN
          EXCHAN=.TRUE.
       ELSE IF(ITEM(1,L+M).EQ.ITEM(1,L)) THEN
          IF(ITEM(2,L+M).LT.ITEM(2,L)) EXCHAN=.TRUE.
       END IF
       IF(EXCHAN) THEN
          ITEMS     =ITEM(1,L+M)       ! exchange first element ...
          ITEM(1,L+M)=ITEM(1,L)
          ITEM(1,L)  =ITEMS
          ITEMS     =ITEM(2,L+M)       ! ... and second element
          ITEM(2,L+M)=ITEM(2,L)
          ITEM(2,L)  =ITEMS
          L=L-M
          IF(L.GT.0) GOTO 30
       END IF
      END DO
      IF(M.GT.1) GOTO 20
      END


      FUNCTION CHINDL(N,ND)
*     return limit in chi^2/ND for N sigmas (N=1, 2 or 3)
      REAL SN(3),TABLE(30,3)
*     REAL PN(3) 
*     DATA PN/0.31731,0.0455002785,2.69985E-3/         ! probabilities
      DATA SN/0.47523,1.690140,2.782170/
      DATA TABLE/
     +  1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630,
     +  1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321,
     +  1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136,
     +  1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040,
     +  4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736,
     +  1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658,
     +  1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346,
     +  1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742,
     +  9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468,
     +  2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635,
     +  2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901,
     +  1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681/
*     ...
      IF(ND.LT.1) THEN
         CHINDL=0.0
      ELSE
         M=MAX(1,MIN(N,3))         ! 1, 2 or 3 sigmas
         IF(ND.LE.30) THEN
            CHINDL=TABLE(ND,M)     ! from table
         ELSE                      ! approximation for ND > 30
            CHINDL=(SN(M)+SQRT(FLOAT(ND+ND-1)))**2/FLOAT(ND+ND)
         END IF
      END IF
      END


 



 

