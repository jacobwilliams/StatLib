      SUBROUTINE DWPVAL (X, IAX, N, M, NCONST, LAG, IALT, E, C, P,
     *                   WS1, WS2, WS3, IFAULT)
C
C     ALGORITHM AS279 APPL. STATIST. (1993) VOL. 42, NO.1
C
C     P-VALUE FOR GENERALIZED/ALTERNATIVE DURBIN-WATSON STATISTIC WITH
C     ARBITRARY LAG BY CHOLESKY TRANSFORMATION METHOD
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N, M, IAX, NCONST, LAG, IFAULT, IALT
      REAL X(IAX,*), E, C, P, WS1(*), WS2(*), WS3(*)
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER J, K, MTOT, MNET, ITAIL
      REAL POINT5, POINT1, ZERO, ONE, E1, E2, E3, F1, F2, U, UJ,
     *     PI, RBD, UMAX, DELTA, ZFLOAT, ZLOG, ZIMAG, ZABS
      COMPLEX CF
C
C       CONSTANTS
C
      DATA ZERO, ONE, POINT5, POINT1 / 0.E0, 1.E0, 5.E-1, 1.E-1 /
      DATA PI /3.14159265358979E0 /
C
C       STANDARD FUNCTIONS
C
      ZFLOAT(J) = FLOAT(J)
      ZLOG(U) = ALOG(U)
      ZIMAG(CF) = AIMAG(CF)
      ZABS(CF) = CABS(CF)
C
C       INITIALIZATIONS
C
      IFAULT = 0
      MTOT = M + NCONST
      IF (N .LE. MTOT) THEN
        IFAULT = -1
        RETURN
      ELSEIF (N .LT. 2 * LAG) THEN
        IFAULT = -3
        RETURN
      ENDIF
      E1 = E * POINT1
      E3 = (E - E1) * POINT1
      E2 = E - E1 - E3
      E3 = E3 * PI
C
C       INTERVAL SIZE AND TRUNCATION LIMIT FOR TRAPEZOIDAL RULE
C
      MNET = MTOT - NCONST
      CALL GKFTRP (N, MNET, NCONST, LAG, IALT, C, E1, E2, DELTA, K,
     *             WS1, WS1(N+1), ITAIL, IFAULT)
      IF (IFAULT .NE. 0) RETURN
      IF (ITAIL .GT. 0) THEN
        P = ONE
        RETURN
      ELSEIF (ITAIL .LT. 0) THEN
        P = ZERO
        RETURN
      ENDIF
C
C       SCALE FACTOR FOR CHARACTERISTIC FUNCTION
C
      CALL DWSCL (X, IAX, N, M, NCONST, F1, F2, WS2, IFAULT)
      IF (IFAULT .NE. 0) RETURN
C
C       NUMERICAL INTEGRATION
C
      UMAX = ZFLOAT(K) + POINT5
      RBD = ONE
      P = ZERO
      J = 0
   50 IF (RBD .GT. E3 .AND. J .LE. K) THEN
        UJ = ZFLOAT(J) + POINT5
        U = UJ * DELTA
        CALL DWCHF(X, IAX, N, M, NCONST, LAG, IALT, U, C, F1, F2,
     *             CF, WS2, WS3, IFAULT)
        IF (IFAULT .NE. 0) RETURN
        P = P + ZIMAG(CF) / UJ
        RBD = ZABS(CF) * ZLOG(UMAX / UJ)
        J = J + 1
        GOTO 50
      ENDIF
      P = POINT5 - P / PI
      RETURN
      END
C
C
      SUBROUTINE DWCHF (X, IAX, N, M, NCONST, LAG, IALT, T, C, F1,
     *                  F2, CF, S, XT, IFAULT)
C
C     CHARACTERISTIC FUNCTION FOR GENERALIZED/ALTERNATIVE DURBIN-
C     WATSON STATISTIC AND EXACT LBI VERSIONS FOR ARBITRARY LAG
C
C       ARGUMENT DECLARATIONS:
C
      INTEGER IAX, IALT, N, M, NCONST, IFAULT, LAG
      REAL X(IAX,*), T, C, F1, F2
      COMPLEX CF, S(*), XT(LAG+1,*)
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER MTOT, I, J, K, MLAG, II, JJ, KK, NN, N1, IROW, JROW
      REAL ZERO, ONE, TWO, D2
      COMPLEX A1, A2, A, B, LB, ZCMPLX, ZSQRT, D1, LD, LLD
C
C       CONSTANTS
C
      DATA ZERO,ONE,TWO / 0.E0,1.E0,2.E0 /
C
C       STANDARD FUNCTIONS
C
      ZCMPLX(F1,F2) = CMPLX(F1,F2)
      ZSQRT(A) = CSQRT(A)
C
C       INITIALIZATIONS
C
      IFAULT = 0
      CF = F1
      MTOT = M + NCONST
      K = 0
      DO 20 I = 1,MTOT
      DO 10 J = 1,I
      K = K + 1
      S(K) = ZERO
 10   CONTINUE
 20   CONTINUE
      A1 = ZCMPLX(ONE,TWO * T * (C - ONE))
      A2 = A1 - ZCMPLX(ZERO, TWO * T)
      IF (IALT .NE. 0) A1 = A2
      B = ZCMPLX(ZERO, TWO * T)
      NN = (N - 1)/LAG + 1
      N1 = N - LAG +1
      MLAG = LAG + 1
      IROW = 0
      KK = 0
C
C       CHOLESKY TRANSFORMATION CYCLE
C
      DO 110 II = 1,NN
      DO 100 JJ = 1,LAG
      KK = KK + 1
      IF (KK .LE. N) THEN
        IROW = IROW + 1
        IF (IROW .GT. MLAG) IROW = 1
        IF (II .EQ. 1) THEN
          IF (JJ .EQ. 1) LD = ONE / ZSQRT(A1)
          DO 30 I = 1,M
          XT(IROW,I) = X(KK,I) * LD
 30       CONTINUE
          IF (NCONST .NE. 0) XT(IROW,MTOT) = LD
        ELSE
          IF (JJ .EQ. 1 .OR. KK .EQ. N1) THEN
            IF (JJ .NE. 1 .OR. KK .GT. N1) THEN
              LB = B * LLD
            ELSE
              LB = B * LD
            ENDIF
            LLD = LD
            A = A2
            IF (KK .GE. N1) A = A1
            LD = ONE / ZSQRT(A - LB * LB)
          ENDIF
	  JROW = IROW - LAG
	  IF (JROW .LE. 0) JROW = JROW + MLAG
          DO 40 I = 1,M
          XT(IROW,I) = LD * (X(KK,I) - LB * XT(JROW,I))
 40       CONTINUE
          IF (NCONST .NE. 0) XT(IROW,MTOT)=LD*(ONE-LB*XT(JROW,MTOT))
        ENDIF
        CF = CF * LD
C
C       UPDATE CROSS-PRODUCT MATRIX
C
        IF (MTOT .GT. 0) THEN
       	  K = 0
	  DO 60 I = 1,MTOT
	  DO 50 J = 1,I
	  K = K + 1
	  S(K) = S(K) + XT(IROW,I) * XT(IROW,J)
 50       CONTINUE
 60       CONTINUE
        ENDIF
      ENDIF
 100  CONTINUE
 110  CONTINUE
C
C       DETERMINANT OF CROSS-PRODUCT MATRIX
C
      IF (MTOT .GT. 0) THEN
        CALL CCLSKY(S,MTOT,D1,D2,IFAULT)
        CF = CF * D1 * TWO ** (D2 + F2)
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE CCLSKY(S,M,D1,D2,IFAULT)
C
C     CHOLESKY FACTORIZATION IN COMPLEX ARITHMETIC
C
C       ARGUMENT DECLARATIONS:
C
      INTEGER M, IFAULT
      REAL D2
      COMPLEX S(*), D1
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER I, J, JJ, K, L, KK, KI, KJ
      REAL EPS, ZERO, ONE, FOUR, SIXTEN, SIXNTH, ZABS
      COMPLEX SUM, ZSQRT
C
C       CONSTANTS
C
      DATA ZERO,ONE,FOUR,SIXTEN,SIXNTH / 0.E0,1.E0,4.E0,1.6E1,.625E-1 /
      DATA EPS / 1.E-6 /
C
C       STANDARD FUNCTIONS
C
       ZSQRT(SUM) = CSQRT(SUM)
       ZABS(SUM) = CABS(SUM)
C
C       INITIALIZATIONS
C
      IFAULT = 0
      D1 = ONE
      D2 = ZERO
C
C       FACTORIZATION
C
      K = 0
      DO 60 I = 1,M
      KK = 0
      KI = K
      DO 50 J = 1,I
      K = K + 1
      SUM = S(K)
      IF (J .GT. 1) THEN
        JJ = J - 1
        KJ = KI
        DO 10 L = 1,JJ
        KK = KK + 1
        KJ = KJ + 1
        SUM = SUM - S(KK) * S(KJ)
 10     CONTINUE
      ENDIF
      KK = KK + 1
      IF (J .LT. I) THEN
        S(K) = SUM * S(KK)
      ELSE
	IF (ZABS(SUM) .LT. ZABS(S(K)) * EPS) THEN
	  IFAULT = -2
          D1 = ZERO
	  RETURN
        ENDIF
	S(K) = ONE / ZSQRT(SUM)
        D1 = D1 * S(K)
 20     IF (ZABS(D1) .LE. ONE) GOTO 30
        D1 = D1 * SIXNTH
        D2 = D2 + FOUR
        GOTO 20
 30     IF (ZABS(D1) .GT. SIXNTH) GOTO 40
        D1 = D1 * SIXTEN
        D2 = D2 - FOUR
        GOTO 30
 40     CONTINUE
      ENDIF
 50   CONTINUE
 60   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE DWSCL (X, IAX, N, M, NCONST, F1, F2, S, IFAULT)
C
C     SCALE FACTOR FOR CHARACTERISTIC FUNCTION
C
C       ARGUMENT DECLARATIONS
C
      INTEGER IAX, N, M, NCONST, IFAULT
      REAL X(IAX,*), F1, F2, S(2,*)
C
C       LOCAL VARIABLE DECLARATIONS:
C
      INTEGER I, J, K, L, MTOT
      REAL EPS, TEMP, ONE, ZERO, ZREAL
      COMPLEX DET
C
C       CONSTANTS
C
      DATA ONE,ZERO / 1.E0,0.E0 /
C
C       STANDARD FUNCTIONS
C
      ZREAL(DET) = REAL(DET)
C
C       INITIALIZATIONS
C
      IFAULT = 0
      K = 0
      MTOT = M + NCONST
      IF (MTOT .EQ. 0) THEN
        F1 = ONE
        F2 = ZERO
        RETURN
      ENDIF
      DO 20 I = 1,MTOT
      DO 10 J = 1,I
      K = K + 1
      S(1,K) = ZERO
      S(2,K) = ZERO
 10   CONTINUE
 20   CONTINUE
C
C       CROSS-PRODUCT MATRIX
C
      DO 50 L = 1,N
      K = 0
      DO 40 I = 1,MTOT
      TEMP = ONE
      IF (I .LE. M) TEMP = X(L,I)
      DO 30 J = 1,I
      K = K + 1
      IF (J .LE. M) THEN
	S(1,K) = S(1,K) + X(L,J) * TEMP
      ELSE
	S(1,K) = S(1,K) + TEMP
      ENDIF
  30  CONTINUE
  40  CONTINUE
  50  CONTINUE
C
C       SCALE FACTOR: DETERMINANT OF CHOLESKY FACTOR
C
      CALL CCLSKY(S,MTOT,DET,F2,IFAULT)
      IF (IFAULT .NE. 0) RETURN
      F1 = ONE / ZREAL(DET)
      F2 = - F2
      RETURN
      END
C
C
      SUBROUTINE GKFTRP (N, M, NCONST, LAG, IALT, C, E1, E2, DELTA, K,
     *                   WS1, WS2, ITAIL, IFAULT)
C
C      INTERVAL SIZE AND TRUNCATION LIMIT FOR TRAPEZOIDAL RULE
C      EVALUATION OF GIL-PALAEZ FORMULA.
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N, M, NCONST, IALT, K, LAG, ITAIL, IFAULT
      REAL C, E1, E2, DELTA, WS1(*), WS2(*)
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER NNET, I, J, MAX, NG, K1
      REAL D2, ZERO, TWO, E, ZMIN1
C
C       CONSTANTS
C
      DATA MAX / 20 /
      DATA ZERO, TWO / 0.E0, 2.E0 /
C
C       STANDARD FUNCTIONS
C
      ZMIN1(E1,E2) = AMIN1(E1,E2)
C
C       INITIALIZATIONS
C
      IFAULT = 0
      IF (N .LE. M + NCONST) THEN
        IFAULT = -1
        RETURN
      ENDIF
      IF (IALT .EQ. 0) THEN
        NNET = N - M - NCONST
        K1 = NCONST
      ELSE
        NNET = N - M
        K1 = 0
      ENDIF
C
C       EIGENVALUES OF MATRIX DEFINING STATISTIC
C
      CALL EVALUE (WS1, N, LAG, IALT, C)
C
C       INITIAL CHECK FOR TAILS
C
      IF (WS1(K1 + 1) .GE. 0) THEN
        ITAIL = -1
        RETURN
      ELSEIF (WS1(N) .LE. 0) THEN
        ITAIL = 1
        RETURN
      ENDIF
C
C       INTERVAL SIZE
C
      E = E1 + E2
      IF (C .GT. TWO) THEN
        ITAIL = 1
        J = M + K1
      ELSE
        ITAIL = -1
        J = K1
      ENDIF
      DO 10 I = 1, NNET
      J = J + 1
      WS2(I) = WS1(J)
   10 CONTINUE
      J = ITAIL
      CALL INTBD (WS2, NNET, E, E1, DELTA, ITAIL)
      IF (ITAIL .NE. 0) RETURN
      ITAIL = -J
      IF (ITAIL .EQ. 1) THEN
        J = M + K1
      ELSE
        J = K1
      ENDIF
      DO 20 I = 1, NNET
      J = J + 1
      WS2(I) = WS1(J)
   20 CONTINUE
      CALL INTBD (WS2, NNET, E, E1, D2, ITAIL)
      IF (ITAIL .NE. 0) RETURN
      DELTA = ZMIN1 (DELTA, D2)
C
C       TRUNCATION POINT FOR TRAPEZOIDAL RULE
C
      CALL TBDWT (WS1, N, M, NCONST, IALT, WS2, NG)
      IF (NG .EQ. 0) THEN
        IFAULT = 2
        RETURN
      ENDIF
      CALL TBND1 (WS2, NG, E2, DELTA, K1, IFAULT)
      IF (IFAULT .NE. 0) K1 = 0
      IFAULT = 0
      CALL TBND2 (WS2, NG, E2, DELTA, K)
      IF (K1 .GT. 0) K = MIN0(K,K1)
      RETURN
      END
C
C
      SUBROUTINE TBDWT (W, N, M, NCONST, IALT, G, NG)
C
C     WEIGHTS (SQUARED) FOR BOUNDING FUNCTION TO WHICH DAVIES' (APPL.
C     STATIST. 29, 1980, 323-333) TRUNCATION RULES ARE TO BE APPLIED.
C
C       ARGUMENT DECLARATIONS:
C
      INTEGER N, M, NCONST, IALT, NG
      REAL W(N), G(N)
C
C       LOCAL VARIABLE DECLARATIONS:
C
      INTEGER I, MTOT, NNET, L1, L2
      REAL ZERO
C
C       CONSTANTS
C
      DATA ZERO / 0.E0 /
C
C       NUMBER OF REGRESSORS, DEGREES OF FREEDOM
C
      MTOT = M + NCONST
      NNET = N - MTOT
C
C       WEIGHTS
C
      IF (IALT .EQ. 0 .AND. NCONST .NE. 0) THEN
        L1 = 1
      ELSE
        L1 = 0
      ENDIF
      L2 = MTOT
      NG = 0
      DO 30 I = 1, NNET
      L1 = L1 + 1
      L2 = L2 + 1
      IF (W(L1) .GT. ZERO) THEN
        NG = NG + 1
        G(NG) = W(L1) ** 2
      ELSEIF (W(L2) .LT. ZERO) THEN
        NG = NG + 1
        G(NG) = W(L2) ** 2
      ENDIF
   30 CONTINUE
      IF (NG .EQ. N) RETURN
      NNET = NG + 1
      DO 40 I = NNET, N
      G(I) = ZERO
   40 CONTINUE
      RETURN
      END
C                                                                       
C                                                                       
      SUBROUTINE EVALUE (W, N, LAG, IALT, C)
C
C     EIGENVALUES OF DEFINING QUADRATIC FORM                       
C
C       ARGUMENT DECLARATIONS:
C
      INTEGER N, LAG, IALT
      REAL W(N), C
C                                                                       
C       LOCAL VARIABLE DECLARATIONS:
C
      INTEGER I, J, K, L, NVAL, M1, M2                                  
      REAL PI, D1, D2, TWO, C2, TEMP, ZFLOAT, ZCOS
C
C       STANDARD FUNCTIONS
C
      ZFLOAT(I) = FLOAT(I)
      ZCOS(D1) = COS(D1)
C                                                                       
C       CONSTANTS                                                         
C                                                                       
      DATA PI / 3.14159265358979E0 /
      DATA TWO / 2.E0 /
      C2 = TWO - C
C                                                                       
C     FIRST-ORDER CASE                                                  
C                                                                       
      IF (LAG .GT. 1) GOTO 20                                           
      L = IALT                                                          
      IF (IALT .EQ. 0) D1 = PI / ZFLOAT(N)                              
      IF (IALT .NE. 0) D1 = PI / ZFLOAT(N + 1)                          
      DO 10 I = 1, N                                                    
      W(I) = C2 - TWO * ZCOS(L * D1)                                    
      L = L + 1                                                         
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C     GENERAL-ORDER CASE                                                
C                                                                       
   20 NVAL = (N - 1) / LAG + 1                                          
      IF (IALT .EQ. 0) D1 = PI / ZFLOAT(NVAL)                           
      IF (IALT .NE. 0) D1 = PI / ZFLOAT(NVAL + 1)                       
      M1 = N - (NVAL - 1) * LAG                                         
      M2 = M1 + 1
      IF (M1 .EQ. LAG) GOTO 30                                          
      IF (IALT .EQ. 0) D2 = PI / ZFLOAT(NVAL - 1)                       
      IF (IALT .NE. 0) D2 = PI / ZFLOAT(NVAL)                           
   30 K = 0                                                             
      L = IALT - 1                                                      
      DO 60 I = 1, NVAL                                                 
      L = L + 1                                                         
      TEMP = C2 - TWO * ZCOS(L * D1)                                    
      DO 40 J = 1, M1                                                   
      K = K + 1                                                         
      W(K) = TEMP                                                       
   40 CONTINUE                                                          
      IF (M1 .EQ. LAG .OR. K .EQ. N)GOTO 60                             
      TEMP = C2 - TWO * ZCOS(L * D2)                                    
      DO 50 J = M2, LAG                                                 
      K = K + 1                                                         
      W(K) = TEMP                                                       
   50 CONTINUE                                                          
   60 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE LMGFEV (U, G, N, PSI, DPSI, D2PSI)
C
C     LOG MOMENT GENERATING FUNCTION FROM EIGENVALUES ADJUSTED FOR
C     ARGUMENT OF DISTRIBUTION FUNCTION, AND DERIVATIVES
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N
      REAL U, G(N), PSI, DPSI, D2PSI
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER I
      REAL U2, F, ZERO, ONE, TWO, ZLOG
C
C       CONSTANTS
C
      DATA ZERO, ONE, TWO / 0.E0, 1.E0, 2.E0/
C
C       STANDARD FUNCTIONS
C
      ZLOG(F) = ALOG(F)
C
C       FUNCTION AND DERIVATIVES
C
      U2 = - TWO * U
      PSI = ZERO
      DPSI = ZERO
      D2PSI = ZERO
      DO 10 I = 1, N
      F = ONE + G(I) * U2
      PSI = PSI + ZLOG(F)
      F = G(I) / F
      DPSI = DPSI + F
      D2PSI = D2PSI + F ** 2
   10 CONTINUE
      PSI = - PSI / TWO
      D2PSI = TWO * D2PSI
      RETURN
      END
C
C
      SUBROUTINE INTBD (G, N, E, ET, DELTA, ITAIL)
C
C     INTERVAL FOR TRAPEZOIDAL RULE INTEGRATION TO LIMIT ERROR TO E
C     USING MOMENT GENERATING FUNCTION BOUND FOR WEIGHTED SUM OF CHI-
C     SQUARED(1) RANDOM VARIABLES
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N, ITAIL
      REAL DELTA, G(N), E, ET
C
C       LOCAL VARIABLE DECLARATIONS:
C
      INTEGER I
      REAL ZERO, TWO, POINT1, POINT5, TWOPI, ULIM, PSI, DPSI, D2PSI,
     *     HB, U, U1, U2, ELOG, EPS, HM, HSD, SGN, ZABS, ZEXP, ZFLOAT,
     *     ZLOG, ZMAX1, ZMIN1, ZSQRT
C
C       SET CONSTANTS (EPS = LOG(1.1))
C
      DATA ZERO, TWO, POINT1, POINT5 / 0.E0, 2.E0, .1E0, .5E0 /
      DATA TWOPI / 6.28318530717959E0 /
      DATA EPS / .09531E0 /
C
C       STANDARD FUNCTIONS
C
      ZABS(U) = ABS(U)
      ZEXP(U) = EXP(U)
      ZFLOAT(I) = FLOAT(I)
      ZLOG(U) = ALOG(U)
      ZMAX1(U,U1) = AMAX1(U,U1)
      ZMIN1(U,U1) = AMIN1(U,U1)
      ZSQRT(U) = SQRT(U)
C
C       TAIL CHECK
C
      SGN = ZFLOAT(ITAIL)
      IF (ITAIL .EQ. 1) THEN
        IF (G(N) .LE. ZERO) RETURN
        ULIM = POINT5 / G(N)
      ELSE
        IF (G(1) .GE. ZERO) RETURN
        ULIM = POINT5 / G(1)
      ENDIF
      HM = ZERO
      HSD = ZERO
      DO 10 I = 1, N
      HM = HM + G(I)
      HSD = HSD + G(I) ** 2
   10 CONTINUE
C
C       NEWTON SEARCH FOR ZERO LOG(MGF) DERIVATIVE
C
      HSD = ZSQRT(HSD)
      U2 = ZERO
      IF ( SGN * HM .LT. ZERO .AND. - SGN * HM / HSD .GT. TWO) THEN
        DPSI = HM
   20   IF (DPSI * SGN .LT. ZERO) THEN
          U1 = U2
          U2 = POINT5 * (U2 + ULIM)
          CALL LMGFEV (U2, G, N, PSI, DPSI, D2PSI)
          GOTO 20
        ENDIF
   25   IF (DPSI ** 2 .GE. POINT1 * D2PSI) THEN
          U2 = U2 - DPSI / D2PSI
          CALL LMGFEV (U2, G, N, PSI, DPSI, D2PSI)
        ENDIF
        IF (ZEXP(PSI) .LE. E) RETURN
      ENDIF
C
C       BOUND EQUALITY
C
      ELOG = - ZLOG(ET)
      HB = POINT5
   30 IF (HB .GT. ZERO) THEN
        U1 = U2
        U2 = POINT5 * (U2 + ULIM)
        CALL LMGFEV (U2, G, N, PSI, DPSI, D2PSI)
        HB = PSI - U2 * DPSI + ELOG
        GOTO 30
      ENDIF
C
C       NEWTON/BINARY SEARCH SOLUTION
C
      U = U2
   40 IF (ZABS(HB) .GE. EPS) THEN
        U = U - HB / (U * D2PSI)
        IF (U .LE. ZMIN1(U1,U2) .OR. U .GT. ZMAX1(U1,U2))
     *      U = POINT5 * (U1 + U2)
        CALL LMGFEV (U, G, N, PSI, DPSI, D2PSI)
        HB = PSI - U * DPSI + ELOG
        IF (HB .GT. ZERO) THEN
          U1 = U
        ELSE
          U2 = U
        ENDIF
        GOTO 40
      ENDIF
      IF (ITAIL .EQ. 1 .AND. DPSI .GT. ZERO) THEN
        DELTA = TWOPI / DPSI
        ITAIL = 0
      ELSEIF (ITAIL .EQ. -1 .AND. DPSI .LT. ZERO) THEN
        DELTA = - TWOPI / DPSI
        ITAIL = 0
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE TBND1 (G, N, E, DELTA, K, IFAULT)
C
C     UPPER BOUND ON NUMBER OF TERMS REQUIRED IN TRAPEZOIDAL RULE
C     INTEGRATION USING DAVIES' BOUND EQN (6), APPL. STATIST. 29 (1980)
C     P. 324.  WEIGHTS G(J) ARE SQUARES OF DAVIES' COEFFICIENTS.
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N, K, IFAULT
      REAL G(N), E, DELTA
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER I, M
      REAL PIBY2, ZERO, ONE, TWO, FOUR, S, F, DF, C, U, ULAST,
     *     FLAST, DFLAST, EPS, POINT5, ZABS, ZFLOAT, ZLOG, ZSQRT
C
C       CONSTANTS (EPS = 4 * LOG(1.1))
C
      DATA PIBY2 / 1.57079632679490E0 /
      DATA ZERO, ONE, TWO, FOUR, POINT5
     *    / 0.E0, 1.E0, 2.E0, 4.E0, 5.E-1 /
      DATA EPS / .38124E0 /
C
C       STANDARD FUNCTIONS
C
      ZABS(C) = ABS(C)
      ZFLOAT(I) = FLOAT(I)
      ZLOG(C) = ALOG(C)
      ZSQRT(C) = SQRT(C)
C
C       INITIALIZATIONS
C
      IFAULT = 0
      M = 0
      C = ZERO
      DO 10 I = 1, N
      IF (G(I) .GT. ONE) THEN
        M = M + 1
        C = C + ZLOG(G(I))
      ENDIF
   10 CONTINUE
      IF (M .EQ. 0) THEN
        IFAULT = 1
        RETURN
      ENDIF
      S = ZFLOAT(M)
      C = C + FOUR * ZLOG(S * PIBY2 * E)
      I = 0
      U = FOUR * FOUR * DELTA ** 2
C
C       INITIAL ESTIMATE OF BOUND FUNCTION SOLUTION
C
   20 CALL FBND1 (G, N, C, S, U, F, DF)
      IF (F .GT. ZERO .AND. I .EQ. -1) THEN
        U = ULAST
        F = FLAST
        DF = DFLAST
      ELSEIF (F .GT. ZERO) THEN
        I = 1
        U = U * POINT5
        GOTO 20
      ELSEIF (F .LE. ZERO .AND. I .NE. 1) THEN
        I = -1
        ULAST = U
        FLAST = F
        DFLAST = DF
        U = TWO * U
        GOTO 20
      ENDIF
C
C       NEWTON'S METHOD SOLUTION FOR ARGUMENT OF BOUND FUNCTION
C
   30 IF (ZABS(F) .GE. EPS) THEN
        U = U - F / DF
        CALL FBND1 (G, N, C, S, U, F, DF)
        GOTO 30
      ENDIF
      K = ZSQRT(U / FOUR) / DELTA + POINT5
      RETURN
      END
C
C
      SUBROUTINE FBND1 (G, N, C, S, U, F, DF)
C
C     EVALUATION OF DAVIES' BOUND EQN (6), APPL. STATIST. 29 (1980)
C     P. 324.  WEIGHTS G(J) ARE SQUARES OF DAVIES' COEFFICIENTS, AND F
C     IS - 4 * (LOG(BOUND) - LOG(MAX ERROR)).
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N
      REAL G(N), C, S, U, F, DF
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER J
      REAL ONE, TEMP, ZLOG
C
C       CONSTANTS
C
      DATA ONE / 1.E0 /
C
C       STANDARD FUNCTIONS
C
      ZLOG(C) = ALOG(C)
C
      F = C + S * ZLOG(U)
      DF = S / U
      DO 10 J = 1, N
      IF (G(J) .LE. ONE) THEN
        TEMP = ONE + U * G(J)
        F = F + ZLOG(TEMP)
        DF = DF + G(J) / TEMP
      ENDIF
   10 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE TBND2 (G, N, E, DELTA, K)
C
C     UPPER BOUND ON NUMBER OF TERMS REQUIRED IN TRAPEZOIDAL RULE
C     INTEGRATION USING DAVIES' BOUND EQN (8), APPL. STATIST. 29 (1980)
C     P. 324.  WEIGHTS G(J) ARE SQUARES OF DAVIES' COEFFICIENTS.
C     LEADING CONSTANT IS 2/PI AS IN AKS (1990).
C
C       ARGUMENT DECLARATIONS:
C
      INTEGER N, K
      REAL G(N), E, DELTA
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER I, M
      REAL C1, ZERO, ONE, TWO, FOUR, F, DF, CONST, EPS,
     *     POINT5, U, ULAST, FLAST, DFLAST, ZABS, ZLOG, ZSQRT
C
C       CONSTANTS (C1 = - 4 * LOG(2/PI), EPS = 4 * LOG(1.1))
C
      DATA C1 / 1.80633E0 /
      DATA ZERO, ONE, TWO, FOUR, POINT5 /
     *     0.E0, 1.E0, 2.E0, 4.E0, 5.E-1 /
      DATA EPS / .38124E0 /
C
C       STANDARD FUNCTIONS
C
      ZABS(U) = ABS(U)
      ZLOG(U) = ALOG(U)
      ZSQRT(U) = SQRT(U)
C
C       INITIAL ESTIMATE OF BOUND FUNCTION SOLUTION
C
      CONST = C1 + FOUR * ZLOG(E)
      I = 0
      U = FOUR * FOUR * DELTA ** 2
   10 CALL FBND2 (G, N, CONST, U, F, DF)
      IF (F .GT. ZERO .AND. I .EQ. -1) THEN
        U = ULAST
        F = FLAST
        DF = DFLAST
      ELSEIF (F .GT. ZERO) THEN
        I = 1
        U = U * POINT5
        GOTO 10
      ELSEIF (F .LE. ZERO .AND. I .NE. 1) THEN
        I = -1
        ULAST = U
        FLAST = F
        DFLAST = DF
        U = TWO * U
        GOTO 10
      ENDIF
C
C       NEWTON'S METHOD SOLUTION FOR ARGUMENT OF BOUND FUNCTION
C
   20 IF (ZABS(F) .GE. EPS) THEN
        U = U - F / DF
        CALL FBND2 (G, N, CONST, U, F, DF)
        GOTO 20
      ENDIF
C
C       NUMBER OF INTERVALS
C
      K = ZSQRT(U / FOUR) / DELTA + POINT5
      RETURN
      END
C
C
      SUBROUTINE FBND2 (G, N, C, U, F, DF)
C
C     EVALUATION OF DAVIES' BOUND EQN (8), APPL. STATIST. 29 (1980)
C     P. 324, WITH CONSTANT 2/PI AS IN AKS (1990).  WEIGHTS G(J) ARE
C     SQUARES OF DAVIES' COEFFICIENTS, AND F IS -4 * (LOG(BOUND) 
C     - LOG(MAX ERROR)).
C
C       ARGUMENT DECLARATIONS
C
      INTEGER N
      REAL G(N), C, U, F, DF
C
C       LOCAL VARIABLE DECLARATIONS
C
      INTEGER J
      REAL ZERO, ONE, TEMP, ZLOG
C
C       CONSTANTS
C
      DATA ZERO, ONE / 0.E0, 1.E0 /
C
C       STANDARD FUNCTIONS
C
      ZLOG(U) = ALOG(U)
C
C       FUNCTION AND DERIVATIVE
C
      F = C
      DF = ZERO
      DO 20 J = 1, N
      TEMP = ONE + G(J) * U
      F = F + ZLOG(TEMP)
      DF = DF + G(J) / TEMP
   20 CONTINUE
      RETURN
      END


