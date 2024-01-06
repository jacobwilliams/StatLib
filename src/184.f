      SUBROUTINE PROB(EPSI, N, NU, C, M, VM, ANS, BOUND, RND, IFAULT)
C
C<<<<<  Acquired in machine-readable form from 'Applied Statistics'
C<<<<<  algorithms editor, January 1983.
C
C        ALGORITHM AS 184  APPL. STATIST. (1982) VOL.31, NO.3
C
C        CALCULATES EXPRESSION (1) OF THE PURPOSE SECTION.
C        USEFUL FOR NON-CENTRAL STUDENTIZED MAXIMUM, MAXIMUM MODULUS,
C        AND MULTIPLE T-TEST PROBABILITIES.
C
      DIMENSION PROD(2), EK(3), PHI(2, 10), C(2, N), AC(2), AM(2), M(N),
     *  PH(200, 10), S(5), G2(3), G1(3), GMX(4), G2S(200, 3), AD(2, 10),
     *  ER(4), PHM(10), RS(10), VM(N), R(200), PR(200)
C
C        ALR2PI=ALOG(2.0*PI), AL2=ALOG(2.0),
C        ER CONTAINS ROUNDING ERROR BOUNDS
C
      DATA ALR2PI /0.91893853320467274178/,AL2 /0.69314718055994530942/,
     *  ER(1), ER(2), ER(3), ER(4) /5.0E-07, 5.0E-07, 5.0E-07, 5.0E-07/,
     *  ZERO, HALF, ONE, TWO, SEVEN, EIGHT /0.0, 0.5, 1.0, 2.0, 7.0,
     *  8.0/, EITH, SEVEIG, FOUR, TWEL, XNIN6 /0.125, 0.875, 4.0, 12.0,
     *  96.0/, MAXN /10/
C
C        FUNCTION STATEMENTS
C
      BDD1(X, Y, Z) = ABS((X - Z) / Y - X * Y)
      BDD2(W, X, Y, Z) = ABS((W - Y) * (W - Z) + (X * X) * ((-Z * W + Y)
     *  * W + W * W * X * X))
      GABS(X) = ABS(X)
      GEXP(X) = EXP(X)
      GMAX1(X, Y) = AMAX1(X, Y)
      TMAX1(X, Y, Z) = AMAX1(X, Y, Z)
      GMIN1(X, Y) = AMIN1(X, Y)
      GSQRT(X) = SQRT(X)
      GFLOAT(IJ) = FLOAT(IJ)
      GLOG(X) = ALOG(X)
C
C        CHECKING FOR FAULTY INPUT DATA
C
      IFAULT = 0
      IF (NU .LE. 0) IFAULT = 1
      IF (EPSI .LE. ER(3)) IFAULT = 2
      IF (N .LE. 0 .OR. N .GT. MAXN) IFAULT = 3
      IF (IFAULT .NE. 0) RETURN
      IZERO = 0
      DO 2 J = 1, N
      IF (C(1, J) .EQ. C(2, J) .AND. M(J) .EQ. 0) IZERO = 1
      IF (C(1, J) .GT. C(2, J) .AND. M(J) .EQ. 0) IFAULT = 100 + J
      IF (M(J) .EQ. 1 .AND. C(1, J) .NE. ZERO) IFAULT = 200 + J
      IF (M(J) .EQ. 2 .AND. C(2, J) .NE. ZERO) IFAULT = 300 + J
      IF (M(J) .EQ. 3 .AND. (C(1, J) .NE. ZERO .OR. C(2, J) .NE. ZERO))
     *  IFAULT = 400 + J
      IF (M(J) .LT. 0 .OR. M(J) .GT. 3) IFAULT = 500 + J
    2 CONTINUE
      IF (IFAULT .NE. 0) RETURN
C
C        IZERO=1 IF ONE OF THE RANGES (C1,C2) IS DEGENERATE
C
      ANS = ZERO
      BOUND = ZERO
      RND = ZERO
      IF (IZERO .EQ. 1) RETURN
C
C        INITIALIZING VALUES TO BE USED IN PROGRAM
C
      NK = N - 1
      ITCH = 2
      INDXL = 3
      DF = GFLOAT(NU)
      HD = HALF * DF
      COEFF = HD * GLOG(HD) - ALGAMA(HD) + AL2
      COFF = GEXP(COEFF)
      EG = ONE
      FE = ZERO
      K = 1
      TD = HALF / DF
      SDFT = GSQRT(EIGHT * DF - SEVEN) * TD
      S(1) = ZERO
      S(2) = GSQRT(ONE - TD - SDFT)
      S(3) = GSQRT((DF - ONE) / DF)
      S(4) = GSQRT(ONE - TD + SDFT)
      S(5) = GMAX1(TWO, GSQRT(ONE + GSQRT(GSQRT(XNIN6 * (FOUR + DF)
     *  / (EPSI * DF ** 3)))))
      EPSL = EPSI * SEVEIG / S(5)
      ES = S(5) - (S(4) - S(2)) * GMIN1(ONE, GFLOAT(NU / 3))
      EK(2) = ZERO
      EK(3) = ZERO
      PROD(1) = ONE
C
C        THE FIRST 4 LINES OF THE DO 6 LOOP DETERMINE WHERE THE FACTORS,
C        PHI( ,JK), OF PROD ARE MONOTONE.
C        THE LOOP COMPUTES PHI AND PROD FOR S=0.0
C
      DO 6 JK = 1, N
      RS(JK) = -ONE
      SUMC = C(1, JK) + C(2, JK)
      IF (C(1, JK) * C(2, JK) .GT. ZERO) RS(JK) = VM(JK) / SUMC +
     *  GSQRT(VM(JK) * VM(JK) + TWO * SUMC * GLOG(C(2, JK) / C(1, JK))
     *  / (C(2, JK) - C(1, JK))) / GABS(SUMC)
      IF (M(JK) .EQ. 0) PHI(1, JK) = ZERO
      IF (M(JK) .EQ. 1) PHI(1, JK) = ALNORM(-VM(JK), .FALSE.)
      IF (M(JK) .EQ. 2) PHI(1, JK) = ONE - ALNORM(-VM(JK), .FALSE.)
      IF (M(JK) .EQ. 3) PHI(1, JK) = ONE
      PROD(1) = PROD(1) * PHI(1, JK)
    6 CONTINUE
C
C        THESE 5 LINES AND THE FIRST 4 LINES OF THE DO 60 LOOP ARE
C        SPECIAL HANDLING FOR NU.LE.2
C
      IXG = MIN0(2, NU - 1)
      XH = GFLOAT((2 - IXG) / 2)
      G1(1) = XH * COFF
      G1(2) = GFLOAT(IXG - 2 * (IXG / 2)) * COFF
      G1(3) = BDD2(DF, ZERO, ONE, TWO) + GFLOAT(1 - (1 + IXG) / 2)
C
C        INTEGRATION BEGINS.  (FE,R(K)) IS CURRENT SUBINTERVAL
C
      DO 60 L = 2, 5
      IF ((NU .LE. 2 .AND. L .EQ. 2) .OR. (NU .EQ. 1 .AND. L .LE. 3))
     *  GOTO 60
      R(1) = S(L)
      ITCH = 2
      IF (NU .LE. 2 .OR. L .EQ. 2 .OR. L .EQ. 5) GOTO 13
      INDXL = 2
      ITCH = 1
      GOTO 13
C
C        SETTING UP THE MESH.   STATEMENTS 9-12 MAKE THE MESH FINER
C
    9 IF (K .GE. 200) GOTO 62
      R(K + 1) = (R(K) + FE) * HALF
      PR(K) = PROD(2)
      DO 11 J = 1, N
   11 PH(K, J) = PHI(2, J)
      DO 12 J = 1, INDXL
   12 G2S(K, J) = G2(J)
      K = K + 1
      ITCH = INDXL - 1
C
C        COMPUTING PHI, PROD, DENSITY OF S, AND ITS DERIVATIVE AT S=R(K)
C
   13 G2(1) = GEXP(COEFF + (DF - ONE) * GLOG(R(K)) - DF * R(K) * R(K)
     *  * HALF)
      G2(2) = G2(1) * BDD1(DF, R(K), ONE)
      PROD(2) = ONE
      DO 14 J = 1, N
      IF (M(J) .EQ. 0) PHI(2, J) = ALNORM(C(2, J) * R(K) - VM(J),
     *  .FALSE.) - ALNORM(C(1, J) * R(K) - VM(J), .FALSE.)
      IF (M(J) .EQ. 1) PHI(2, J) = ALNORM(C(2, J) * R(K) - VM(J),
     *  .FALSE.)
      IF (M(J) .EQ. 2) PHI(2, J) = ONE - ALNORM(C(1, J) * R(K) - VM(J),
     *  .FALSE.)
      IF (M(J) .EQ. 3) PHI(2, J) = ONE
      PROD(2) = PROD(2) * PHI(2, J)
   14 CONTINUE
C
C        STATEMENTS 18-22 DECIDE WHETHER TO USE FORM A (INDXL=3)
C        OR FORM B (INDXL=2).
C        ITCH DETERMINES WHETHER OR NOT G2(3) NEEDS TO BE COMPUTED
C
   18 EPS = R(K) * EPSL
      GMIN = GMIN1(G1(1), G2(1))
      IF (NU .LE. 2) GOTO 22
      IF (GMIN .LT. ONE) GOTO 20
      INDXL = 2
      GOTO 26
   20 IF (((INDXL .EQ. 3) .OR. (L .EQ. 3)) .OR. (L .EQ. 4)) GOTO 22
      G1(3) = BDD2(DF, FE, ONE, TWO) / (FE * FE)
      INDXL = 3
      ITCH = 2
   22 IF (ITCH .EQ. 2) G2(3) = BDD2(DF, R(K), ONE, TWO) / (R(K) * R(K))
C
C        THE DO 32 LOOP COMPUTES BOUNDS FOR FIRST (AD(1,J)) AND
C        SECOND (AD(2,J)) DERIVATIVES OF PHI( ,J)
C
   26 PRM = ONE
      DO 32 J = 1, N
      PHM(J) = GMAX1(PHI(1, J), PHI(2, J))
      IF (FE .LT. RS(J) .AND. R(K) .GT. RS(J)) PHM(J) = ALNORM(C(2, J)
     *  * RS(J) - VM(J), .FALSE.) - ALNORM(C(1, J) * RS(J) - VM(J),
     *  .FALSE.)
      PRM = PRM * PHM(J)
      DO 30 I = 1, 2
      AC1 = GABS(GABS(C(I, J) * FE) - GABS(VM(J)))
      AC2 = GABS(GABS(C(I, J) * R(K)) - GABS(VM(J)))
      AC(I) = GMIN1(AC1, AC2)
      IF (C(I, J) .EQ. ZERO) GOTO 28
      RATC = VM(J) / C(I, J)
      IF (FE .LT. RATC .AND. R(K) .GT. RATC) AC(I) = ZERO
   28 AM(I) = GMAX1(AC1, AC2)
      AC(I) = GABS(C(I, J)) * GEXP(-HALF * AC(I) * AC(I) - ALR2PI)
   30 CONTINUE
      AD(1, J) = AC(1) + AC(2)
      AD(2, J) = GABS(C(2, J)) * AM(2) * AC(2) + GABS(C(1, J)) * AM(1)
     *  * AC(1)
   32 CONTINUE
C
C        THE DO 33 AND DO 34 LOOPS BOUND FIRST (BOU1) AND
C        SECOND (BOU2) DERIVATIVES OF PROD
C
      BOU1 = ZERO
      BOU2 = ZERO
      DO 33 J = 1, N
      BOU2 = BOU2 + PRM / PHM(J) * AD(2, J)
   33 BOU1 = BOU1 + PRM / PHM(J) * AD(1, J)
      IF (N .EQ. 1) GOTO 36
      DO 34 J = 1, NK
      JJ = J + 1
      DO 34 KM = JJ, N
      BOU2 = BOU2 + TWO * PRM / (PHM(J) * PHM(KM)) * AD(1, J) *
     *  AD(1, KM)
   34 CONTINUE
   36 IF (INDXL .EQ. 2) GOTO 40
C
C        CHECKING THE BOUND FOR THE FORM A CASES
C
      DEL = R(K) - FE
      DO 38 J = 1, 3
   38 GMX(J) = GMAX1(G1(J), G2(J))
      GMX(4) = GMX(1) * GMX(3)
      IF (NU .EQ. 3 .AND. FE .EQ. ZERO) GMX(4) = GMAX1(GMX(4), TWO *
     *  COFF)
      CHECK = (DEL ** 3) * (PRM * GMX(4) + BOU2 * GMX(1) + BOU1 * TWO *
     *  GMX(2)) / TWEL
      GA = G1(1)
      GB = G2(1)
      GOTO 43
C
C        CHECKING THE BOUND FOR THE FORM B CASES
C
   40 DEL = GAMAIN(DF, R(K) * R(K) * DF) - GAMAIN(DF, FE * FE * DF)
      BOU3 = GMAX1(BDD1(DF, FE, ONE), BDD1(DF, R(K), ONE))
      CHECK = (DEL ** 3) * (BOU2 + BOU1 * BOU3) / (GMIN * GMIN * TWEL)
      GA = ONE
      GB = ONE
   43 BOU4 = BOUND + CHECK
      IF (BOU4 .GT. EPS) GOTO 9
      ANS = ANS + HALF * (PROD(1) * GA + PROD(2) * GB) * DEL
      EK(INDXL) = EK(INDXL) + ONE
C
C        ADVANCING THE INTEGRATION IF THE BOUND IS SMALL ENOUGH
C
      BOUND = BOU4
      IF (INDXL .EQ. 3) EG = TMAX1(EG, G1(1), G2(1))
      FE = R(K)
      PROD(1) = PROD(2)
      DO 46 J = 1, INDXL
   46 G1(J) = G2(J)
      DO 50 J = 1, N
   50 PHI(1, J) = PHI(2, J)
      IF (K .LE. 1) GOTO 60
      K = K - 1
      ITCH = 1
      DO 54 J = 1, N
   54 PHI(2, J) = PH(K, J)
      PROD(2) = PR(K)
      DO 56 J = 1, INDXL
   56 G2(J) = G2S(K, J)
      GOTO 18
   60 CONTINUE
      BOUND = BOUND + EPSI * EITH
C
C        COMPUTING THE BOUND ON ROUNDING ERROR
C
      AN = GFLOAT(N)
      EA = GABS(GEXP(ER(4)) - ONE)
      RND = (EG * TWO * AN * ER(1) + EG * EA / (ONE + EA) + ER(3))
     *  * (ES + EK(3) * ER(3)) + ER(3) * EK(3) * EG + TWO * AN * ER(1)
     *  + (ONE + TWO * AN * ER(1)) * EK(2) * ER(2)
      RETURN
   62 IFAULT = 4
      RETURN
      END
