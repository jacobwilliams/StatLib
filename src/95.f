      SUBROUTINE CURVE (P, X, N0, N, EPS, MAXITR, MU, SIGMA, ITER,
     1   SEMU, SESIG, COVAR, E0, EX, CHISQ, IFAULT)
C 
C       ALGORITHM AS 95 APPL. STATIST. (1976) VOL.25, NO.1
C 
C       ESTIMATES MU AND SIGMA OF DISTRIBUTION FUNCTION
C       F( (X-MU)/SIGMA ) FROM A GROUPED SAMPLE OF X VALUES.
C       NOTE ON ARRAY SIZES
C       THE ARRAYS IN THE SECOND DIMENSION STATEMENT MUST HAVE
C       MINIMUM SIZE P.  IF P IS TO EXCEED 20, A SUITABLE SIZE
C       MUST BE SET FOR THEM, AND THE IF STATEMENT WHICH CHECKS
C       THE VALUE OF P MUST BE AMENDED.
C
C     Auxiliary routines required: FUNC & DEVIAT (both user-supplied)
C
      INTEGER P
      REAL NN, NI, NP, MU, ONE, ZERO
      DIMENSION X(P), N(P), EX(P)
      DIMENSION F(20), F1(20), XN(20)
      DATA RR/1.0E-10/, ONE/1.0/, ZERO/0.0/
C 
C       ERROR EXIT IF P TOO SMALL OR TOO LARGE
C 
      IF (P.LT.2 .OR. P.GT.20) GO TO 80
      IFAULT = 0
C 
C       SET FREQUENCIES IN FLOATING POINT
C 
      XN0 = N0
      NSUM = N0
      DO 10 I = 1, P
        XN(I) = N(I)
        NSUM = NSUM + N(I)
 10   CONTINUE
      K = P - 1
      XNSUM = FLOAT(NSUM)
      NP = XN(P)
C 
C       ITERATIVE APPROXIMATION
C 
      DO 40 ITER = 1, MAXITR
C 
C       COMPUTE VALUES OF DISTRIBUTION AND DENSITY FUNCTIONS,
C       USING CURRENT VALUES OF MU, SIGMA
C 
        DO 20 I = 1, P
 20       CALL FUNC ((X(I) - MU)/SIGMA, F(I), F1(I))
        DM = ONE - F(P)
C 
C 
C       TEST FOR SMALL DIVISOR TO AVOID OVERFLOW
C 
        IF (ABS(DM).LT.RR) GO TO 90
        F1P = F1(P)
        IF (ABS(F(1)).LT.RR) GO TO 90
        XI1 = X(1) - MU
        XP = X(P) - MU
        R = F1(1)/F(1)
        S = F1P/DM
        T = -XN0*R
        U = NP*S
        A = T + U
        B = XI1*T + XP*U
        R = F1(1)*R
        S = F1P*S
        C = R + S
        R = XI1*S
        S = XP*S
        D = R + S
        E = XI1*R + XP*S
        DO 30 I = 1, K
          FI = F(I)
          FI1 = F(I + 1)
          F1I1 = F1(I + 1)
          F1I = F1(I)
          XI = XI1
          XI1 = X(I + 1) - MU
          NI = XN(I)
          R = FI1 - FI
          IF (ABS(R).LT.RR) GO TO 90
          S = F1I1 - F1I
          U = XI1*F1I1 - XI*F1I
          SR = S/R
          UR = U/R
          A = A - NI*SR
          B = B - NI*UR
          C = C + S*SR
          D = D + S*UR
          E = E + U*UR
 30     CONTINUE
        DENOM = (C*E - D*D)*XNSUM
C 
C       COMPUTE ADJUSTMENTS TO MU, SIGMA
C 
        SIGDEN = SIGMA/DENOM
        DMU = (E*A - D*B)*SIGDEN
        DSIGMA = (C*B - D*A)*SIGMA*SIGDEN
        MU = MU + DMU
        SIGMA = SIGMA + DSIGMA
        ERR = ABS(DMU) + ABS(DSIGMA)
C 
C       TEST FOR CONVERGENCE
C 
        IF (ERR.LT.EPS) GO TO 50
 40   CONTINUE
C 
C       SET FAULT IF LIMIT FOR NUMBER OF ITERATIONS IS
C       REACHED, THEN PROCEED
C 
      IFAULT = 4
      ITER = MAXITR
 50   DO 60 I = 1, P
 60     CALL FUNC ((X(I) - MU)/SIGMA, F(I), DUM)
C 
C       COMPUTE VARIANCES AND COVARIANCE OF ESTIMATES
C 
      SIGDEN = SIGMA*SIGMA/DENOM
      VARMU = E*SIGDEN
      SIGDEN = SIGMA*SIGDEN
      COVAR = -D*SIGDEN
      VARSIG = C*SIGMA*SIGDEN
      IF (VARMU.LT.ZERO .OR. VARSIG.LT.ZERO) GO TO 100
      SEMU = SQRT(VARMU)
      SESIG = SQRT(VARSIG)
C 
C       COMPUTE EXPECTED FREQUENCIES AND CHI SQUARE
C 
      E0 = XNSUM*F(1)
      EP = XNSUM*(ONE - F(P))
      EX(P) = EP
      CHISQ = ((XN0 - E0)**2)/E0 + ((NP - EP)**2)/EP
      DO 70 I = 1, K
        NN = XNSUM*(F(I+1) - F(I))
        CHISQ = CHISQ + ((NN - XN(I))**2)/NN
        EX(I) = NN
 70   CONTINUE
      RETURN
C 
C       ERROR EXITS
C 
 80   IFAULT = 1
      RETURN
 90   IFAULT = 2
      RETURN
 100  IFAULT = 3
      RETURN
      END
C
C
      SUBROUTINE INITL (P, X, N0, N, MU, SIGMA, IFAULT)
C 
C       ALGORITHM AS 95.1 APPL. STATIST. (1976) VOL.25, NO.1
C 
C       COMPUTES ROUGH LEAST SQUARES ESTIMATES OF MU AND SIGMA
C       ( FOR DEFINITION OF MU AND SIGMA SEE SUBROUTINE CURVE ).
C 
      INTEGER P
      REAL MU, ONE, ZERO
      DIMENSION X(P), N(P)
      DATA ONE/1.0/, ZERO/0.0/
C 
C       ERROR EXIT IF P TOO SMALL
C 
      IFAULT = 1
      IF (P.LT.2) RETURN
      IFAULT = 0
C 
C       COMPUTE AND FLOAT SUM OF FREQUENCIES
C 
      NSUM = N0
      DO 10 I = 1, P
 10     NSUM = NSUM + N(I)
      XNSUM = FLOAT(NSUM)
C 
C       ZERO ACCUMULATORS
C 
      NPAR = N0
      XBAR = ZERO
      YBAR = ZERO
      SXX = ZERO
      SXY = ZERO
      SW = ZERO
C 
C       COMPUTE WEIGHTED MEANS XAR, YBAR, AND CORRECTED SUMS
C       OF X*X AND X*Y.
C 
      DO 30 I = 1, P
C 
C       NULL FREQUENCIES AT EITHER END OF THE RANGE ARE
C       ZERO WEIGHTED
C 
        IF (NPAR.EQ.0 .OR. NPAR.EQ.NSUM) GO TO 20
        PROB = FLOAT(NPAR)/XNSUM
        Y = DEVIAT(PROB)
        CALL FUNC (Y, DUMMY, DFY)
        DX = X(I) - XBAR
        DY = Y - YBAR
        W = DFY*DFY/(PROB*(ONE - PROB))
        SW = SW + W
        FAC = W/SW
        XBAR = XBAR + FAC*DX
        YBAR = YBAR + FAC*DY
        FAC = W*DX*(ONE - FAC)
        SXX = SXX + FAC*DX
        SXY = SXY + FAC*DY
 20     NPAR = NPAR + N(I)
 30   CONTINUE
      SIGMA = SXX/SXY
      MU = XBAR - SIGMA*YBAR
      RETURN
      END

