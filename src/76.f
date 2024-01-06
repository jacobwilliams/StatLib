C Two versions of algorithm AS 76 are given here; the original with one
C correction incorporated, and AS R55, also amended.   AS R55 requires
C AS 76.   N.B. The accuracy of AS 76 could be increased by using more
C Gaussian quadrature points, or better, by using Hermite integration.
C
      FUNCTION TFN(X, FX)
C
C     ALGORITHM AS 76  APPL. STATIST. (1974) VOL.23, NO.3
C
C     Calculates the T-function of Owen, using Gaussian quadrature.
C     Incorporates correction AS R30 (vol.28, no.1, 1979)
C
      REAL U(5), R(5)
C
      DATA U /0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533/
      DATA R /0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357/
      DATA NG,    TP,    TV1,     TV2,     TV3,     TV4
     *   /  5, 0.159155, 1.E-35,  15.0,    15.0,   1.E-5 /
      DATA ZERO, QUART, HALF, ONE,  TWO
     *   / 0.0,  0.25,  0.5,  1.0,  2.0 /
C
C     Test for X near zero
C
      IF (ABS(X) .GE. TV1) GO TO 5
      TFN = TP * ATAN(FX)
      RETURN
C
C     Test for large values of abs(X)
C
    5 IF (ABS(X) .GT. TV2) GO TO 10
C
C     Test for FX near zero
C
      IF (ABS(FX) .GE. TV1) GO TO 15
   10 TFN = ZERO
      RETURN
C
C     Test whether abs(FX) is so large that it must be truncated
C
   15 XS = -HALF * X * X
      X2 = FX
      FXS = FX * FX
      IF (LOG(ONE + FXS) - XS * FXS .LT. TV3) GO TO 25
C
C     Computation of truncation point by Newton iteration
C
      X1 = HALF * FX
      FXS = QUART * FXS
   20 RT = FXS + ONE
      X2 = X1 + (XS * FXS + TV3 - LOG(RT)) / (TWO * X1 * (ONE/RT - XS))
      FXS = X2 * X2
      IF (ABS(X2 - X1) .LT. TV4) GO TO 25
      X1 = X2
      GO TO 20
C
C     Gaussian quadrature
C
   25 RT = ZERO
      DO 30 I = 1, NG
	R1 = ONE + FXS * (HALF + U(I))**2
	R2 = ONE + FXS * (HALF - U(I))**2
	RT = RT + R(I) * (EXP(XS * R1) / R1 + EXP(XS * R2) / R2)
   30 CONTINUE
      TFN = RT * X2 * TP
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
      REAL FUNCTION THA(H1, H2, A1, A2)
C
C     AS R55  APPL. STATIST. (1985) VOL.34, NO.1
C
C     A remark on AS 76
C     Incorporating improvements in AS R80 (Appl. Statist. (1989) 
C     vol.38, no.3), and AS R89 (Appl. Statist. (1992) vol.41, no.2).
C
C     Computes T(H1/H2, A1/A2) for any real numbers H1, H2, A1 and A2
C
C     Auxiliary function required: ALNORM (= AS 66) and AS 76
C
      REAL A, ALNORM, A1, A2, G, H, H1, H2, TFN, ABSA, AH, GH, GAH,
     *  TWOPI, LAM, EX, C1, C2,
     *  AH, ZERO, ONE, TWO, PT3, SEVEN, HALF, SIX, QUART
C
      DATA TWOPI /6.2831853/, ZERO /0.0/, ONE /1.0/, TWO /2.0/,
     *   PT3 /0.3/, SEVEN /7.0/, HALF /0.5/, SIX /6.0/, QUART /0.25/
C
      IF (H2 .NE. ZERO) GO TO 1
      THA = ZERO
      RETURN
C
    1 H = H1 / H2
      IF (A2 .EQ. ZERO) GO TO 2
      A = A1 / A2
      IF ((ABS(H) .LT. PT3) .AND. (ABS(A) .GT. SEVEN)) GO TO 6
C
C     Correction AS R89
C
      ABSA = ABS(A)
      IF (ABSA .GT. ONE) GO TO 7
      THA = TFN(H, A)
      RETURN
    7 AH = ABSA * H
      GH = ALNORM(H, .FALSE.)
      GAH = ALNORM(AH, .FALSE.)
      THA = HALF * (GH + GAH) - GH * GAH - TFN(AH, ONE/ABSA)
      IF (A .LT. ZERO) THA = - THA
      RETURN
C
    2 G = ALNORM(H, .FALSE.)
      IF (H .GE. ZERO) GO TO 3
      THA = G / TWO
      GO TO 4
    3 THA = (ONE - G) / TWO
    4 IF (A1 .GE. ZERO) RETURN
      THA = -THA
      RETURN
C
    6 LAM = ABS(A * H)
      EX = EXP(-LAM * LAM / TWO)
      G = ALNORM(LAM, .FALSE.)
      C1 = (EX/LAM + SQRT(TWOPI) * (G - HALF)) / (TWOPI)
      C2 = ((LAM * LAM + TWO) * EX/LAM**3 + SQRT(TWOPI) * (G - HALF))
     *       / (SIX * TWOPI)
      AH = ABS(H)
      THA = QUART - C1 * AH + C2 * AH**3
      THA = SIGN(THA, A)
      RETURN
      END

