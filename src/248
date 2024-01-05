      SUBROUTINE EDFGOF(X, N, ITEST, Z, XBAR, S2, D, V, W2, W2P, U2,
     +   U2P, A2, A2P, IFAULT)
C
C     ALGORITHM AS 248.1  APPL. STATIST. (1989), VOL.38, NO. 3
C
C     Tests a sample for uniformity, normality or exponentiality
C     using goodness-of-fit statistics based on the empirical
C     distribution function
C
C     Auxiliary routines called: ALNORM from AS 66
C
      INTEGER I, IFAULT, ITEST, J, N
      REAL  X(N), Z(N), XBAR, S2, D, V, W2, W2P, U2, U2P, A2, A2P, RN,
     +      ROOTN, RN2, S, ZERO, PT01, PT05, PT1, PT11, PT12, PT155,
     +      PT16, PT2, PT24, PT26, PT3, PT35, PT4, PT5, PT6, PT75, PT8,
     +      PT82, PT85, ONE, TWOP25, WCUT2(3), WCOEF2(4, 3), WCUT3(3),
     +      WCOEF3(4, 3), UCUT2(3), UCOEF2(4, 3), UCUT3(3),
     +      UCOEF3(4, 3), ACUT2(3), ACOEF2(4, 3), ACUT3(3), ACOEF3(4, 3)
      REAL ALNORM, PVALUE
      DATA ZERO, PT01, PT05, PT1, PT11, PT12, PT155, PT16, PT2, PT24
     +     /0.0, 0.01, 0.05, 0.1, 0.11, 0.12, 0.155, 0.16, 0.2, 0.24/
      DATA PT26, PT3, PT35, PT4, PT5, PT6, PT75, PT8, PT82, PT85, ONE
     +     /0.26, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75, 0.8, 0.82, 0.85, 1.0/
      DATA TWOP25 /2.25/
      DATA WCUT2, WCUT3, UCUT2, UCUT3, ACUT2, ACUT3 /0.0275, 0.051,
     +     0.092, 0.035, 0.074, 0.160, 0.0262, 0.048, 0.094, 0.029,
     +     0.062, 0.120, 0.200, 0.340, 0.600, 0.260, 0.510, 0.950/
      DATA ((WCOEF2(I,J), J=1,3), I=1,4) / -13.953, 775.5, -12542.61,
     +     -5.903, 179.546, -1515.29, 0.886, -31.62, 10.897, 1.111,
     +     -34.242, 12.832/
      DATA ((WCOEF3(I,J), J=1,3), I=1,4) / -11.334, 459.098, -5652.1,
     +     -5.779, 132.89, -866.58, 0.586, -17.87, 7.417, 0.447,
     +     -16.592, 4.849/
      DATA ((UCOEF2(I,J), J=1,3), I=1,4) / -13.642, 766.31, -12432.74,
     +     -6.3328, 214.57, -2022.28, 0.8510, -32.006, -3.45, 1.325,
     +     -38.918, 16.45/
      DATA ((UCOEF3(I,J), J=1,3), I=1,4) / -11.703, 542.5, -7574.59,
     +     -6.3288, 178.1, -1399.49, 0.8071, -25.166, 8.44, 0.7663,
     +     -24.359, 4.539/
      DATA ((ACOEF2(I,J), J=1,3), I=1,4) / -13.436, 101.14, -223.73,
     +     -8.318, 42.796, -59.938, 0.9177, -4.279, -1.38, 1.2937,
     +     -5.709, 0.0186/
      DATA ((ACOEF3(I,J), J=1,3), I=1,4) / -12.2204, 67.459, -110.3,
     +     -6.1327, 20.218, -18.663, 0.9209, -3.353, 0.3, 0.731, -3.009,
     +     0.15/
C
      IFAULT = 1
      IF (ITEST .LT. 1 .OR. ITEST .GT. 3) RETURN
      IFAULT = 2
      IF (N .LT. 2 .OR. (ITEST .EQ. 2 .AND. N .EQ. 2)) RETURN
      IFAULT = 3
      DO 10 I = 2, N
	IF (X(I) .LT. X(I-1)) RETURN
   10 CONTINUE
      RN = N
      ROOTN = SQRT(RN)
      RN2 = RN * RN
      IF (ITEST .EQ. 1) THEN
C
C     Test for uniformity (F(x) is completely specified)
C
	IFAULT = 0
	CALL STATS(X, N, D, V, W2, U2, A2, IFAULT)
	IF (IFAULT .NE. 0) RETURN
C
C     Modifications when F(x) is completely specified
C
	D = D * (ROOTN + PT12 + PT11 / ROOTN)
	V = V * (ROOTN + PT155 + PT24 / ROOTN)
	W2 = (W2 - PT4 / RN + PT6 / RN2) * (ONE + ONE / RN)
	U2 = (U2 - PT1 / RN + PT1 / RN2) * (ONE + PT8 / RN)
	RETURN
      ELSE
C
C     Estimate the mean XBAR
C
	XBAR = ZERO
	DO 20 I = 1, N
	  XBAR = XBAR + X(I)
   20   CONTINUE
	XBAR = XBAR / RN
	IF (ITEST .EQ. 2) THEN
C
C     Test for normality (MU and SIGMA**2 unspecified).
C     First estimate the variance S**2
C
	  IFAULT = 5
	  S2 = ZERO
	  DO 30 I = 1, N
	    S2 = S2 + (X(I) - XBAR)**2
   30     CONTINUE
	  S2 = S2 / (RN - ONE)
	  IF (S2 .LE. ZERO) RETURN
C
C     Compute Z(I) = F(XI - XBAR)/S)
C
	  S = SQRT(S2)
	  DO 40 I = 1, N
	    Z(I) = ALNORM((X(I) - XBAR)/S, .FALSE.)
   40     CONTINUE
	  CALL STATS(Z, N, D, V, W2, U2, A2, IFAULT)
	  IF (IFAULT .NE. 0) RETURN
C
C     Modifications when F(x) is the normal distribution
C
	  D = D * (ROOTN - PT01 + PT85 / ROOTN)
	  V = V * (ROOTN + PT05 + PT82 / ROOTN)
	  W2 = W2 * (ONE + PT5 / RN)
	  W2P = PVALUE(W2, WCUT2, WCOEF2)
	  U2 = U2 * (ONE + PT5 / RN)
	  U2P = PVALUE(U2, UCUT2, UCOEF2)
	  A2 = A2 * (ONE + PT75 / RN + TWOP25 / RN2)
	  A2P = PVALUE(A2, ACUT2, ACOEF2)
	  RETURN
	ELSE
C
C     Test for exponentiality (scale parameter unspecified)
C
	  IFAULT = 6
	  IF (XBAR .LE. ZERO) RETURN
	  DO 50 I = 1, N
	    Z(I) = ONE - EXP(-X(I) / XBAR)
   50     CONTINUE
	  CALL STATS(Z, N, D, V, W2, U2, A2, IFAULT)
	  IF (IFAULT .NE. 0) RETURN
C
C     Modifications when F(x) is the exponential distribution
C
	  D = (D - PT2 / RN) * (ROOTN + PT26 + PT5 / ROOTN)
	  V = (V - PT2 / RN) * (ROOTN + PT24 + PT35 / ROOTN)
	  W2 = W2 * (ONE + PT16 / RN)
	  W2P = PVALUE(W2, WCUT3, WCOEF3)
	  U2 = U2 * (ONE + PT16 / RN)
	  U2P = PVALUE(U2, UCUT3, UCOEF3)
	  A2 = A2 * (ONE + PT3 / RN)
	  A2P = PVALUE(A2, ACUT3, ACOEF3)
	  RETURN
	END IF
      END IF
      END
C
C
      SUBROUTINE STATS(Z, N, D, V, W2, U2, A2, IFAULT)
C
C     ALGORITHM AS 248.2  APPL. STATIST. (1989), VOL.38, NO. 3
C
C     Computes the goodness-of-fit statistics D, V, W**2, U**2 and
C     A**2 from the transformed Z values
C
      INTEGER I, IFAULT, N, NI
      REAL Z(N), D, V, W2, U2, A2, RI, RN, DPLUS, DMINUS, D1, D2, SUMZ,
     +   TWON, ZM1, A2SUM, ZERO, SMALL, HALF, ONE, TWO, TWELVE
C
C     Initialize constants
C
      DATA ZERO/0.0/, SMALL/1.E-37/, HALF/0.5/, ONE/1.0/,
     +     TWO/2.0/, TWELVE/12.0/
C
      RN = N
      TWON = TWO * RN
C
C     Calculating the Kolmorogov statistics DPLUS, DMINUS, D and
C     the Kuiper statistic V.
C
      DPLUS = ZERO
      DMINUS = ZERO
      IFAULT = 4
      DO 10 I = 1, N
	IF (Z(I) .LE. ZERO .OR. Z(I) .GE. ONE) RETURN
	RI = I
	D1 = RI / RN - Z(I)
	IF (D1 .GT. DPLUS) DPLUS = D1
	D2 = Z(I) - (RI - ONE) / RN
	IF (D2 .GT. DMINUS) DMINUS = D2
   10 CONTINUE
      IFAULT = 0
      D = DPLUS
      IF (DMINUS .GT. DPLUS) D = DMINUS
      V = DPLUS + DMINUS
C
C     Calculating the Cramer-von Mises statistic W2 and the Watson
C     statistic U2
C
      W2 = ONE / (TWELVE * RN)
      SUMZ = ZERO
      DO 20 I = 1, N
	W2 = W2 + (Z(I) - (TWO * I - ONE) / TWON)**2
	SUMZ = SUMZ + Z(I)
   20 CONTINUE
      U2 = W2 - RN * (SUMZ / RN - HALF)**2
C
C     Calculating the Anderson-Darling statistic A2
C
      A2SUM = ZERO
      NI = N
      DO 30 I = 1, N
	IF (Z(I) .LT. SMALL) Z(I) = SMALL
	ZM1 = ONE - Z(NI)
	IF (ZM1 .LT. SMALL) ZM1 = SMALL
	A2SUM = A2SUM + (TWO * I - ONE) * LOG(Z(I) * ZM1)
	NI = NI - 1
   30 CONTINUE
      A2 = -A2SUM / RN - RN
C
      RETURN
      END
C
C
      REAL FUNCTION PVALUE(T, CUT, COEF)
C
C     ALGORITHM AS 248.3  APPL. STATIST. (1989), VOL.38, NO. 3
C
C     Computes the approximate significance level for W**2, U**2 or
C     A**2 when testing for normality or exponentiality
C
      INTEGER I, J
      REAL T, CUT(3), COEF(4, 3), ONE
      DATA ONE /1.0/
C
      I = 4
      DO 10 J = 1, 3
	IF (T .LT. CUT(J)) I = I - 1
   10 CONTINUE
      PVALUE = EXP(COEF(I, 1) + T * (COEF(I, 2) + T * COEF(I, 3)))
      IF (I .LT. 3) PVALUE = ONE - PVALUE
C
      RETURN
      END
