      SUBROUTINE ORPS1(ALPHA, LAMDA, DELTA, IP, SIGMA, KNEW, FM10, FM1K,
     *		FM20, FM2K, ITER, IFAULT)
C
C     ALGORITHM AS223  APPL. STATIST. (1987) VOL. 36, NO. 1
C
C     Calculates the optimum ridge parameter under the criterion of
C     minimizing the mean squared error of parameter estimation.
C
      INTEGER FM20, IMAX, IP, ITER, IFAULT
      REAL ALPHA(IP), LAMDA(IP), DELTA, K, KNEW, KOLD, DK, FM10, FM1K,
     *	FM2K, DFM1K, DDFM1K, S, SIGMA
      REAL ZERO, TWO, THREE
C
      FM1(X, Y) = (Y * S**2 + (X * K)**2) / ((Y + K)**2 * S**2)
      FM2(X, Y) = Y * FM1(X, Y)
      DFM1(X, Y) = Y * (X**2 * K - S**2) / (Y + K)**3
      DDFM1(X, Y) = Y * (Y * X**2 - TWO * X**2 * K + THREE * S**2) /
     *			(Y + K)**4
      DATA ZERO/0.0/, TWO/2.0/, THREE/3.0/
C
      IFAULT = 1
      IF (DELTA .LE. ZERO) RETURN
      IFAULT = 2
      IF (SIGMA .LE. ZERO) RETURN
      IFAULT = 3
      DO 10 J = 1, IP
	IF (LAMDA(J) .LE. ZERO) RETURN
   10 CONTINUE
      IFAULT = 0
      K = ZERO
      FM10 = ZERO
      S = SIGMA
      DO 20 J = 1, IP
	FM10 = FM10 + FM1(ALPHA(J), LAMDA(J))
   20 CONTINUE
      FM20 = IP
      ITER = 0
      IMAX = 20
      KNEW = ZERO
      DO 40 I = 1, IMAX
	ITER = ITER + 1
	FM1K = ZERO
	FM2K = ZERO
	DFM1K = ZERO
	DDFM1K = ZERO
	DO 30 J = 1, IP
	  FM1K = FM1K + FM1(ALPHA(J), LAMDA(J))
	  FM2K = FM2K + FM2(ALPHA(J), LAMDA(J))
	  DFM1K = DFM1K + DFM1(ALPHA(J), LAMDA(J))
	  DDFM1K = DDFM1K + DDFM1(ALPHA(J), LAMDA(J))
   30   CONTINUE
	KOLD = KNEW
	KNEW = KOLD - DFM1K / DDFM1K
	K = KNEW
	DK = ABS(KNEW - KOLD)
	IF (DK .LE. DELTA) GO TO 50
   40 CONTINUE
   50 RETURN
      END
C
C
C
      SUBROUTINE ORPS2(ALPHA, LAMDA, DELTA, IP, SIGMA, KNEW, FM10, FM1K,
     *		FM20, FM2K, ITER, IFAULT)
C
C     ALGORITHM AS223  APPL. STATIST. (1987) VOL. 36, NO. 1
C
C     Calculates the optimum ridge parameter under the criterion of
C     minimizing the mean squared error of prediction.
C
      INTEGER FM20, IMAX, IP, ITER, IFAULT
      REAL ALPHA(IP), LAMDA(IP), DELTA, K, KNEW, KOLD, DK, FM10, FM1K,
     *	FM2K, DFM2K, DDFM2K, S, SIGMA
      REAL ZERO, TWO, THREE
C
      FM1(X, Y) = (Y * S**2 + (X * K)**2) / ((Y + K)**2 * S**2)
      FM2(X, Y) = Y * FM1(X, Y)
      DFM1(X, Y) = Y * (X**2 * K - S**2) / (Y + K)**3
      DDFM1(X, Y) = Y * (Y * X**2 - TWO * X**2 * K + THREE * S**2) /
     *			(Y + K)**4
      DDFM2(X, Y) = Y * DDFM1(X, Y)
      DATA ZERO/0.0/, TWO/2.0/, THREE/3.0/
C
      IFAULT = 1
      IF (DELTA .LE. ZERO) RETURN
      IFAULT = 2
      IF (SIGMA .LE. ZERO) RETURN
      IFAULT = 3
      DO 10 J = 1, IP
	IF (LAMDA(J) .LE. ZERO) RETURN
   10 CONTINUE
      IFAULT = 0
      K = ZERO
      FM10 = ZERO
      S = SIGMA
      DO 20 J = 1, IP
	FM10 = FM10 + FM1(ALPHA(J), LAMDA(J))
   20 CONTINUE
      FM20 = IP
      ITER = 0
      IMAX = 20
      KNEW = ZERO
      DO 40 I = 1, IMAX
	ITER = ITER + 1
	FM1K = ZERO
	FM2K = ZERO
	DFM2K = ZERO
	DDFM2K = ZERO
	DO 30 J = 1, IP
	  FM1K = FM1K + FM1(ALPHA(J), LAMDA(J))
	  FM2K = FM2K + FM2(ALPHA(J), LAMDA(J))
	  DFM2K = DFM2K + DFM2(ALPHA(J), LAMDA(J))
	  DDFM2K = DDFM2K + DDFM2(ALPHA(J), LAMDA(J))
   30   CONTINUE
	KOLD = KNEW
	KNEW = KOLD - DFM2K / DDFM2K
	K = KNEW
	DK = ABS(KNEW - KOLD)
	IF (DK .LE. DELTA) GO TO 50
   40 CONTINUE
   50 RETURN
      END
