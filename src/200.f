      SUBROUTINE SUMSQ(N, APPROX, BL, BU, IFAULT)
C
C     ALGORITHM AS 200  APPL. STATIST. (1984) VOL.33, NO.2
C
C     Approximates the sum of squares of normal scores and gives bounds
C     for the truncation error involved in the approximation.
C
C     Coefficients SA(K) are the squares of the Fourier coefficients
C     A(2*K-1).
C
      REAL APPROX, BL, BU, SA(14), ZERO, ONE, PROD, RI, RK, RN, SUM1,
     *     SUM2
C
      DATA SA/
     *  0.95492 96583, 0.03349 20160, 0.00667 47227, 0.00227 80925,
     *  0.00101 63627, 0.00053 26443, 0.00031 08437, 0.00019 58978,
     *  0.00013 08187, 0.00009 13941, 0.00006 63703, 0.00004 96799,
     *  0.00004 11471, 0.00003 49349/
      DATA ZERO/0.0/, ONE/1.0/, NMAX/27/
C
C     Value of zero is arbitrarily given in case of failure
C
      APPROX = ZERO
      BL = ZERO
      BU = ZERO
C
C     Check consistency of input parameter
C
      IFAULT = 1
      IF (N .LE. 1) RETURN
      IFAULT = 0
C
      RN = N
      NM = NMAX - 2
      IF (N .GT. NMAX) GO TO 1
      M = N / 2
      NM = 2 * M - 1
    1 SUM1 = ZERO
      SUM2 = ZERO
      KK = 0
      DO 3 K = 1, NM, 2
	RK = K
	PROD = ONE
	DO 2 I = 1, K
	  RI = I
	  PROD = PROD * (RN - RI) / (RN + RK - RI + ONE)
    2   CONTINUE
	KK = KK + 1
	PROD = PROD * RN
	SUM1 = SUM1 + SA(KK) * PROD
	SUM2 = SUM2 + SA(KK)
    3 CONTINUE
C
C     Approximate value is computed
C
      APPROX = SUM1
      IF (N .LE. NMAX) RETURN
      PROD = ONE
      DO 4 I = 1, NMAX
	RI = I
	PROD = PROD * (RN - RI) / (RN + RI)
    4 CONTINUE
      PROD = PROD * RN
C
C     Lower & upper bounds for the truncation error are computed
C
      BL = PROD * SA(14)
      BU = PROD * (ONE - SUM2)
C
      RETURN
      END
