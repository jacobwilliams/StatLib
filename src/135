      SUBROUTINE LFNORM(N, M, NDIM, MDIM, X, Y, BETA, Z, KY, IFAULT)
C
C     ALGORITHM AS 135  APPL. STATIST. (1979) VOL.28, NO. 1
C
C     Min-Max estimates for a linear multiple regression problem
C
      DOUBLE PRECISION X(NDIM, MDIM), Y(NDIM), BETA(MDIM), Z, HILO(20),
     +	XRXF(20), XSXF(20), LU(20, 20), ACU, BIG, ZERO, ONE, DEV1,
     +	SIGR, SUMXR, DEVIAT, YEST, YDEV, SUMXS, RATIO, TEST, DELTA,
     +	DIV, SWING, SAVE, SIGS, TOP, TWO
      INTEGER IBASE(20), SSS, RRR, INDX(20)
      LOGICAL INTL
C
      DATA ACU/1.D-15/, BIG/1.D30/, ZERO/0.D0/, ONE/1.D0/, TWO/2.D0/
C
      IFAULT = 0
      KY = 0
      Z = ZERO
      M1 = M - 1
C
C     Set up initial LU decomposition
C
      DO 10 I = 1, M
   10 INDX(I) = I
      INTL = .TRUE.
      KKK = 1
      CALL UPDATE(KKK, X, LU, IBASE, INDX, INTL, N, M, NDIM, MDIM,
     +	IFAULT)
      IF (IFAULT .NE. 0) RETURN
      INTL = .FALSE.
      IROW = KKK
C
C     Calculate beta value
C
      K = INDX(1)
      K1 = IBASE(1)
      BETA(K) = Y(K1) / LU(K, 1)
      DO 30 II = 2, M
	K = INDX(II)
	K1 = IBASE(II)
	BETA(K) = Y(K1)
	II1 = II - 1
	DO 20 I = 1, II1
	  KK = INDX(I)
	  BETA(K) = BETA(K) - LU(KK, II) * BETA(KK)
   20   CONTINUE
	BETA(K) = BETA(K) / LU(K, II)
   30 CONTINUE
      DO 40 II = 1, M1
	K1 = M - II
	K = INDX(K1)
	DO 40 I = 1, II
	  KK = M - I + 1
	  K2 = INDX(KK)
	  BETA(K) = BETA(K) - LU(K2, K1) * BETA(K2)
   40 CONTINUE
C
C     Search for and set first violated constraint as R-th constraint
C
   50 IROW = IROW + 1
      IF (IROW .GT. N) RETURN
      DEV1 = ZERO
      DO 60 I = 1, M
   60 DEV1 = DEV1 + X(IROW, I) * BETA(I)
      DEV1 = DEV1 - Y(IROW)
      IF (ABS(DEV1) .LT. ACU) GO TO 50
      SIGR = SIGN(ONE, DEV1)
      RRR = IROW
C
C     Adjust for the R-th constraint
C
      K = INDX(1)
      XRXF(1) = X(RRR, K)
      DO 80 II = 2, M
	K = INDX(II)
	XRXF(II) = X(RRR, K)
	II1 = II - 1
	DO 70 I = 1, II1
   70   XRXF(II) = XRXF(II) - LU(K, I) * XRXF(I)
   80 CONTINUE
      K = INDX(M)
      XRXF(M) = XRXF(M) / LU(K, M)
      HILO(M) = SIGN(ONE, -SIGR * XRXF(M))
      SUMXR = SIGR - HILO(M) * XRXF(M)
      DO 100 II = 1, M1
	K1 = M - II
	K = INDX(K1)
	DO 90 I = 1, II
	  K2 = M - I + 1
	  XRXF(K1) = XRXF(K1) - LU(K, K2) * XRXF(K2)
   90   CONTINUE
	XRXF(K1) = XRXF(K1) / LU(K, K1)
	HILO(K1) = SIGN(ONE, -SIGR * XRXF(K1))
	SUMXR = SUMXR - HILO(K1) * XRXF(K1)
  100 CONTINUE
      Z = ABS(DEV1 / SUMXR)
C
C     Start of main iterative loop.
C     Search for the most violated S-th constraint
C
  110 SSS = 0
      DEVIAT = ACU
C
C     Calculate beta value
C
      K = INDX(1)
      K1 = IBASE(1)
      BETA(K) = (Y(K1) + Z * HILO(1)) / LU(K, 1)
      DO 130 II = 2, M
	K = INDX(II)
	K1 = IBASE(II)
	BETA(K) = Y(K1) + Z * HILO(II)
	II1 = II - 1
	DO 120 I = 1, II1
	  KK = INDX(I)
	  BETA(K) = BETA(K) - LU(KK, II) * BETA(KK)
  120   CONTINUE
	BETA(K) = BETA(K) / LU(K, II)
  130 CONTINUE
      DO 140 II = 1, M1
	K1 = M - II
	K = INDX(K1)
	DO 140 I = 1, II
	  KK = M - I + 1
	  K2 = INDX(KK)
	  BETA(K) = BETA(K) - LU(K2, K1) * BETA(K2)
  140 CONTINUE
C
C     Calculate residuals
C
      DO 160 I = 1, N
	YEST = ZERO
	DO 150 J = 1, M
  150   YEST = YEST + X(I, J) * BETA(J)
	DEV1 = ABS(Y(I) - YEST) - Z
	IF (DEV1 .LE. DEVIAT) GO TO 160
	YDEV = YEST - Y(I)
	DEVIAT = DEV1
	SSS = I
  160 CONTINUE
C
C     Check if at optimum
C
      IF (SSS .EQ. 0) RETURN
C
C     Set up information on the S-th constraint
C
      SIGS = SIGN(ONE, YDEV)
      K = INDX(1)
      XSXF(1) = X(SSS, K)
      DO 180 II = 2, M
	K = INDX(II)
	XSXF(II) = X(SSS, K)
	II1 = II - 1
	DO 170 I = 1, II1
  170   XSXF(II) = XSXF(II) - LU(K, I) * XSXF(I)
  180 CONTINUE
      K = INDX(M)
      XSXF(M) = XSXF(M) / LU(K, M)
      SUMXS = -SIGS + HILO(M) * XSXF(M)
      DO 200 II = 1, M1
	K1 = M - II
	K = INDX(K1)
	DO 190 I = 1, II
	  K2 = M - I + 1
	  XSXF(K1) = XSXF(K1) - LU(K, K2) * XSXF(K2)
  190   CONTINUE
	XSXF(K1) = XSXF(K1) / LU(K, K1)
	SUMXS = SUMXS + HILO(K1) * XSXF(K1)
  200 CONTINUE
C
C     Search for minimum ratio
C
  210 KKK = 0
      RATIO = BIG
      DO 220 I = 1, M
	IF (SIGS * SIGN(ONE, XSXF(I)) .NE. HILO(I) .OR.
     +		ABS(XSXF(I)) .LT. ACU) GO TO 220
	TEST = ABS(XRXF(I) / XSXF(I))
	IF (TEST .GE. RATIO) GO TO 220
	RATIO = TEST
	KKK = I
  220 CONTINUE
C
C     Check if R-th constraint moves interior
C
      IF (KKK .NE. 0) GO TO 260
C
C     Process the movement of the R-th constraint
C
      DELTA = ABS(DEVIAT / SUMXS)
C
C     Calculate the largest tolerable delta
C
      DIV = ABS(SUMXR) - TWO
      IF (DIV .LT. ACU) GO TO 240
      SWING = TWO * Z / DIV
      IF (SWING .GE. DELTA) GO TO 240
C
C     Switch R and S constraint indicators
C
      SAVE = SUMXS
      SUMXS = -SUMXR + SIGR + SIGR
      SUMXR = -SAVE
      SAVE = SIGR
      SIGR = SIGS
      SIGS = -SAVE
      DEVIAT = ABS(SUMXS * DELTA) - TWO * Z
      Z = Z + DELTA
      DO 230 I = 1, M
	SAVE = XSXF(I)
	XSXF(I) = XRXF(I)
	XRXF(I) = SAVE
  230 CONTINUE
      I = RRR
      RRR = SSS
      SSS = I
      GO TO 210
C
C     Replace the R-th constraint with the S-th constraint
C
  240 SIGR = SIGS
      DO 250 I = 1, M
  250 XRXF(I) = XSXF(I)
      SUMXR = -SUMXS
      Z = Z + DELTA
      RRR = SSS
      GO TO 110
C
C     Process the movement of the K-th constraint
C
  260 DELTA = ABS(XRXF(KKK) * DEVIAT /
     +		(XRXF(KKK) * SUMXS + XSXF(KKK) * SUMXR))
      TOP = -TWO * Z * XRXF(KKK)
      DIV = XRXF(KKK) * XRXF(KKK) + HILO(KKK) * SUMXR
      IF (SIGN(ONE, TOP) .NE. SIGN(ONE, DIV)) GO TO 270
      IF (ABS(DIV) .LT. ACU) GO TO 270
      SWING = TOP / DIV
C
C     Check to see if the K-th constraint swings across
C
      IF (SWING .GE. DELTA) GO TO 270
      Z = Z + SWING
      DEVIAT = DEVIAT - SWING *
     +		ABS(SUMXS + XSXF(KKK) * SUMXR / XRXF(KKK))
      SUMXR = SUMXR + TWO * HILO(KKK) * XRXF(KKK)
      SUMXS = SUMXS - TWO * HILO(KKK) * XSXF(KKK)
      HILO(KKK) = -HILO(KKK)
      GO TO 210
C
C     Update XRXF and the LU of the current basis
C
  270 HILO(KKK) = SIGS
      SUMXR = SIGR
      XRXF(KKK) = XRXF(KKK) / XSXF(KKK)
      SUMXR = SUMXR - HILO(KKK) * XRXF(KKK)
      DO 280 I = 1, M
	IF (I .EQ. KKK) GO TO 280
	XRXF(I) = XRXF(I) - XSXF(I) * XRXF(KKK)
	SUMXR = SUMXR - HILO(I) * XRXF(I)
  280 CONTINUE
      IBASE(KKK) = SSS
C
C     Update LU decomposition
C
      CALL UPDATE(KKK, X, LU, IBASE, INDX, INTL, N, M, NDIM, MDIM,
     +	IFAULT)
      IF (IFAULT .NE. 0) RETURN
      Z = Z + DELTA
      KY = KY + 1
      GO TO 110
      END
C
      SUBROUTINE UPDATE(KKK, X, LU, IBASE, INDX, INTL, N, M, NDIM,
     +	MDIM, IFAULT)
C
C     ALGORITHM AS 135.1  APPL. STATIST. (1979) VOL.28, NO. 1
C
C     Update LU decomposition matrix
C
      DOUBLE PRECISION X(NDIM, MDIM), LU(20, 20), ACU, SUBT, PIVOT
      INTEGER IBASE(20), INDX(20)
      LOGICAL INTL
      DATA ACU/1.D-15/
C
      IROW = 0
      DO 90 II = KKK, M
	IF (INTL) GO TO 10
	IROW = IBASE(II)
	GO TO 20
   10   IROW = IROW + 1
	IF (IROW .LE. N) GO TO 20
	IFAULT = 1
	RETURN
   20   DO 30 I = 1, M
   30   LU(I, II) = X(IROW, I)
C
C     Set up representation of incoming row
C
	IF (II .EQ. 1) GO TO 60
	II1 = II - 1
	DO 50 ICOL = 1, II1
	  K = INDX(ICOL)
	  SUBT = LU(K, II)
	  J = ICOL + 1
	  DO 40 I = J, M
	    K = INDX(I)
	    LU(K, II) = LU(K, II) - SUBT * LU(K, ICOL)
   40     CONTINUE
   50   CONTINUE
C
C     Find maximum entry
C
   60   PIVOT = ACU
        KK = 0
	DO 70 I = II, M
	  K = INDX(I)
	  IF (ABS(LU(K, II)) .LE. PIVOT) GO TO 70
	  PIVOT = ABS(LU(K, II))
	  KK = I
   70   CONTINUE
	IF (KK .EQ. 0) GO TO 10
C
C     Switch order
C
	ISAVE = INDX(KK)
	INDX(KK) = INDX(II)
	INDX(II) = ISAVE
C
C     Put into columns of LU one at a time
C
	IF (INTL) IBASE(II) = IROW
	IF (II .EQ. M) GO TO 90
	J = II + 1
	DO 80 I = J, M
	  K = INDX(I)
	  LU(K, II) = LU(K, II) / LU(ISAVE, II)
   80   CONTINUE
   90 CONTINUE
      KKK = IROW
      RETURN
      END

