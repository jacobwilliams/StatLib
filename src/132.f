      SUBROUTINE SIMLP(N, X, Y, SAD, ALPHA, BETA, D, ITER, INEXT,
     *   IFAULT)
C
C     ALGORITHM AS 132  APPL. STATIST. (1978) VOL.27, NO.3
C
C     SIMPL:   Fit  Y = ALPHA + BETA.X + error
C
      REAL X(N), Y(N), D(N)
      INTEGER INEXT(N)
      DATA ACU/1.0E-06/, BIG/1.0E19/, HALF/0.5/, ZERO/0.0/, ONE/1.0/,
     *     TWO/2.0/
C
C     Initial settings
C
      IFAULT = 0
      ITER = 0
      AHALF = HALF + ACU
      AONE = AHALF + AHALF
C
C     Determine initial basis
C
      D(1) = ZERO
      Y1 = Y(1)
      IBAS1 = 1
      A1 = X(1)
      DO 10 I = 2, N
	IF (ABS(A1 - X(I)) .LT. ACU) GO TO 10
	A2 = X(I)
	IBAS2 = I
	Y2 = Y(I)
	GO TO 20
   10 CONTINUE
      IFAULT = 1
      RETURN
C
C     Calculate initial beta value
C
   20 DET = ONE / (A2 - A1)
      AAAA = (A2 * Y1 - A1 * Y2) * DET
      BBBB = (Y2 - Y1) * DET
C
C     Calculate initial D-vector
C
      DO 30 I = 2, N
	DDD = Y(I) - (AAAA + BBBB * X(I))
	D(I) = SIGN(ONE, DDD)
   30 CONTINUE
      TOT1 = ONE
      TOT2 = X(IBAS2)
      D(IBAS2) = - ONE
      DO 40 I = 2, N
	TOT1 = TOT1 + D(I)
	TOT2 = TOT2 + D(I) * X(I)
   40 CONTINUE
      T = (A2 * TOT1 - TOT2) * DET
      IF (ABS(T) .LT. AONE) GO TO 50
      DET = - DET
      GO TO 70
C
C     Main iterative loop begins
C
   50 T = (TOT2 - A1 * TOT1) * DET
      IF (ABS(T) .LT. AONE) GO TO 160
      IFLAG = 2
      IOUT = IBAS2
      X(IOUT) = A1
      AAA = A1
      BBB = A2
      GO TO 80
   60 T = (TOT2 - A2 * TOT1) * DET
      IF (ABS(T) .LT. AONE) GO TO 160
   70 IFLAG = 1
      BBB = A1
      AAA = A2
      IOUT = IBAS1
   80 RHO = SIGN(ONE, T)
      T = HALF * ABS(T)
      DET = DET * RHO
C
C     Perform partial sort of ratios
C
      INEXT(IBAS1) = IBAS2
      RATIO = BIG
      SUM = AHALF
      DO 120 I = 1, N
	DDD = (X(I) - AAA) * DET
	IF (DDD * D(I) .LE. ACU) GO TO 120
	TEST = (Y(I) - AAAA - BBBB * X(I)) / DDD
	IF (TEST .GE. RATIO) GO TO 120
	J = IBAS1
	SUM = SUM + ABS(DDD)
   90   ISAVE = ABS(INEXT(J))
	IF (TEST .GE. D(ISAVE)) GO TO 110
	IF (SUM .LT. T) GO TO 100
	SUBT = ABS((X(ISAVE) - AAA) * DET)
	IF (SUM - SUBT .LT. T) GO TO 100
	SUM = SUM - SUBT
	D(ISAVE) = SIGN(1, INEXT(J))
	INEXT(J) = INEXT(ISAVE)
	GO TO 90
  100   J = ISAVE
	ISAVE = ABS(INEXT(J))
	IF (TEST .LT. D(ISAVE)) GO TO 100
  110   INEXT(I) = INEXT(J)
	INEXT(J) = SIGN(I, INT(D(I)))
	D(I) = TEST
	IF (SUM .LT. T) GO TO 120
	IIN = ABS(INEXT(IBAS1))
	RATIO = D(IIN)
  120 CONTINUE
C
C     Update basic indicators
C
      IIN = ABS(INEXT(IBAS1))
      J = IIN
  130 ISAVE = ABS(INEXT(J))
      IF (ISAVE .EQ. IBAS2) GO TO 140
      ZZZ = SIGN(1, INEXT(J))
      TOT1 = TOT1 - ZZZ - ZZZ
      TOT2 = TOT2 - TWO * ZZZ * X(ISAVE)
      D(ISAVE) = - ZZZ
      J = ISAVE
      GO TO 130
  140 ZZZ = SIGN(1, INEXT(IBAS1))
      TOT1 = TOT1 - RHO - ZZZ
      TOT2 = TOT2 - RHO * BBB - ZZZ * X(IIN)
      D(IOUT) = - RHO
      ITER = ITER + 1
      IF (IFLAG .EQ. 1) GO TO 150
      X(IBAS2) = A2
      IBAS2 = IIN
      D(IBAS2) = - ONE
      A2 = X(IIN)
      Y2 = Y(IIN)
      DET = ONE / (A1 - A2)
      AAAA = (A1 * Y2 - A2 * Y1) * DET
      BBBB = (Y1 - Y2) * DET
      GO TO 60
  150 IBAS1 = IIN
      A1 = X(IIN)
      D(IBAS1) = ZERO
      Y1 = Y(IIN)
      DET = ONE / (A2 - A1)
      AAAA = (A2 * Y1 - A1 * Y2) * DET
      BBBB = (Y2 - Y1) * DET
      GO TO 50
C
C     Calculate optimal sum of absolute deviations
C
  160 SAD = ZERO
      DO 170 I = 1, N
	D(I) = Y(I) - (AAAA + BBBB * X(I))
	SAD = SAD + ABS(D(I))
  170 CONTINUE
      ALPHA = AAAA
      BETA = BBBB
C
      RETURN
      END
