      SUBROUTINE QDISX(A, B, M, T, MAX, CZERO, MAX1, P0, P, Y, C0, C,
     +		MEAN, SD, SUM, IFAULT)
C
C     ALGORITHM AS230  APPL. STATIST. (1986) VOL. 36, NO. 3
C
C     Tabulates approximate distribution of customers in multi-server
C     queues having random arrivals and Erlangian service times.
C     Hokstad's method is used.
C
      REAL P(MAX), Y(MAX), C(MAX)
      REAL A, B, C0, CZERO, FI, MEAN, ONE, PM1, PNNN, P0, R, R1, SD,
     +		SUM, TERM, FT, TWO, ZERO, WT
      INTEGER T
      DATA ZERO/0.0/, ONE/1.0/, TWO/2.0/, PNNN/0.999/
C
      IFAULT = 0
      FT = FLOAT(T)
      MAX0 = MAX
      IF (MAX .GT. 4000) MAX0 = 4000
      DO 5 I = 1, MAX0
	P(I) = ZERO
	Y(I) = ZERO
	C(I) = ZERO
    5 CONTINUE
C
C     Calculate P0 to P(M-1)
C
      R = A * B / FLOAT(M)
      IF (R .GT. PNNN) GO TO 10
      P0 = ONE - R
      PM1 = P0
      IF (M .EQ. 1) GO TO 15
      P0 = ONE
      TERM = ONE
      M1 = M - 1
      DO 20 I = 1, M1
	TERM = TERM * A * B / FLOAT(I)
	P0 = P0 + TERM
	P(I) = TERM
   20 CONTINUE
      TERM = TERM * A * B / ((ONE - R) * FLOAT(M))
      P0 = ONE / (P0 + TERM)
      DO 25 I = 1, M1
   25 P(I) = P(I) * P0
      PM1 = P(M1)
   15 IF (T .GT. 100) GO TO 30
C
C     Calculate Y-values for Erlangian service
C
      C0 = (ONE + R / FT) ** FT
      P(M) = (C0 - ONE) * PM1
      Y(1) = (ONE + R / FT) ** FT - R / (ONE + R / FT)
      Y(2) = -(FT * (FT + ONE) / TWO) / (FT / R + ONE) ** 2 
      DO 35 I = 3, MAX0
   35 Y(I) = ZERO
      DO 40 I = 3, MAX0
	FI = FLOAT(I)
	IF (ABS(Y(I-1)) .LT. CZERO) GO TO 45
	Y(I) = (Y(I-1) * (FT + FI - ONE)) / (FI * (FT / R + ONE))
   40 CONTINUE
      GO TO 45
C
C     Calculate Y-values for constant service times
C
   30 C0 = EXP(R)
      P(M) = (C0 - ONE) * PM1
      Y(1) = C0 - R
      Y(2) = -R * R / TWO
      DO 50 I = 3, MAX0
	Y(I) = ZERO
   50 CONTINUE
      DO 55 I = 3, MAX0
	IF (ABS(Y(I-1)) .LT. CZERO) GO TO 45
	Y(I) = Y(I-1) * R / FLOAT(I)
   55 CONTINUE
C
C     Calculate P(M+1) to P(MAX0)
C
   45 R1 = ONE / (ONE - R)
      C(1) = C0 * Y(1)
      DO 60 I = 2, MAX0
	C(I) = ZERO
	I1 = I - 1
	DO 65 J = 1, I1
	  C(I) = C(I) + C(I-J) * Y(J)
   65   CONTINUE
	C(I) = C(I) + C0 * Y(I)
	IF (C(I) .GE. R1 .AND. C(I-1) .GE. R1) GO TO 95
   60 CONTINUE
      I = MAX0
   95 MAX1 = I
      M2 = MAX1 - M
      P(M+1) = (C(1) - C0) * PM1
      DO 70 I = 2, M2
   70 P(M+I) = ZERO
      DO 75 I = 2, M2
        IF (ABS(P(M+I-1)) .LT. CZERO .AND. P(M+I-1) .LT. P(M+I-2))
     +		GO TO 80
        P(M+I) = (C(I) - C(I-1)) * PM1
        IF (P(M+I) .GE. ZERO) GO TO 75
        P(M+I) = ZERO
        GO TO 80
   75 CONTINUE
C
C     Calculate mean, standard deviation and sum of terms.
C     The algorithm for calculating the mean & standard deviation used
C     here is better than that in the journal.
C
   80 SUM = P0
      MEAN = ZERO
      SD = ZERO
      DO 90 I = 1, MAX1
	FI = FLOAT(I)
	SUM = SUM + P(I)
	WT = P(I) / SUM
	TERM = FI - MEAN
	MEAN = MEAN + WT * TERM
        SD = SD + P(I) * TERM * (FI - MEAN)
   90 CONTINUE
      SD = SQRT(SD / SUM)
C
      RETURN
C
   10 IFAULT = -1
      RETURN
      END
