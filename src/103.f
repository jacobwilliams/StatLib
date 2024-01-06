      REAL FUNCTION DIGAMA(X, IFAULT)
C
C     ALGORITHM AS 103  APPL. STATIST. (1976) VOL.25, NO.3
C
C     Calculates DIGAMMA(X) = D( LOG( GAMMA(X))) / DX
C
      REAL ZERO, HALF, ONE
C
C     Set constants, SN = Nth Stirling coefficient, D1 = DIGAMMA(1.0)
C
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/
      DATA S, C, S3, S4, S5, D1 /1.E-05, 8.5, 8.333333333E-02,
     *    8.3333333333E-03, 3.96825 3968E-03, -0.57721 56649/
C
C     Check argument is positive
C
      DIGAMA = ZERO
      Y = X
      IFAULT = 1
      IF (Y .LE. ZERO) RETURN
      IFAULT = 0
C
C     Use approximation if argument <= S
C
      IF (Y .LE. S) THEN
	DIGAMA = D1 - ONE / Y
	RETURN
      END IF
C
C     Reduce to DIGAMA(X + N) where (X + N) >= C
C
    1 IF (Y .GE. C) GO TO 2
      DIGAMA = DIGAMA - ONE/Y
      Y = Y + ONE
      GO TO 1
C
C     Use Stirling's (actually de Moivre's) expansion if argument > C
C
    2 R = ONE / Y
      DIGAMA = DIGAMA + LOG(Y) - HALF*R
      R = R * R
      DIGAMA = DIGAMA - R*(S3 - R*(S4 - R*S5))
      RETURN
      END

