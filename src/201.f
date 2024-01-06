      SUBROUTINE COMICS(ALP, NP, P, TOL, PI, CPI, CP, SP, IFAULT)
C
C     ALGORITHM AS 201  APPL. STATIST. (1984) VOL.33, NO.2
C
C     Subroutine to combine predictions about a statistic based on the
C     orderings of a set of means with an F-test to test for differences
C     between the means.
C
      INTEGER NP, IFAULT, I, J
      REAL ALP, P(NP), TOL, PI(NP), CPI(NP), CP(NP), SP(NP), ZERO, ONE,
     *     CUM, PIJ, X, XALP, ABS
C
      DATA ZERO /0.0/, ONE /1.0/
C
C     Test ALP
C
      IFAULT = 1
      IF (ALP .LE. ZERO .OR. ALP .GE. ONE) RETURN
C
C     Compute cumulative distribution of statistic
C
      IFAULT = 2
      CP(1) = P(1)
      IF (P(1) .LE. ZERO) RETURN
      DO 10 I = 2, NP
	IF (P(I) .LE. ZERO) RETURN
	CP(I) = CP(I-1) + P(I)
   10 CONTINUE
C
C     Check that distribution sums to one
C
      IFAULT = 3
      IF (ABS(CP(NP) - ONE) .GT. TOL) RETURN
      IFAULT = 0
C
C     Compute distribution of combined likelihoods (PI) and cumulative
C     combined distribution (CPI)
C
      CUM = ZERO
      DO 30 J = 1, NP
	PIJ = ZERO
	DO 20 I = J, NP
   20   PIJ = PIJ + P(I) / CP(I)
	PI(J) = P(J) * PIJ
	CUM = CUM + PI(J)
	CPI(J) = CUM
   30 CONTINUE
C
C     Determine criteria for significance on combined dist X
C
      DO 40 I = 1, NP
	IF (CPI(I) .GT. ALP) GO TO 50
   40 CONTINUE
      I = NP
   50 XALP = CPI(I) - ALP
      X = -XALP * P(I) / PI(I) + CP(I)
C
C     Determine significance level required of F
C
      DO 60 I = 1, NP
	SP(I) = X / CP(I)
	IF (SP(I) .GT. ONE) SP(I) = ONE
   60 CONTINUE
C
      RETURN
      END
