      SUBROUTINE CF(NORD, X, AC, A, D, H, P, DEL, IFAULT)
C
C     ALGORITHM AS 269 APPL.STATIST. (1992), VOL.41, NO.1
C
C     Calculates the Cornish-Fisher adjustment to the normal deviate.
C
      INTEGER NORD, IFAULT
      REAL X, AC(NORD), A(NORD), D(NORD), H(3 * NORD),
     *     P(3 * NORD * (NORD+1)/2), DEL(NORD)
C
C     Local variables
C
      INTEGER J, JA, JAL, JB, JBL, K, L
      REAL AA, BC, CC, DD, FAC, LIMIT, ONE, ZERO
      DATA LIMIT, ONE, ZERO / 3.719017274, 1.0, 0.0 /
C
C     Check input arguments
C
      IFAULT = 0
      IF (NORD .GT. 18) THEN
	IFAULT = 1
      ELSE IF (X .LT. -LIMIT .OR. X .GT. LIMIT) THEN
	IFAULT = 2
      END IF
      IF (IFAULT .NE. 0) RETURN
C
C     Compute the adjusted cumulants
C
      CC = -ONE
      DO 10 J = 1, NORD
	A(J) = CC * AC(J) / ((J+1) * (J+2))
	CC = -CC
   10 CONTINUE
C
C     Compute the Hermite polynomial values.
C
      H(1) = -X
      H(2) = X * X - ONE
      DO 20 J = 3, 3 * NORD
	H(J) = - (X * H(J-1) + (J-1) * H(J-2))
   20 CONTINUE
C
C     Clear the polynomial array.
C
      DO 30 J = 1, 3 * NORD * (NORD+1)/2
	P(J) = ZERO
   30 CONTINUE
      D(1) = - A(1) * H(2)
      DEL(1) = D(1)
      P(1) = D(1)
      P(3) = A(1)
      JA = ZERO
      FAC = ONE
C
C     Main loop
C
      DO 70 J = 2, NORD
C
C     Initialize.
C
	FAC = FAC * J
	JA = JA + 3 * (J-1)
	JB = JA
	BC = ONE
C
C     Calculate coefficients of Hermite polynomials.
C
	DO 50 K = 1, J-1
	  DD = BC * D(K)
	  AA = BC * A(K)
	  JB = JB - 3 * (J - K)
	  DO 40 L = 1, 3 * (J - K)
	    JBL = JB + L
	    JAL = JA + L
	    P(JAL+1) = P(JAL+1) + DD * P(JBL)
	    P(JAL+K+2) = P(JAL+K+2) + AA * P(JBL)
   40     CONTINUE
	  BC = BC * (J - K) / K
   50   CONTINUE
	P(JA+J+2) = P(JA+J+2) + A(J)
C
C     Calculate the adjustments.
C
	D(J) = ZERO
	DO 60 L = 2, 3 * J
	  D(J) = D(J) - P(JA+L) * H(L-1)
   60   CONTINUE
	P(JA+1) = D(J)
	DEL(J) = D(J) / FAC
   70 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CUM(CU, MU, NORD)
C
C     Calculates the cumulants from the moment series by a
C     recursion relation.
C
C     This code was published in Appl. Statist. vol.42 (1993),
C     pp. 268-269.
C
      INTEGER NORD
      REAL CU(NORD), MU(NORD)
C
C     Local variables
C
      INTEGER J, JJ
      REAL ONE, C
      DATA ONE /1.0/
C
      CU(1) = MU(1)
      DO 1 J = 1, NORD - 1
	C = ONE
	CU(J+1) = MU(J+1)
	DO 1 JJ = 0, J-1
	  CU(J+1) = CU(J+1) - C * CU(JJ+1) * MU(J-JJ)
	  C = C * (J-JJ) / (JJ+1)
    1 CONTINUE
      RETURN
      END

