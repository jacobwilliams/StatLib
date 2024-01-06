      DOUBLE PRECISION FUNCTION ROY TST(NN1, NN2, NNP, X, IFAULT)
C
C     ALGORITHM AS 77  APPL. STATIST. (1974) VOL.23, NO.3
C
C     A function for the exact null distribution of the largest root
C     of a beta matrix.
C
C     Restrictions:
C       The value of NN2 - NNP - 1 must be an even number.
C       Local storage is sufficient for min(NN1, NNP) at most 10.
C       The number of terms is limited to at most 150,000.
C
      INTEGER ND(10), NP(10), K(56), R, P, PMAX, TERMS
      DOUBLE PRECISION T(10), S(10), A, B, C, CX, CXS, HALF, HCX, ONE,
     +	P1, P2, TL, TWO, X, XM, XMM, XMN2, XN1, XN11, ZERO
      LOGICAL C1, C2
      DATA ZERO, HALF, ONE, TWO /0.D0, 0.5D0, 1.D0, 2.D0/
      DATA PMAX, TERMS /10, 150000/
C
C     Check for errors in the parameters.
C     Commence initialisation.
C
      N1 = MAX(NN1, NNP)
      N2 = MIN(NN2, NN1 + NN2 - NNP)
      P = MIN(NN1, NNP)
      IFAULT = 0
      ROY TST = ZERO
      IF (N1 .GT. 0 .AND. N2 .GT. 0 .AND. P .GT. 0 .AND. X .GE. ZERO
     +	.AND. X .LE. ONE) GO TO 10
      IFAULT = 1
      RETURN
   10 M = (N2 - P - 1) / 2
      IF (M .GE. 0 .AND. M + M .EQ. N2 - P - 1) GO TO 20
      IFAULT = 2
      RETURN
   20 P1 = ONE
      P2 = ZERO
      IF (M .EQ. 0) GO TO 160
      A = ONE
      B = M + P
      R = MIN(M, P)
      C = R
      DO 30 J = 1, R
	A = A * B / C
	B = B - ONE
	C = C - ONE
   30 CONTINUE
      R = A
      IF (P .LE. PMAX .AND. R .LE. TERMS) GO TO 40
      IFAULT = 3
      RETURN
C
C     Error checks complete; initialisation continues.
C
   40 XM = P
      XN11 = N1 - 1
      CX = ONE - X
      P2 = ONE
      DO 50 I = 1, P
	A = -I
	DO 50 J = 1, M
	  A = A + TWO
	  P2 = P2 * CX * (XN11 + A) / (XM + A)
   50 CONTINUE
      S(1) = P2
      MP = M * P
      MN = MP - 2
      IF (MN .LT. 0) GO TO 160
      R = 0
      L = 0
      N = 1
      LE = 1
      K(1) = M
      XMM = P + P + 1 - N1 - N2
      XMN2 = P - N2
      XN1 = N1
      HCX = HALF * CX
      CXS = ONE / (CX * CX)
      T(1) = ONE
      GO TO 80
C
C     Initialisation finished.   Main loop starts here.
C
C     Check if the current partition cannot have a unit part appended.
C
   60 C1 = (R .EQ. P)
C
C     Check if the current partition cannot have its last part
C     increased.
C
      C2 = (K(LE) .EQ. K(LE-1))
      IF (C1 .AND. C2) GO TO 150
      IF (C1) GO TO 100
      IF (C2) GO TO 90
      LB = LE - R + 1
      DO 70 I = LB, LE
	J = I + R
	K(J) = K(I)
   70 CONTINUE
      S(L+1) = S(L)
      T(L+1) = T(L)
C
C     Signal that a second partition is waiting in the list.
C
      ASSIGN 80 TO NEXT
      GO TO 100
   80 ASSIGN 140 TO NEXT
      L = L + 1
      LE = LE + R
   90 LE = LE + 1
      R = R + 1
      NP(L) = R
      K(LE) = 0
  100 K(LE) = K(LE) + 1
      ND(L) = N
C
C     Use the recurrence formulae for the terms of the current partition
C
      B = 2 * K(LE) - R - 1
      A = K(LE)
      TL = T(L)
      TL = TL * HCX * (XN1 + B) * (XM + B) / (A * (A + A - ONE))
      JP = R - 1
      IF (JP .EQ. 0) GO TO 130
      LB = LE - JP
      C = R + 2
      DO 120 J = 1, JP
	C = C - ONE
	TL = TL * (ONE - TWO / (TWO * (K(LB) - A) + C))
	LB = LB + 1
  120 CONTINUE
  130 P1 = P1 + TL
      T(L) = TL
C
C     Check if conjugate partition terms are necessary.
C
      IF (MN .EQ. 0) GO TO NEXT, (80, 140)
      B = B + ONE
      S(L) = S(L) * CXS * (XMN2 + B) * (XM + B) / ((XMM + B) *
     +			(XN11 + B))
      P2 = P2 + TL * S(L)
C
C     Check if a second partition has been generated.
C
      GO TO NEXT, (80, 140)
  140 N = ND(L) + 1
      R = NP(L)
C
C     Check that the degree of the partition is within limits.
C
      MN = MP - N - N
      IF (MN .GE. 0) GO TO 60
  150 L = L - 1
      LE = LE - R
C
C     If L = 0 then list is empty and calculation complete.
C
      IF (L .GT. 0) GO TO 140
C
C     Multiply by the external factor, and return
C
  160 ROY TST = (P1 + P2) * X ** (HALF * (N1 * P))
      RETURN
      END

