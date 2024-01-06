      SUBROUTINE PROBS(K, W, P, IFAULT)
C
C          ALGORITHM AS 158 APPL. STATIST. (1981) VOL.30, NO.1
C
C          CALCULATION OF THE PROBABILITIES P(L,K) FOR
C          THE CASE OF SIMPLE ORDER,
C          FOR EQUAL WEIGHTS, K .GE. 3 AND .LE. 10 .
C          FOR UNEQUAL WEIGHTS, K .GE. 3 AND .LE. 6
C
      DIMENSION W(10), P(10), Q(10, 10)
      REAL ZERO, HALF, ONE, THREE, SIX
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/, THREE/3.0/,
     *     SIX/6.0/, C1 /1.0E-6/
C
C        CHECK THAT WEIGHTS ARE POSITIVE
C
      IFAULT = 0
      DO 1 I = 1, K
	IF (W(I) .LE. ZERO) GOTO 101
    1 CONTINUE
C
C        CHECK THAT K .GE. 3 AND .LE. 10
C
      IF (K .LT. 2 .OR. K .GT. 10) GOTO 102
      WW = W(1)
      DO 2 I = 2, K
	IF (ABS(WW - W(I)) .GT. C1) GOTO 7
    2 CONTINUE
C
C        EQUAL WEIGHTS
C
      Q(1, 3) = ONE / THREE
      Q(2, 3) = HALF
      Q(3, 3) = ONE / SIX
      IF (K .EQ. 3) GOTO 5
      DO 4 J = 4, K
	AJ = J
	A1 = ONE /AJ
	A2 = (AJ - ONE) * A1
	Q(1, J) = A1
	J1 = J - 1
	DO 3 L = 2, J1
	  L1 = L- 1
    3   Q(L, J) = A1 * Q(L1, J1) + A2 * Q(L, J1)
	Q(J, J) = ONE / FACT(J, IFAULT)
	IF (IFAULT .NE. 0) RETURN
    4 CONTINUE
    5 CONTINUE
      DO 6 J = 1, K
    6 P(J) = Q(J, K)
      RETURN
C
C        UNEQUAL WEIGHTS - CHECK THAT K .LE. 6
C
    7 IF (K .GT. 6) GOTO 102
      K2 = K - 2
      GOTO (8, 9, 10, 11), K2
    8 P(1) = PR1(1, 3, W)
      P(2) = HALF
      P(3) = HALF - P(1)
      RETURN
    9 P(1) = PR1(1, 4, W)
      P(4) = PR1(4, 4, W)
      P(2) = HALF - P(4)
      P(3) = HALF - P(1)
      RETURN
   10 P(5) = PR1(5, 5, W)
      P(4) = PR1(4, 5, W)
      P(2) = HALF - P(4)
      P(1) = PR1(1, 5, W)
      P(3) = HALF - P(1) - P(5)
      RETURN
   11 SUM = ZERO
      DO 12 I = 2, 6
	PP = PR2(I, 2)
	P(I) = PP
	SUM = SUM + PP
   12 CONTINUE
      P(1) = ONE - SUM
      IF (P(1) .LT. ZERO) P(1) = ZERO
      RETURN
  101 IFAULT = 1
      RETURN
  102 IFAULT = 2
      RETURN
      END
C
      FUNCTION PR1(I, J, W)
C
C        ALGORITHM AS 158.1 APPL. STATIST. (1981) VOL.30, NO.1
C
C        EXPLICIT CALCULATION OF PROBABILITIES FOR K .LE. 5
C        ALSO CALLED BY FUNCTION F2
C
      DIMENSION W(10)
      REAL ZERO, HALF, QTR, EIGHTH, PT0625, PT375, PII
      DATA ZERO/0.0/, HALF/0.5/, QTR/0.25/, EIGHTH/0.125/,
     *     PT0625/0.0625/, PT375/0.375/, PII/0.318309886/
C
      IF (J .NE. 3) GOTO 40
      C = HALF * PII * F1(W(1), W(2), W(3))
      IF (I .EQ. 3) GOTO 30
      PR1 = QTR - C
      RETURN
   30 PR1 = QTR + C
      RETURN
   40 W1 = W(1)
      W2 = W(2)
      W3 = W(3)
      W4 = W(4)
      W12 = W1 + W2
      W23 = W2 + W3
      W34 = W3 + W4
      S12 = F1(W1, W2, W3)
      S23 = F1(W2, W3, W4)
      IF (J .EQ. 5) GOTO 50
      IF (I .EQ. 4) GOTO 41
      C1 = QTR * PII *
     *   (F1(W1, W2, W34) + F1(W1, W23, W4) + F1(W12, W3, W4))
      PR1 = EIGHTH - C1
      RETURN
   41 C2 = QTR *PII * (S12 + S23)
      PR1 = EIGHTH + C2
      RETURN
   50 W5 = W(5)
      W45 = W4 + W5
      W123 = W12 + W3
      W234 = W23 + W4
      W345 = W34 + W5
      S34 = F1(W3, W4, W5)
      IF (I .EQ. 4) GOTO 52
      C5 = PT0625 + EIGHTH * PII * (S12 + S23 + S34) +
     *   QTR * PII * PII * S12 * S34
      IF (I .EQ. 1) GOTO 51
      PR1 = C5
      IF (PR1 .LT. ZERO) PR1 = ZERO
      RETURN
   51 S113 = F1(W1, W2, W345)
      S131 = F1(W1, W234, W5)
      S311 = F1(W123, W4, W5)
      C3 = PT375 + EIGHTH * PII * (S113 + S131 +S311 + F1(W1, W23, W45)
     *  + F1(W12, W3, W45) + F1(W12, W34, W5) - S12 - S23 - S34)
     *  - QTR * PII * PII * (S12 * S311 + S23 * S131 + S34 * S113)
      PR1 = HALF - C3 - C5
      RETURN
   52 C2 = EIGHTH * PII * (S12 + S34 + F1(W1, W2, W34) + F1(W12, W3, W4)
     *  + F1(W1, W23, W4) + F1( W2, W3, W45) + F1(W23, W4, W5) +
     *  F1(W2, W34, W5))
      PR1 = QTR + C2
      RETURN
      END
C
C
      FUNCTION F1(V1, V2, V3)
C
C      ALGORITHM AS 158.2 APPL. STATIST. (1981) VOL.30, NO.1
C
      RHO = -SQRT(V1 * V3 / ((V1 + V2) * (V2 + V3)))
      F1 = ASIN(RHO)
      RETURN
      END
C
      FUNCTION PR2(I, W)
C
C        ALGORITHM AS 158.3 APPL. STATIST. (1981) VOL.30, NO.1
C
C        CALCULATION OF PROBABILITIES FOR K .EQ. 6 USING
C        RECURRENCE RELATION
C
      DIMENSION W(10)
      REAL ZERO, HALF, EIGHTH, PT0625, P03125, QTR, PII
      DATA ZERO/0.0/, HALF/0.5/, EIGHTH/0.125/, PT0625/0.0625/,
     *     P03125/0.03125/, QTR/0.25/
      DATA PII /0.318309886/
C
      W1 = W(1)
      W2 = W(2)
      W3 = W(3)
      W4 = W(4)
      W5 = W(5)
      W6 = W(6)
      W12 = W1 + W2
      W23 = W2 + W3
      W34 = W3 + W4
      W45 = W4 + W5
      W56 = W5 + W6
      W123 = W12 + W3
      W234 = W23 + W4
      W345 = W34 + W5
      W456 = W45 + W6
      W1234 = W123 + W4
      W2345 = W234 + W5
      W3456 = W345 + W6
      D4 = ZERO
      D5 = ZERO
      I1 = I - 1
      GOTO (2, 3, 4, 5, 6), I1
    2 PR2 = HALF * (F2(1, 5, W2, W3, W4, W5, W6) +
     *  F2(1, 5, W1, W2, W3, W4, W5) + F2(1, 3, W1, W2, W3, D4, D5) *
     *  F2(1, 3, W4, W5, W6, D4, D5)) + QTR *
     *  (F2(1, 4, W3, W4, W5, W6, D5) + F2(1, 4, W1, W2, W3, W4, D5))
      RETURN
    3 PR2 = F2(3, 3, W1, W2, W3456, D4, D5) *
     *  F2(1, 4, W3, W4, W5, W6, D5) + F2(3, 3, W1, W2345, W6, D4, D5) *
     *  F2(1, 4, W2, W3, W4, W5, D5) + F2(3, 3, W1234, W5, W6, D4, D5) *
     *  F2(1, 4, W1, W2, W3, W4, D5) + HALF *
     *  (F2(3, 3, W1, W23, W456, D4, D5) * F2(1, 3, W4, W5, W6, D4, D5)
     * + F2(3, 3, W1, W234, W56, D4, D5) * F2(1, 3, W2, W3, W4, D4, D5)
     * + F2(3, 3, W12, W3, W456, D4, D5) * F2(1, 3, W4, W5, W6, D4, D5)
     * + F2(3, 3, W12, W345, W6, D4, D5) * F2(1, 3, W3, W4, W5, D4, D5)
     * + F2(3, 3, W123, W4, W56, D4, D5) * F2(1, 3, W1, W2, W3, D4, D5)
     * + F2(3, 3, W123, W45, W6, D4, D5) * F2(1, 3, W1, W2, W3, D4, D5))
     * + EIGHTH * F2(3, 3, W12, W34, W56, D4, D5)
      RETURN
    4 PR2 = F2(4, 4, W1, W2, W3, W456, D5) *
     *  F2(1, 3, W4, W5, W6, D4, D5) + F2(4, 4, W1, W2, W345, W6, D5) *
     *  F2(1, 3, W3, W4, W5, D4, D5) + F2(4, 4, W1, W234, W5, W6, D5) *
     *  F2(1, 3, W2, W3, W4, D4, D5) + F2(4, 4, W123, W4, W5, W6, D5) *
     *  F2(1, 3, W1, W2, W3, D4, D5) + QTR *
     *  (F2(4, 4, W1, W2, W34, W56, D5) + F2(4, 4, W1, W23, W4, W56, D5)
     *  + F2(4, 4, W1, W23, W45, W6, D5) +
     *  F2(4, 4, W12, W3, W4, W56, D5) + F2(4, 4, W12, W3, W45, W6, D5)
     *  + F2(4, 4, W12, W34, W5, W6, D5))
      RETURN
    5 PR2 = HALF * (F2(5, 5, W1, W2, W3, W4, W56) +
     *  F2(5, 5, W1, W2, W3, W45, W6) + F2(5, 5, W1, W2, W34, W5, W6) +
     *  F2(5, 5, W1, W23, W4, W5, W6) + F2(5, 5, W12, W3, W4, W5, W6))
      RETURN
    6 S12 = F1(W1, W2, W3)
      S23 = F1(W2, W3, W4)
      S34 = F1(W3, W4, W5)
      S45 = F1(W4, W5, W6)
      PR2 = P03125 + PT0625 * PII * (S12 + S23 + S34 + S45) +
     *  EIGHTH * PII * PII * (S12 * S34 + S12 * S45 + S23 * S45)
      IF (PR2 .LT. ZERO) PR2 = ZERO
      RETURN
      END
C
      FUNCTION F2(I, J, V1, V2, V3, V4, V5)
C
C        ALGORITHM AS 158.4 APPL. STATIST. (1981) VOL.30, NO.1
C
      DIMENSION VV(10)
      VV(1) = V1
      VV(2) = V2
      VV(3) = V3
      VV(4) = V4
      VV(5) = V5
      F2 = PR1(I, J, VV)
      RETURN
      END
C
      FUNCTION FACT(M, IFAULT)
C
C        ALGORITHM AS 158.5 APPL. STATIST. (1981) VOL.30, NO.1
C
C        CALCULATION OF M FACTORIAL
C
      REAL ONE
      DATA ONE/1.0/, MAXM /50/
C
      IFAULT = 3
      IF (M .LT. 0 .OR. M .GT. MAXM) RETURN
      IFAULT = 0
      FACT = ONE
      IF (M .LE. 1) RETURN
      A1 = ONE
      DO 1 I = 2, M
	AI = I
	A1 = A1 * AI
    1 CONTINUE
      FACT = A1
      RETURN
      END
