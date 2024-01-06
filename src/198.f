      SUBROUTINE PROBS(K, W, P, IFAULT)
C
C        ALGORITHM AS 198 APPL. STATIST. (1984) VOL.33,NO.1
C
C        CALCULATION OF THE PROBABILITIES P(L,K) FOR THE CASE OF
C        SIMPLE ORDER.
C
C     Auxiliary functions required: PR1 & F1 from AS 158
C
      DIMENSION W(20), P(20), Q(20, 20), CH(20, 20)
      REAL ZERO, HALF
      DATA C1 /1.0E-6/, ZERO/0.0/, HALF/0.5/
C
C        CHECK THAT K .GE. 3 AND .LE. 20
C
      IFAULT = 2
      IF (K .LT. 3 .OR. K .GT. 20) RETURN
C
C        CHECK THAT THE WEIGHTS ARE POSITIVE
C
      IFAULT = 1
      DO 1 I = 1, K
      IF (W(I) .LE. ZERO) RETURN
    1 CONTINUE
C
C        CHECK WEIGHTS FOR EQUALITY
C
      IFAULT = 0
      WW = W(1)
      DO 2 I = 2, K
        IF (ABS(WW - W(I)) .GT. C1) GOTO 7
    2 CONTINUE
C
C        EQUAL WEIGHTS
C
      CALL CHASE(K, Q, .TRUE.)
      DO 3 J = 1, K
    3 P(J) = A(J, K)
      RETURN
C
C        UNEQUAL WEIGHTS
C
    7 IF (K .GT. 5) GOTO 11
      K2 = K - 2
      GOTO (8, 9, 10), K2
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
   11 CALL CHASE(K, CH, .FALSE.)
      CALL CHASE(K, Q, .TRUE.)
      CALL PAPRX(K, W, CH, Q, P)
      RETURN
      END
C
      SUBROUTINE PAPRX(K, W, CH, P, PA)
C
C        ALGORITHM AS 198.1 APPL. STATIST. (1984) VOL.33, NO.1
C
C        THIS SUBROUNTINE COMPUTES THE APPROXIMATE PLK-S.
C
      REAL W(20), CH(20, 20), P(20, 20), PA(20), PB(20), PBP(20),
     *  SUMS(2), ZERO, ONE, THREE, PT65, PT35
      INTEGER INDEX(20), A, B
      DATA ZERO/0.0/, ONE/1.0/, THREE/3.0/, PT65/0.65/, PT35/0.35/
C
      ALPHA = ONE / THREE
C
C        INITIALISE PA,PB,PBP
C
      DO 10 I = 1, 20
        PA(I) = ZERO
        PB(I) = ZERO
        PBP(I) = ZERO
   10 CONTINUE
C
C         DETERMINE MASIMUM AND MINIMUM OF WEIGHTS
C
      WMAX = W(1)
      WMIN = W(1)
      DO 20 I = 2, K
        WW = W(I)
        IF (WW .LT. WMIN) WMIN = WW
        IF (WW .LT. WMAX) WMAX = WW
   20 CONTINUE
      CUT = PT65 * WMIN + PT35 * WMAX
C
C        DETERMINE THE INDICES OF THE WEIGHTS AND THE NUMBER OF
C        LARGE WEIGHTS (M)
C
      M = 0
      DO 40 I = 1, K
        IF (W(I) .LT. CUT) GOTO 30
        INDEX(I) = 1
        M = M + 1
        GOTO 40
   30   INDEX(I) = 0
   40 CONTINUE
C
C        IF A=0 AND B=0 SET PB=PLM
C
      N = INDEX(1) + INDEX(K)
      IF (N .NE. 2) GO TO 90
      DO 80 L = 1, M
   80 PB(L) = P(L, M)
      GOTO 240
C
C        DETERMINE A,B
C
   90 A = 1
      B = 1
      DO 100 I = 1, K
        IF (INDEX(I) .NE. 0) GOTO 110
        A = A + 1
  100 CONTINUE
  110 DO 120 I = 1, K
        J = K - I + 1
        IF (INDEX(J) .NE. 0) GOTO 130
        B = B + 1
  120 CONTINUE
  130 IF (N .EQ. 0) GOTO 190
      N = MAX0(A, B)
      DO 180 L = 1, K
        SUM = ZERO
        DO 170 I = 1, L
          J = L - I + 1
          SUM = SUM + P(I, M) * CH(J, N)
  170   CONTINUE
        PB(L) = SUM
  180 CONTINUE
      GOTO 240
  190 DO 210 L = 1, K
        SUM = ZERO
        DO 200 I = 1, L
          J = L - I + 1
          SUM = SUM + CH(I, A) * CH(J, B)
  200   CONTINUE
        PBP(L) = SUM
  210 CONTINUE
      DO 230 L = 1, K
        SUM = ZERO
        DO 220 I = 1, L
          J = L - I + 1
          SUM = SUM + PBP(I) * P(J, M)
  220   CONTINUE
        PB(L) = SUM
  230 CONTINUE
C
C        DETERMINE R
C
  240 SUMS(1) = ZERO
      SUMS(2) = ZERO
      DO 250 I = 1, K
        J = INDEX(I) + 1
        SUMS(J) = SUMS(J) + W(I)
  250 CONTINUE
      X = K - M
      Y = M
      R = SUMS(2) * X / (SUMS(1) * Y)
      RI = ONE / R ** ALPHA
      DO 260 L = 1, K
  260 PA(L) = (ONE - RI) * PB(L) + RI * P(L, K)
      RETURN
      END
C
      SUBROUTINE CHASE(K, CH, EQUAL)
C
C        ALGORITHM AS 198.2 APPL. STATIST. (1984) VOL.33, NO.1
C
C        THIS SUBROUTINE COMPUTES CHASE PLK-S IF EQUAL=.FALSE.
C        OR EQUAL WEIGHT PROBABILITIES IF EQUAL=.TRUE.
C
      LOGICAL EQUAL
      REAL ZERO, HALF, ONE, CH(20, 20)
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/
C
      DO 20 J = 1, K
        DO 10 I = 1, K
   10   CH(I, J) = ZERO
   20 CONTINUE
      CH(1, 1) = ONE
      CH(1, 2) = HALF
      CH(2, 2) = HALF
      KM = K - 1
      DO 60 J = 1, KM
        J1 = J + 1
        IF (EQUAL) GOTO 30
        X = J + J - 1
        Y = J + J
        GOTO 40
   30   X = J
        Y = J1
   40   Y1 = ONE / Y
        CH(1, J1) = X * Y1 * CH(1, J)
        CH(J1, J1) = Y1 * CH(J, J)
        DO 50 I = 2, J
          CH(I, J1)  = Y1 * CH(I - 1, J) + X * Y1 * CH(I, J)
   50   CONTINUE
   60 CONTINUE
      RETURN
      END
