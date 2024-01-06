C UKC NETLIB DISTRIBUTION COPYRIGHT 1990 RSS
C
      FUNCTION BMIX(X, N, L, Q, ITYPE, IFAULT)
C
C        Algorithm AS 123 (Applied Statistics, 27 (1978), pp. 104-109);
C
C        CALCULATES C.D.F.S OF MIXTURES OF BETA DISTRIBUTIONS.
C
      DIMENSION Q(101)
      LOGICAL SEC, ALT
      PARAMETER (C1=0.636619772368, C2=76.7625)
C
C        C1 = 2.0/PI.  C2 MUST NOT EXCEED LOG (BASE 10) OF LARGEST
C        POSSIBLE FLOATING POINT NUMBER.
C
      BMIX = 0.0
C
C        CHECK FOR ERRORS AND INITIALISE.
C
      IFAULT = 1
      I = L + ITYPE
      IF (N .LE. 0 .OR. N .GE. 101 .OR. I .LE. 1) RETURN
      IFAULT = 2
      XF1 = 1.0 - X
      IF (X .LT. 0.0 .OR. XF1 .LE. 0.0) RETURN
      RX = SQRT(X)
      RXF1 = SQRT(XF1)
      IFAULT = 3
      IF (ITYPE .EQ. 2) GO TO 40
      IF (ITYPE .NE. 1) RETURN
      M = L - 2
      D = M
      E = 1.0
      F = L
      C = (F - 1.0)*XF1
      Z = X
      GO TO 50
   40 I = L + N
      M = I - 3
      D = I
      E = -1.0
      F = I - 2
      C = F
      Z = X/XF1
   50 IFAULT = 4
      NN = N
      ALT = .FALSE.
      K = -2
      J = N + 1
C
C        CHECK FOR POSSIBILITY OF OVERFLOW.
C
   60 A = -F*ALOG10(RXF1)
      IF (A .LT. C2) GO TO 90
      IF (ITYPE .EQ. 1 .OR. ALT) RETURN
      IF (I .LT. 101) GO TO 70
      IFAULT = 5
      RETURN
C
C        IF ITYPE = 2 TRY ALTERNATIVE APPROACH.
C
   70 ALT = .TRUE.
      NN = I
      Z = XF1/X
      XF1 = X
      A = RX
      RX = RXF1
      RXF1 = A
      IF (L .EQ. 0) GO TO 85
      J = N + 2
      K = I + 1
      DO 80 I = J, K
   80 Q(I) = 0.0
   85 J = 1
      K = 2
      GO TO 60
C
C        MAIN SECTION OF ALGORITHM BEGINS HERE.
C
   90 IFAULT = 0
      SEC = .FALSE.
      I = NN
  100 A = 0.0
      B = 0.0
C
C        INNER LOOP.
C
  110 B = B + Q(J)
      IF (I .LE. 1) GO TO 120
      FI = I
      A = A*Z*(D/FI + E) - B
      I = I - 2
      J = J + K
      GO TO 110
  120 IF (.NOT. ALT) GO TO 130
      BMIX = BMIX + B
      A = -A
      B = -B
  130 IF (I .EQ. 1) GO TO 140
C
C        COMPLETION OF EVEN CALCULATIONS.
C
      BMIX = BMIX + B + A*RXF1**F
      GO TO 180
  140 IF (M .GE. 0) A = B + C*A
      I = M
C
C        OUTER LOOP.
C
  150 IF (I .LE. 1) GO TO 160
      FI = I
      A = (A - A/FI)*XF1 + B
      I = I - 2
      GO TO 150
C
C        COMPLETION OF ODD CALCULATIONS.
C
  160 A = A*RX
      IF (I .EQ. 0 .OR. I .EQ. -2) GO TO 170
      BMIX = BMIX + C1*(ATAN(RX/RXF1)*B + RXF1*A)
      GO TO 180
  170 BMIX = BMIX + A
  180 IF (SEC) RETURN
      SEC = .TRUE.
      I = NN - 1
      J = N
      IF (ALT) J = 2
      GO TO 100
      END
