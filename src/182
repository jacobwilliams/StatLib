      SUBROUTINE FORKAL(IP, IQ, IR, NP, IRD, IRZ, ID, IL, N, NRBAR,
     *   PHI, THETA, DELTA, W, Y, AMSE, A, P, V, RESID, E, XNEXT, XROW,
     *   RBAR, THETAB, STORE, IFAULT)
C
C     ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2
C
C     Finite sample prediction from ARIMA processes.
C
C     Auxiliary routines required: KARMA & STARMA from AS 154 and
C     routines called by them: INCLU2 from ASR 17 (a slight variant on
C     AS 75, and REGRES from AS 75.
C
      REAL PHI(IR), THETA(IR), DELTA(ID), W(N), Y(IL), AMSE(IL),
     *   A(IRD), P(IRZ), V(NP), RESID(N), E(IR), XNEXT(NP), XROW(NP),
     *   RBAR(NRBAR), THETAB(NP), STORE(IRD)
      REAL ZERO, ONE, TWO
      DATA ZERO/0.0/, ONE/1.0/, TWO/2.0/
C
C     Invoking this routine will calculate the finite sample predictions
C     and their conditional mean square errors for any ARIMA process.
C
C     Check for input faults.
C
      IFAULT = 0
      IF (IP .LT. 0) IFAULT = 1
      IF (IQ .LT. 0) IFAULT = IFAULT + 2
      IF (IP * IP + IQ * IQ .EQ. 0) IFAULT = 4
      K = IQ + 1
      IF (K .LT. IP) K = IP
      IF (IR .NE. K) IFAULT = 5
      IF (NP .NE. IR * (IR + 1) / 2) IFAULT = 6
      IF (NRBAR .NE. NP * (NP - 1) / 2) IFAULT = 7
      IF (ID .LT. 0) IFAULT = 8
      IF (IRD .NE. IR + ID) IFAULT = 9
      IF (IRZ .NE. IRD * (IRD + 1) / 2) IFAULT = 10
      IF (IL .LT. 1) IFAULT = 11
      IF (IFAULT .NE. 0) RETURN
C
C     Calculate initial conditions for Kalman filter
C
      A(1) = ZERO
      V(1) = ONE
      IF (NP .EQ. 1) GO TO 130
      DO 100 I = 2, NP
  100 V(I) = ZERO
      IF (IQ .EQ. 0) GO TO 130
      IQ1 = IQ + 1
      DO 110 I = 2, IQ1
  110 V(I) = THETA(I-1)
      DO 120 J = 1, IQ
	LL = J * (2*IR + 1 - J) / 2
	DO 120 I = J, IQ
	  LLI = LL + I
	  V(LLI) = THETA(I) * THETA(J)
  120 CONTINUE
C
C     Find initial likelihood conditions.
C     IFAULT not tested on exit from STARMA as all possible errors
C     have been checked above.
C
  130 IF (IR .EQ. 1) P(1) = ONE / (ONE - PHI(1) * PHI(1))
      IF (IR .NE. 1) CALL STARMA(IP, IQ, IR, NP, PHI, THETA, A, P, V,
     *  THETAB, XNEXT, XROW, RBAR, NRBAR, IFAULT)
C
C     Calculate data transformations
C
      NT = N - ID
      IF (ID .EQ. 0) GO TO 170
      DO 140 J = 1, ID
	NJ = N - J
	STORE(J) = W(NJ)
  140 CONTINUE
      DO 160 I = 1, NT
	AA = ZERO
	DO 150 K = 1, ID
	  IDK = ID + I - K
	  AA = AA - DELTA(K) * W(IDK)
  150   CONTINUE
	IID = I + ID
	W(I) = W(IID) + AA
  160 CONTINUE
C
C     Evaluate likelihood to obtain final KF conditions
C
  170 SUMLOG = ZERO
      SSQ = ZERO
      IUPD = 1
      DEL = - ONE
      NIT = 0
      CALL KARMA(IP, IQ, IR, NP, PHI, THETA, A, P, V, NT, W, RESID,
     *    SUMLOG, SSQ, IUPD, DEL, E, NIT)
C
C     Calculate M.L.E. of sigma squared
C
      SIGMA = ZERO
      DO 200 J = 1, NT
  200 SIGMA = SIGMA + RESID(J)**2
      SIGMA = SIGMA / NT
C
C     Reset the initial A and P when differencing occurs
C
      IF (ID .EQ. 0) GO TO 250
      DO 210 I = 1, NP
  210 XROW(I) = P(I)
      DO 220 I = 1, IRZ
  220 P(I) = ZERO
      IND = 0
      DO 230 J = 1, IR
	K = (J-1) * (ID + IR + 1) - (J-1) * J / 2
	DO 230 I = J, IR
	  IND = IND + 1
	  K = K + 1
	  P(K) = XROW(IND)
  230 CONTINUE
      DO 240 J = 1, ID
	IRJ = IR + J
	A(IRJ) = STORE(J)
  240 CONTINUE
C
C     Set up constants
C
  250 IR2 = IR + 1
      IR1 = IR - 1
      ID1 = ID - 1
      ID2R = 2 * IRD
      ID2R1 = ID2R - 1
      IDD1 = 2 * ID + 1
      IDD2 = IDD1 + 1
      I45 = ID2R + 1
      IDRR1 = IRD + 1
      IDDR = 2 * ID + IR
      JKL = IR * (IDDR + 1) / 2
      JKL1 = JKL + 1
      ID2R2 = ID2R + 2
      IBC = IR * (I45 - IR) / 2
      DO 560 L = 1, IL
C
C     Predict A
C
	A1 = A(1)
	IF (IR .EQ. 1) GO TO 310
	DO 300 I = 1, IR1
  300   A(I) = A(I+1)
  310   A(IR) = ZERO
	IF (IP .EQ. 0) GO TO 330
	DO 320 J = 1, IP
  320   A(J) = A(J) + PHI(J) * A1
  330   IF (ID .EQ. 0) GO TO 360
	DO 340 J = 1, ID
	  IRJ = IR + J
	  A1 = A1 + DELTA(J) * A(IRJ)
  340   CONTINUE
	IF (ID .LT. 2) GO TO 360
	DO 350 I = 1, ID1
	  IRI1 = IRD - I
	  A(IRI1 + 1) = A(IRI1)
  350   CONTINUE
  360   A(IR2) = A1
C
C     Predict P
C
	IF (ID .EQ. 0) GO TO 480
	DO 370 I = 1, ID
	  STORE(I) = ZERO
	  DO 370 J = 1, ID
	    LL = MAX(I,J)
	    K = MIN(I,J)
	    JJ = JKL + (LL - K) + 1 + (K-1) * (IDD2 - K) / 2
	    STORE(I) = STORE(I) + DELTA(J) * P(JJ)
  370   CONTINUE
	IF (ID .EQ. 1) GO TO 400
	DO 380 J = 1, ID1
	  JJ = ID - J
	  LK = (JJ-1) * (IDD2 - JJ) / 2 + JKL
	  LK1 = JJ * (IDD1 - JJ) / 2 + JKL
	  DO 380 I = 1, J
	    LK = LK + 1
	    LK1 = LK1 + 1
	    P(LK1) = P(LK)
  380   CONTINUE
        DO 390 J = 1, ID1
	  JKLJ = JKL1 + J
	  IRJ = IR + J
	  P(JKLJ) = STORE(J) + P(IRJ)
  390   CONTINUE
  400   P(JKL1) = P(1)
	DO 410 I = 1, ID
	  IRI = IR + I
	  P(JKL1) = P(JKL1) + DELTA(I) * (STORE(I) + TWO * P(IRI))
  410   CONTINUE
	DO 420 I = 1, ID
	  IRI = IR + I
	  STORE(I) = P(IRI)
  420   CONTINUE
	DO 430 J = 1, IR
	  KK1 = J * (ID2R1 - J) / 2 + IR
	  K1 = (J-1) * (ID2R - J) / 2 + IR
	  DO 430 I = 1, ID
	    KK = KK1 + I
	    K = K1 + I
	    P(K) = PHI(J) * STORE(I)
	    IF (J .NE. IR) P(K) = P(K) + P(KK)
  430   CONTINUE
C
	DO 440 J = 1, IR
	  STORE(J) = ZERO
	  KKK = J * (I45 - J) / 2 - ID
	  DO 440 I = 1, ID
	    KKK = KKK + 1
	    STORE(J) = STORE(J) + DELTA(I) * P(KKK)
  440   CONTINUE
	IF (ID .EQ. 1) GO TO 460
	DO 450 J = 1, IR
	  K = J * IDRR1 - J * (J+1) / 2 + 1
	  DO 450 I = 1, ID1
	    K = K - 1
	    P(K) = P(K-1)
  450   CONTINUE
  460   DO 470 J = 1, IR
          K = (J-1) * (ID2R - J) / 2 + IR + 1
	  P(K) = STORE(J) + PHI(J) * P(1)
	  IF (J .LT. IR) P(K) = P(K) + P(J+1)
  470   CONTINUE
  480   DO 490 I = 1, IR
  490   STORE(I) = P(I)
C
	IND = 0
	DT = P(1)
	DO 500 J = 1, IR
	  PHIJ = PHI(J)
	  PHIJDT = PHIJ * DT
	  IND2 = (J-1) * (ID2R2 - J) / 2
	  IND1 = J * (I45 - J) / 2
	  DO 500 I = J, IR
	    IND = IND + 1
	    IND2 = IND2 + 1
	    PHII = PHI(I)
	    P(IND2) = V(IND) + PHII * PHIJDT
	    IF (J .LT. IR) P(IND2) = P(IND2) + STORE(J+1) * PHII
	    IF (I .EQ. IR) GO TO 500
	    IND1 = IND1 + 1
	    P(IND2) = P(IND2) + STORE(I+1) * PHIJ + P(IND1)
  500   CONTINUE
C
C     Predict Y
C
        Y(L) = A(1)
	IF (ID .EQ. 0) GO TO 520
	DO 510 J = 1, ID
	  IRJ = IR + J
	  Y(L) = Y(L) + A(IRJ) * DELTA(J)
  510   CONTINUE
C
C     Calculate M.S.E. of Y
C
  520   AMS = P(1)
	IF (ID .EQ. 0) GO TO 550
	DO 530 J = 1, ID
	  JRJ = IBC + (J-1) * (IDD2 - J) / 2
	  IRJ = IR + J
	  AMS = AMS + TWO * DELTA(J) * P(IRJ) + P(JRJ+1) * DELTA(J)**2
  530   CONTINUE
	IF (ID .EQ. 1) GO TO 550
	DO 540 J = 1, ID1
	  J1 = J + 1
	  JRK = IBC + 1 + (J-1) * (IDD2 - J) / 2
	  DO 540 I = J1, ID
	    JRK = JRK + 1
	    AMS = AMS + TWO * DELTA(I) * DELTA(J) * P(JRK)
  540   CONTINUE
  550   AMSE(L) = AMS * SIGMA
  560 CONTINUE
C
      RETURN
      END
