      SUBROUTINE MODLIK(NN, NNP2, Y, DEL, ALPHA, BETA, DELTA1, DELTA2,
     *  NGR, K, INK, NK, NFAIL, H, A, VML, CUMCEN, DAT, IFAULT)
C
C     ALGORITHM AS 202  APPL. STATIST. (1984) VOL.33, NO.2
C
C     MODLIK uses cross-validation to select parameters for the
C     3-parameter nonparametric hazard estimator.
C
      INTEGER NN, NNP2, DEL(NN), NGR, K, INK, NK, NFAIL, CUMCEN(NN),
     *  IFAULT, I, IM, ITEST, J, JM, JM1, JP1, KTH, KTH2, LKTH, MAXIT,
     *  NJ1, NT
      REAL Y(NN), ALPHA, BETA, DELTA1, DELTA2, H(NN,NN), A(NN,NN), VML,
     *  DAT(NNP2), ZERO, EPS2, BIG, ALEPH1, ALEPH2, ALEPH3, BET1, BET2,
     *  BET3, DHAA, FMAB, FMAX, SS, TEST
C
      DATA ZERO, EPS2, BIG /0.0, 1.0E-04, 1.0E25/, MAXIT/10/
C
      NT = NN - 1
C
C     Check for input faults
C
      NFAIL = 0
      IFAULT = 1
      IF (NN .LE. 0 .OR. ALPHA .LE. ZERO .OR. BETA .LE. ZERO .OR.
     *  DELTA1 .LE. ZERO .OR. DELTA2 .LE. ZERO .OR. NGR .LE. 0 .OR.
     *  K .LE. 0 .OR. INK .LE. 0 .OR. NK .LE. 0 .OR. NNP2 .NE. NN+2)
     *  RETURN
      IFAULT = 2
      DO 25 I = 1, NN
	J = DEL(I)**2
	IF (Y(I) .LT. ZERO .OR. J .NE. DEL(I)) RETURN
	NFAIL = NFAIL + DEL(I)
   25 CONTINUE
      DHAA = NFAIL
      IFAULT = 3
      DO 30 I = 1, NT
	IF (Y(I) .GT. Y(I+1)) RETURN
   30 CONTINUE
      IFAULT = 4
      LKTH = K + (NK-1)*INK
      ITEST = NFAIL - LKTH - 1
      IF (ITEST .LE. 0) RETURN
      IM = 1
      JM = 2
      DO 55 I = 1, NT
	KTH = 0
   35   IF (Y(IM) .LT. Y(JM) .OR. DEL(IM) .LT. DEL(JM)) GO TO 45
	IF (DEL(IM) .GT. DEL(JM)) RETURN
	IF (DEL(JM) .LT. 1) GO TO 45
	KTH = KTH + 1
	IF (K .LE. KTH) RETURN
	IF (JM .GE. NN) GO TO 60
	JM = JM + 1
	GO TO 35
   45   IM = JM
	IF (JM .GE. NN) GO TO 60
	JM = JM + 1
   55 CONTINUE
C
C     Set up the A matrix
C
   60 SS = DEL(1) / FLOAT(NN-1)
      DO 65 J = 2, NN
	NJ1 = NN - J + 1
	A(1,J) = SS
	A(NN,NJ1) = DEL(NN)
   65 CONTINUE
      DO 95 J = 2, NT
	SS = DEL(J) / FLOAT(NN-J)
	JP1 = J + 1
	DO 75 I = JP1, NN
   75   A(J,I) = SS
	JM1 = J - 1
	SS = DEL(J) / FLOAT(NN-JM1)
	DO 85 I = 1, JM1
   85   A(J,I) = SS
   95 CONTINUE
      DO 99 J = 1, NN
   99 A(J,J) = ZERO
C
C     Locate the maximizer of the log-likelihood function.
C     For a given value of K, evaluate the log-modified likelihood
C     function over a grid of alpha and beta points, then use the
C     maximizer as the starting point for the N-R algorithm.
C
      FMAB = - BIG
      DO 250 KTH = K, LKTH, INK
	CALL KTHNN(NN, NNP2, DEL, Y, KTH, NFAIL, 1, H, CUMCEN, DAT)
	FMAX = - BIG
	DO 150 IM = 1, NGR
	  ALEPH1 = ALPHA + (IM - 1) * DELTA1
          DO 100 JM = 1, NGR
	    BET1 = BETA + (JM - 1) * DELTA2
	    CALL EVDER(NN, DEL, Y, A, H, 0, ALEPH1, BET1, DHAA, VML,
     *                 TEST, IFAULT)
	    IF (VML .LT. FMAX .OR. IFAULT .GT. 0) GO TO 100
	    FMAX = VML
	    ALEPH2 = ALEPH1
	    BET2 = BET1
  100     CONTINUE
  150   CONTINUE
C
C     This loop controls the maximum number (MAXIT) of N-R iterations.
C
	DO 175 IM = 1, MAXIT
	  CALL EVDER(NN, DEL, Y, A, H, 1, ALEPH2, BET2, DHAA, VML, TEST,
     *               IFAULT)
C
C     The following IF statement tests for non-positive values of alpha
C     and beta.   If the condition is true, user must specify another
C     grid search for this value of K.   Irrespective of the outcome of
C     this test, the program proceeds to a larger (valid) value of K.
C
	  IF (IFAULT .EQ. 5) GO TO 250
C
C     Test if the algorithm has converged.
C
    	  IF (TEST .LT. EPS2) GO TO 200
  175   CONTINUE
C
C     Convergence achieved or maximum number of iterations performed.
C     Compare with previous maximizer vector.   Update vector if nec.
C
  200   IF (VML .LT. FMAB) GO TO 250
	FMAB = VML
	KTH2 = KTH
	ALEPH3 = ALEPH2
	BET3 = BET2
  250 CONTINUE
      K = KTH2
      ALPHA = ALEPH3
      BETA = BET3
      VML = FMAB
C
      RETURN
      END
C
C
      SUBROUTINE KTHNN(NN, NNP2, DEL, Y, K, NFAIL, IKTH, H, CUMCEN, DAT)
C
C     ALGORITHM AS 202.1  APPL. STATIST. (1984) VOL.33, NO.2
C
C     KTHNN computes the matrix H of K-th nearest neighbour distances.
C
      INTEGER NN, NNP2, DEL(NN), K, NFAIL, IKTH, CUMCEN(NN), CENTRE,
     *   FSD, I, IDUM, J, K1, KP1, LOW, NJ1, NR
      REAL Y(NN), H(NN,NN), DAT(NNP2), ZERO, BIG, DIFF, KNN, KP1NN
C
      DATA ZERO, BIG /0.0, 1.0E25/
C
      CENTRE = 2
      LOW = 2
      DAT(1) = - BIG
      DAT(NFAIL+2) = BIG
      KP1 = K + 1
      NR = NN
      IF (IKTH .LT. 1) NR = 1
C
C     Collapse the observation vector, but record the cumulative number
C     of censored observations between failures.
C
      J = 0
      K1 = 0
      DO 150 I = 1, NFAIL
   50   J = J + 1
	IF (DEL(J) .GT. 0) GO TO 100
	K1 = K1 + 1
	GO TO 50
  100   DAT(I+1) = Y(J)
	CUMCEN(I) = K1
  150 CONTINUE
C
C     This loop moves the point of interest to the right.
C
      DO 750 I = 1, NFAIL
C
C     Move the window to the right.
C
  180   IDUM = LOW + K + 1
        DIFF = (DAT(IDUM) - DAT(CENTRE)) - (DAT(CENTRE) - DAT(LOW))
        IF (DIFF .GT. ZERO) GO TO 250
        LOW = LOW + 1
        GO TO 180
C
C     Determine which edge is the K-th nearest neighbour.
C
  250   IDUM = LOW + K
	DIFF = (DAT(IDUM) - DAT(CENTRE)) - (DAT(CENTRE) - DAT(LOW))
	KNN = DAT(IDUM)
	IF (DIFF .LE. ZERO) KNN = DAT(LOW)
	IF (IKTH .GT. 0) GO TO 400
	H(I,1) = ABS(KNN - DAT(CENTRE))
	GO TO 725
  400   FSD = LOW + CUMCEN(LOW-1) - 1
C
C     Find the K+1st nearest neighbour.
C
	IDUM = LOW + K + 1
	DIFF = (DAT(IDUM) - DAT(CENTRE)) - (DAT(CENTRE) - DAT(LOW-1))
	KP1NN = DAT(IDUM)
	IF (DIFF .GT. ZERO) KP1NN = DAT(LOW-1)
C
C     Compute the K-th nearest neighbour matrix H.   A typical row
C     consists of three segments - pre-window, window and post-window.
C
	DIFF = ABS(KNN - DAT(CENTRE))
	IF (FSD .LE. 1) GO TO 625
	K1 = FSD - 1
	DO 600 J = 1, K1
  600   H(I,J) = DIFF
  625   K1 = 1
	DO 650 J = FSD, NN
	  IF (DEL(J) .GT. 0) GO TO 635
	  H(I,J) = DIFF
	  GO TO 650
  635     H(I,J) = ABS(KP1NN - DAT(CENTRE))
	  K1 = K1 + 1
	  IF (K1 .GT. KP1) GO TO 675
  650   CONTINUE
  675   IF (J .GE. NN) GO TO 725
	K1 = J + 1
        DO 700 J = K1, NN
  700   H(I,J) = DIFF
  725   IF (CENTRE .LE. NFAIL) CENTRE = CENTRE + 1
  750 CONTINUE
C
C     Expand the H matrix.
C
      J = 0
      DO 850 K1 = 1, NFAIL
  775   J = J + 1
	NJ1 = NN - J + 1
	IF (DEL(NJ1) .LT. 1) GO TO 775
	IDUM = NFAIL - K1 + 1
	DO 800 I = 1, NR
  800   H(NJ1,I) = H(IDUM,I)
  850 CONTINUE
      DO 950 I = 1, NN
	IF (DEL(I) .GT. 0) GO TO 950
	DO 900 J = 1, NR
  900   H(I,J) = ZERO
  950 CONTINUE
      IF (IKTH .LT. 1) GO TO 990
      DO 975 I = 1, NN
  975 H(I,I) = ZERO
  990 RETURN
      END
C
C
      SUBROUTINE EVDER(NN, DEL, Y, A, H, IFLAG, ALPHA, BETA, DHAA, VML,
     *  TEST, IFAULT)
C
C     ALGORITHM AS 202.2  APPL. STATIST. (1984) VOL.33, NO.2
C
C     EVDER evaluates the log-modified-likelihood function (ML) and its
C     gradient vector and Hessian matrix.   If FLAG = 1, an N-R step is
C     performed.
C
      INTEGER NN, DEL(NN), IFLAG, IFAULT, I, J, KM
      REAL Y(NN), A(NN,NN), H(NN,NN), ALPHA, BETA, DHAA, VML, TEST,
     *  ZERO, EPS1, ONFI, ONTH, TWTH, FS, ONE, TWO, THREE, FF, FOUR,
     *  FIVE, TBIG, BIG, A1, B1, CT, D1, D0K, D1K, D2K, DHB, D2HB, DSA,
     *  G1, G2, SK, SDK, SD2K, T, T1, T2, X1, X3, X5
C
      DATA ZERO /0.0/, EPS1 /1.0E-4/, ONFI /0.2/, ONTH /0.3333333333/,
     *  TWTH /0.6666666667/, FS /0.9375/, ONE /1.0/, TWO /2.0/, THREE
     *  /3.0/, FF /3.75/, FOUR /4.0/, FIVE /5.0/, TBIG /350.0/, BIG
     *  /1.0E25/
C
      IFAULT = 5
      D0K = ZERO
      D1K = ZERO
      D2K = ZERO
      DHB = ZERO
      D2HB = ZERO
      VML = ZERO
C
C     Delete each observation, one at a time.
C
      DO 150 KM = 1, NN
	SK = ZERO
	SDK = ZERO
	SD2K = ZERO
C
C     Compute the contribution of the deleted point to ML and its
C     first two derivatives for fixed values of alpha and beta.
C
	DO 100 I = 1, NN
	  IF (A(I,KM) .LE. ZERO) GO TO 100
	  T = BETA * H(I,KM)
	  T2 = MIN(Y(KM), Y(I) + T)
	  T1 = MAX(ZERO, Y(I) - T)
	  IF (T1 .GT. T2) GO TO 100
	  A1 = (T2 - Y(I)) / T
	  B1 = (T1 - Y(I)) / T
	  X1 = A1 - B1
	  X3 = A1**3 - B1**3
	  X5 = A1**5 - B1**5
	  D0K = D0K + (X1 - X3*TWTH + X5*ONFI) * A(I,KM)
	  IF (IFLAG .LT. 1) GO TO 50
	  D1K = D1K + (X3*ONTH - X5*ONFI)* A(I,KM)
	  D2K = D2K + (X5 - X3) * A(I,KM)
   50     IF (DEL(KM) .LT. 1) GO TO 100
	  A1 = (Y(KM) - Y(I)) / T
	  A1 = A1**2
	  IF (A1 .GT. ONE) GO TO 100
	  SK = SK + (ONE - A1)**2 / H(I,KM) * A(I,KM)
	  IF (IFLAG .LT. 1) GO TO 100
	  SDK = SDK + (ONE - A1) * A1 / H(I,KM) * A(I,KM)
	  SD2K = SD2K + (FIVE*A1 - THREE) * A1 / H(I,KM) * A(I,KM)
  100   CONTINUE
	IF (DEL(KM) .LT. 1) GO TO 150
C
C     If SK <= 0 restart the N-R algorithm.
C
	IF (SK .LE. ZERO) RETURN
	VML = VML + LOG(SK * FS / ALPHA)
	IF (IFLAG .LT. 1) GO TO 150
	DHB = DHB + SDK / SK
	D2HB = D2HB + SD2K/SK - (SDK/SK)**2 * FOUR
  150 CONTINUE
      DSA = D0K * FS * BETA / ALPHA**2
      VML = VML - DSA * ALPHA
      IF (IFLAG .LT. 1) GO TO 300
C
C     Compute the gradient vector (G1, G2), the Hessian matrix
C     (H11=A1, H12=H21=B1, H22=D1), and update alpha and beta.
C
      G1 = -DHAA/ALPHA + DSA
      G2 = -D1K*FF/ALPHA + DHB*FOUR/BETA
      A1 = DHAA/ALPHA**2 - D0K*FS*BETA*TWO/ALPHA**3
      B1 = D1K*FF/ALPHA**2
      D1 = D2HB*FOUR/BETA**2 - D2K*FF/(ALPHA*BETA)
      CT = A1*D1 - B1**2
      ALPHA = ALPHA - (G1*D1 - G2*B1)/CT
      BETA = BETA + (G1*B1 - G2*A1)/CT
      TEST = SQRT(G1**2 + G2**2)
C
C     Test for negative alpha and beta values.
C     If condition is true, the algorithm is diverging.
C     Also test for unbounded increase of beta.
C
      IF (BETA .LT. EPS1 .OR. ALPHA .LT. EPS1) RETURN
      IF (BETA .LT. TBIG) GO TO 300
C
C     If beta is large, compute alpha and VML analytically.
C
      SD2K = ZERO
      VML = ZERO
      DO 250 I = 1, NN
	SDK = ZERO
	DO 200 J = 1, NN
  200   IF (A(J,I) .GT. ZERO) SDK = SDK + A(J,I)/H(J,I)
	IF (DEL(I) .GT. 0) VML = VML + LOG(SDK * FS)
	SD2K = SD2K + SDK*Y(I)
  250 CONTINUE
      ALPHA = FS * SD2K / DHAA
      BETA = BIG
      VML = -DHAA - DHAA*LOG(ALPHA) + VML
      TEST = ZERO
C
  300 IFAULT = 0
      RETURN
      END
