      SUBROUTINE CCLASS (DD, DI, N1, N2, NP, KMAX, A1, COND, RAW, DATA,
     *                   SX, SY, TMV, CUT, ICC, PQ, PMIN, LKA, Z, UPMQM,
     *                   ICR, CGFE, IFAULT)
C
C        ALGORITHM AS 276.1 APPL.STATIST. (1992), VOL.41, NO.2
C
C        Obtains the optimal combinatoric classification procedure for
C        the two normal problem
C
      INTEGER DD, DI, N1, N2, NP, KMAX, COND, RAW, ICC(NP), LKA(NP + 2),
     *        ICR(2, NP), IFAULT
      REAL A1, DATA(DD, NP), SX(DI, NP), SY(DI, NP), TMV(DI, 2),
     *     CUT(NP), PQ(2, NP), PMIN, Z(DI, NP + 1), UPMQM(3, NP *
     *     (NP + 1)/2), CGFE(DI, 3 * NP + 1)
C
      EXTERNAL DIAG2, ERPROB, MCE, SWAP
C
      INTEGER HOLD, I, IFLT, IZ, J, JP1, K, MX, MY, NDEN1, NDEN2, N1N2,
     *        N1P1
      REAL ONE, SSQ, STORE, SUM, XYMEAN, ZERO
C
      PARAMETER (ZERO = 0.0E0, ONE = 1.0E0)
C
C        Initialise local variables
C
      IZ = NP + 1
      N1N2 = N1 + N2
      MX = N1N2 + 1
      MY = MX + 1
      N1P1 = N1 + 1
C
C        Check for input failures
C

      IFAULT = 0
      IF (N1 .LT. 0 .OR. N2 .LT. 0 .OR. MY .GT. DD)
     *   IFAULT = 1
      IF (NP .LT. 2 .OR. NP .GT. DI) IFAULT = IFAULT + 2
      IF (KMAX .LT. 1 .OR. KMAX .GT. NP) IFAULT = IFAULT + 4
      IF (A1 .LE. ZERO .OR. A1 .GE. ONE) IFAULT = IFAULT + 8
      IF (COND .LT. 1 .OR. COND .GT. 4) IFAULT = IFAULT + 16
      IF (RAW .LT. 0 .OR. RAW .GT. 1) IFAULT = IFAULT + 32
      IF (RAW .NE. 0 .AND. COND .EQ. 1) IFAULT = IFAULT + 64
      IF (IFAULT .GT. 0) RETURN
C
C        Check for known population parameter
C
      IF (RAW .EQ. 1) THEN
	 NDEN1 = N1
	 NDEN2 = N2
	 IF (COND.NE.2) THEN
C
C        Estimation of individual mean vector
C
	      DO 30 J = 1, NP
		    SUM = ZERO
		 DO 10 I = 1, N1
		    SUM = SUM + DATA(I, J)
  10             CONTINUE
		 DATA(MX, J) = SUM / N1
		 SUM = ZERO
		 DO 20 I = N1P1, N1N2
		    SUM = SUM + DATA(I, J)
  20             CONTINUE
		 DATA(MY, J) = SUM / N2
  30          CONTINUE
C
C        Estimation of covariance matrices
C
	      NDEN1 = N1 - 1
	      NDEN2 = N2 - 1
	   ENDIF
	   DO 60 J = 1, NP
	      XYMEAN = DATA(MX, J)
	      DO 40 I = 1, N1
		 DATA(I, J) = DATA(I, J) - XYMEAN
  40          CONTINUE
	      XYMEAN = DATA(MY, J)
	      DO 50 I = N1P1, N1N2
		 DATA(I, J) = DATA(I, J) - XYMEAN
  50          CONTINUE
  60       CONTINUE
	   DO 100 J = 1, NP
	      DO 90 K = 1, J
		 SSQ = ZERO
		 DO 70 I = 1, N1
		    SSQ = SSQ + DATA(I, J) * DATA(I, K)
  70             CONTINUE
		 SX(J, K) = SSQ / NDEN1
		 SSQ = ZERO
		 DO 80 I = N1P1, N1N2
		    SSQ = SSQ + DATA(I, J) * DATA(I, K)
  80             CONTINUE
		 SY(J, K) = SSQ / NDEN2
  90          CONTINUE
 100       CONTINUE
	   DO 120 J = 2, NP
	      DO 110 K = 1, J - 1
		 SX(K, J) = SX(J, K)
		 SY(K, J) = SY(J, K)
 110          CONTINUE
 120       CONTINUE
	   IF (COND .EQ. 3) THEN
C
C        Estimation of common mean vector
C
	    DO 130 J = 1, NP
	       DATA(MX, J) = (N1 * DATA(MX, J) + N2 * DATA(MY, J)) /
     *                       N1N2
	       DATA(MY, J) = DATA(MX, J)
  130       CONTINUE
	 ENDIF
      ENDIF
C
C         Simultaneous diagonalisation of the covariance matrices
C
      CALL DIAG2 (DI, NP, SX, SY, TMV(1, 2), Z, CGFE(1, 1),
     *    CGFE(1, NP+1), CGFE(1, 2*NP+1), CGFE(1, 3*NP+1), IFLT)
      IF (IFLT .GT. 0) THEN
	  IFAULT = 128
	  RETURN
      ENDIF
C
C         Compute population 2 mean vector (transformed)
C
      IF (COND .EQ. 3) THEN
	 DO 140 I = 1, NP
	    TMV(I, 1) = ZERO
  140    CONTINUE
      ELSE
	 DO 160 I = 1, NP
	    STORE = ZERO
	    DO 150 J = 1, NP
	       STORE = Z(J, I) * (DATA(MY, J) - DATA(MX, J)) + STORE
  150       CONTINUE
	    TMV(I, 1) = STORE
  160    CONTINUE
      ENDIF
C
C         Compute Z' * XBAR
C
      DO 180 I = 1, NP
	 STORE = ZERO
	 DO 170 J = 1, NP
	    STORE = STORE + Z(J, I) * DATA(MX, J)
  170    CONTINUE
	 Z(I, IZ) = STORE
  180 CONTINUE
C
C         Determination of marginal misclassification probabilities
C
      CALL ERPROB (DI, NP, TMV, A1, PQ, ICC, CUT)
C
C          Re-indexing the variables such that the values of p are
C          increasing
C
      DO 210 I = 1, NP - 1
	 DO 200 J = 1, NP - I
	    JP1 = J + 1
	    IF ( (PQ(1, J) .EQ. PQ(1, JP1) .AND. PQ(2, J) .GT.
     *            PQ(2, JP1)) .OR. PQ(1, J) .GT. PQ(1, JP1)) THEN
	       CALL SWAP (PQ(1, J), PQ(1, JP1))
	       CALL SWAP (PQ(2, J), PQ(2, JP1))
	       CALL SWAP (CUT(J), CUT(JP1))
	       HOLD = ICC(J)
	       ICC(J) = ICC(JP1)
	       ICC(JP1) = HOLD
	       CALL SWAP (TMV(J, 2), TMV(JP1, 2))
	       CALL SWAP (TMV(J, 1), TMV(JP1, 1))
	       DO 190 K = 1, NP
		  CALL SWAP (Z(K, J), Z(K, JP1))
  190          CONTINUE
	       CALL SWAP (Z(J, IZ), Z(JP1, IZ))
	    ENDIF
  200    CONTINUE
  210 CONTINUE
C
C         Determination of the optimal joint classification rule
C
      CALL MCE (NP, KMAX, A1, PQ, PMIN, LKA, UPMQM, ICR, IFAULT)
      END

      SUBROUTINE DIAG2 (DI, NP, SX, SY, D, Z, C, G, F, E, IFAULT)
C
C        ALGORITHM AS 276.2 APPL.STATIST. (1992), VOL.41, NO.2
C
C        Simultaneously diagonalizes two covariance matrices such that
C        Z'(SX)Z = 1 and Z'(SY)Z = D where D is a diagonal matrix
C
      INTEGER DI, NP, IFAULT
      REAL SX(DI, NP), SY(DI, NP), D(NP), Z(DI, NP), C(DI, NP),
     *     G(DI, NP), F(DI, NP), E(NP)
C
      EXTERNAL TQL2, TRED2
C
      INTRINSIC SQRT
C
      INTEGER I, IFLT, J, K
      REAL ROOTD, SUM, ZERO
      PARAMETER (ZERO = 0.0E0)
C
      IFAULT = 0
C
C        Find G such that G' (SX) G = I
C
      CALL TRED2 (DI, NP, SX, D, E, G)
      CALL TQL2 (DI, NP, D, E, G, IFLT)
      IF (IFLT .NE. 0) THEN
	 IFAULT = 1
	 RETURN
      ENDIF
      DO 10 I = 1, NP
	 IF (D(I) .LE. ZERO) THEN
	    IFAULT = 2
	    RETURN
	 ENDIF
   10 CONTINUE
      DO 30 J = 1, NP
	 ROOTD = SQRT(D(J))
	 DO 20 I = 1, NP
	    G(I, J) = G(I, J) / ROOTD
   20    CONTINUE
   30 CONTINUE
C
C        Create matrix F = G' (SY) G
C
      DO 60 I = 1, NP
	 DO 50 J = 1, NP
	    SUM = ZERO
	    DO 40 K = 1, NP
	       SUM = SY(I, K) * G(K, J) + SUM
   40       CONTINUE
	    Z(I, J) = SUM
   50    CONTINUE
   60 CONTINUE
      DO 90 I = 1, NP
	DO 80 J = 1, NP
	   SUM = ZERO
	   DO 70 K = 1, NP
	      SUM = Z(K, I) * G(K, J) + SUM
   70       CONTINUE
	   F(I, J) = SUM
   80    CONTINUE
   90 CONTINUE
C
C        Find C such that C' F C = D
C
      CALL TRED2 (DI, NP, F, D, E, C)
      CALL TQL2 (DI, NP, D, E, C, IFLT)
      IF (IFLT .NE. 0) THEN
	 IFAULT = 3
	 RETURN

      ENDIF
      DO 100 I = 1, NP
	 IF (D(I) .LE. ZERO) THEN
	    IFAULT = 4
	    RETURN
	 ENDIF
  100 CONTINUE
C
C        Create matrix Z = GC
C
      DO 130 I = 1, NP
	 DO 120 J = 1, NP
	    SUM = ZERO
	    DO 110 K = 1, NP
	       SUM = G(I, K) * C(K, J) + SUM
  110       CONTINUE
	    Z(I, J) = SUM
  120    CONTINUE
  130 CONTINUE
      RETURN
      END

      SUBROUTINE ERPROB (DI, N, TMV, A1, PQ, ICC, CUT)
C
C        ALGORITHM AS 276.3 APPL.STATIST. (1992), VOL.41, NO.2
C
C        Gives the classification rules and misclassification
C        probabilities for each of the marginal classification
C        procedures in normal population combinatoric classification
C
      INTEGER DI, N, ICC(N)
      REAL TMV(DI, 2), A1, PQ(2, N), CUT(N)
C
      REAL ALNORM, CHISQ1
C
      EXTERNAL ALNORM, CHISQ1
C
      INTRINSIC LOG, SQRT
C
      INTEGER I
      REAL B, G, H, HALF, ONE, S, T, TWO, VU, ZERO
C
      PARAMETER (ONE = 1.0E0, ZERO = 0.0E0, TWO = 2.0E0, HALF = 0.5E0)
C
C        Misclassification probabilities for cases 1 and 2
C
      B = ((ONE - A1) / A1) ** 2
      DO 10 I = 1, N
	 S = TMV(I, 2)
	 T = TMV(I, 1)
	 CUT(I) = ZERO
	 IF (S .NE. ONE) THEN
	    G = S * LOG(B / S) / (ONE - S) + (T / (ONE - S)) ** 2 * S
	    IF (G .GT. ZERO) THEN
	       VU = T / (S - ONE)
	       CUT(I) = G
	       ICC(I) = 1
	       PQ(1, I) = ONE - CHISQ1(G, VU)
	       G = G / S
	       VU = (VU + T) / SQRT(S)
	       PQ(2, I) = CHISQ1(G, VU)
	       IF (S .LT. ONE) THEN
		  ICC(I) = 2
		  PQ(1, I) = ONE - PQ(1, I)
		  PQ(2, I) = ONE - PQ(2, I)
	       ENDIF
C
C        Misclassification probabilities for cases 3 and 4
C
	    ELSE IF (S .GT. ONE) THEN
	       ICC(I) = 3
	       PQ(1, I) = ONE
	       PQ(2, I) = ZERO
	    ELSE
	       ICC(I) = 4
	       PQ(1, I) = ZERO
	       PQ(2, I) = ONE
	    ENDIF
C
C        Misclassification probabilities for cases 5 and 6
C
	    ELSE IF (T .EQ. ZERO) THEN
	       IF (A1 .LT. HALF) THEN
		  ICC(I) = 5
		  PQ(1, I) = ONE
		  PQ(2, I) = ZERO
	       ELSE
		  ICC(I) = 6
		  PQ(1, I) = ZERO
		  PQ(2, I) = ONE
	       ENDIF
C
C        Misclassification probabilities for cases 7 and 8
C
	    ELSE
	       H = T / TWO - LOG(B) / TWO / T
	       CUT(I) = H
	       ICC(I) = 7
	       PQ(1, I) = ALNORM(H, .TRUE.)
	       H = H - T
	       PQ(2, I) = ALNORM(H, .FALSE.)
	       IF (T .LT. ZERO) THEN
		  ICC(I) = 8
		  PQ(1, I) = ONE - PQ(1, I)
		  PQ(2, I) = ONE - PQ(2, I)
	       ENDIF
	 ENDIF
   10 CONTINUE
      RETURN
      END

      SUBROUTINE SWAP (X, Y)
C
C        ALGORITHM AS 276.4 APPL.STATIST. (1992), VOL.41, NO.2
C
C        Interchanges the values of two real variables : X and Y
C
      REAL X, Y
C
      REAL STORE
C
      STORE = X
      X = Y
      Y = STORE
      RETURN
      END

      REAL FUNCTION CHISQ1(G, VU)
C
C        ALGORITHM AS 276.5 APPL.STATIST. (1992), VOL.41, NO.2
C
C        Calculates the probability that a noncentral chi-square
C        variable is less than G (with DF=1 and noncentrality parameter
C        VU**2)
C
      REAL G, VU
C
      REAL ALNORM
      EXTERNAL ALNORM
      INTRINSIC SQRT
C
      REAL X, Y, ZERO
      PARAMETER (ZERO = 0.0E0)
C
      IF (G .GT. ZERO) THEN
	 X = SQRT(G) - VU
	 Y = -SQRT(G) - VU
	 CHISQ1 = ALNORM(X, .FALSE.) - ALNORM(Y, .FALSE.)
	 RETURN
      ELSE

	 CHISQ1 = ZERO
	 RETURN
      ENDIF
      END

C-------------------------------------------------------------------------


