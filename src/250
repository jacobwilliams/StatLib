      SUBROUTINE DMAT(KSET, NP, ND, NM, NG, X, XSIGMA, MDF, CHI, DETMUL,
     *               W2, IFAULT)
C
C        ALGORITHM AS 250 APPL. STATIST. (1989), VOL. 38, NO. 3
C
C        Test of equality of dispersion matrices using
C        Sen & Puri's non-parametric and Anderson's parametric methods.
C
C     Auxiliary routine required: CHOL from AS 6
C
      INTEGER IA, IB, IC, ID, IE
      PARAMETER (IA=9, IB=10, IC=IA * (IA + 1)/2, ID=1000,
     *	        IE=IC * (IC + 1)/2)
C
C        Arguments
      INTEGER NG( * ), IFAULT, KSET, MDF, ND, NM, NP
      REAL X(IA, * ), XSIGMA(IC, * ), CHI, DETMUL, W2
C
C        Local variables
      INTEGER I, I1, II, IM, J, J1, JJ, K, L, M, MP, MR, MS, N, NPT, NT,
     *        NULLTY
      REAL XMEAN(IA, IB), XMALL(IA), DET(IB), XSALL(IC), Y(ID), VRS(IE),
     *     W(IE), DE, DETALL, F48, F6, HALF, ONE, R1, R2, RHO, S12,
     *     TEMP, TWO
      DATA HALF, ONE, TWO, F6, S12, F48 / 0.5, 1.0, 2.0, 6.0,
     *     3.4641016151377, 48.0 /
C
      NT = 0
      DO 10 K = 1, KSET
	 NT = NT + NG(K)
   10 CONTINUE
      MP = NP * (NP + 1) / 2
      MDF = 0
      CHI = 0.0
      W2 = 0.0
      DETMUL = 0.0
C
C        Check input data for failure.
C
      IFAULT = 0
      IF (KSET .LT. 1 .OR. NP .LT. 1 .OR.
     *    (ND .NE. 1 .AND. ND .NE. 2) .OR.
     *    (NM .NE. 1 .AND. NM .NE. 2)) IFAULT = 1
      IF (ND .EQ. 2 .AND. NM .NE. 2) IFAULT = 2
      DO 20 K = 1, KSET
	 IF (NG(K) .LT. 1) IFAULT = 3
   20 CONTINUE
C
C        Partial check on dimensions.
C
      IF (IA .LT. NP .OR. IB .LT. KSET .OR. ID .LT. NT) IFAULT = 4
      IF (IFAULT .NE. 0) RETURN
C
C        Check the `type' of data input
C
      IF (ND .EQ. 2) GO TO 290
C
C        Manipulate raw data according to flag given by NM and
C        calculate MEAN and SIGMA.  Note that `raw' data is `ordered'
C        if Sen & Puri's method(NM=1) is used.
C
C        Rank order the raw data if NM=1 and load it in X(L,N).
C        Ordered data for each variable are temporarily stored in Y(N).
C        Map the ranked data into its general score such that variate
C        means = 0 and variate std = 1 across all samples.
C
      IF (NM .EQ. 1) THEN
	 DO 70 L = 1, NP
	    DO 30 N = 1, NT
	       Y(N) = ONE
   30       CONTINUE
	    DO 50 I = 1, NT
	       DO 40 N = 1, I - 1
C
		  IF (X(L, I) .LT. X(L, N)) THEN
		     Y(N) = Y(N) + ONE
                  ELSE IF (X(L, I) .EQ. X(L, N)) THEN
		     Y(N) = Y(N) + HALF
		     Y(I) = Y(I) + HALF
                  ELSE
		     Y(I) = Y(I) + ONE
                  END IF
C
   40          CONTINUE
   50       CONTINUE
C
	    DO 60 N = 1, NT
	       X(L, N) = S12 * (AINT(Y(N)) / (NT + 1) - HALF)
   60       CONTINUE
   70    CONTINUE
      END IF
C
C        The `raw' or `ordered' data are now in X(L,N).  Compute mean
C        of each variate for each sample and store in XMEAN(L,K).
C
C        Compute mean of each variate over all samples, store in
C        XMALL(L).  Note that when NM=1, this step is known and is
C        equal to zero since the ranks are scaled.  When NM=2, XMALL(L)
C        is not needed.
C
      DO 100 L = 1, NP
         NPT = 0
         XMALL(L) = 0.0
         DO 90 K = 1, KSET
	    XMEAN(L, K) = 0.0
	    DO 80 N = 1, NG(K)
	       NPT = NPT + 1
	       XMEAN(L, K) = XMEAN(L, K) + X(L, NPT)
   80       CONTINUE
	    XMEAN(L, K) = XMEAN(L, K) / NG(K)
   90    CONTINUE
  100 CONTINUE
C
C        Compute sigma matrix for each sample XSIGMA(M,K) where M takes
C        values 1 through NP*(NP+1)/2.  Compute sigma of all samples
C        XSALL(M).  Calculation of these matrices is controlled by NM.
C
      M = 0
      DO 140 I = 1, NP
	 DO 130 J = 1, NP
	    M = M + 1
	    NPT = 0
	    XSALL(M) = 0.0
	    DO 120 K = 1, KSET
	      XSIGMA(M, K) = 0.0
	      DO 110 N = 1, NG(K)
		 NPT = NPT + 1
		 XSIGMA(M, K) = XSIGMA(M, K) +
     *                          (X(I, NPT) - XMEAN(I, K)) *
     *                          (X(J, NPT) - XMEAN(J, K))
  110            CONTINUE
C
		 IF (NM .EQ. 1) THEN
		    XSALL(M) = XSALL(M) + XSIGMA(M, K) +
     *                         XMEAN(I, K) * XMEAN(J, K) * NG(K)
		    XSIGMA(M, K) = XSIGMA(M, K) / (NG(K) - 1)
                 ELSE
		    XSALL(M) = XSALL(M) + XSIGMA(M, K)
                 END IF
C
  120         CONTINUE
	      IF (NM .EQ. 1) XSALL(M) = XSALL(M) / (NT - 1)
  130      CONTINUE
  140   CONTINUE
C
C          Check the method chosen: Sen and Puri (NM=1) or Anderson (NM=2)
C
        IF (NM .EQ. 2) GO TO 330
C
C          Sen and Puri's method
C
C          Compute the lower triangle of 4th moment matrix VRS of order
C          (MP,MP), where MP=MP*(NP+1)/2, in the one dimensional array
C          VRS.
C
      L = 0
      MR = 0
      DO 200 I = 1, NP
	 DO 190 J = I, NP
	    MR = MR + 1
	    MS = 0
	    DO 180 I1 = 1, I
	       DO 170 J1 = I1, NP
		  MS = MS + 1
 	          IF (MR .GE. MS) THEN
		     L = L + 1
		     NPT = 0
		     VRS(L) = 0.0
		     DO 160 K = 1, KSET
			DO 150 N = 1, NG(K)
			   NPT = NPT + 1
			   VRS(L) = VRS(L) + X(I, NPT) * X(J, NPT) *
     *                              X(I1, NPT) * X(J1, NPT)
  150                   CONTINUE
  160                CONTINUE
		     VRS(L) = VRS(L) / NT - XSALL(MR) * XSALL(MS)
                  END IF
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C        Compute XSIGMA (M,K)-XSALL(M) for all M,K.  Sotre the results
C        in XSIGMA(M,K).
C
      DO 220 K = 1, KSET
	 DO 210 M = 1, MP
	    XSIGMA(M, K) = XSIGMA(M, K) - XSALL(M)
  210    CONTINUE
  220 CONTINUE
C
C        Compure the quadratic form:
C        CHI = (XSIGMA transpose * VRS inverse * XSIGMA)
C        by first getting the square root W of VRS (see sub. CHOL), then
C        finding W inverse * XSIGMA, and finally CHI.
C
      CALL CHOL(VRS, MP, MP * (MP + 1) / 2, W, NULLTY, IFAULT)
C
C        Invert root of VRS (=W).  Note that W is an upper diagonal 
C        matrix and is inverted `in place'.
C
      IF (NULLTY .GT. 0 .OR. IFAULT .GT. 0) THEN
	 IFAULT = 5
	 RETURN
      END IF
      DO 250 J = MP, 1, -1
	 JJ = J * (J + 1) / 2
	 W(JJ) = 1.0 / W(JJ)
	 DO 240 I = J - 1, 1, -1
	    II = I * (I + 1) / 2
	    IM = II
	    TEMP = 0.0
	    DO 230 M = I + 1, J
	       IM = IM + M - 1
	       TEMP = TEMP + W(IM) * W(JJ + M - J)
  230       CONTINUE
	    W(IM) = -TEMP / W(II)
  240    CONTINUE
  250 CONTINUE
C
C        Compute CHI
C
      CHI = 0.0
      DO 280 K = 1, KSET
         DO 270 I = 1, MP
	    N = I * (I - 1) / 2
	    TEMP = 0.0
	    DO 260 J = 1, I
	       TEMP = TEMP + W(N + J) * XSIGMA(J, K)
  260       CONTINUE
	    CHI = CHI + TEMP ** 2
  270    CONTINUE
  280 CONTINUE
      MDF = (KSET - 1) * MP
C
      RETURN
C
C        Anderson's method, XSIGMA input
C
C        Calculate the total matrix for all samples and store it in
C        XSALL(M).
C
  290 CONTINUE
      DO 300 M = 1, MP
	 XSALL(M) = 0.0
  300 CONTINUE
      DO 320 K = 1, KSET
	 DO 310 M = 1, MP
	    XSALL(M) = XSALL(M) + XSIGMA(M, K)
  310    CONTINUE
  320 CONTINUE
C
C        Anderson's method, raw data or XSIGMA inpur.
C        From the diagonal XSIGMA matrix and XSALL generate the lower
C        triangle of respective covariance matrices.  Store these lower
C        triangles in the one dimensional array VRS.  Get the
C        log(determinant) for each matrix by calling subroutine `CHOL'.
C        store results in DET(K) and DETALL.
C
  330 DO 370 K = 1, KSET
	 M = 0
	 DO 350 I = 1, NP
	    N = I * (I - 1) /2 + 1
	    DO 340 J = I, NP
	       M = M + 1
	       N = N + J - 1
	       VRS(N) = XSIGMA(M, K) / (NG(K) - 1)
  340       CONTINUE
  350    CONTINUE
	 CALL CHOL(VRS, NP, NP * (NP + 1) / 2, W, NULLTY, IFAULT)
	 IF (NULLTY .GT. 0 .OR. IFAULT .GT. 0) THEN
	    IFAULT = 5
	    RETURN
         END IF
	 DE = 1.0
	 N = 0
	 DO 360 I = 1, NP
	    N = N + I
	    DE = DE * W(N)
  360    CONTINUE
	 DET(K) = TWO * LOG(DE)
  370 CONTINUE
C
      M = 0
      DO 390 I = 1, NP
	 N = I * (I - 1) / 2 + 1
	 DO 380 J = I, NP
	    M = M + 1
	    N = N + J - 1
	    VRS(N) = XSALL(M) / (NT - KSET)
  380    CONTINUE
  390 CONTINUE
      CALL CHOL(VRS, NP, NP * (NP + 1) / 2, W, NULLTY, IFAULT)
      IF (NULLTY .GT. 0 .OR. IFAULT .GT. 0) THEN
	 IFAULT = 5
	 RETURN
      END IF
      DE = 1.0
      N = 0
      DO 400 I = 1, NP
	 N = N + I
	 DE = DE * W(N)
  400 CONTINUE
      DETALL = TWO * LOG(DE)
C
C        Compute log(V), store in DETMUL.  Compute CHI.
C
      DETMUL = 0.0
      DO 410 K = 1, KSET
	 DETMUL = DETMUL + (DET(K) - DETALL) * (NG(K) - 1) / TWO
  410 CONTINUE
C
C        Calculate RHO and W2
C
      M = NT - KSET
      R1 = 0.0
      R2 = 0.0
      DO 420 K = 1, KSET
         R1 = R1 + ONE / (NG(K) - 1)
         R2 = R2 + ONE / (NG(K) - 1) ** 2
  420 CONTINUE
      RHO = ONE - (R1 - ONE / M) * (2 * (NP ** 2) + 3 * NP - 1) /
     *      (F6 * (NP + 1) * (KSET - 1))
      W2 = NP * (NP + 1) * ((NP - 1) * (NP + 2) * (R2 - ONE / M ** 2) -
     *     F6 * (KSET - 1) * (1 - RHO) * * 2) / (F48 * RHO ** 2)
C
      MDF = (KSET - 1) * MP
      CHI = -TWO * RHO * DETMUL
      RETURN 
      END
