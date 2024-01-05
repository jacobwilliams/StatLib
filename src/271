      SUBROUTINE MCE (N, KMAX, A1, PQ, PMIN, LKA, UPMQM, ICR, IFAULT)
C
C        ALGORITHM AS 271 APPL.STATIST. (1992), VOL.41, NO.1
C
C        Obtains the optimal (smallest misclassification probability)
C        joint classification procedure combining at most KMAX of N
C        independent marginal classification rules.
C
      INTEGER N, KMAX, LKA(N  + 2), ICR(2, N), IFAULT
      REAL A1, PQ(2, N), PMIN, UPMQM(3, N * (N + 1) / 2)
C
      INTRINSIC ABS, MIN
C
      INTEGER HELP, HOLD, I, I0, I2, II, J, J0, JI, JJ, K, L, LL, M,
     *        MM, N1, N2, UPPER, R
      REAL A2, ERP, ERQ, NEGONE, ONE, PSTAR, SIGN, STORE, ZERO
C
      PARAMETER (ZERO = 0.0E0, ONE = 1.0E0 , NEGONE = -1.0E0)
C
C        Check for errors in the input variables
C
      IFAULT = 0
      IF (N .LT. 2) IFAULT = 1
      IF (A1 .LE. ZERO .OR. A1 .GE. ONE) IFAULT = IFAULT + 2
      DO 10 I = 2, N
	 IF (PQ(1, I - 1) .GT. PQ(1, I)) THEN
	    IFAULT = IFAULT + 4
	    GOTO 20
	 ENDIF
   10 CONTINUE
   20 DO 30 I = 1, N
	 IF (PQ(2, I) .LT. ZERO .OR. PQ(2, I) .GT. ONE) THEN
	    IFAULT = IFAULT + 8
	    GOTO 40
	 ENDIF
   30 CONTINUE
   40 IF (PQ(1, 1) .LT. ZERO .OR. PQ(1, N) .GT. ONE)
     *   IFAULT = IFAULT + 16
      IF (KMAX .LT. 1 .OR. KMAX .GT. N) IFAULT = IFAULT + 32
      IF (IFAULT .GT. 0) RETURN
C
C         Initialization step
C
      A2 = ONE - A1
      PMIN = ONE
      R = 1
      ICR(2, 1) = 1
      DO 50 I = 2, N
	 IF (PQ(2, I - 1) .GT. PQ(2, I)) R = R + 1
	 ICR(2, I) = R
   50 CONTINUE
      UPMQM(1, 1) = ONE
      UPMQM(1, 2) = NEGONE
      UPMQM(1, 3) = ONE
      I0 = 2
      II = 4
      DO 70 J = 3, N
	 UPMQM(1, II) = -UPMQM(1, I0)
	 SIGN = UPMQM(1, I0)
	 II = II + 1
	 DO 60 I = 3, J
	    UPMQM(1, II) = SIGN * (ABS(UPMQM(1, I0)) +
     *                     ABS(UPMQM(1, I0 + 1)))
	    SIGN = -SIGN
	    I0 = I0 + 1
	    II = II + 1
   60    CONTINUE
      UPMQM(1, II) = ONE
      I0 = I0 + 1
      II = II + 1
   70 CONTINUE
C
C         Start with the full candidate problem (1, 2, ..., N)
C
      K = N
      MM = N
      DO 80 I = 1, N
	 ICR(1, I) = I
   80 CONTINUE
      N1 = N + 1
      N2 = N + 2
C
C         Optimality check
C
   90 M = K - MM
      HOLD = M + 1
      HELP = M + MM
      UPPER = MIN(HELP, KMAX)
      IF (HOLD .EQ. 1) THEN
	 UPMQM(2, 1) = PQ(1, ICR(1, 1))
	 UPMQM(3, 1) = PQ(2, ICR(1, 1))
	 PSTAR  = A1 * UPMQM(2, 1) + A2 * UPMQM(3, 1)
	 IF (PSTAR .LT. PMIN) THEN
	    PMIN = PSTAR
	    LKA(N1) = 1
	    LKA(N2) = 1
	    LKA(1) = ICR(1, 1)
	 ENDIF
	 HOLD = 2
      ENDIF
      IF (HELP .NE. 1) THEN
	 II = HOLD * (HOLD - 1) / 2 + 1
	 I0 = II - HOLD + 1
	 DO 110 J = HOLD, UPPER
	    L = ICR(1, J)
	    UPMQM(2, II) = UPMQM(2, I0) + PQ(1, L)
	    UPMQM(3, II) = UPMQM(3, I0) + PQ(2, L)
	    II = II + 1
	    DO 100 I = 3, J
	       UPMQM(2, II) = UPMQM(2, I0 + 1) +
     *                        UPMQM(2, I0) * PQ(1, L)
	       UPMQM(3, II) = UPMQM(3, I0 + 1) +
     *                        UPMQM(3, I0) * PQ(2, L)
	       II = II + 1
	       I0 = I0 + 1
  100       CONTINUE
	    UPMQM(2, II) = UPMQM(2, I0) * PQ(1, L)
	    UPMQM(3, II) = UPMQM(3, I0) * PQ(2, L)
	    II = II + 1
	    I0 = I0 + 1
  110   CONTINUE
      ENDIF
      J0 = HOLD * (HOLD - 1) / 2
      DO 160 J = HOLD, UPPER
	 I0 = 1
	 DO 150 I = 1, J
	    ERP = ZERO
	    ERQ = ZERO
	    II = I0
	    JI = J0 + I
	    DO 120 L = I, J
	       ERQ = ERQ + UPMQM(1, II) * UPMQM(3, JI)
	       II = II + L
	       JI = JI + 1
  120       CONTINUE
	    LL = J - I + 1
	    II = LL * (LL + 1) / 2
	    JI = J0 + LL
	    DO 130 L = LL, J
	       ERP = ERP + UPMQM(1, II) * UPMQM(2, JI)
	       II = II + L
	       JI = JI + 1
  130       CONTINUE
	    I0 = I0 + I + 1
	    PSTAR = A1 * ERP + A2 * ERQ
	    IF (PSTAR .LT. PMIN .OR.
     *         (PSTAR .EQ. PMIN .AND. LKA(N2) .GT. J)) THEN
	       PMIN = PSTAR
	       LKA(N2) = J
	       LKA(N1) = I
	       DO 140 L = 1, J
		  LKA(L) = ICR(1, L)
  140          CONTINUE
	   ENDIF
  150 CONTINUE
      J0 = J0 + J
  160 CONTINUE
      MM = 0
C
C        Backward step
C
  170 IF (ICR(2, ICR(1, 1)) .EQ. R) RETURN
      L = ICR(2, ICR(1, K))
  180 K = K - 1
      IF (K .GT. 0) THEN
	 IF (ICR(2, ICR(1, K)) .GE. L) GOTO 180
      ENDIF
      IF (K .GT. 0 .AND. L .EQ. R) THEN
	L = ICR(2, ICR(1, K))
	K = K - 1
      ENDIF
      L = L + 1
      K = K + 1
      HELP = ICR(1, K) + 1
      DO 190 I = HELP, N
	 IF (ICR(2, I) .EQ. L) GOTO 200
  190 CONTINUE
  200 ICR(1, K) = I
C
C        Fathoming step
C
  210 JJ = ICR(1, K) - 1
      J = 1
      STORE = PQ(2, ICR(1, K))
      DO 230 I = 1, JJ
	 IF (ICR(1, J) .GT. I) THEN
	    IF (PQ(2, I) .LE. STORE) THEN
	      IF (MM .EQ. 0) THEN
		 GOTO 170
	      ELSE
C
C        Side step
C
		 HOLD = ICR(2, ICR(1, K))
		 K = K - 1
		 IF (HOLD .EQ. R) GOTO 90
		 K = K + 1
		 HELP = ICR(1, K) + 1
		 DO 220 I2 = HELP, N
		    IF (ICR(2, I2) .NE. HOLD) THEN
		       ICR(1, K) = I2
		       GOTO 210
		    ENDIF
 220             CONTINUE
	      ENDIF
	   ENDIF
	   ELSE
	      J = J + 1
	ENDIF
  230 CONTINUE
C
C        Forward step
C
      MM = MM + 1
      IF (ICR(1, K) .EQ. N) GOTO 90
      IF (ICR(1, K) .LT. N) THEN
	  K = K + 1
	  ICR(1, K) = ICR(1, K - 1) + 1
	  GOTO 210
      ENDIF
      RETURN
      END
