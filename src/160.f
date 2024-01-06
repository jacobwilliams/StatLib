      SUBROUTINE SCREEN(NVAR, MP, MM, NTAB, TABLE, DIM, GSQ, DGFR,
     *   PART, MARG, DFS, IP, IM, ISET, JSET, CONFIG, FIT, SIZE,
     *   COORD, X, Y, IFAULT)
C
C     ALGORITHM AS 160 APPL. STATIST. (1981) VOL.30, NO.1
C
C     Screen all efects for partial and marginal association.
C
      INTEGER NVAR, MP, MM, NTAB, IP(NVAR,MP), IM(NVAR,MM), DGFR(NVAR),
     *   DFS(NVAR,MP), ISET(NVAR), JSET(NVAR), CONFIG(NVAR,MP), 
     *   DIM(NVAR), DF, SIZE(NVAR), COORD(NVAR)
      REAL GSQ(NVAR), PART(NVAR,MP), MARG(NVAR,MP), TABLE(NTAB),
     *   FIT(NTAB), X(NTAB), Y(NTAB), ZERO
      DATA ZERO /0.0/
C
C     Check for input errors
C
      IFAULT = 1
      IF (NVAR .LE. 1) RETURN
      ISZ = 1
      DO 100 I = 1, NVAR
	IF (DIM(I) .LE. 1) IFAULT = 2
	ISZ = ISZ * DIM(I)
  100 CONTINUE
      IF (ISZ .NE. NTAB) IFAULT = 2
      MAX = 1
      LIM = NVAR / 2
      DO 110 I = 1, LIM
  110 MAX = MAX * (NVAR - I + 1) / I
      IF (MP .LT. MAX) IFAULT = 3
      MAX = 1
      LIM = (NVAR - 1) / 2
      DO 120 I = 1, LIM
  120 MAX = MAX * (NVAR - I) / I
      MAX = MAX * NVAR
      IF (MM .LT. MAX) IFAULT = 4
      IF (IFAULT .GT. 1) RETURN
C
C     Fit the no effect model
C
      DGFR(NVAR) = NTAB - 1
      AVG = ZERO
      IFAULT = 5
      DO 130 I = 1, NTAB
	IF (TABLE(I) .LT. ZERO) RETURN
	AVG = AVG + TABLE(I)
  130 CONTINUE
      IFAULT = 0
      AVG = AVG / NTAB
      CALL RESET(FIT, NTAB, AVG)
      CALL LIKE(GSQ, FIT, TABLE, NTAB)
C
C     Begin fitting effects
C
      NV1 = NVAR - 1
      DO 200 M = 1, NV1
C
C     Set up the arrays IP and IM
C
        CALL CONF(NVAR, M, MP, MM, ISET, JSET, IP, IM, NP)
C
C     Fit the saturated model
C
        CALL RESET(FIT, NTAB, AVG)
        CALL EVAL(IP, NP, M, 1, NVAR, MP, CONFIG, DIM, DGFR(M))
        CALL LOGFIT(NVAR, NTAB, NP, DIM, CONFIG, TABLE, FIT, SIZE,
     *     COORD, X, Y)
        CALL LIKE(GSQ(M+1), FIT, TABLE, NTAB)
C
C     Move the first column of IP to the last
C
        DO 150 I = 1, M
	  ITP = IP(I,1)
	  NP1 = NP - 1
	  DO 140 J = 1, NP1
  140     IP(I,J) = IP(I,J+1)
	  IP(I,NP) = ITP
  150   CONTINUE
        L3 = -M + 1
        DO 190 J = 1, NP
C
C     Fit the effects in IP ignoring the last column
C
	  CALL RESET(FIT, NTAB, AVG)
	  CALL EVAL(IP, NP-1, M, 1, NVAR, MP, CONFIG, DIM, DF)
	  CALL LOGFIT(NVAR, NTAB, NP-1, DIM, CONFIG, TABLE, FIT, SIZE,
     *         COORD, X, Y)
	  CALL LIKE(G21, GIT, TABLE, NTAB)
	  DFS(M,J) = DGFR(M) - DF
	  PART(M,J) = G21 - GSQ(M+1)
C
C     For M = 1, partials and marginals are equal
C
  	  IF (M .GT. 1) GO TO 160
	  MARG(1,J) = PART(1,J)
	  GO TO 170
C
C     Fit the last column alone
C
  160     CALL RESET(FIT, NTAB, AVG)
	  CALL EVAL(IP, 1, M, NP, NVAR, MP, CONFIG, DIM, DF)
	  CALL LOGFIT(NVAR, NTAB, 1, DIM, CONFIG, TABLE, FIT, SIZE,
     *         COORD, X, Y)
	  CALL LIKE(G22, FIT, TABLE, NTAB)
C
C     Locate the appropriate columns of IM and fit them
C
	L3 = L3 + M
	  CALL RESET(FIT, NTAB, AVG)
	  CALL EVAL(IM, M, M-1, L3, NVAR, MM, CONFIG, DIM, DF)
	  CALL LOGFIT(NVAR, NTAB, M, DIM, CONFIG, TABLE, FIT, SIZE, 
     *         COORD, X, Y)
	  CALL LIKE(G23, FIT, TABLE, NTAB)
	  MARG(M,J) = G23 - G22
C
C     Move the next effect to be ignored to the last in IP
C
  170     DO 180 I = 1, M
	    ITP = IP(I,NP)
	    IP(I,NP) = IP(I,J)
	    IP(I,J) = ITP
  180     CONTINUE
  190   CONTINUE
C
        DGFR(NVAR) = DGFR(NVAR) - DGFR(M)
        GSQ(M) = GSQ(M) - GSQ(M+1)
  200 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE CONF(N, M, MP, MM, ISET, JSET, IP, IM, NP)
C
C     ALGORITHM AS 160.1 APPL. STATIST. (1981) VOL.30, NO.1
C
C     Set up the arrays IP and IM for a given N and M.   Essentially
C     IP contains all possible combinations of (N choose M).   For each
C     combination found IM contains all combinations of degree M-1.
C
      INTEGER ISET(N), JSET(N), IP(N,MP), IM(N,MM)
      LOGICAL ILAST, JLAST
C
      ILAST = .TRUE.
      NP = 0
      NM = 0
C
C     Get IP
C
  100 CALL COMBO(ISET, N, M, ILAST)
      IF (ILAST) RETURN
      NP = NP + 1
      DO 110 I = 1, M
  110 IP(I,NP) = ISET(I)
      IF (M .EQ. 1) GO TO 100
C
C     Get IM
C
      JLAST = .TRUE.
      L = M - 1
  120 CALL COMBO(JSET, M, L, JLAST)
      IF (JLAST) GO TO 100
      NM = NM + 1
      DO 130 I = 1, L
	JS = JSET(I)
	IM(I,NM) = ISET(JS)
  130 CONTINUE
      GO TO 120
C
      END
C
C
      SUBROUTINE COMBO(ISET, N, M, LAST)
C
C     ALGORITHM AS 160.2  APPL. STATIST. (1981) VOL.30, NO.1
C
C     Subroutine to generate all possible combinations of M of the
C     integers from 1 to N in a stepwise fashion.   Prior to the first
C     call, LAST should be set to .FALSE.   Thereafter, as long as LAST
C     is returned .FALSE., a new valid combination has been generated.
C     When LAST goes .TRUE., there are no more combinations.
C
      LOGICAL LAST
      INTEGER N, M, ISET(M)
C
      IF (LAST) GO TO 110
C
C     Get next element to increment
C
      K = M
100   L = ISET(K) + 1
      IF (L + M - K .LE. N) GO TO 150
      K = K - 1
C
C     See if we are done
C
      IF (K .LE. 0) GO TO 130
      GO TO 100
C
C     Initialize first combination
C
110   DO 120 I = 1, M
120   ISET(I) = I
130   LAST = .NOT. LAST
      RETURN
C
C     Fill in remainder of combination.
C
150   DO 160 I = K, M
	ISET(I) = L
	L = L + 1
160   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE EVAL(IAR, NC, NV, IBEG, NVAR, MAX, CONFIG, DIM, DF)
C
C     ALGORITHM AS 160.3 APPL. STATIST. (1981) VOL.30, NO.1
C
C     IAR  = array containing the effects to be fitted
C     NC   = number of columns of IAR to be used
C     NV   = number of variables in each effect
C     IBEG = gebinning column
C     DF   = degrees of freedom
C
C     CONFIG is in a format compatible with algorithm AS 51
C
      INTEGER IAR(NVAR,MAX), CONFIG(NVAR,NC), DIM(NVAR), DF
C
      DF = 0
      DO 110 J = 1, NC
	KK = 1
	DO 100 I = 1, NV
	  L = IBEG + J - 1
	  K = IAR(I,L)
	  KK = KK * (DIM(K) - 1)
	  CONFIG(I,J) = K
  100   CONTINUE
	CONFIG(NV+1,J) = 0
	DF = DF + KK
  110 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE RESET(FIT, NTAB, AVG)
C
C     ALGORITHM AS 160.4 APPL. STATIST. (1981) VOL.30, NO.1
C
C     Initialize the fitted values to the average entry
C
      REAL FIT(NTAB)
C
      DO 100 I = 1, NTAB
	FIT(I) = AVG
  100 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE LIKE(GSQ, FIT, TABLE, NTAB)
C
C     ALGORITHM AS 160.5 APPL. STATIST. (1981) VOL.30, NO.1
C
C     Compute the likelihood-ration chi-square
C
      REAL FIT(NTAB), TABLE(NTAB), ZERO, TWO
      DATA ZERO /0.0/, TWO /2.0/
C
      GSQ = ZERO
      DO 100 I = 1, NTAB
	IF (FIT(I) .EQ. ZERO .OR. TABLE(I) .EQ. ZERO) GO TO 100
	GSQ = GSQ + TABLE(I) * LOG(TABLE(I) / FIT(I))
  100 CONTINUE
      GSQ = TWO * GSQ
      RETURN
      END
C
C
      SUBROUTINE LOGFIT(NVAR, NTAB, NCON, DIM, CONFIG, TABLE, FIT, SIZE,
     *     COORD, X, Y)
C
C     ALGORITHM AS 160.6 APPL. STATIST. (1981) VOL.30, NO.1
C
C     Iterative proportional fitting of the marginals of a contingency
C     table.   Relevant code from AS 51 is used.
C
      REAL TABLE(NTAB), FIT(NTAB), MAXDEV, X(NTAB), Y(NTAB), ZERO
      INTEGER CONFIG(NVAR,NCON), DIM(NVAR), SIZE(NVAR), COORD(NVAR)
      LOGICAL OPTION
      DATA MAXDEV /0.25/, MAXIT /25/, ZERO /0.0/
C
      DO 230 KK = 1, MAXIT
C
C     XMAX is the maximum deviation between fitted and true marginal
C
	XMAX = ZERO
	DO 220 II = 1, NCON
	  OPTION = .TRUE.
C
C     Initialize arrays
C
	  SIZE(1) = 1
	  NV1 = NVAR - 1
	  DO 100 K = 1, NV1
	    L = CONFIG(K,II)
	    IF (L .EQ. 0) GO TO 110
	    SIZE(K+1) = SIZE(K) * DIM(L)
  100     CONTINUE
	  K = NVAR
  110     N = K - 1
	  ISZ = SIZE(K)
	  DO 120 J = 1, ISZ
	    X(J) = ZERO
	    Y(J) = ZERO
  120     CONTINUE
C
C     Initialize co-ordinates
C
  130     DO 140 K = 1, NVAR
  140     COORD(K) = 0
C
C     Find locations in tables
C
	  I = 1
  150     J = 1
	  DO 160 K = 1, N
	    L = CONFIG(K,II)
	    J = J + COORD(L) * SIZE(K)
  160     CONTINUE
          IF (.NOT. OPTION) GO TO 170
C
C     Compute marginals
C
	  X(J) = X(J) + TABLE(I)
	  Y(J) = Y(J) + FIT(I)
	  GO TO 180
C
C     Make adjustments
C
  170     IF (Y(J) .LE. ZERO) FIT(I) = ZERO
	  IF (Y(J) .GT. ZERO) FIT(I) = FIT(I) * X(J) / Y(J)
C
C     Update co-ordinates
C
  180     I = I + 1
	  DO 190 K = 1, NVAR
	    COORD(K) = COORD(K) + 1
	    IF (COORD(K) .LT. DIM(K)) GO TO 150
	    COORD(K) = 0
  190     CONTINUE
	  IF (.NOT. OPTION) GO TO 200
	  OPTION = .FALSE.
	  GO TO 130
C
C     Find the largest deviation
C
  200     DO 210 I = 1, ISZ
	    E = ABS(X(I) - Y(I))
	    IF (E .GT. XMAX) XMAX = E
  210     CONTINUE
  220   CONTINUE
C
C     Test convergence
C
	IF (XMAX .LT. MAXDEV) RETURN
  230 CONTINUE
C
      RETURN
      END

