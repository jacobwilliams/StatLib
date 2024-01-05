      SUBROUTINE RCINT(N, RC, L, X, WT, Y, M, WORK, NS, WORK2)
C
C     ALGORITHM AS 281.1 APPL. STATIST. (1993) VOL.42, NO.1
C
C     Scaling and rounding regression coefficients to integers.
C
C     N.B. No checks of input parameters are made in this routine.
C     N.B. The published routine is in double precision contrary to
C          the rules of the journal.
C
      INTEGER N, L, M, NS
      DOUBLE PRECISION RC(N), X(N,L), WT(L), Y(L), WORK(3,M), WORK2(N)
C
C     Local variables
C
      INTEGER I, J, NRC, JMIN, JMAX, K, NSC, KK
      LOGICAL OK
      DOUBLE PRECISION ZERO, EPS, HALF, BIG, SCMIN, RCSUM, RCT, SCMAX,
     +    WTSUM, YSS, T, SC, RCSS, RCSP, SDP, SD, SDPP, SCP
C     EXTERNAL RCINT2     ! This is a subroutine, not a function!
C
      DATA ZERO,  EPS,    HALF,  BIG
     +    /0.D0, 1.D-08, 0.5D0, 1.D08/
C
C     Identify sums and range of coefficients
C
      NS = 0
      NRC = 0
      SCMIN = BIG
      RCSUM = ZERO
      DO 10 I = 1, N
	RCT = ABS( RC(I) )
	IF ( RCT .EQ. ZERO ) GO TO 10
	NRC = NRC + 1
	RCSUM = RCSUM + RCT
	IF ( RCT .LT. SCMIN ) SCMIN = RCT
   10 CONTINUE
      IF ( NRC .EQ. 0 ) RETURN
C
C     Derive range for scaling factor K.
C
      SCMIN = HALF / SCMIN
      SCMAX = SCMIN + M / RCSUM
C
C     Derive linear predictor, sum of squares and sum of weights.
C
      WTSUM = ZERO
      YSS = ZERO
      DO 30 I = 1, L
	IF (WT(I) .EQ. ZERO ) GO TO 30
	T = ZERO
	DO 20 J = 1, N
	  T = T + RC(J) * X(J,I)
   20   CONTINUE
	Y(I) = T
	WTSUM = WTSUM + WT(I)
	YSS = YSS + T * T * WT(I)
   30 CONTINUE
      IF ( WTSUM .LE. ZERO .OR. YSS .LE. ZERO ) RETURN
      NSC = 0
      DO 40 I = 1, M
	WORK(1,I) = ZERO
   40 CONTINUE
C
C     For each coefficient in turn ...
C
      DO 130 I = 1, N
	RCT = ABS( RC(I) )
	IF ( RCT .EQ. ZERO ) GO TO 130
C
C     ... derive scaling factors at integer boundaries ...
C
	JMIN = NINT( RCT * SCMIN + EPS )
	JMAX = NINT( RCT * SCMAX + EPS )
	DO 120 J = JMIN, JMAX
	  SC = ( J - HALF ) / RCT
	  IF ( SC .LT. SCMIN .OR. SC .GT. SCMAX ) GO TO 120
C
C     ... derive sum of squares and sum of products for each
C     scaling factor.
C
	  RCSS = ZERO
	  RCSP = ZERO
	  DO 50 K = 1, N
	    WORK2(K) = NINT( RC(K) * SC + SIGN( EPS, RC(K) ) )
   50     CONTINUE
C
	  DO 70 K = 1, L
	    IF ( WT(K) .EQ. ZERO ) GO TO 70
	    T = ZERO
	    DO 60 KK = 1, N
	      T = T + WORK2(KK) * X(KK,K)
   60       CONTINUE
	    RCSS = RCSS + T * T * WT(K)
	    RCSP = RCSP + T * Y(K) * WT(K)
   70     CONTINUE
C
C     Sort by scaling factor.
C
	  NRC = 0
   80     NRC = NRC + 1
	  IF ( NRC .GT. M ) GO TO 120
	  IF ( ABS( SC - WORK(1,NRC) ) .LT. SC * EPS ) GO TO 120
	  IF ( NRC .LE. NSC .AND. SC .GT. WORK(1,NRC) ) GO TO 80
	  IF ( NSC .LT. M ) NSC = NSC + 1
	  IF ( NRC .EQ. NSC ) GO TO 110
	  DO 100 K = NSC - 1, NRC, -1
	    DO 90 KK = 1, 3
	      WORK(KK, K+1) = WORK(KK,K)
   90       CONTINUE
  100     CONTINUE
  110     WORK(1,NRC) = SC
	  WORK(2,NRC) = RCSS
	  WORK(3,NRC) = RCSP
  120   CONTINUE
  130 CONTINUE
C
C     Identify optimal values for scaling factor.
C
      IF ( NSC .LE. 1 ) RETURN
      SDP = BIG
      DO 150 K = 1, NSC - 1
	RCSS = WORK(2,K)
	RCSP = WORK(3,K)
	IF ( RCSP .EQ. ZERO ) RETURN
C
C     Adjust scaling factor and derive residual sum of squares.
C
	SC = RCSS / RCSP
	IF ( SC .LT. WORK(1,K) ) THEN
	  SC = WORK(1,K)
	ELSE IF ( SC .GT. WORK(1,K+1) ) THEN
	  SC = WORK(1,K+1)
	END IF
	SD = YSS - 2 * RCSP / SC + RCSS / ( SC * SC )
	IF ( SD .LT. EPS * YSS ) SD = ZERO
	IF ( K .EQ. 1 .OR. SDP .GT. SD .OR. .NOT. OK
     +      .OR. ( NS .GT. 0 .AND. SDP .GE. SDPP ) ) GO TO 140
	NS = NS + 1
	WORK(1,NS) = SCP
	WORK(2,NS) = SQRT( SDP / WTSUM )
	SDPP = SDP
  140   OK = SD .LT. SDP
	SDP = SD
	SCP = SC
  150 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE RCINT2(N, RC, INT, SCALE)
C
C     ALGORITHM AS 281.2 APPL. STATIST. (1993) VOL.42, NO.1
C
      INTEGER N, INT(N)
      DOUBLE PRECISION RC(N), SCALE
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION EPS
      DATA EPS / 1.D-08 /
C
      DO 10 I = 1, N
	INT(I) = NINT( RC(I) * SCALE + SIGN( EPS, RC(I) ) )
   10 CONTINUE
C
      RETURN
      END
