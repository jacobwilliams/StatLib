      SUBROUTINE MEDIAN(X, Y, N, IP, IT, IFAULT)
C
C     ALGORITHM AS 78  APPL. STATIST. (1974) VOL.23, NO.3
C
C     The mediancentre, generalising the median, is the point with
C     minimum total distance from a set of multivariate samples
C
      REAL X(N, IP), Y(IP)
      REAL Z(50), LAMBDA, LEPSI, LEPSR, LEPSD
C
      DOUBLE PRECISION C(50), COMP, D, DD, DELTA, EPSR, EPSD, SLAM,
     +	ZERO, ONE
C
      DATA LEPSD /0.0001/, LEPSR / 0.00001/, LEPSI /0.000001/
      DATA ICOUNT /100/, LCOUNT /50/
      DATA ZERO /0.D0/, ONE /1.D0/
C
C     Initial settings
C
      IFAULT = 0
      LL = 0
      II = 1
      IF (N .EQ. 1) GO TO 25
      IF (N .LE. 0 .OR. IP .LE. 0) GO TO 5
C
C     Calculate the diameter
C
      DIAM = 0.0
      DO 2 I = 2, N
	DO 2 J = 1, I-1
	  S = 0.0
	  DO 1 K = 1, IP
    1     S = S + (X(I,K) - X(J,K))**2
	  DIAM = MAX(S, DIAM)
    2 CONTINUE
      DIAM = SQRT(DIAM)
      EPSR = LEPSR * DIAM
      EPSI = LEPSI * DIAM
      EPSD = LEPSD * DIAM
C
C     Initial median centre = the centroid
C
      U1 = 1.0 / FLOAT(N)
      DO 4 J = 1, IP
	S = 0.0
	DO 3 I = 1, N
    3   S = S + X(I,J)
	Y(J) = S * U1
    4 CONTINUE
      IT = ICOUNT
      IF (IP .LE. 50) GO TO 6
    5 IFAULT = 1
      IT = 0
      RETURN
C
C     Main iterative loop
C
    6 DO 23 L = 1, ICOUNT
C
C     Direction cosines and resultant
C
	CORNER = 0.0
	DO 7 J = 1, IP
    7   C(J) = ZERO
	LL = L
	DO 11 I = 1, N
	  D = ZERO
	  DO 8 J = 1, IP
    8     D = D + (X(I,J) - Y(J))**2
          DD = SQRT(D)
	  IF (DD .GT. EPSD) GO TO 9
	  CORNER = CORNER + 1.0
	  II = I
	  GO TO 11
    9     D = ONE / DD
	  DO 10 J = 1, IP
   10     C(J) = C(J) + (X(I,J) - Y(J)) * D
   11   CONTINUE
	D = ZERO
	DO 12 J = 1, IP
   12   D = D + C(J)**2
	D = SQRT(D)
	DD = D
C
C     Tests for zero resultant or degenerate solution
C
        IF (CORNER .EQ. 0.0) GO TO 13
        IF (D .LE. CORNER) GO TO 25
        D = D - CORNER
   13   IF (D .LE. EPSR) GO TO 24
	DD = ONE / DD
        DO 14 J = 1, IP
   14   C(J) = C(J) * DD
C
C     Step by bisection to give zero component at lambda
C
	U1 = 0.0
	U2 = DIAM
        DO 20 LC = 1, LCOUNT
	  COMP = ZERO
	  LAMBDA = 0.5 * (U1 + U2)
	  SLAM = LAMBDA * LAMBDA
	  DO 15 J = 1, IP
   15     Z(J) = Y(J) + LAMBDA * C(J)
	  DO 17 I = 1, N
	    DELTA = ZERO
	    D = SLAM
	    DO 16 J = 1, IP
	      XX = X(I,J)
	      D = D - (XX - Y(J))**2
	      DELTA = DELTA + (XX - Z(J))**2
   16       CONTINUE
	    DD = SQRT(DELTA)
	    IF (DD .LT. EPSD) GO TO 21
	    COMP = COMP - (D + DELTA) / DD
   17     CONTINUE
	  IF (COMP .GT. ZERO) GO TO 18
	  U2 = LAMBDA
	  GO TO 19
   18     U1 = LAMBDA
   19     IF ((U2 - U1) .LE. EPSI) GO TO 21
   20   CONTINUE
C
   21   DO 22 J = 1, IP
   22   Y(J) = Y(J) + C(J) * LAMBDA
   23 CONTINUE
C
   24 IT = LL
      RETURN
   25 IT = -LL
      DO 26 J = 1, IP
   26 Y(J) = X(II,J)
      RETURN
      END

