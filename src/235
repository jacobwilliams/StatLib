	SUBROUTINE RTALLY(NEWVAL, EMAX, MVCODE, MMAX, N, VALUES, COUNT,
     +		K, IFAULT)
C
C	ALGORITHM AS235  APPL. STATIST. (1988) VOL. 37, NO. 2
C
	INTEGER KCUT
	PARAMETER (KCUT = 25)
	INTEGER COUNT(*)
	INTEGER HI, I, J, K, KMAX, LO, M, N, IFAULT, MMAX
	REAL VALUES(*), NEWVAL, MVCODE(*)
C
C	Tallies (tabulates) input values NEWVAL with their frequencies.
C	The K distinct values are stored in array VALUES in ascending
C	order, and their frequencies are stored in COUNT.
C	N is the current number of observations.
C
	IFAULT = 0
	IF (MMAX .GT. 0) THEN
	  DO 99 I = 1, MMAX
	    IF (NEWVAL .EQ. MVCODE(I)) RETURN
   99	  CONTINUE
	END IF
	N = N + 1
C
C	Initialize
C
	IF (N .LE. 1) THEN
	  N = 1
	  K = 1
	  J = 1
	  GO TO 5
	END IF
	J = K + 1
	IF (NEWVAL .GT. VALUES(K)) GO TO 3
C
	IF (K .LT. KCUT) THEN
C
C	Simple searching of VALUES from first element, for small K
C
	  DO 1 I = 1, K
	    J = I
	    IF (NEWVAL - VALUES(J)) 3, 6, 1
    1	  CONTINUE
	ELSE
C
C	Bisection search
C
	  J = 1
	  IF (NEWVAL .LT. VALUES(1)) GO TO 3
	  LO = 0
	  HI = K + 1
    2	  J = (HI + LO) / 2
	  IF (NEWVAL .EQ. VALUES(J)) GO TO 6
C
C	J = LO means that NEWVAL is bracketed between VALUES(LO) and
C	VALUES(HI).   Insert NEWVAL at position HI and move the rest of
C	VALUES and COUNT up.
C
	  IF (J .EQ. LO) THEN
	    J = HI
	    GO TO 3
	  END IF
	  IF (NEWVAL .LT. VALUES(J)) THEN
	    HI = J
	  ELSE
	    LO = J
	  END IF
	  GO TO 2
	END IF
C
    3	IF (K .EQ. KMAX) THEN
	  IFAULT = 1
	  N = N - 1
	  RETURN
	END IF
	K = K + 1
	IF (K .GT. J) THEN
	  DO 4 M = K, J+1, -1
	    VALUES(M) = VALUES(M-1)
	    COUNT(M) = COUNT(M-1)
    4	  CONTINUE
	END IF
    5	VALUES(J) = NEWVAL
    	COUNT(J) = 1
	RETURN
C
    6	COUNT(J) = COUNT(J) + 1
	RETURN
C
	END
