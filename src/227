	SUBROUTINE GCOUNT(N, APPLY, IFAULT)
C
C	ALGORITHM AS227  APPL. STATIST. (1987) VOL. 36, NO. 2
C
C	Generates all possible N-bit binary codes, and applies a users
C	procedure for each code generated.
C	N must be <= NMAX, which is set in the PARAMETER statement below.
C	IFAULT = 2 is returned otherwise.
C
C	Translated from Algol 60.
C
	INTEGER N, IFAULT
	EXTERNAL APPLY
C
C	Local variables
C
	INTEGER NMAX
	PARAMETER (NMAX = 100)
	INTEGER CHANGE, I, TPOINT(NMAX)
	LOGICAL STATUS(NMAX)
C
	IF (N .LT. 1) THEN
	  IFAULT = 1
	  RETURN
	END IF
	IF (N .GT. NMAX) THEN
	  IFAULT = 2
	  RETURN
	END IF
	IFAULT = 0
C
C	Initialize and make first call to user's routine.
C
	DO 10 I = 1, N
	  STATUS(I) = .FALSE.
	  TPOINT(I) = I + 1
   10	CONTINUE
	CALL APPLY(N, N, STATUS)
C
C	Generate a new code.   The user's routine is called twice each
C	cycle; the first time the bit which changes is bit 1.
C
   20	IF (STATUS(1)) THEN
	  STATUS(1) = .FALSE.
	  CHANGE = TPOINT(2)
	ELSE
	  STATUS(1) = .TRUE.
	  CHANGE = 2
	END IF
	CALL APPLY(N, 1, STATUS)
C
C	Check if count exhausted.
C
	IF (CHANGE .GT. N) RETURN
C
	IF (STATUS(CHANGE)) THEN
	  STATUS(CHANGE) = .FALSE.
	  TPOINT(CHANGE) = TPOINT(CHANGE + 1)
	ELSE
	  STATUS(CHANGE) = .TRUE.
	  TPOINT(CHANGE) = CHANGE + 1
	END IF
	CALL APPLY(N, CHANGE, STATUS)
C
	GO TO 20
	END
