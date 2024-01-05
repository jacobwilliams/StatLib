      SUBROUTINE BAB(N, M, SI, T, S, A, MAX, CRIT, BEST, NEVAL, IFAULT)
C
C     ALGORITHM AS  199  APPL. STATIST. (1984) VOL.33, NO.2
C
C     Efficient branch and bound algorithm for determining the optimal
C     feature subset of given size.
C
C     Incorporates a corrected version of the corrections in ASR 67.
C     The second of those corrections caused a jump into a DO loop.
C
      REAL SI(N,N), T(N), S(N,N,N), CRIT, ZERO, D
      INTEGER A(N), MAX(M), BEST(M), NEVAL(N), AI8, AJ8, OMIT
      LOGICAL SKIP
C
      DATA ZERO /0.0/
C
C     Check for invalid values of M and N
C
      IFAULT = 1
      IF (M .LT. 2 .OR. N .LE. M) RETURN
      IFAULT = 0
C
C     Initialize arrays - A holds current M-subset, MAX holds upper
C     limit for possible member features in M-subset.
C     Initial optimal criterion is CRIT = 0
C     SKIP = .FALSE. indicates that current position on solution tree
C     follows consecutively in enumeration sequence rather than a jump
C     from a rejected branch.
C
      DO 4 I = 1, N
	DO 3 J = 1, N
    3   S(I,J,N) = SI(I,J)
	NEVAL(I) = 0
    4 CONTINUE
      SKIP = .FALSE.
      CRIT = ZERO
      DO 5 I = 1, M
	A(I) = I
	MAX(I) = N - M + I
	BEST(I) = 0
    5 CONTINUE
   10 K1 = A(M)
      K2 = MAX(M)
      DO 30 I = K1, K2
C
C     Enumerate members of next subset
C
	DO 13 J = I, N
	  JJ = M - I + J
	  A(JJ) = J
   13   CONTINUE
	MM = M + N - I
C
C     If subset of size N, skip to calculation of criterion for
C     subsequent descendent subsets.
C     If consecutive enumeration (SKIP = .FALSE.) next feature omitted
C     is A(M).
C     If following jump in enumeration sequence (SKIP = .TRUE.) next
C     feature omitted is A(M-L).
C
	IF (MM .EQ. N) GO TO 22
	OMIT = M
	IF (.NOT. SKIP) GO TO 15
	OMIT = M - L
	SKIP = .FALSE.
C
C     Evaluate criterion function for current subset from that of
C     previous subset, excluding feature A(OMIT).
C
   15   M9 = MM + 1
	D = ZERO
	DO 17 I8 = 1, MM
	  I9 = I8
	  AI8 = A(I8)
	  IF (I8 .GE. OMIT) I9 = I8 + 1
	  DO 16 J8 = 1, MM
	    J9 = J8
	    AJ8 = A(J8)
	    IF (J8 .GE. OMIT) J9 = J8 + 1
	    S(I8,J8,MM) = S(I9,J9,M9) - S(OMIT,I9,M9) * S(OMIT,J9,M9) /
     *                     S(OMIT,OMIT,M9)
	    D = D + T(AI8) * T(AJ8) * S(I8,J8,MM)
   16     CONTINUE
   17   CONTINUE
	NEVAL(MM) = NEVAL(MM) + 1
C
C     Check criteraion of current subset, D, relative to optimal value
C     CRIT.
C
	IF (D .GT. CRIT .AND. MM .EQ. M) GO TO 28
	IF (D .GT. CRIT) GO TO 22
C
C     Optimal value not exceeded, so abandon downward search of branch.
C     Prepare for return to next lower level branching mode not yet
C     explored.   SKIP = .TRUE. indicates jump to lower level before
C     termination of branch.
C
   60   M1 = M - 1
	DO 18 L = 1, M1
	  ML = M - L
	  IF (A(ML) .LT. A(ML+1) - 1) GO TO 19
   18   CONTINUE
	L = M
   19   DO 20 LL = 1, L
	ML = M - LL + 1
	A(ML) = MAX(ML)
   20   CONTINUE
	SKIP = .TRUE.
	GO TO 35
C
C     Optimal value CRIT exceeded by lower level subset.   Reduce by
C     one feature at a time until either criterion falls below CRIT or
C     subset size reaches M.
C
   22   MX = MM - M
	DO 27 MO = 1, MX
	  D = ZERO
	  M8 = MM - MO
	  M9 = M8 + 1
	  DO 24 I8 = 1, M8
	    I9 = I8
	    IF (I8 .GT. M) I9 = I8 + 1
	    A(I8) = A(I9)
	    DO 23 J8 = 1, M8
	      J9 = J8
	      IF (J8 .GT. M) J9 = J8 + 1
	      S(I8,J8,M8) = S(I9,J9,M9) - S(M+1,I9,M9) * S(M+1,J9,M9) /
     *                        S(M+1,M+1,M9)
   23       CONTINUE
   24     CONTINUE
	  DO 26 I8 = 1, M8
	    AI8 = A(I8)
	    DO 25 J8 = 1, M8
	      AJ8 = A(J8)
	      D = D + T(AI8) * T(AJ8) * S(I8,J8,M8)
   25       CONTINUE
   26     CONTINUE
	  NEVAL(M8) = NEVAL(M8) + 1
	  IF (D .LE. CRIT) GO TO 30
   27   CONTINUE
C
C     Update optimal criterion function, CRIT, and corresponding
C     M-subset BEST(.).
C
   28   CRIT = D
	DO 29 L1 = 1, M
   29   BEST(L1) = A(L1)
	IF (I .EQ. K2) GO TO 60
   30 CONTINUE
C
C     Update feature subset required for return to lower level branching
C     mode of tree.   Search complete if first member A(1) has taken all
C     its possible feature values.
C
   35 DO 50 I = 2, M
	II = M - I + 1
	IF (A(II) .EQ. MAX(II)) GO TO 50
	A(II) = A(II) + 1
	I1 = II + 1
	DO 40 J = I1, M
   40   A(J) = A(J-1) + 1
	GO TO 10
   50 CONTINUE
C
C     Successful completion of search
C
      RETURN
      END
