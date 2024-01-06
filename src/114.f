      SUBROUTINE PANDQ(S, F, R, C, ROWSUM, COLSUM, N, IFAULT)
C
C     ALGORITHM AS 114  APPL. STATIST. (1977) VOL. 26, NO. 2
C
C     Finds the numerator S required by several measures of association.
C
      INTEGER S, R, C, F(R,C), ROWSUM(R), COLSUM(C), N, IFAULT
C
C     Local variables
C
      INTEGER P, Q, FI1, FIJ, I, J, NC, ND
C
      IFAULT = 1
      IF (R .LT. 2 .OR. C .LT. 2) RETURN
      IFAULT = 0
C
      ROWSUM(1) = 0
      DO 10 J = 1, C
	COLSUM(J) = F(1,J)
	ROWSUM(1) = ROWSUM(1) + F(1,J)
   10 CONTINUE
      NC = 0
      P = 0
      Q = 0
      DO 20 I = 2, R
	FI1 = F(I,1)
	ND = NC + ROWSUM(I-1) - COLSUM(1)
	Q = Q + FI1 * ND
	NC = COLSUM(1)
	COLSUM(1) = COLSUM(1) + FI1
	ROWSUM(I) = FI1
	DO 15 J = 2, C
	  FIJ = F(I,J)
	  ND = ND - COLSUM(J)
	  Q = Q + FIJ * ND
	  P = P + FIJ * NC
	  NC = NC + COLSUM(J)
	  COLSUM(J) = COLSUM(J) + FIJ
	  ROWSUM(I) = ROWSUM(I) + FIJ
   15   CONTINUE
   20 CONTINUE
      N = NC + ROWSUM(R)
      S = P - Q
      RETURN
      END
