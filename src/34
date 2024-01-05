      SUBROUTINE BANINV(N, K, SIGMA, S, G)
C
C     ALGORITHM AS 34 APPL. STATIST. (1970) VOL.19, NO.3
C
C     Update inverse of symmetric banded matrices
C
      INTEGER N, K
      DOUBLE PRECISION SIGMA(N), S(N), G(N)
C
C     Local variables
C
      INTEGER I, J, INDX, M, JJ, MM, NN
      DOUBLE PRECISION B, C, CG, U, V, ZERO, ONE
      DATA ZERO /0.D0/, ONE /1.D0/
C
C     Form C
C
      C = S(1)
      I = N
      DO 1 J = 2, K
	C = C - S(J) * G(I)
	I = I - 1
    1 CONTINUE
      C = ONE / C
C
C     Form SIGMA at N+1
C
      INDX = 0
      M = N * (N+1) / 2
      DO 2 J = 1, N
	CG = C * G(J)
	M = M + 1
	SIGMA(M) = -CG
	DO 2 JJ = 1, J
	  INDX = INDX + 1
	  SIGMA(INDX) = SIGMA(INDX) + CG * G(JJ)
    2 CONTINUE
      SIGMA(M+1) = C
C
C     Form G at N+1
C
      B = ZERO
      DO 4 J = 2, K
    4 B = B + S(J) * G(J-1)
      NN = N / 2
      M = N - 2 * NN
      B = B * C
      U = G(1)
      G(1) = -B
      DO 5 J = 1, NN
	MM = N - J + 1
	G(MM+1) = G(MM) + B * U
	V = G(J+1)
	G(J+1) = U + B * G(MM)
	U = V
    5 CONTINUE
      IF (M .EQ. 0) GO TO 7
      G(NN+2) = U + B * U
    7 RETURN
      END
