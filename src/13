CSTART OF AS 13
      SUBROUTINE PRTREE(N, M, A, DLARGE, D, B, C, IFAULT)
C
C        ALGORITHM AS 13  APPL. STATIST. (1969) VOL.18, P.103
C
C        COMPUTES THE MINIMUM SPANNING TREE OF A DISTANCE MATRIX
C
      REAL D(M), C(N), AM, DIST, DLARGE
      INTEGER B(N)
      LOGICAL A(N)
      IFAULT = 1
      IF (N .LT. 2) RETURN
      IFAULT = 2
      IF (N * (N - 1) / 2 .NE. M) RETURN
      IFAULT = 0
C
C        A(I) IS .FALSE. IF I IS ALREADY ASSIGNED TO
C        THE TREE (INITIALLY CONSISTING OF NO. 1 ONLY),
C        OR .TRUE. OTHERWISE
C
      DO 10 I = 2, N
      A(I) = .TRUE.
      B(I) = 0
      C(I) = DLARGE
   10 CONTINUE
      J = 1
      DO 40 I = 2, N
      AM = DLARGE
      DO 30 K = 2, N
      IF (.NOT. A(K)) GOTO 30
      IF (J .GT. K) L = (J - 1) * (J - 2) / 2 + K
      IF (J .LE. K) L = (K - 1) * (K - 2) / 2 + J
      DIST = D(L)
      IF (DIST .GE. C(K)) GOTO 20
      C(K) = DIST
      B(K) = J
   20 IF (AM .LE. C(K)) GOTO 30
      AM = C(K)
      NEXT = K
   30 CONTINUE
      J = NEXT
      A(J) = .FALSE.
   40 CONTINUE
      RETURN
      END
CEND OF AS 13
