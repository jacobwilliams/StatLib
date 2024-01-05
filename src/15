CSTART OF AS 15
      SUBROUTINE SLINK(N, N1, N20, DELTA, B, C, DLARGE,
     $  G, H, X, W1, W2, GROUPP, TOPP, PRINTX, IFAULT)
C
C        ALGORITHM AS 15  APPL. STATIST. (1969) VOL.18, P.106
C
C        PERFORMS SINGLE LINKAGE CLUSTERING. INFORMATION IS
C        SUPPLIED BY SUBROUTINE PRTREE, IN ARRAYS B AND C.
C        POINTS ARE LISTED IN SORTED ORDER IN ARRAY G, AND THE
C        CORRESPONDING ARRAY H MARKS THE LAST MEMBER OF EACH
C        CLUSTER WITH A 1, OTHERWISE THE ENTRY IS 0. THE ARRAY
C        X STORES CODE NUMBERS FOR OUTPUT OF THE DENDROGRAM
C
      INTEGER P, Q, R, S, T, U, V, W, B(N), G(N), H(N),
     $  X(N20), W1(N1), W2(N1)
      REAL DELTA, C(N), DLARGE, D, LEVEL, ONE, ZINT
C
      DATA ONE /1.0/
C
      ZINT(D) = AINT(D)
C
      IFAULT = 1
      IF (N1 .NE. N + 1) RETURN
      IF (N20 .NE. 20 * N) RETURN
      IFAULT = 0
C
C        CLUSTERING STARTS AT THE FIRST INTEGRAL MULTIPLE OF
C        DELTA WHICH IS GREATER THAN D, THE SHORTEST LINK OF
C        THE MINIMUM SPANNING TREE
C
      D = C(2)
      IF (N .LT. 3) GOTO 15
      DO 10 I = 3, N
   10 IF (D .GT. C(I)) D = C(I)
   15 DO 20 I = 1, N
      G(I) = I
      H(I) = 1
      X(I) = 3
   20 CONTINUE
      P = 0
      LEVEL = DELTA * (ONE + ZINT(D / DELTA))
C
C        FOR EACH LINK IN ARRAY C THAT IS SHORTER THAN LEVEL,
C        TWO CLUSTERS ARE AMALGAMATED. THE AMALGAMATION
C        INVOLVES A REORDERING OF ARRAYS G AND H AND REMOVAL
C        OF THE END-OF-CLUSTER MARKER FROM THE EARLIER
C        CLUSTER. LINKS ONCE USED ARE INCREASED BY DLARGE TO
C        PREVENT RE-USE
C
   30 P = P + 1
      DO 150 I = 2, N
      IF (C(I) .GE. LEVEL) GOTO 150
      J = B(I)
      C(I) = C(I) + DLARGE
      K = I
      DO 40 M = 1, N
      IF (G(M) .EQ. J) Q = M
      IF (G(M) .EQ. K) R = M
   40 CONTINUE
      IF (Q .LE. R) GOTO 50
      M = R
      R = Q
      Q = M
   50 DO 60 S = Q, N
      IF (H(S) .NE. 0) GOTO 70
   60 CONTINUE
   70 T = R
   80 T = T - 1
      IF (T .LT. 1) GOTO 90
      IF (H(T) .EQ. 0) GOTO 80
   90 T = T + 1
      H(S) = 0
      L = 0
      DO 100 R = T, N
      L = L + 1
      W1(L) = G(R)
      W2(L) = H(R)
      IF (H(R) .NE. 0) GOTO 110
  100 CONTINUE
  110 W = S + 1
      M = T
      L = R + 1
  120 M = M - 1
      IF (M .LT. W) GOTO 130
      L = L - 1
      G(L) = G(M)
      H(L) = H(M)
      GOTO 120
  130 U = R - T + 1
      L = W - 1
      DO 140 M = 1, U
      L = L + 1
      G(L) = W1(M)
      H(L) = W2(M)
  140 CONTINUE
  150 CONTINUE
      CALL GROUPP(LEVEL, N, G, H)
      W = N * P
      U = 0
      V = 0
      K = 0
C
C        DENDROGRAM INDICATORS ARE NOW COMPILED AND STORED
C        IN X. POINTS ARE EXAMINED IN THE ORDER DEFINED BY
C        THE ARRAY G. S IS THE CORRESPONDING ENTRY ON THE
C        PREVIOUS ITERATION, U = 0 FOR THE FIRST MEMBER OF
C        A CLUSTER, V = 1 WHEN AMALGAMATIONS OCCUR, FROM
C        THE LAST MEMBER OF THE FIRST COMPONENT CLUSTER
C        UNTIL THE LAST MEMBER OF THE AMALGAMATED CLUSTER.
C        K IS THE TOTAL NUMBER OF CLUSTERS
C
      DO 160 I = 2, N
  160 K = K + H(I)
      IF (P .GT. 19) GOTO 270
      DO 260 I = 1, N
      J = G(I)
      L = J + W - N
      S = X(L)
      IF (U .NE. 0) GOTO 190
      IF (H(I) .NE. 1) GOTO 170
      T = 3
      GOTO 250
  170 IF (S .NE. 3) GOTO 180
      T = 1
      U = 1
      V = 1
      GOTO 250
  180 T = 0
      U = 1
      GOTO 250
  190 IF (H(I) .NE. 1) GOTO 210
      IF (V .NE. 0) GOTO 200
      T = 3
      U = 0
      GOTO 250
  200 T = 2
      U = 0
      V = 0
      GOTO 250
  210 IF (S .LT. 2 .OR. S .GT. 3) GOTO 230
      IF (V .NE. 0) GOTO 220
      T = 1
      U = 1
      V = 1
      GOTO 250
  220 T = 5
      U = 1
      GOTO 250
  230 IF (V .NE. 0) GOTO 240
      T = 0
      U = 1
      GOTO 250
  240 T = 4
  250 L = J + W
      X(L) = T
  260 CONTINUE
  270 LEVEL = LEVEL + DELTA
      IF (K .NE. 1) GOTO 30
      CALL TOPP
      DO 280 I = 1, N
  280 CALL PRINTX(I, G, N, P, X, N20)
      RETURN
      END
CEND OF AS 15
