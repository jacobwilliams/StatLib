      SUBROUTINE AMALGM(K, XO, WO, X, W, XA, IFAULT)
C
C<<<<<  Acquired in machine-readable form from 'Applied Statistics'
C<<<<<  algorithms editor, January 1983.
C
C
C        ALGORITHM AS 149 APPL. STATIST. (1980) VOL.29, NO.2
C
C        AMALGAMATION  OF  MEANS  BY  THE  UP-AND-DOWN  BLOCKS  ALGORITHM
C        OF  KRUSKAL  (BARLOW  ET  AL. ,1972, P.72 )
C
      DIMENSION XO(K),WO(K),X(K),W(K),XA(K)
      DATA TOL /1.0E-6/
      IFAULT = 1
C
C        CHECK  THAT  K .GE. 2
C
      IF (K .LT. 2) RETURN
      IFAULT = 2
C
C        CHECK  THAT  THE  WEIGHTS  ARE  POSITIVE
C
      DO 1 I = 1,K
      IF (WO(I) .LE. 0.0) RETURN
    1 CONTINUE
      IFAULT = 0
      DO 2 I = 1,K
      X(I) = XO(I)
      W(I) = WO(I)
    2 CONTINUE
      M = K
      I = 1
    3 IF (I .EQ. M) GOTO 4
      IF (X(I) .GT. X(I+1)) GOTO 9
      IF (I .EQ. 1) GOTO 13
    4 IF (X(I-1) .GT. X(I)) GOTO 6
      IF (I .LT. M) GOTO 13
      GOTO 14
C
C        POOL  THE  ACTIVE  BLOCK  WITH  THE  NEXT  LOWER  BLOCK
C
    6 IM1 = I - 1
      WW = W(IM1) + W(I)
      X(IM1) = (W(IM1)*X(IM1) + W(I)*X(I))/WW
      W(IM1) = WW
      MM1 = M  - 1
      IF (I .EQ. M) GOTO 8
      DO 7 J = I,MM1
      J1 = J + 1
      X(J) = X(J1)
      W(J) = W(J1)
    7 CONTINUE
    8 I = IM1
      M = MM1
      IF (M .EQ. 1) GOTO 14
      GOTO 3
C
C        POOL  THE  ACTIVE  BLOCK  WITH  THE  NEXT  HIGHER  BLOCK
C
    9 I1 = I + 1
      WW = W(I) + W(I1)
      X(I) = (W(I)*X(I) + W(I1)*X(I1))/WW
      W(I) = WW
      MM1 = M - 1
      IF (I1 .EQ. M) GOTO 11
      DO 10 J = I1,MM1
      J1 = J + 1
      X(J) = X(J1)
      W(J) = W(J1)
   10 CONTINUE
   11 M = MM1
      IF (M .EQ. 1) GOTO 14
      IF (I .EQ. 1) GOTO 12
      IF (X(I-1) .GT. X(I)) GOTO 6
   12 IF (I .EQ. M) GOTO 14
      IF (X(I) .GT. X(I+1)) GOTO 9
   13 I = I + 1
      GO TO 12
C
C        OBTAIN  THE  AMALGAMATED  MEANS  XA(K) FROM  THE  WORKING  ARRAY  X(M)
C
   14 I1 = 1
      DO  17 I = 1,M
      S = 0.0
      DO 15 J = I1,K
      S = S + WO(J)
      XA(J) = X(I)
      IF (ABS(S - W(I)) .LT. TOL) GOTO 16
   15 CONTINUE
   16 I1 = J + 1
   17 CONTINUE
      RETURN
      END
