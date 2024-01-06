CSTART OF AS 3
      REAL FUNCTION PROBST(T, IDF, IFAULT)
C
C        ALGORITHM AS 3  APPL. STATIST. (1968) VOL.17, P.189
C
C        STUDENT T PROBABILITY (LOWER TAIL)
C
      REAL A, B, C, F, G1, S, FK, T, ZERO, ONE, TWO, HALF, ZSQRT, ZATAN
C
C        G1 IS RECIPROCAL OF PI
C
      DATA ZERO, ONE, TWO, HALF, G1
     $     /0.0, 1.0, 2.0,  0.5, 0.3183098862/
C
      ZSQRT(A) = SQRT(A)
      ZATAN(A) = ATAN(A)
C
      IFAULT = 1
      PROBST = ZERO
      IF (IDF .LT. 1) RETURN
      IFAULT = 0
      F = IDF
      A = T / ZSQRT(F)
      B = F / (F + T ** 2)
      IM2 = IDF - 2
      IOE = MOD(IDF, 2)
      S = ONE
      C = ONE
      F = ONE
      KS = 2 + IOE
      FK = KS
      IF (IM2 .LT. 2) GOTO 20
      DO 10 K = KS, IM2, 2
      C = C * B * (FK - ONE) / FK
      S = S + C
      IF (S .EQ. F) GOTO 20
      F = S
      FK = FK + TWO
   10 CONTINUE
   20 IF (IOE .EQ. 1) GOTO 30
      PROBST = HALF + HALF * A * ZSQRT(B) * S
      RETURN
   30 IF (IDF .EQ. 1) S = ZERO
      PROBST = HALF + (A * B * S + ZATAN(A)) * G1
      RETURN
      END
CEND OF AS 3
