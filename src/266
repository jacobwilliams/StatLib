      SUBROUTINE DIRICH(K, N, X, IX, INIT, ALPHA, RLOGL, V, G, NITER, S,
     *                  EPS, WORK, IFAULT)
C
C        ALGORITHM AS 266  APPL.STATIST. (1991), VOL.40, NO.2
C
C     Auxiliary routines required: ALOGAM (CACM algorithm 291 or AS 245),
C     DIGAMA (AS 103), GAMMAD (AS 239), PPCHI2 (AS 91), TRIGAM (AS 121).
C
      INTEGER K, N, IX, INIT, NITER, IFAULT
      REAL X(IX, K), ALPHA(K), RLOGL, V(K * (K + 1) / 2), G(K), S, EPS,
     *     WORK(2 * K)
      INTEGER I, J, IF1, KK, I2, ITN, MAXIT
      REAL ALOGAM, AN, BETA, CHI2, DIGAMA, GAMMA, GAMMAD, GG, ONE,
     *     PPCHI2, RK, SUM, SUM1, TEMP, TRIGAM, TWO, VARP1, X11, X12,
     *     ZERO
      PARAMETER (GAMMA=0.0001, ZERO=0.0, ONE=1.0, TWO=2.0, MAXIT=100)
      INTRINSIC ABS, LOG, MIN, REAL
      EXTERNAL ALOGAM, DIGAMA, GAMMAD, PPCHI2, TRIGAM
C
C        Check input arguments
C
      AN = REAL(N)
      RK = REAL(K)
      IF (K .LT. 2) THEN
         IFAULT = 1
         RETURN
      ELSE IF (N .LE. K) THEN
         IFAULT = 2
         RETURN
      ELSE IF (IX .LT. N) THEN
         IFAULT = 3
         RETURN
      END IF
C
      IFAULT = 4
      DO 20 I = 1, N
         SUM = ZERO
         NITER = I
         DO 10 J = 1, K
            IF (X(I, J) .LE. ZERO) RETURN
            SUM = SUM + X(I, J)
   10    CONTINUE
         IF (ABS(SUM - ONE) .GE. GAMMA) RETURN
   20 CONTINUE
      IFAULT = 0
C
      IF (INIT .EQ. 1) THEN
C
C        Calculate initial estimates using the method of moments
C
         SUM = ZERO
         DO 40 J = 1, K - 1
            X12 = ZERO
            DO 30 I = 1, N
               X12 = X12 + X(I, J)
   30       CONTINUE
            ALPHA(J) = X12 / AN
            SUM = SUM + ALPHA(J)
   40    CONTINUE
C
         X12 = ZERO
         DO 50 I = 1, N
            X12 = X12 + X(I, 1) ** 2
   50    CONTINUE
C
         ALPHA(K) = ONE - SUM
         X12 = X12 / AN
         VARP1 = X12 - ALPHA(1) ** 2
C
         X11 = (ALPHA(1) - X12) / VARP1
         DO 60 J = 1, K
            ALPHA(J) = X11 * ALPHA(J)
   60    CONTINUE
C
      END IF
C
      IF (INIT .EQ. 2) THEN
C
C        Calculate initial estimates using Ronning's suggestion
C
         SUM = X(1, 1)
         DO 80 J = 1, K
            DO 70 I = 1, N
               SUM = MIN(SUM, X(I, J))
   70       CONTINUE
   80    CONTINUE
         DO 90 J = 1, K
            ALPHA(J) = SUM
   90    CONTINUE
C
      END IF
C
C        Check whether any ALPHAs are negative or zero
C
      NITER = 0
      DO 100 J = 1, K
         IF (ALPHA(J) .LE. ZERO) THEN
            IFAULT = 6
            RETURN
         END IF
  100 CONTINUE
C
C     Calculate n * log(G(j)) for j = 1,2,...,k and store in WORK array
C
      DO 120 J = 1, K
         SUM = ZERO
         DO 110 I = 1, N
            SUM = SUM + LOG(X(I, J))
  110    CONTINUE
         WORK(J) = SUM
  120 CONTINUE
C
C        Note that ALOGAM cannot fail
C
      GG = ALOGAM(RK / TWO, IFAULT)
C
C        Call Algorithm AS 91 to compute chi-squared value
C
      CHI2 = PPCHI2(GAMMA, RK, GG, IFAULT)
      IF (IFAULT .GT. 0) THEN
         IFAULT = 5
         RETURN
      END IF
C
      DO 220 ITN = 1, MAXIT
C
         SUM = ZERO
         DO 130 J = 1, K
            SUM = SUM + ALPHA(J)
  130    CONTINUE
C
C        Note that TRIGAM and DIGAMA cannot fail if the first argument
C        is positive
C
         BETA = TRIGAM(SUM, IFAULT)
         TEMP = DIGAMA(SUM, IFAULT)
         SUM1 = ZERO
         DO 140 J = 1, K
            WORK(K + J) = TRIGAM(ALPHA(J), IFAULT)
            SUM1 = SUM1 + ONE / WORK(K + J)
            G(J) = AN * (TEMP - DIGAMA(ALPHA(J), IFAULT)) + WORK(J)
  140    CONTINUE
         BETA = AN * BETA / (ONE - BETA * SUM1)
C
C        Calculate the lower triangle of the Variance-Covariance matrix
C        (V)
C
         SUM = BETA / (AN * AN)
         DO 160 I = 1, K
            DO 150 J = 1, I
               KK = I * (I - 1) / 2 + J
               V(KK) = SUM / (WORK(K + I) * WORK(K + J))
               IF (I .EQ. J) V(KK) = V(KK) + ONE / (AN * WORK(K + J))
  150       CONTINUE
  160    CONTINUE
C
C        Postmultiply the Variance-Covariance matrix (V) by G and store
C        in the last k elements of WORK
C
         DO 190 I = 1, K
            SUM = ZERO
            I2 = I * (I - 1) / 2
            DO 170 J = 1, I - 1
               SUM = SUM + V(I2 + J) * G(J)
  170       CONTINUE
            DO 180 J = I + 1, K
               SUM = SUM + V(J * (J - 1) / 2 + I) * G(J)
  180       CONTINUE
            WORK(K + I) = SUM + V(I * (I + 1) / 2) * G(I)
  190    CONTINUE
C
C        Update the ALPHAs
C
         NITER = ITN
         DO 200 J = 1, K
            ALPHA(J) = ALPHA(J) + WORK(K + J)
            IF (ALPHA(J) .LE. ZERO) THEN
               IFAULT = 6
               RETURN
            END IF
  200    CONTINUE
C
C        Test for convergence
C
         S = ZERO
         DO 210 J = 1, K
            S = S + G(J) * WORK(K + J)
  210    CONTINUE
C
         IF (S .LT. CHI2) GO TO 230
C
  220 CONTINUE
C
      IFAULT = 7
C
C        Note that GAMMAD cannot fail
C
  230 EPS = GAMMAD(S / TWO, RK / TWO, IF1)
      RLOGL = ZERO
      SUM = ZERO
      DO 240 J = 1, K
         SUM = SUM + ALPHA(J)
         RLOGL = RLOGL + (ALPHA(J) - ONE) * WORK(J) -
     *           AN * ALOGAM(ALPHA(J), IF1)
  240 CONTINUE
      RLOGL = RLOGL + AN * ALOGAM(SUM, IF1)
C
      END

