      SUBROUTINE PPRANK(TA, UA, ET, SD, K3, K4, N, M, ISEL, IFAULT)
C
C        ALGORITHM AS234  APPL. STATIST. (1988) VOL. 37, NO. 2
C
C        Declaration formal parameters
C
      INTEGER N, M, ISEL, IFAULT
      REAL TA, UA, ET, SD, K3, K4
C
C        Auxiliary algorithm
C
      REAL SCORE
      EXTERNAL SCORE
C
C        Constants
C
      REAL ZERO, ONE, TWO, THREE, SIX, HALF, TWOHLF
      REAL CMT1, CMT2, CMT3, CTA1, CTA2, CTA3
      DATA ZERO, ONE, TWO, THREE /0.0E0, 1.0E0, 2.0E0, 3.0E0/
      DATA SIX, HALF, TWOHLF /6.0E0, 0.5E0, 2.5E0/
      DATA CMT1, CMT2, CMT3 /0.0833333E0, 2.3333333E0, 0.0125E0/
      DATA CTA1, CTA2, CTA3 /0.4166667E0, 0.0555555E0, 0.1666667E0/
C
C        Declaration internal variables
C
      INTEGER I
      REAL S1A, S2A, S3A, S4A, S1C, S2C, S3C, S4C, W, H1, H2, H3, H4,
     *  MAXU
      DATA MAXU /1.0E+10/
C
C        Checking input
C
      IFAULT = 1
      IF ((ISEL .LE. 0) .OR. (ISEL .GT. 4)) RETURN
      IF (ISEL .EQ. 4) GOTO 400
      IFAULT = 2
      IF (N .LE. 3) RETURN
C
C        Initializing auxiliary variables
C
      H1 = N - THREE
      H2 = N - TWO
      H3 = N - ONE
      H4 = N + ONE
C
C        Compute S1a, S2A, S3a, S4a
C
      S1A = ZERO
      DO 10 I = 1, N
      S1A = S1A + SCORE(I / H4)
   10 CONTINUE
      S1A = S1A / N
      S2A = ZERO
      S3A = ZERO
      S4A = ZERO
      DO 20 I = 1, N
      W = SCORE(I / H4) - S1A
      S2A = S2A + W ** 2
      S3A = S3A + W ** 3
      S4A = S4A + W ** 4
   20 CONTINUE
      IFAULT = 3
      IF (S2A .LE. ZERO) RETURN
      S2A = S2A / H3
      S3A = S3A / H3
      S4A = S4A / H3
      GOTO (100, 200, 300), ISEL
C
C        Two-sample
C
  100 IFAULT = 4
      IF (N .LT. 2 * M) RETURN
      S1C = (ONE * M) / N
      W = ONE - S1C
      S2C = (M * W ** 2 + (N - M) * S1C ** 2) / H3
      S3C = (M * W ** 3 + (N - M) * S1C ** 3) / H3
      S4C = (M * W ** 4 + (N - M) * S1C ** 4) / H3
      GOTO 500
C
C        Monotone trend
C
  200 S1C = H4 * HALF
      S2C = H4 * N * CMT1
      S3C = ZERO
      S4C = (N * N - CMT2) * H4 * N * CMT3
      GOTO 500
C
C        Independence
C
  300 S1C = S1A
      S2C = S2A
      S3C = S3A
      S4C = S4A
      GOTO 500
C
C        Compute ET, SD, K3, K4
C
  500 W = (N / H3) * (H4 / H3)
      ET = N * S1A * S1C
      SD = SQRT(H3 * S2A * S2C)
      K3 = (N / H2) * (S3A / S2A) * (S3C / S2C) / SD
      K4 = (H3 / H1) * (W * S4A / (S2A ** 2) - THREE) * (W * S4C /
     *  (S2C ** 2) - THREE) / (W * H2) - SIX / H4
C
C        Compute percentage point
C
  400 IFAULT = 5
      IF (ABS(UA) .GT. MAXU) RETURN
      W = UA ** 2
      TA = (ONE + (W - THREE) * K4 * CTA1 - (W - TWOHLF) * (K3 ** 2)
     *  * CTA2) * UA + (W - ONE) * K3 * CTA3
      TA = TA * SD + ET
      IFAULT = 0
      RETURN
      END
