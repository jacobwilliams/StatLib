C        .. ULCC Toolpack/1 2.1
      DOUBLE PRECISION FUNCTION BINOM(N, LP1, LP2, LMP1, LMP2, X, R)
C
C        FAST BIVARIATE BINOMIAL DENSITY
C
      DOUBLE PRECISION ALNGAM, LP1, LP2, LMP1, LMP2, ONE, TWO
      INTEGER N, X, R, IER
      DATA ONE / 1.0 / , TWO / 2.0 /
C
      BINOM = DEXP(TWO * ALNGAM(N + ONE, IER) -
     *        (ALNGAM(X + ONE, IER) + ALNGAM(N - X + ONE,
     *        IER) + ALNGAM(R - X + ONE, IER) + ALNGAM(N - (R - X) +
     *        ONE, IER)) + X * LP1 + (R - X) * LP2 + (N - X) * LMP1 +
     *        (N - (R - X)) * LMP2)
      RETURN
      END
 
 
      INTEGER FUNCTION XFMAX(N, R, THETA)
C
C        RETURNS X SUCH THAT F(X,R) IS A MAXIMUM
C
      INTEGER N, R
      DOUBLE PRECISION THETA, A, B, C, Q, ONE, HALF, FOUR
      DATA ONE / 1.0 / , TWO / 2.0 / , HALF / 0.5 / , FOUR / 4.0 /
C
      IF (THETA .NE. ONE) THEN
         A = - (ONE - THETA)
         B = - ((R + N + TWO) * THETA + N - R)
         C = (N + ONE) * (R + ONE) * THETA
         Q = -HALF * (B - DSQRT(B * B - (FOUR * A * C)))
         XFMAX = INT(C / Q)
      ELSE
         XFMAX = (R + 1) / 2
      END IF
      RETURN
      END
 
 
      INTEGER FUNCTION CV(N, R, ALPHA)
C
C        RETURNS LARGEST X SUCH THAT PROB(X<=X) < ALPHA
C
      INTEGER N, R, XMIN, X, XSTART, IER
      DOUBLE PRECISION ALNGAM, ALPHA, HYPSUM, LNUM, LDEN, PROB, ONE,
     *                 HALF, TWO
      DATA ONE / 1.0 / , HALF / 0.5 / , TWO / 2.0 /
C
      LDEN = ALNGAM(2 * N + ONE, IER) - TWO * ALNGAM(N + ONE, IER) -
     *       (ALNGAM(R + ONE, IER) + ALNGAM(2 * N - R + ONE, IER))
      XSTART = (R - 1) / 2
      X = XSTART + 1
      LNUM = - (ALNGAM(R - X + ONE, IER) +
     *       ALNGAM(N - (R - X) + ONE, IER) + ALNGAM(X + ONE, IER) +
     *       ALNGAM(N - X + ONE, IER))
      PROB = DEXP(LNUM - LDEN)
      HYPSUM = HALF
      IF (MOD(R, 2) .EQ. 0) HYPSUM = HYPSUM + HALF * PROB
      XMIN = 0
      IF (R .GT. N) XMIN = R - N
C
C        HYPERGEOMETRIC ADJACENCY
C
      DO 10 X = XSTART, XMIN, -1
         PROB = PROB * (X + ONE) * (N - R + X + ONE) /
     *          ((R - X) * ONE * (N - X))
         HYPSUM = HYPSUM + PROB
         IF (HYPSUM .GT. ONE - ALPHA) GO TO 20
   10 CONTINUE
      X = XMIN
C
C        IF NO CRITICAL VALUE, RETURN -1
C
   20 CV = -1
      IF (X - 1 .GE. XMIN) CV = X - 1
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION FSUM(FROM, TO, INC, F, FLIM, N, THETA,
     *                               R, RCNT)
C
C        SUMS BINOMIAL PROBABILITIES ALONG FIXED R
C
      DOUBLE PRECISION FP, F, FLIM, THETA, ONE, ZERO
      INTEGER X, N, FROM, TO, R, INC, RCNT
      DATA ONE / 1.0 / , ZERO / 0.0 /
C
      FP = ZERO
      RCNT = 0
      X = FROM
   10 IF (F .LT. FLIM) GO TO 20
      RCNT = RCNT + 1
      FP = FP + F
C
C        BINOMIAL ADJACENCY
C
      IF (INC .GT. 0) F = F * THETA * (N - X) * (R - X) /
     *                    ((X + ONE) * (N - R + X + 1))
      IF (INC .LE. 0) F = F * (ONE / THETA) * X * (N - R + X) /
     *                    ((N - X + ONE) * (R - X + 1))
      X = X + INC
      IF (X .NE. TO + INC) GO TO 10
   20 FSUM = FP
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION FXPOWER(N, P1, P2, ALPHA, EPS, ERROR,
     *                 COUNT, TOTAL, IFAULT)
C
C        POWER FOR FISHER'S EXACT TEST
C
C        RETURN THE POWER OF A LEVEL ALPHA TEST FOR A GIVEN
C        COMMON SAMPLE SIZE AND TRUE PROBABILITIES OF EVENTS P1 AND P2.
C
      DOUBLE PRECISION P1, P2, ALPHA, EPS, ERROR, LP1, LP2, LMP1, LMP2,
     *                 F, FLIM, POWSUM, THETA, FMAX, FSUM, BINOM, ONE,
     *                 ZERO, TWO
      INTEGER XM, N, TOTAL, COUNT, RMAX, XMAX, YMAX, XMIN, XCV, XFM, R,
     *        I, J, RCNT, CV, XFMAX
      DATA ONE / 1.0 / , ZERO / 0.0 / , TWO / 2.0 /
C
      IFAULT = 0
      ERROR = ZERO
      POWSUM = ZERO
      COUNT = 0
      TOTAL = 0
C
C        CHECK INPUT PARMS
C
      IF (N .LE. 0 .OR. P1 .LE. ZERO .OR. P1 .GE. ONE .OR.
     *    P2 .LE. ZERO .OR. P2 .GE. ONE .OR. ALPHA .LE. ZERO .OR.
     *    ALPHA .GE. ONE .OR. EPS .LE. ZERO .OR. EPS .GE. ONE) THEN
         IFAULT = 1
         RETURN
      END IF
C
      LP1 = DLOG(P1)
      LP2 = DLOG(P2)
      LMP1 = DLOG(ONE - P1)
      LMP2 = DLOG(ONE - P2)
      XMAX = (N + 1) * P1
      YMAX = (N + 1) * P2
      RMAX = XMAX + YMAX
      FMAX = BINOM(N, LP1, LP2, LMP1, LMP2, XMAX, RMAX)
      FLIM = FMAX * EPS
      THETA = P1 * (ONE - P2) / (P2 * (ONE - P1))
      XM = 0
      IF (RMAX .GT. N) XM = RMAX - N
      TOTAL = (N + 1 - RMAX + TWO * (CV(N, RMAX, ALPHA) - XM)) *
     *        (N + 2 - RMAX + TWO * (CV(N, RMAX, ALPHA) - XM)) / 2
      R = RMAX
      I = 0
      J = 1
C
C        MAIN LOOP ON R
C
   10 XMIN = 0
      IF (R .GT. N) XMIN = R - N
      XCV = CV(N, R, ALPHA)
      XFM = XFMAX(N, R, THETA)
      IF (XCV .LE. 0) GO TO 30
      IF (XCV .GT. XFM) GO TO 20
C
C        SUM FROM CRITICAL VALUE DOWN TO MINIMUM
C
      F = BINOM(N, LP1, LP2, LMP1, LMP2, XCV, R)
      POWSUM = POWSUM + FSUM(XCV, XMIN, -1, F, FLIM, N, THETA, R, RCNT)
      IF (RCNT .EQ. 0) GO TO 40
      COUNT = COUNT + RCNT
      GO TO 30
C
C        SUM FROM DENSITY MAX DOWN TO MINIMUM
C
   20 F = BINOM(N, LP1, LP2, LMP1, LMP2, XFM, R)
      POWSUM = POWSUM + FSUM(XFM, XMIN, -1, F, FLIM, N, THETA, R, RCNT)
      IF (RCNT .EQ. 0) GO TO 40
      COUNT = COUNT + RCNT
C
C        SUM FROM DENSITY MAX UP TO CRITICAL VALUE
C
      F = BINOM(N, LP1, LP2, LMP1, LMP2, XFM + 1, R)
      POWSUM = POWSUM + FSUM(XFM + 1, XCV, 1, F, FLIM, N, THETA, R,
     *         RCNT)
      COUNT = COUNT + RCNT
C
C        WORK OUT FROM RMAX ALTERNATING BETWEEN SMALLER AND LARGER R
C
   30 I = I + 1
      J = -1 * J
      R = R + I * J
      IF (R .GE. 0 .AND. R .LE. 2 * N) GO TO 10
C
C        CALCULATE ERROR BOUND AND RETURN
C
   40 ERROR = ONE - POWSUM
      IF (ONE - POWSUM .GT. (TOTAL - COUNT) * FLIM) ERROR = FLIM *
     *    (TOTAL - COUNT)
      FXPOWER = POWSUM
      RETURN
      END
