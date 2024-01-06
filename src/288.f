      SUBROUTINE GSMIRN (NX, NY, KIND, M, DSTAT, P, Q, IFAULT)
C
C  ALGORITHM AS 288 APPL.STATIST. (1994), VOL.43, NO.1
C
C  P-value calculation for the generalized two-sample 
C  Smirnov tests.
C  The tests are conditional on ties in the pooled sample.
C
      DIMENSION P(*), M(*)
      REAL P, DSTAT, Q, SLOPE, DELTA, DEVIAT, DL, ONE, ZERO, EPS, FLO
     *, SMALL, SCL, SMALLN, T, FL, ALN2, CHKNUM
      LOGICAL NEWRCT
      DATA ONE /1.0/ ZERO /0.0/ EPS /1E-6/
      DATA SMALL /1E-35/ SMALLN /-80.5904782547916/
      DATA ALN2 /0.69314718056/ CHKNUM /1E32/ ITERUP /116/
      FLO(L) = FLOAT(L)
C
      N = NX + NY
      IFAULT = 1
      IF (NX .LE. 0 .OR. NY .LE. 0) RETURN
      IFAULT = 2
      IF (KIND .LT. 1 .OR. KIND .GT. 3) RETURN
      IFAULT = 4
      ICAT = 0
      L = 0
    1 L = L + 1
      IF (M(L) .LE. 0) RETURN
      ICAT = ICAT + M(L)
      IF (ICAT - N) 1, 3, 2
    2 RETURN
C
    3 IFAULT = 0
      Q = ONE
      DELTA = DSTAT - EPS
      IF (DELTA .LE. ZERO) RETURN
      P(1) = ONE
C
C  Parameters to define a set for trajectories to lie within it
C
      SLOPE = NX / FLO(N)
      DEVIAT = SLOPE * DELTA * NY
C*    DELTA = DEVIAT
      NEWRCT = .TRUE.
      ICAT = 1
      NTIES = M(1)
      ILE = 0
      IRI = 0
C
C  Variables to prevent from overflows in P
C
      IC = ITERUP
      NOFDIV = 0
      SCL = ONE
C
      DO 100 L = 1, N - 1
        IF (NTIES .EQ. 1) THEN
C
C  Calculate boundaries for current L (if Lth value in the
C  pooled sample is unique or `last' in a series of tied values)
C
          DL = L * SLOPE
C*        T = FLOAT (L) / N
C*        DEVIAT = DELTA * SQRT (T * (ONE - T))
          IRI = MIN0 (INT (DL + DEVIAT), L, NX)
          ILE = MAX0 (INT (DL - DEVIAT + ONE), L - NY, 0)
          ICAT = ICAT + 1
          NTIES = M(ICAT)
          NEWRCT = .TRUE.
        ELSE
C
C      Calculations for tied observations
C
          NTIES = NTIES - 1
C
C  If we have the first observation with the new value (that
C  is not unique), then determine the ICATth `rectangle'
C
          IF (NEWRCT) THEN
            NEWRCT = .FALSE.
            L2 = L + NTIES
C*          IF (L2 .EQ. N) THEN
C*            IRI2 = NX
C*            ILE2 = NX
C*          ELSE
              DL = L2 * SLOPE
C*            T = FLOAT (L2) / N
C*            DEVIAT = DELTA * SQRT (T * (ONE - T))
C
C  X axis boundaries of the subset of `line' L(l+1) within 
C  ICATth rectangle
C
              IRI2 = MIN0 (INT (DL + DEVIAT), L2, NX)
              ILE2 = MAX0 (INT (DL - DEVIAT + ONE),  L2 - NY, 0)
C*          ENDIF
C
C  Four sides of the rectangle on the X and Y axes
C
            ILEFT = ILE
            IRIGHT = IRI2
            JUPP = L2 - ILE2
            JLOW = L - IRI - 1
          ENDIF
C
C  Calculate boundaries for current L (Lth value is tied)
C
          ILE = MAX0 (ILEFT, L - JUPP)
          IRI = MIN0 (IRIGHT, L - JLOW)
        ENDIF
C
C  Set the left (right) boundary for the one-sided test
C
        GOTO (30, 20, 10) KIND
   10   IRI = MIN0 (NX, L)
        GOTO 30
   20   ILE = MAX0 (0, L - NY)
C
C  Calculate the number of trajectories p(i,j) for current L
C
   30   ILES = MAX0 (1, ILE)
        IRIS = MIN0 (L - 1, IRI)
 
        DO 50 I = IRIS, ILES, - 1
   50   P(I + 1) = P(I + 1) + P(I)
C
C  Check whether elements of P are large enough to multiply 
C  them by SMALL
C
        IC = IC - 1
        IF (IC .LE. 0) THEN
          DL = ZERO
          DO 60 I = ILES + 1, IRIS + 1
   60     DL = MAX (P(I), DL)
          IF (DL .EQ. ZERO) RETURN
          IF (DL .GT. CHKNUM) THEN
            DO 65 I = ILES + 1, IRIS + 1
   65       P(I) = P(I) * SMALL
            IC = ITERUP
            NOFDIV = NOFDIV + 1
            SCL = SCL * SMALL
          ELSE
C
C  Estimate the number of iterations for DL=Pmax to became 
C  of order 1/SMALL
C
            IC = (- SMALLN - LOG (DL)) / ALN2
          END IF
        END IF
C
C  Define, whether boundaries lie on the left and lower sides
C  of the rectangle R and define boundary values for the next 
C  iterations
C
        IF (ILE .EQ. 0) THEN
          P(ILES) = SCL
        ELSE
          P(ILES) = ZERO
        ENDIF
 
        IF (IRI .EQ. L) THEN
          P(IRIS + 2) = SCL
        ELSE
          P(IRIS + 2) = ZERO
        ENDIF
C
  100 CONTINUE
C
      DL = P(NX + 1) + P(NX)
      IF (DL .EQ. ZERO) RETURN
C
C  The P-value
C
      Q= ONE - EXP (FL(NX) + FL(NY) + LOG(DL) - NOFDIV * SMALLN - FL(N))
C
C  Q=1 is allowable, Q<=0 is not (accuracy loss due to 
C  rounding errors)
C
      IF (Q .LE. ZERO) IFAULT = 3
      END
