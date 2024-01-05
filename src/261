      REAL FUNCTION SQMCQ(CDF, IP, N, RHO2, IFAULT)
C
C        ALGORITHM AS 261  APPL. STATIST. (1991) VOL. 40, NO. 1
C
C        Returns the quantile of the distribution of the square
C        of the sample multiple correlation coefficient for
C        given values of CDF, IP, N, and RHO2
C
C        A modification of the secant method, called the
C        Illinois method, is used
C
C        The auxiliary algorithm SQMCOR, used to compute the
C        C.D.F. of the square of the sample multiple correlation
C        coefficient, is required
C
      INTEGER IP, N, IFAULT
      REAL CDF, RHO2
      REAL DIFF, EPS, F0, F1, F2, X0, X1, X2, ZERO, ONE, TWO
      REAL SQMCOR
      EXTERNAL SQMCOR
      DATA ZERO, ONE, TWO / 0.0, 1.0, 2.0 /
      DATA EPS / 1.0E-6 /
C
      SQMCQ = CDF
      IFAULT = 2
C
C        Perform domain check
C
      IF (RHO2 .LT. ZERO .OR. RHO2 .GT. ONE .OR. IP .LT. 2 .OR.
     *    N .LE. IP) RETURN
      IFAULT = 3
      IF (CDF .LT. ZERO .OR. CDF .GT. ONE) RETURN
      IFAULT = 0
      IF (CDF .EQ. ZERO .OR. CDF .EQ. ONE) RETURN
C
C        Use ONE and ZERO as two starting points for the
C        Illinois method
C
      X0 = ONE
      F0 = ONE - CDF
      X1 = ZERO
      F1 = -CDF
C
C        Continue iterations until convergence is achieved
C
   10 DIFF = F1 * (X1 - X0) / (F1 - F0)
      X2 = X1 - DIFF
      F2 = SQMCOR(X2, IP, N, RHO2, IFAULT) - CDF
      IF (IFAULT .NE. 0) RETURN
C
C        Check for convergence
C
      IF (ABS(F2) .LE. EPS * CDF) GO TO 20
      IF (ABS(DIFF) .LE. EPS * X2) GO TO 20
      IF (F2 * F1 .GE. ZERO) THEN
         F0 = F0 / TWO
      ELSE
         X0 = X1
         F0 = F1
      END IF
      X1 = X2
      F1 = F2
      GO TO 10
C
   20 SQMCQ = X2
      RETURN
      END
===========================
The author of algorithms AS 260 & 261 has found an improved algorithm
which has been published in the journal `Computational Statistics &
Data Analysis'.   Here it is.
C------------------------------------------------------------
C
C The following code for calculating quantiles of R^2 (qr2.for)
C is by the author of AS 261 and uses an improved algorithm.

      REAL FUNCTION QR2(M, SIZE, RHO2, P, IFAULT)
C
C  Computes the quantile of the distribution of the square of the
C  sample multiple correlation coefficient for given number of
C  random variables M, sample size SIZE, square of the population
C  multiple correlation coefficient RHO2, and lower tail area P
C
C  Reference:
C  Ding, C.G. (1996) `On the computation of the distribution of
C  the square of the sample multiple correlation coefficient',
C  Comput. Statist. & Data Analysis, vol. 22, 345-350.
C
C  IFAULT is a fault indicator:
C  = 1 if there is no convergence after 10 Newton's iterations
C  = 2 if any of the input values is illegal
C  = 0 otherwise
C
C  No auxiliary algorithm is required
C
      INTEGER SIZE
      REAL N
      DATA ZERO, HALF, ONE, TWO / 0.0, 0.5, 1.0, 2.0 /
      DATA EPS, DELTA / 1.0E-6, 1.0E-4 /
      DATA ITRMAX / 10 /
      DATA RP / 1.772453850905516028 /

      QR2 = P
C
C        Test for admissibility of arguments
C
      IFAULT = 2
      IF (M .LE. 1 .OR. SIZE .LE. M .OR. RHO2 .LT. ZERO .OR.
     $    RHO2 .GT. ONE .OR. P .LT. ZERO .OR. P .GT. ONE) RETURN
      IFAULT = 0
      IF (P .EQ. ZERO .OR. P .EQ. ONE) RETURN
C
C        Calculate the constants needed for each Newton's iteration
C
      A = (M - 1) / TWO
      B = (SIZE - M) / TWO
      AB = (SIZE - 1) / TWO
      IF (MOD(M + 1, 2) .EQ. 0) THEN
         NA = A + HALF
         GA = ONE
         DO 10 I = 1, NA
   10    GA = GA * I
      ELSE
         NA = A + ONE
         GA = RP
         DO 20 I = 1, NA
   20    GA = GA * (I - HALF)
      ENDIF
      IF (MOD(SIZE - M, 2) .EQ. 0) THEN
         NB = B - HALF
         GB = ONE
         IF (NB .EQ. 0) GO TO 50
         DO 30 I = 1, NB
   30    GB = GB * I
      ELSE
         NB = B
         GB = RP
         IF (NB .EQ. 0) GO TO 50
         DO 40 I = 1, NB
   40    GB = GB * (I - HALF)
      ENDIF
   50 IF (MOD(SIZE - 1, 2) .EQ. 0) THEN
         NAB = AB - HALF
         GAB = ONE
         IF (NAB .EQ. 0) GO TO 75
         DO 60 I = 1, NAB
   60    GAB = GAB * I
      ELSE
         NAB = AB
         GAB = RP
         DO 70 I = 1, NAB
   70    GAB = GAB * (I - HALF)
      ENDIF
   75 Q0 = (ONE - RHO2) ** AB
      COEFF = GAB / GA / GB
C
C        Use 0.5 as a starting value for Newton's iterations
C
      Y = HALF
C
C        Perform Newton's iterations
C
      DO 120 ITER = 1, ITRMAX
C
C        Evaluate the first terms of the series for CDF (distribution
C        function) and PDF (density)
C
      N = ONE
      YP = ONE - Y
      T = COEFF * Y ** A * YP ** B
      S = A * T / Y / YP
      Q = Q0
      V = Q
      CDF = V * T
      PDF = Q * S
C
C        Check if a + n > (a + b + n)y
C
   80 IF (A + N .GT. (A + B + N) * Y) GO TO 90
C
C        Evaluate the next terms of two series and then the
C        partial sums
C
      Q = Q * (A + B + N - ONE) * RHO2 / N
      V = V + Q
      S = T * (A + B + N - ONE) / YP
      T = T * Y * (A + B + N - ONE) / (A + N)
      CDF = CDF + V * T
      PDF = PDF + Q * S
      N = N + ONE
      GO TO 80
C
C        Find the error bounds and check for convergence for both
C        series
C
   90 BNDCDF = T * Y * (A + B + N - ONE) / (A + N - (A + B + N) * Y)
      BNDPDF = T * (A + B + N - ONE) * (ONE - V) / YP
  100 IF (BNDCDF .LE. EPS .AND. BNDPDF .LE. EPS) GO TO 110
C
C        Continue to update the terms and then accumulate
C
      Q = Q * (A + B + N - ONE) * RHO2 / N
      V = V + Q
      IF (BNDCDF .LE. EPS) THEN
         S = S * Y * (A + B + N - ONE) / (A + N - ONE)
         PDF = PDF + Q * S
         N = N + ONE
         BNDPDF = S * Y * (A + B + N - ONE) * (ONE - V) / (A + N - ONE)
         GO TO 100
      ELSE IF (BNDPDF .LE. EPS) THEN
              T = T * Y * (A + B + N - ONE) / (A + N)
              CDF = CDF + V * T
              N = N + ONE
              BNDCDF = T * Y * (A+B+N-ONE) / (A+N - (A+B+N)*Y)
              GO TO 100
           ELSE
              S = T * (A + B + N - ONE) / YP
              T = T * Y * (A + B + N - ONE) / (A + N)
              CDF = CDF + V * T
              PDF = PDF + Q * S
              N = N + ONE
              GO TO 90
      ENDIF
C
C        Obtain a new Y and make changes if it is illegal
C
  110 DIFF = (CDF - P) / PDF
      YNEW = Y - DIFF
      IF (YNEW .LE. ZERO) THEN
         Y = Y / TWO
      ELSE IF (YNEW .GE. ONE) THEN
              Y = (Y + ONE) / TWO
           ELSE
              Y = YNEW
      ENDIF
C
C        Check for convergence of Newton's iterations
C
      IF (ABS(DIFF) .LE. DELTA * Y) THEN
         QR2 = Y
         RETURN
      ENDIF
  120 CONTINUE
      IFAULT = 1
      RETURN
      END
C
      PROGRAM MAIN
C
C        This is a driver program that calls QR2 and produces output
C
      INTEGER SIZE

   10 WRITE (*,11)
   11 FORMAT (/1X,'ENTER M (>1), N (>M), RHO2 (BETWEEN 0 AND 1), and',
     $        /1X,'P (BETWEEN 0 AND 1) FOR QR2 ==> ')
      READ (*,*) M, SIZE, RHO2, P
      Y = QR2(M, SIZE, RHO2, P, IFAULT)
      ICODE = IFAULT + 1
      GO TO (20, 30, 40), ICODE
   20 WRITE (*,21) M, SIZE, RHO2, P, Y
   21 FORMAT (/1X, 'QR2(', I3, ',', I3, ',', E10.4, ',', E10.4, ') =',
     $        E11.5)
      GO TO 60
   30 WRITE (*,31)
   31 FORMAT (/1X,43HNO CONVERGENCE AFTER 10 NEWTON'S ITERATIONS)
      GO TO 60
   40 WRITE (*,41)
   41 FORMAT (/1X,'THE INPUT VALUE IS ILLEGAL!')
   60 WRITE (*,61)
   61 FORMAT (/1X, 'ENTER 1 TO CONTINUE OR 0 TO QUIT ==> ')
      READ (*,*) K
      IF (K .EQ. 1) GO TO 10
      STOP
      END
