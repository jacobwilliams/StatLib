      REAL FUNCTION SQMCOR(X, IP, N, RHO2, IFAULT)
C
C        ALGORITHM AS 260  APPL. STATIST. (1991) VOL. 40, NO. 1
C
C        Computes the C.D.F. for the distribution of the
C        square of the multiple correlation coefficient
C        with parameters IP, N, and RHO2
C
C        The following auxiliary algorithms are required:
C        ALOGAM - Log-gamma function (CACM 291)
C         (or  ALNGAM    AS 245)
C        BETAIN - Incomplete beta function (AS 63)
C
      REAL X, RHO2
      INTEGER IP, N, IFAULT
      REAL A, B, BETA, ERRBD, ERRMAX, GX, Q, SUMQ, TEMP, TERM, XJ,
     *     ZERO, HALF, ONE
      INTEGER ITRMAX
      REAL ALOGAM, BETAIN
      EXTERNAL ALOGAM, BETAIN
C
      DATA ERRMAX, ITRMAX / 1.0E-6, 100 /
      DATA ZERO, HALF, ONE / 0.0, 0.5, 1.0 /
C
      SQMCOR = X
      IFAULT = 2
      IF (RHO2 .LT. ZERO .OR. RHO2 .GT. ONE .OR. IP .LT. 2 .OR.
     *    N .LE. IP) RETURN
      IFAULT = 3
      IF (X .LT. ZERO .OR. X .GT. ONE) RETURN
      IFAULT = 0
      IF (X .EQ. ZERO .OR. X .EQ. ONE) RETURN
C
      A = HALF * (IP - 1)
      B = HALF * (N - IP)
C
C        Initialize the series
C
      BETA = EXP(ALOGAM(A, IFAULT) + ALOGAM(B, IFAULT) -
     *       ALOGAM(A + B, IFAULT))
      TEMP = BETAIN(X, A, B, BETA, IFAULT)
C
C        There is no need to test IFAULT since all of the
C        parameter values have already been checked
C
      GX = EXP(A * LOG(X) + B * LOG(ONE - X) - LOG(A)) / BETA
      Q = (ONE - RHO2) ** (A + B)
      XJ = ZERO
      TERM = Q * TEMP
      SUMQ = ONE - Q
      SQMCOR = TERM
C
C        Perform recurrence until convergence is achieved
C
   10 XJ = XJ + ONE
      TEMP = TEMP - GX
      GX = GX * (A + B + XJ - ONE) * X / (A + XJ)
      Q = Q * (A + B + XJ - ONE) * RHO2 / XJ
      SUMQ = SUMQ - Q
      TERM = TEMP * Q
      SQMCOR = SQMCOR + TERM
C
C        Check for convergence and act accordingly
C
      ERRBD = (TEMP - GX) * SUMQ
      IF ((INT(XJ) .LT. ITRMAX) .AND. (ERRBD .GT. ERRMAX)) GO TO 10
      IF (ERRBD .GT. ERRMAX) IFAULT = 1
      RETURN
      END
=============================
The author of algorithms AS 260 & 261 has found an improved algorithm
which has been published in the journal `Computational Statistics &
Data Analysis'.   Here it is.
C------------------------------------------------------------
C
C The following code for the distribution function of R^2 (cdfr2.for)
C and for its probability density (pdfr2.for) is by the author of AS 260
C and uses an improved algorithm.

      REAL FUNCTION CDFR2(Y, M, SIZE, RHO2, IFAULT)
C
C  Computes the distribution function of the square of the sample
C  multiple correlation coefficient for given abscissa Y, number
C  of random variables M, sample size SIZE, and square of the population
C  multiple correlation coefficient RHO2
C
C  Reference:
C  Ding, C.G. (1996) `On the computation of the distribution of
C  the square of the sample multiple correlation coefficient',
C  Comput. Statist. & Data Analysis, vol. 22, 345-350.
C
C  IFAULT is a fault indicator:
C  = 1 if any of the input values is illegal
C  = 0 otherwise
C
C  No auxiliary algorithm is required
C
      INTEGER SIZE
      REAL N
      DATA EPS / 1.0E-6 /
      DATA HALF, ZERO, ONE, TWO / 0.5, 0.0, 1.0, 2.0 /
      DATA RP / 1.772453850905516028 /

      CDFR2 = ZERO
      IFAULT = 1
      IF (M .LE. 1 .OR. SIZE .LE. M .OR. RHO2 .LT. ZERO .OR.
     $    RHO2 .GT. ONE) RETURN
      IFAULT = 0
      IF (Y .LE. ZERO) RETURN
      CDFR2 = ONE
      IF (Y .GE. ONE) RETURN
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
         IF (NAB .EQ. 0) GO TO 80
         DO 60 I = 1, NAB
   60    GAB = GAB * I
      ELSE
         NAB = AB
         GAB = RP
         DO 70 I = 1, NAB
   70    GAB = GAB * (I - HALF)
      ENDIF

C
C        Evaluate the first term
C
   80 N = ONE
      Q = (ONE - RHO2) ** AB
      V = Q
      T = Y ** A * (ONE - Y) ** B * GAB / GA / GB
      TERM = V * T
      CDFR2 = TERM
C
C        Check if a + n > (a + b + n)y
C
   90 IF (A + N .GT. (A + B + N) * Y) GO TO 100
C
C        Evaluate the next term of the expansion and then the
C        partial sum
C
      Q = Q * (A + B + N - ONE) * RHO2 / N
      V = V + Q
      T = T * Y * (A + B + N - ONE) / (A + N)
      TERM = V * T
      CDFR2 = CDFR2 + TERM
      N = N + ONE
      GO TO 90
C
C        Find the error bound and check for convergence
C
  100 BOUND = T * Y * (A + B + N - ONE) / ((A + N) - (A + B + N) * Y)
      IF (BOUND .LE. EPS) RETURN
C
C        Evaluate the next term of the expansion and then the
C        partial sum
C
      Q = Q * (A + B + N - ONE) * RHO2 / N
      V = V + Q
      T = T * Y * (A + B + N - ONE) / (A + N)
      TERM = V * T
      CDFR2 = CDFR2 + TERM
      N = N + ONE
      GO TO 100
      END
C
      PROGRAM MAIN
C
C        is a driver program that calls CDFR2 and produces output
C
      INTEGER SIZE
   10 WRITE (*,11)
   11 FORMAT (/1X,'ENTER Y, M (>1), N (>M), RHO2 (BETWEEN 0 AND 1)',
     $        /1X,'FOR CDFR2 ==>')
      READ (*,*) Y, M, SIZE, RHO2
      CDF = CDFR2(Y, M, SIZE, RHO2, IFAULT)
      IF (IFAULT .EQ. 0) THEN
         WRITE (*,21) Y, M, SIZE, RHO2, CDF
   21    FORMAT (/1X, 'CDFR2(', E10.4, ',', 2(I3, ','), E10.4, ') =',
     $           E12.6)
      ELSE
         WRITE (*,31)
   31    FORMAT (/1X,'THE INPUT VALUE IS ILLEGAL!')
      ENDIF
      WRITE(*,41)
   41 FORMAT (/1X, 'ENTER 1 TO CONTINUE OR 0 TO QUIT ==> ')
      READ (*,*)K
      IF (K .EQ. 1) GOTO 10
      STOP
      END


C------------------------------------------------------------
C


      REAL FUNCTION PDFR2(Y, M, SIZE, RHO2, IFAULT)
C
C  Computes the density (pdf) of the square of the sample multiple
C  correlation coefficient for given abscissa Y, number of random
C  variables M, sample size SIZE, and square of the population
C  multiple correlation coefficient RHO2
C
C  Reference:
C  Ding, C.G. (1996) `On the computation of the distribution of
C  the square of the sample multiple correlation coefficient',
C  Comput. Statist. & Data Analysis, vol. 22, 345-350.
C
C  IFAULT is a fault indicator:
C  = 1 if any of the input values is illegal
C  = 0 otherwise
C
C  No auxiliary algorithm is required
C
      INTEGER SIZE
      REAL N
      DATA EPS / 1.0E-6 /
      DATA HALF, ZERO, ONE, TWO / 0.5, 0.0, 1.0, 2.0 /
      DATA RP / 1.772453850905516028 /

      PDFR2 = ZERO
      IFAULT = 1
      IF (M .LE. 1 .OR. SIZE .LE. M .OR. RHO2 .LT. ZERO .OR.
     $    RHO2 .GT. ONE) RETURN
      IFAULT = 0
      IF (Y .LE. ZERO .OR. Y .GE. ONE) RETURN
      A = (M - 1) / TWO
      B = (SIZE - M) / TWO
      AB = (SIZE - 1) / TWO
      IF (MOD(M - 1, 2) .EQ. 0) THEN
         NA = A - HALF
         GA = ONE
         IF (NA .EQ. 0) GO TO 25
         DO 10 I = 1, NA
   10    GA = GA * I
      ELSE
         NA = A
         GA = RP
         IF (NA .EQ. 0) GO TO 25
         DO 20 I = 1, NA
   20    GA = GA * (I - HALF)
      ENDIF
   25 IF (MOD(SIZE - M, 2) .EQ. 0) THEN
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
         IF (NAB .EQ. 0) GO TO 80
         DO 60 I = 1, NAB
   60    GAB = GAB * I
      ELSE
         NAB = AB
         GAB = RP
         DO 70 I = 1, NAB
   70    GAB = GAB * (I - HALF)
      ENDIF

C
C        Evaluate the first term
C
   80 N = ONE
      Q = (ONE - RHO2) ** AB
      V = Q
      S = Y ** (A - ONE) * (ONE - Y) ** (B - ONE) * GAB / GA / GB
      TERM = Q * S
      PDFR2 = TERM
C
C        Check if a + n > (a + b + n)y
C
   90 IF (A + N .GT. (A + B + N) * Y) GO TO 100
C
C        Evaluate the next term of the expansion and then the
C        partial sum
C
      Q = Q * (A + B + N - ONE) * RHO2 / N
      V = V + Q
      S = S * Y * (A + B + N - ONE) / (A + N - ONE)
      TERM = Q * S
      PDFR2 = PDFR2 + TERM
      N = N + ONE
      GO TO 90
C
C        Find the error bound and check for convergence
C
  100 BOUND = S * Y * (A + B + N - ONE) * (ONE - V) / (A + N - ONE)
      IF (BOUND .LE. EPS) RETURN
C
C        Evaluate the next term of the expansion and then the
C        partial sum
C
      Q = Q * (A + B + N - ONE) * RHO2 / N
      V = V + Q
      S = S * Y * (A + B + N - ONE) / (A + N - ONE)
      TERM = Q * S
      PDFR2 = PDFR2 + TERM
      N = N + ONE
      GO TO 100
      END
C
      PROGRAM MAIN
C
C        is a driver program that calls PDFR2 and produces output
C
      INTEGER SIZE

   10 WRITE (*,11)
   11 FORMAT (/1X,'ENTER Y, M (>1), N (>M), RHO2 (BETWEEN 0 AND 1)',
     $        /1X,'FOR PDFR2 ==> ')
      READ (*,*) Y, M, SIZE, RHO2
      PDF = PDFR2(Y, M, SIZE, RHO2, IFAULT)
      IF (IFAULT .EQ. 0) THEN
         WRITE (*,21) Y, M, SIZE, RHO2, PDF
   21    FORMAT (/1X, 'PDFR2(', E10.4, ',', 2(I3, ','), E10.4, ') =',
     $           E12.6)
      ELSE
         WRITE (*,31)
   31    FORMAT (/1X, 'THE INPUT VALUE IS ILLEGAL!')
      ENDIF
      WRITE(*,41)
   41 FORMAT (/1X, 'ENTER 1 TO CONTINUE OR 0 TO QUIT ==> ')
      READ (*,*)K
      IF (K .EQ. 1) GOTO 10
      STOP
      END
