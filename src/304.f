      SUBROUTINE FISHER (X, M, Y, N, TOTAL, POSSIB, P, IFAULT)
C
C        ALGORITHM AS 304.1 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Fisher's non-parametric randomization test for two small
C        independent random samples
C
      INTEGER M, N, TOTAL, POSSIB, IFAULT
      REAL X(*), Y(*), P
C
      INTEGER MAXSAM, MAXSIZ
      PARAMETER (MAXSAM = 14, MAXSIZ = 3432)
C
C        Important : set MAXSIZ >= COMB(MAXSAM, MAXSAM / 2)
C
      INTEGER K, SIZE1, SIZE2, WS3(MAXSAM), COUNT
      REAL SUMX, SUMY, WS1(MAXSIZ), WS2(MAXSIZ), SUMM
C
      INTEGER COMB, TRADES
      REAL MEAN, SUM
      EXTERNAL COMB, SUM, TRADES
C
      EXTERNAL CMPLMT, EXCHNG, KTRADE
C
      MEAN(SUMM, COUNT) = SUMM / REAL(COUNT)
C
      IF (COMB(MAXSAM, MAXSAM/2) .GT. MAXSIZ) THEN
         IFAULT = 1
      ELSE IF (M .GT. MAXSAM .OR. N .GT. MAXSAM) THEN
              IFAULT = 2
           ELSE
              IFAULT = 0
              SUMX = SUM(X, M)
              SUMY = SUM(Y, N)
              IF (MEAN(SUMX, M) .GT. MEAN(SUMY, N)) THEN
                 CALL EXCHNG(X, M, Y, N, SUMX, SUMY)
           END IF
           TOTAL = 0
           IF (M .EQ. N) THEN
              DO 10 K = 1, (M-1)/2
                 CALL KTRADE(X, M, WS1, SIZE1, WS3, K)
                 CALL KTRADE(Y, N, WS2, SIZE2, WS3, K)
                 TOTAL = TOTAL + TRADES(WS1, SIZE1, WS2, SIZE2)
                 CALL CMPLMT(WS1, SIZE1, SUMX)
                 CALL CMPLMT(WS2, SIZE2, SUMY)
                 TOTAL = TOTAL + TRADES(WS1, SIZE1, WS2, SIZE2)
   10         CONTINUE
              IF (MOD(M, 2) .EQ. 0) THEN
                 CALL KTRADE(X, M, WS1, SIZE1, WS3, K)
                 CALL KTRADE(Y, N, WS2, SIZE2, WS3, K)
                 TOTAL = TOTAL + TRADES(WS1, SIZE1, WS2, SIZE2)
              END IF
         ELSE
            DO 20 K = 1, MIN(M, N)
               CALL KTRADE(X, M, WS1, SIZE1, WS3, K)
               CALL KTRADE(Y, N, WS2, SIZE2, WS3, K)
               TOTAL = TOTAL + TRADES(WS1, SIZE1, WS2, SIZE2)
   20       CONTINUE
         END IF
         POSSIB = COMB(M+N, M)
         P = REAL(TOTAL + 1) / REAL(POSSIB)
      END IF
C
      RETURN
      END
C
      SUBROUTINE EXCHNG (X, M, Y, N, SX, SY)
C
C        ALGORITHM AS 304.2 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Exchanges the sample data.  Assumes both X and Y have been
C        previously dimensioned to at least max(M, N) elements
C
      INTEGER M, N
      REAL X(*), Y(*), SX, SY
C
      INTEGER C, K
      REAL TEMP
C
      TEMP = SX
      SX = SY
      SY = TEMP
C
      C = MIN(M, N)
      DO 10 K = 1, C
         TEMP = X(K)
         X(K) = Y(K)
         Y(K) = TEMP
   10 CONTINUE
      IF (M .GT. N) THEN
         DO 20 K = C+1, M
            Y(K) = X(K)
   20    CONTINUE
         N = M
         M = C
      ELSE IF (M .LT. N) THEN
         DO 30 K = C+1, N
            X(K) = Y(K)
   30    CONTINUE
         M = N
         N = C
      END IF
C
      RETURN
      END
C
      SUBROUTINE KTRADE (W, K, WPRIME, KPRIME, WS, R)
C
C        ALGORITHM AS 304.3 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Generates and sorts the sums of the R-combinations of the
C        elements of W
C
      INTEGER K, KPRIME, WS(*), R
      REAL W(*), WPRIME(*)
C
      INTEGER COMB
      REAL SUM
      EXTERNAL COMB, SUM
C
      EXTERNAL CMPLMT, GENER, SORT
C
      KPRIME = COMB(K, R)
      IF (R .LE. K - R .OR. R .EQ. K) THEN
         CALL GENER(W, K, WPRIME, KPRIME, WS, R)
         CALL SORT(WPRIME, KPRIME)
      ELSE
         CALL GENER(W, K, WPRIME, KPRIME, WS, K - R)
         CALL SORT(WPRIME, KPRIME)
         CALL CMPLMT(WPRIME, KPRIME, SUM(W, K))
      ENDIF
C
      RETURN
      END
C
      INTEGER FUNCTION TRADES (XPRIME, MPRIME, YPRIME, NPRIME)
C
C        ALGORITHM AS 304.4 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Returns the number of 1-for-1 trades that refutes the null
C        hypothesis.  Assumes that XPRIME has the smaller mean and
C        that both arrays are sorted in ascending order.
C
      INTEGER MPRIME, NPRIME
      REAL XPRIME(*), YPRIME(*)
C
      INTEGER I, J
C
      TRADES = 0
      I = 1
      J = 1
   10 IF (J .GT. NPRIME) GOTO 40
   20 IF (XPRIME(I) .GE. YPRIME(J)) GOTO 30
      I = I + 1
      IF (I .LE. MPRIME) GOTO 20
   30 TRADES = TRADES + (MPRIME - I + 1)
      J = J + 1
      IF (I .LE. MPRIME) GOTO 10
C
   40 RETURN
      END
C
      SUBROUTINE CMPLMT (WPRIME, KPRIME, SUM)
C
C        ALGORITHM AS 304.5 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Reverse and complement the data in WPRIME
C
      INTEGER KPRIME
      REAL WPRIME(*), SUM
C
      INTEGER I, J
      REAL TEMP
C
      J = KPRIME
      DO 10 I = 1, KPRIME / 2 + MOD(KPRIME, 2)
         TEMP = WPRIME(I)
         WPRIME(I) = REAL(DBLE(SUM) - DBLE(WPRIME(J)))
         WPRIME(J) = REAL(DBLE(SUM) - DBLE(TEMP))
         J = J - 1
   10 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE GENER (W, N, WPRIME, NPRIME, INDEX, R)
C
C        ALGORITHM AS 304.6 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Computes an array of sums of the various R-combinations of
C        the elements of W
C
      INTEGER N, NPRIME, R, INDEX(R)
      REAL W(N), WPRIME(NPRIME)
C
      INTEGER I, J
      DOUBLE PRECISION SUM
      LOGICAL INIT
C
      EXTERNAL NEXT
C
      INIT = .TRUE.
C
      DO 20 I = 1, NPRIME
         CALL NEXT(INDEX, R, N, INIT)
         SUM = 0.0D0
         DO 10 J = 1, R
            SUM = SUM + DBLE(W(INDEX(J)))
   10    CONTINUE
         WPRIME(I) = REAL(SUM)
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE NEXT (RCOMBO, R, N, INIT)
C
C        ALGORITHM AS 304.7 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Accepts some R-combination of the first N integers and then
C        computes the next R-combination in the lexicographic
C        ordering of the N! / (R! * (N - R)!) such R-combinations.  
C        Returns the first R-combination if the initialization 
C        indicator is .true. and then resets the indicator.
C
      INTEGER R, N, RCOMBO(R)
      LOGICAL INIT
C
      INTEGER I, J, D
C
      IF (INIT) THEN
         DO 10 I = 1, R
            RCOMBO(I) = I
   10    CONTINUE
         INIT = .FALSE.
      ELSE
         D = N - R
         J = R
C
C        The counter J is not prevented from going out of bounds
C        which will happen if there is no next R-combination
C
   20    IF (RCOMBO(J) .LT. D + J) GOTO 30
         J = J - 1
         GOTO 20
   30    RCOMBO(J) = RCOMBO(J) + 1
         DO 40 I = J + 1, R
            RCOMBO(I) = RCOMBO(I - 1) + 1
   40    CONTINUE
      END IF
C
      RETURN
      END
C
C        General purpose subroutines
C
      SUBROUTINE SORT (X, N)
C
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Sorts the N values stored in array X in ascending order
C
      INTEGER N
      REAL X(N)
C
      INTEGER I, J, INCR
      REAL TEMP
C
      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. N) GOTO 10

C
C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
   40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50 X(J) = TEMP
      I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20
C
      RETURN
      END
C
      INTEGER FUNCTION COMB (N, K)
C
C        ALGORITHM AS 304.9 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Returns the number of combinations of N things taken K at a
C        time or 0 if the parameters are incompatible
C
      INTEGER N, K
C
      INTEGER M, I
      REAL NUMER, DENOM, FACT
      EXTERNAL FACT
C
      IF (K .LT. 0 .OR. K .GT. N) THEN
         COMB = 0
      ELSE
         M = N - K
         NUMER = 1.0
         DO 10 I = N, 1 + MAX(K, M), -1
            NUMER = NUMER * REAL(I)
   10    CONTINUE
         DENOM = FACT(MIN(K, M))
         COMB = NINT(NUMER / DENOM)
      END IF
C
      RETURN
      END
C
      REAL FUNCTION FACT (N)
C
C        ALGORITHM AS 304.10 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Returns the factorial of N or 0.0 if the parameter is
C        negative
C
      INTEGER N
C
      INTEGER I
C
      IF (N .LT. 0) THEN
         FACT = 0.0
      ELSE
        FACT = 1.0
        DO 10 I = 2, N
           FACT = FACT * REAL(I)
   10   CONTINUE
      END IF
C
      RETURN
      END
C
      REAL FUNCTION SUM (X, N)
C
C        ALGORITHM AS 304.11 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Returns the sum of the values stored in X
C
      INTEGER N
      REAL X(N)
C
      INTEGER I
      DOUBLE PRECISION ACCUM
C
      ACCUM = 0.0D0
      DO 10 I = 1, N
         ACCUM = ACCUM + DBLE(X(I))
   10 CONTINUE
      SUM = REAL(ACCUM)
C
      RETURN
      END
