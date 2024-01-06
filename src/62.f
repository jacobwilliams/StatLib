c AS 62 generates the frequencies for the Mann-Whitney U-statistic.
c Users are much more likely to need the distribution function.
c Code to return the distribution function has been added at the end
c of AS 62 by Alan Miller.   Remove the C's in column 1 to activate it.
c
      SUBROUTINE UDIST(M, N, FRQNCY, LFR, WORK, LWRK, IFAULT)
C
C     ALGORITHM AS 62  APPL. STATIST. (1973) VOL.22, NO.2
C
C     The distribution of the Mann-Whitney U-statistic is generated for
C     the two given sample sizes
C
      INTEGER M, N, LFR, LWRK, IFAULT
      REAL FRQNCY(LFR), WORK(LWRK)
C
C     Local variables
C
      INTEGER MINMN, MN1, MAXMN, N1, I, IN, L, K, J
      REAL ZERO, ONE, SUM
      DATA ZERO /0.0/, ONE /1.0/
C
C     Check smaller sample size
C
      IFAULT = 1
      MINMN = MIN(M, N)
      IF (MINMN .LT. 1) RETURN
C
C     Check size of results array
C
      IFAULT = 2
      MN1 = M * N + 1
      IF (LFR .LT. MN1) RETURN
C
C     Set up results for 1st cycle and return if MINMN = 1
C
      MAXMN = MAX(M, N)
      N1 = MAXMN + 1
      DO 1 I = 1, N1
    1 FRQNCY(I) = ONE
      IF (MINMN .EQ. 1) GO TO 4
C
C     Check length of work array
C
      IFAULT = 3
      IF (LWRK .LT. (MN1 + 1) / 2 + MINMN) RETURN
C
C     Clear rest of FREQNCY
C
      N1 = N1 + 1
      DO 2 I = N1, MN1
    2 FRQNCY(I) = ZERO
C
C     Generate successively higher order distributions
C
      WORK(1) = ZERO
      IN = MAXMN
      DO 3 I = 2, MINMN
        WORK(I) = ZERO
        IN = IN + MAXMN
        N1 = IN + 2
        L = 1 + IN / 2
        K = I
C
C     Generate complete distribution from outside inwards
C
        DO 3 J = 1, L
          K = K + 1
          N1 = N1 - 1
          SUM = FRQNCY(J) + WORK(J)
          FRQNCY(J) = SUM
          WORK(K) = SUM - FRQNCY(N1)
          FRQNCY(N1) = SUM
    3 CONTINUE
C
    4 IFAULT = 0
C
C     Code to overwrite the frequency function with the distribution
C     function.   N.B. The frequency in FRQNCY(1) is for U = 0, and
C     that in FRQNCY(I) is for U = I - 1.
C
C     SUM = ZERO
C     DO 10 I = 1, MN1
C       SUM = SUM + FRQNCY(I)
C       FRQNCY(I) = SUM
C  10 CONTINUE
C     DO 20 I = 1, MN1
C  20 FRQNCY(I) = FRQNCY(I) / SUM
C
      RETURN
      END

