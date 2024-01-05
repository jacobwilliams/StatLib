      SUBROUTINE OPART(N, K, NI, MCOUNT, IFAULT)
C
C        ALGORITHM AS 303.1 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Generation of ordered multinomial frequencies.  Given integers
C        N (N >= 0) and K (K >= 1), generates all possible integer arrays
C        NI of size K, whose elements are non-negative, non-decreasing
C        and sum to N.
C
      INTEGER N, K, NI(K), MCOUNT, IFAULT
C
      INTEGER I, I0, I1, J, K1, KSUM, NSUM
C
      EXTERNAL JOB
C
C        Perform domain checks.  Initialize MCOUNT, NI, K1 and KSUM.
C        Set NI(I0) = 0; for I0 = 1 to K-1.
C
      IFAULT = 1
      MCOUNT = 0
      IF(K .LT. 1) RETURN
      IF(N .LT. 0) RETURN
      IFAULT = 0
      K1 = K - 1
      IF(K .EQ. 1) K1 = 1
      KSUM = 0
C
      DO 10 I = 1, K1
         NI(I) = 0
   10 CONTINUE
C
C        Begin with pivot element in position K-1 of array NI
C
      I0 = K1
   20 IF(I0 .EQ. K1) GO TO 40
C
C        Set NI(J); for J = I0+1 to K-1, equal to pivot cell
C        frequency.  Update KSUM.
C
      I1 = I0 + 1
C
      DO 30 J = I1, K1
         NI(J) = NI(I0)
         KSUM = KSUM + NI(I0)
   30 CONTINUE
C
C        Assign NI(K)
C
   40 MCOUNT = MCOUNT + 1
      NSUM = NI(K1)
      NI(K) = N - KSUM
C
      CALL JOB(K, NI, MCOUNT)
C
C        Reset pivot cell to K-1 and update NSUM.  Check upper limit
C        for NI(I0).
C
      I0 = K1
      NSUM = NSUM + NI(K)
   50 IF( NI(I0) .LT. NSUM / (K - I0 + 1) ) GO TO 60
C
C        Move pivot I0 to next lower cell.  Update KSUM and NSUM.
C
      I0 = I0 - 1
      IF(I0 .LE. 0) RETURN
      KSUM = KSUM - NI(I0+1)
      NSUM = NSUM + NI(I0)
      GO TO 50
C
C        Increment NI(I0) by 1 and update KSUM.
C
   60 NI(I0) = NI(I0) + 1
      KSUM = KSUM + 1
      GO TO 20
      END
C
      SUBROUTINE JOB(K, NI, MCOUNT)
C
C        ALGORITHM AS 303.2 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Evaluates C(X) for X = NI and writes out ordered K-tuple
C
      INTEGER K, NI(K), MCOUNT
C
      REAL SIZE, TIES
      INTEGER J, K1, KTEMP, NTEMP, RK
C
      K1 = K-1
      NTEMP = NI(1)
      SIZE = K
      KTEMP = K
      TIES = 1.0D0
C
      DO 30 J = 1, K1
         KTEMP = KTEMP - 1
         RK = KTEMP
         IF(NI(J + 1) .EQ. NTEMP) GO TO 10
         NTEMP = NI(J + 1)
         TIES = 1.0D0
         GO TO 20
   10    TIES = TIES + 1.0D0
   20    SIZE = SIZE * RK / TIES
   30 CONTINUE
C
      WRITE(*, 40) MCOUNT, SIZE, (NI(J), J = 1, K)
   40 FORMAT(1X,I5,1X,F15.0,7X,10I4/
     *       1X,21X,10I4)
      RETURN
      END
