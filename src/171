      SUBROUTINE TPSHV(N, T, NTMAX, SY2, PROB, KGP, KPR, A, KF, DLF,
     *     IFAULT)
C
C     ALGORITHM AS 171  APPL. STATIST. (1982) VOL.31, NO.3
C
C     Fisher's exact variance test for the Poisson distribution
C
      INTEGER T, SY2, TS, A(T), KF(T)
      REAL DD, DC, DM, DMIN, DP, D1, DZ, DPRB, DLF(NTMAX)
      DATA MAXP /1000000/, MAXT /80/
C
      IFAULT = 0
      DPRB = DZ
      DM = D1 - DMIN
      IF (N .LT. 2 .OR. N .GT. NTMAX) IFAULT = 2
      IF (T .LT. 2 .OR. T .GT. NTMAX .OR. T .GT. MAXT) IFAULT = 3
      IF (IFAULT .NE. 0) GO TO 200
      KGP = 0
      KPR = 0
C
C     Generate log factorials and constant
C
      DLF(1) = DZ
      DO 10 I = 2, NTMAX
   10 DLF(I) = DLF(I-1) + LOG(FLOAT(I))
      DC = - T * LOG(FLOAT(N)) + DLF(N) + DLF(T)
C
C     Generate partitions of sample total (T) starting with M the
C     number of non-zero values in the sample.
C
      M = MIN(T, N)
      DO 15 J = 1, M
   15 KF(J) = 0
   20 DO 25 J = 1, M
   25 A(J) = 1
      TS = M - 1
   30 A(M) = T
      DO 35 J = 2, M
   35 A(M) = A(M) - A(J-1)
      TS = TS + A(M)**2
      KGP = KGP + 1
C
C     The current M-part partition of T:
C         A(1) + A(2) + ... + A(M) = T
C     with sums of squares  TS = SUM A(J)**2
C     is admissible if TS < SY2
C
      IF (TS .GE. SY2) GO TO 80
C
C     Convert sample into frequency distribution.
C     KZ is the number of zero values in the sample.
C     KF(J) is the number of A()'s = J.
C
      KZ = N - M
      MX = A(M)
      DO 60 J = 1, M
	IAK = A(J)
	KF(IAK) = KF(IAK) + 1
   60 CONTINUE
C
C     Compute DP = probability for this sample
C
      KPR = KPR + 1
      DD = DZ
      IF (KZ .GT. 0) DD = DLF(KZ)
      IAK = KF(1)
      IF (KF(1) .GT. 0) DD = DD + DLF(IAK)
      KF(1) = 0
      IF (MX .LT. 2) GO TO 70
      DO 65 J = 2, MX
	IF (KF(J) .EQ. 0) GO TO 65
	IAK = KF(J)
	DD = DD + DLF(IAK) + KF(J) * DLF(J)
	KF(J) = 0
   65 CONTINUE
   70 DP = EXP(DC - DD)
      DPRB = DPRB + DP
      IF (DPRB .LT. DM) GO TO 80
      IFAULT = 1
      DPRB = DM
      GO TO 200
   80 CONTINUE
      IT = M - 1
   90 IF (IT .EQ. 0) GO TO 120
C
C     Determine next M-partition & update sum of squares.
C
      IF ((A(M) - A(IT)) .GT. 1) GO TO 100
      IT = IT - 1
      GO TO 90
  100 IAT = A(IT) + 1
      M1 = M - 1
      DO 110 J = IT, M1
	TS = TS - A(J)**2 + IAT**2
	A(J) = IAT
  110 CONTINUE
      TS = TS - A(M)**2
C
C     New values for A(1), ..., A(M-1) are determined.
C
      IF (KGP .LT. MAXP) GO TO 30
      IFAULT = 1
      GO TO 200
  120 M = M - 1
C
C     End evaluation of M-part partitions.
C     Decrease M by 1 and continue if SY2 > sum of squares for
C     the last partition.
C
      IF (TS .GE. SY2) GO TO 200
      IF (M .GT. 1) GO TO 20
C
  200 PROB = D1 - DPRB
      RETURN
      END
