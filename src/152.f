      REAL FUNCTION CHYPER(POINT, KK, LL, MM, NN, IFAULT)
C
C     ALGORITHM AS R77  APPL. STATIST. (1989), VOL.38, NO.1
C     Replaces AS 59 and AS 152
C     Incorporates AS R86 from vol.40(2)
C
C     Auxiliary routines required: ALNFAC (AS 245), ALNORM (AS 66)
C
      INTEGER          KK, LL, MM, NN, IFAULT, K, L, M, N, I, J, NL, KL,
     *                 MNKL, MVBIG, MBIG
      REAL             ZERO, ONE, P, PT, HALF, ALNFAC, ELIMIT, MEAN, 
     *                 SIG, ALNORM, SXTEEN, SCALE, ROOTPI, ARG, HUNDRD
      LOGICAL          POINT, DIR
      PARAMETER (ZERO = 0.0, HALF = 0.5, ONE = 1.0, MVBIG = 1000,
     *           MBIG = 600, ELIMIT = -88.0, SXTEEN = 16.0,
     *           SCALE = 1.0E35, ROOTPI = 2.50662 82746 31001,
     *           HUNDRD = 100.0)
C
      K = KK + 1
      L = LL + 1
      M = MM + 1
      N = NN + 1
      DIR = .TRUE.
C
C     Check arguments are within permitted limits
C
      IFAULT = 1
      CHYPER = ZERO
      IF (N .LT. 1 .OR. M .LT. N .OR. K .LT. 1 .OR. K .GT. M) RETURN
C
      IFAULT = 2
      IF (L .LT. 1 .OR. K-L .GT. M-N) RETURN
      IF (.NOT. POINT) CHYPER = ONE
      IF (L .GT. N .OR. L .GT. K) RETURN
      IFAULT = 0
      CHYPER = ONE
      IF (K .EQ. 1 .OR. K .EQ. M .OR. N .EQ. 1 .OR. N .EQ. M) RETURN
      IF (.NOT. POINT .AND. LL .EQ. MIN(KK, NN)) RETURN
C
      P = REAL(NN) / REAL(MM - NN)
      IF (REAL(MIN(KK, MM-KK)) .GT. SXTEEN * MAX(P, ONE/P) .AND.
     *   MM .GT. MVBIG .AND. ELIMIT .GT. -HUNDRD) THEN
C
C     Use a normal approximation
C
      MEAN = REAL(KK) * REAL(NN) / REAL(MM)
      SIG = SQRT(MEAN * (REAL(MM-NN) / REAL(MM)) * (REAL(MM-KK) /
     *             (REAL(MM-1))))
      IF (POINT) THEN
        ARG = -HALF * (((REAL(LL) - MEAN) / SIG)**2)
        CHYPER = ZERO
        IF (ARG .GE. ELIMIT) CHYPER = EXP(ARG) / (SIG * ROOTPI)
      ELSE
        CHYPER = ALNORM((REAL(LL) + HALF - MEAN) / SIG, .FALSE.)
      END IF
C
      ELSE
C
C     Calculate exact hypergeometric probabilities.
C     Interchange K and N if this saves calculations.
C
      IF (MIN(K-1, M-K) .GT. MIN(N-1, M-N)) THEN
        I = K
        K = N
        N = I
      END IF
      IF (M-K .LT. K-1) THEN
        DIR = .NOT. DIR
        L = N - L + 1
        K = M - K + 1
      END IF
C
      IF (MM .GT. MBIG) THEN
C
C     Take logarithms of factorials.
C
        P = ALNFAC(NN) - ALNFAC(MM) + ALNFAC(MM-KK) + ALNFAC(KK) +
     *        ALNFAC(MM-NN) - ALNFAC(LL) - ALNFAC(NN-LL) - ALNFAC(KK-LL)
     *      - ALNFAC(MM-NN-KK+LL)
        CHYPER = ZERO
        IF (P .GE. ELIMIT) CHYPER = EXP(P)
        ELSE
C
C     Use Freeman/Lund algorithm
C
        DO 3 I = 1, L-1
          CHYPER = CHYPER * REAL(K-I) * REAL(N-I) / (REAL(L-I) *
     *                      REAL(M-I))
    3     CONTINUE
        IF (L .NE. K) THEN
          J = M - N + L
          DO 5 I = L, K-1
            CHYPER = CHYPER * REAL(J-I) / REAL(M-I)
    5          CONTINUE
          END IF
C
        END IF
C
        IF (POINT) RETURN
        IF (CHYPER .EQ. ZERO) THEN
C
C     We must recompute the point probability since it has underflowed.
C
        IF (MM .LE. MBIG) P = ALNFAC(NN) - ALNFAC(MM) + ALNFAC(KK) +
     *      ALNFAC(MM-NN) - ALNFAC(LL) - ALNFAC(NN-LL) - ALNFAC(KK-LL) -
     *      ALNFAC(MM-NN-KK+LL) + ALNFAC(MM-KK)
        P = P + LOG(SCALE)
        IF (P .LT. ELIMIT) THEN
          IFAULT = 3
          IF (LL .GT. REAL(NN*KK + NN + KK +1)/(MM +2)) CHYPER = ONE
          RETURN
        ELSE
          P = EXP(P)
        END IF
        ELSE
C
C     Scale up at this point.
C
        P = CHYPER * SCALE
        END IF
C
        PT = ZERO
        NL = N - L
        KL = K - L
        MNKL = M - N - KL + 1
        IF (L .LE. KL) THEN
        DO 7 I = 1, L-1
          P = P * REAL(L-I) * REAL(MNKL-I) /
     *                         (REAL(NL+I) * REAL(KL+I))
          PT = PT + P
    7     CONTINUE
        IF (P .EQ. ZERO) IFAULT = 3
        ELSE
        DIR = .NOT. DIR
        DO 9 J = 0, KL-1
          P = P * REAL(NL-J) * REAL(KL-J) / (REAL(L+J) * REAL(MNKL+J))
            PT = PT + P
    9     CONTINUE
        IF (P .EQ. ZERO) IFAULT = 3
        END IF
C
        IF (DIR) THEN
        CHYPER = CHYPER + (PT / SCALE)
        ELSE
        CHYPER = ONE - (PT / SCALE)
        END IF
C
      END IF
C
      END

