      SUBROUTINE WL(NR, N, NN, G, X, S, D, SIGMA, SIGINV, WLT, NRUNIV,
     *              CHIOMB, NRSTOC, METHOD, KERROR, IFAULT)
C
C        ALGORITHM AS 262  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Wei-Lachin multivariate generalization of the generalized
C        Wilcoxon test of Gehan and the log-rank test for incomplete
C        multivariate observations with censoring.
C
      INTEGER NREP, NOBS
      PARAMETER (NREP=10, NOBS=500)
C
C        NREP     the maximum number of repeat times
C        NOBS     the maximum number of observations
C
      INTEGER NR, N(2), NN, G(NOBS), S(NOBS, NREP), METHOD, KERROR,
     *        IFAULT
      REAL X(NOBS, NREP), D(2, NREP), SIGMA(NN), SIGINV(NN), WLT(NREP),
     *     NRUNIV(NREP), CHIOMB, NRSTOC
      INTEGER Y(2, NOBS), I, IFAIL, IPOINT, NULLTY, J, J2, K, K1,
     *        K2, NT
      REAL EE(2, NREP), MU(2, NOBS, NREP), PSI(2, NOBS, NREP), QE(2),  
     *     SIG(2, NREP, NREP), WORK(NREP), Q, SIGSUM, SUM, TSUM, Z,
     *     ZERO, ONE, TWO
C
      DATA ZERO, ONE, TWO / 0.0, 1.0, 2.0 /
C
      ZSQRT(I) = SQRT(Z)
      ZREAL(Z) = REAL(I)
C
C        Initialization
C
      IFAULT = 4
      IF (NR .GT. NREP) RETURN
      IFAULT = 5
      NT = N(1) + N(2)
      IF (NT .GT. NOBS) RETURN
      IFAULT = 0
      DO 30 K = 1, NR
         DO 20 I = 1, 2
            EE(I, K) = ZERO
            D(I, K) = ZERO
            DO 10 K2 = 1, NR
               SIG(I, K, K2) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C        Y    number at risk of failure
C        EE   partial sum used in equation 1
C        QE   Q*E (see equation 4)
C        Evaluate MU (equation 4)
C
      DO 110 K = 1, NR
         DO 70 J = 1, NT
            DO 40 I = 1, 2
               Y(I, J) = 0
               MU(I, J, K) = ZERO
   40       CONTINUE
            IF (S(J, K) .EQ. 1) THEN
               DO 50 J2 = 1, NT
                  IF (X(J2, K) .GE. X(J, K))
     *                Y(G(J2), J) = Y(G(J2), J) + 1
   50          CONTINUE
               SUM = Y(1, J) + Y(2, J)
               IF (SUM .EQ. ZERO) THEN
                  IFAULT = 1
                  RETURN
               END IF
               IF (METHOD .EQ. 1) Q = SUM / NT
               IF (METHOD .EQ. 2) Q = ONE
               DO 60 I = 1, 2
                  QE(I) = Q * Y(I, J) / SUM
   60          CONTINUE
               MU(2, J, K) = QE(1)
               MU(1, J, K) = QE(2)
               EE(G(J), K) = EE(G(J), K) + MU(G(J), J, K)
               D(G(J), K) = D(G(J), K) + ONE
            ELSE
               MU(1, J, K) = ZERO
               MU(2, J, K) = ZERO
            END IF
   70    CONTINUE
C
C        Evaluate PSI (equation 5)
C
         DO 100 J = 1, NT
            DO 90 I = 1, 2
               PSI(I, J, K) = ZERO
               DO 80 J2 = 1, NT
                  IF (G(J2) .EQ. I .AND. X(J2, K) .LE. X(J, K)) THEN
                     IF (S(J2, K) .EQ. 1) THEN
                        PSI(I, J, K) = PSI(I, J, K) +
     *                                 MU(I, J2, K) / Y(I, J2)
                     END IF
                  END IF
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
C
C        Compute entries in covariance matrix (see equations 2 and 3)
C
      DO 140 K1 = 1, NR
         DO 130 K2 = 1, K1
            DO 120 J = 1, NT
               I = G(J)
               IF (N(I) .EQ. 0) THEN
                  IFAULT = 2
                  RETURN
               END IF
               SIG(I, K1, K2) = SIG(I, K1, K2) +
     *                          (MU(I, J, K1) * S(J, K1) -
     *                          PSI(I, J, K1)) * (MU(I, J, K2) *
     *                          S(J, K2) - PSI(I, J, K2)) / N(I)
  120       CONTINUE
            IPOINT = K1 * (K1 - 1) / 2 + K2
            SIGMA(IPOINT) = (N(1) * SIG(1, K1, K2) +
     *                      N(2) * SIG(2, K1, K2)) / NT
  130    CONTINUE
  140 CONTINUE
C
C        Compute Wei-Lachin univariate test NRUNIV  (see eqns 1 and 6)
C
      DO 150 K = 1, NR
         WLT(K) = (EE(1, K) - EE(2, K)) / ZSQRT(ZREAL(NT))
         IPOINT = K * (K + 1) / 2
         IF (SIGMA(IPOINT) .EQ. ZERO) THEN
            IFAULT = 3
            KERROR = K
            RETURN
         END IF
         NRUNIV(K) = WLT(K) / ZSQRT(SIGMA(IPOINT))
  150 CONTINUE
C
C        Compute inverse of covariance matrix and Wei-Lachin
C        multivariate statistics CHIOMB (for omnibus test) and
C        NRSTOC (for test of stochastic ordering)  (see eqn 7)
C
      CALL SYMINV(SIGMA, NR, NN, SIGINV, WORK, NULLTY, IFAIL)
      CHIOMB = ZERO
      TSUM = ZERO
      SIGSUM = ZERO
      DO 170 J = 1, NR
         TSUM = TSUM + WLT(J)
         IPOINT = J * (J + 1) / 2
         CHIOMB = CHIOMB + WLT(J) * WLT(J) * SIGINV(IPOINT)
         SIGSUM = SIGSUM + SIGMA(IPOINT)
         DO 160 J2 = 1, J - 1
            IPOINT = J * (J - 1) / 2 + J2
            CHIOMB = CHIOMB + TWO * WLT(J) * WLT(J2) * SIGINV(IPOINT)
            SIGSUM = SIGSUM + TWO * SIGMA(IPOINT)
  160    CONTINUE
  170 CONTINUE
      NRSTOC = TSUM / ZSQRT(SIGSUM)
      RETURN
      END
