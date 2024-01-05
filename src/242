      SUBROUTINE MARMA(K, N, P, Q, W, IK, PHI, THETA, MEAN, MEANS,
     *                 QQ, XTOL, RLOGL, V, WORK, LWORK, IWORK, LIWORK,
     *                 IFAULT)
C
C     ALGORITHM AS 242.1  APPL. STATIST. (1989), VOL.38, NO. 1
C
C     Auxiliary routines required: DECOMP & SOLVE from CACM algorithm 423
C
      INTEGER          P, Q, R, K, N, LWORK, IFAULT, KR, KMAT, KC, LW1,
     *                 LW2, LW3, LW4, LW5, LW6, LW7, LW8, LW9, LW10,
     *                 LW11, LW12, LW13, IK, I, J, LIWORK, IWORK(LIWORK)
      LOGICAL          MEAN
      REAL             MEANS(K), PHI(IK,P*K+1), QQ(IK,K), REAL,
     *                 V(IK,N), W(IK,N), WORK(LWORK), THETA(IK,Q*K+1),
     *                 TWO, LOG2PI, CONST, RLOGL, XTOL, TSIG
      PARAMETER        (TWO = 2.0, LOG2PI = 1.8378771)
      INTRINSIC        MAX, REAL
      EXTERNAL         VARMA, CHOL
C
C     test for errors in the input data
C
      IFAULT = 0
      IF(K .LT. 1) IFAULT = 1
      IF(IFAULT .EQ. 1) RETURN
      IF(N .LT. 1) IFAULT = 2
      IF(IFAULT .EQ. 2) RETURN
      IF(P .LT. 0) IFAULT = 3
      IF(IFAULT .EQ. 3) RETURN
      IF(Q .LT. 0) IFAULT = 4
      IF(IFAULT .EQ. 4) RETURN
      IF(P .EQ. 0 .AND. Q .EQ. 0) IFAULT = 5
      IF(IFAULT .EQ. 5) RETURN
      IF(IK .LT. K) IFAULT = 6
      IF(IFAULT .EQ. 6) RETURN
C
C     check that workspace arrays are big enough
C
      R = MAX(P,Q)
      KR = K * R
      KMAT = K * MAX(R, K * (P + 1))
      CONST = (REAL(N * K) / TWO) * LOG2PI
C
      KC = K * K
      LW1 = 1
      LW2 = K * KR + LW1
      LW3 = K * KR + LW2
      LW4 = K * KR + LW3
      LW5 = KC + LW4
      LW6 = KC + LW5
      LW7 = KC + LW6
      LW8 = KC + LW7
      LW9 = KC + LW8
      LW10 = KC + LW9
      LW11 = KC * (P + 1) + LW10
      LW12 = KC * (R + 1) + LW11
      LW13 = KMAT * KMAT + LW12
C
      I = MAX(IK * K, LW13 + KC - 1)
      IF(LWORK .LT. I) IFAULT = 7
      IF(IFAULT .EQ. 7) RETURN
      IF(LIWORK .LT. KMAT) IFAULT = 8
      IF(IFAULT .EQ. 8) RETURN
C
C     set upper triangle of QQ to lower triangle of QQ
C
      DO 30 I = 2, K
         DO 31 J = 1, I - 1
            QQ(J,I) = QQ(I,J)
   31    CONTINUE
   30 CONTINUE
C
C     check whether QQ is positive-definite and evaluate its determinant
C     (storing it in TSIG)
C
      CALL CHOL(QQ, IK, K, WORK, IFAULT)
      IF(IFAULT .GT. 0) THEN
         IFAULT = 9
         RETURN
      ELSE
         TSIG = WORK(1)
         DO 20 I = 2, K
            TSIG = TSIG * WORK((I - 1) * IK + I)
   20    CONTINUE
         TSIG = TSIG * TSIG
      ENDIF
C
C     now call VARMA to calculate the log likelihood function
C
      CALL VARMA(K, N, P, Q, R, W, IK, PHI, THETA, MEAN, MEANS, QQ,
     *           TSIG, CONST, V, XTOL, RLOGL, KR, KMAT, WORK(LW1),
     *           WORK(LW2), WORK(LW3), WORK(LW4), WORK(LW5),
     *           WORK(LW6), WORK(LW7), WORK(LW8), WORK(LW9),
     *           WORK(LW10), WORK(LW11), WORK(LW12), WORK(LW13),
     *           IWORK, IFAULT)
C
      END
      SUBROUTINE VARMA(K, N, P, Q, R, W, IK, PHI, THETA, MEAN, MEANS,
     *                 QQ, TSIG, CONST, V, XTOL, RLOGL, KR, KMAT,
     *                 GAMMA, TEMP, TEMPK, TEMPL, F, INVF, TEMPM, A,
     *                 MT, B, Z, MAT, WA, IWORK, IFAULT)
C
C     ALGORITHM AS 242.2  APPL. STATIST. (1989), VOL.38, NO. 1
C
      INTEGER          IFAULT, K, KMAT, KR, N, P, Q, R, IK, I, J,
     *                 J7, K1, K2, L, L8, M, T, IWORK(KMAT)
      LOGICAL          MEAN, ANOTT, DELTA, PEARL
      REAL             CONST, RLOGL, TSIG, XTOL, ZERO, ONE, TWO,
     *                 A(K,K), B(K*K*(P+1)), F(K,K), GAMMA(K,KR),
     *                 INVF(K,K), MAT(KMAT,KMAT), MEANS(K), MT(K,K),
     *                 PHI(IK,P*K+1), QQ(IK,K), TEMP(K,KR),
     *                 TEMPL(K,K), TEMPM(K,K), THETA(IK,Q*K+1),
     *                 W(IK,N), Z(K*K*(R+1)), WA(K,K), TEMPK(K,KR),
     *                 DETP, SIG, SM, SSQ, SUM, REAL, V(IK,N)
      PARAMETER        (ZERO = 0.0, ONE = 1.0, TWO = 2.0)
      INTRINSIC        LOG, REAL, ABS, MAX
      EXTERNAL         COVARS, CHOL, BKSB
C
C     This subroutine computes the log likelihood function of a vector
C     ARMA model
C
      ANOTT = .FALSE.
      IF(P .GT. Q) ANOTT = .TRUE.
C
C     call COVARS to calculate the autocovariances of W(t) and the
C     cross-covariances between W(t) and E(t)
C
      CALL COVARS(K, P, Q, R, IK, PHI, THETA, QQ, GAMMA, Z, KMAT,
     *            MAT, B, IWORK, IFAULT)
C
      IF(IFAULT .GT. 0) IFAULT = 10
      IF(IFAULT .EQ. 10) RETURN
C
C     calculate the first k columns of P(1/0) and store as TEMPK
C
      DO 140 K1 = 1, R
         DO 120 I = 1, K
            DO 100 J = 1, K
               SUM = ZERO
               DO 40 K2 = 1, K
                  DO 20 M = K1, P
                     SUM = SUM + PHI(I, (M - 1) * K + K2) *
     *                     Z((M - K1 + 1) * K * K + (K2 - 1) * K + J)
   20             CONTINUE
   40          CONTINUE
               DO 80 M = K1, Q
                  DO 60 K2 = 1, K
                     SUM = SUM - THETA(I, (M - 1) * K + K2) *
     *                     GAMMA(J, (M - K1) * K + K2)
   60             CONTINUE
   80          CONTINUE
C
               TEMPK(I, (K1 - 1) * K + J) = SUM
  100       CONTINUE
  120    CONTINUE
  140 CONTINUE
C
C     initialize A(1/0), V(1), SSQ, F(1) and DETP
C
      DO 160 I = 1, KR
         Z(I) = ZERO
  160 CONTINUE
C
      DO 200 I = 1, K
         DO 180 J = 1, K
            F(I,J) = TEMPK(I,J) + QQ(I,J)
  180    CONTINUE
         V(I,1) = W(I,1)
         IF(MEAN) V(I,1) = V(I,1) - MEANS(I)
         B(I) = V(I,1)
  200 CONTINUE
C
C     factorise F(1) as LL' and store L as A
C
      CALL CHOL(F, K, K, A, IFAULT)
      IF(IFAULT .NE. 0) THEN
         IFAULT = 11
         RETURN
      ENDIF
C
C     calculate determinant (DETP) of F(1)
C
      DETP = A(1,1)
      DO 220 I = 2, K
         DETP = DETP * A(I,I)
  220 CONTINUE
      DETP = DETP * DETP
C
C     set MT = F(1) inverse
C
      DO 260 I = 1, K
         DO 240 J = 1, K
            MT(I,J) = ZERO
  240    CONTINUE
         MT(I,I) = ONE
  260 CONTINUE
C
C     note that BKSB cannot fail
C
      CALL BKSB(A, K, K, K, .FALSE., MT, IFAULT)
C
C     set upper triangle of INVF to A'
C
      DO 300 J = 1, K
         DO 280 I = 1, J
            INVF(I,J) = A(J,I)
  280    CONTINUE
  300 CONTINUE
      CALL BKSB(INVF, K, K, K, .TRUE., MT, IFAULT)
C
      CALL BKSB(A, K, K, 1, .FALSE., B, IFAULT)
      SSQ = ZERO
      DO 320 I = 1, K
         SSQ = SSQ + B(I) * B(I)
  320 CONTINUE
C
      DETP = LOG(DETP)
C
C     calculate L(1) and K(1)
C
      DO 360 I = 1, K
         DO 340 J = 1, KR
            MAT(I,J) = ZERO
            GAMMA(I,J) = ZERO
  340    CONTINUE
  360 CONTINUE
C
      DO 480 L = 1, R
         DO 460 I = 1, K
            DO 440 J = 1, K
               SUM = ZERO
               IF(L .LE. P) THEN
                  DO 380 K2 = 1, K
                     SUM = SUM + PHI(I, (L - 1) * K + K2) *
     *                     TEMPK(K2,J)
  380             CONTINUE
               ENDIF
               IF(L .LT. R) SUM = SUM + TEMPK(I, L * K + J)
               DO 420 K2 = 1, K
                  SM = ZERO
                  IF(L .LE. P) SM = PHI(I, (L - 1) * K + K2)
                  IF(L .LE. Q) SM = SM - THETA(I, (L - 1) * K + K2)
                  SUM = SUM + SM * QQ(K2,J)
  420          CONTINUE
               MAT(I, (L - 1) * K + J) = SUM
C
  440       CONTINUE
  460    CONTINUE
  480 CONTINUE
C
      IF(ANOTT) THEN
C
         DO 580 L = 1, R
            DO 560 I = 1, K
               DO 540 J = 1, K
                  SUM = ZERO
                  IF(L .LE. Q) THEN
                     DO 500 K2 = 1, K
                        SUM = SUM + THETA(I, (L - 1) * K + K2) *
     *                        TEMPK(K2,J)
  500                CONTINUE
                  ENDIF
                  IF(L .LT. R) SUM = SUM + TEMPK(I, L * K + J)
                  GAMMA(I, (L - 1) * K + J) = SUM
  540          CONTINUE
  560       CONTINUE
  580    CONTINUE
C
      ELSE
C
         DO 640 L = 1, R
            DO 620 I = 1, K
               DO 600 J = 1, K
                  GAMMA(I, (L - 1) * K + J) = MAT(I, (L - 1) * K + J)
  600          CONTINUE
  620       CONTINUE
  640    CONTINUE
C
      ENDIF
C
C     start the recursions
C
      DELTA = .FALSE.
      PEARL = .FALSE.
      J7 = N
C
      DO 2240 T = 2, N
C
         IF(ANOTT  .AND.  T .GT. P - Q) PEARL = .TRUE.
C
C        calculate TEMPK
C
         IF(.NOT.DELTA) THEN
C
C           find L' (inverse) and store as INVF
C
            DO 680 J = 1, K
               DO 660 I = 1, K
                  TEMPK(I,J) = A(J,I)
                  INVF(I,J) = ZERO
  660          CONTINUE
               INVF(J,J) = ONE
  680       CONTINUE
            CALL BKSB(TEMPK, K, K, K, .TRUE., INVF, IFAULT)
C
            DO 760 L = 1, R
               DO 740 I = 1, K
                  DO 720 J = 1, K
                     SUM = ZERO
                     DO 700 K2 = 1, K
                        SUM = SUM + GAMMA(I, (L - 1) * K + K2) *
     *                        INVF(K2,J)
  700                CONTINUE
                     TEMPK(I, (L - 1) * K + J) = SUM
  720             CONTINUE
  740          CONTINUE
  760       CONTINUE
C
C           calculate TEMP
C
         ENDIF
         DO 800 L = 1, KR
            TEMP(1,L) = ZERO
  800    CONTINUE
C
         IF(ANOTT) GOTO 900
C
C        TEMP = T * A(t-1/t-2)
C
         IF(T .EQ. 2) GOTO 1080
         DO 880 L = 1, R
            DO 860 I = 1, K
               SUM = ZERO
               IF(L .LE. P) THEN
                  DO 820 K2 = 1, K
                     SUM = SUM + PHI(I, (L - 1) * K + K2) * Z(K2)
  820             CONTINUE
               ENDIF
               IF(L .LT. R) SUM = SUM + Z(L * K + I)
               TEMP(1, (L - 1) * K + I) = SUM
  860       CONTINUE
  880    CONTINUE
         GOTO 1080
  900    IF(T .GT. 2) THEN
C
C           TEMP = A * A(t-1/t-2) + R * W(t-1)
C
            DO 980 L = 1, R
               DO 960 I = 1, K
                  SUM = ZERO
                  IF(L .LE. Q) THEN
                     DO 920 K2 = 1, K
                        SUM = SUM + THETA(I, (L - 1) * K + K2) * Z(K2)
  920                CONTINUE
                  ENDIF
                  IF(L .LT. R) SUM = SUM + Z(L * K + I)
                  TEMP(1, (L - 1) * K + I) = SUM
  960          CONTINUE
  980       CONTINUE
C
         ENDIF
         DO 1060 L = 1, R
            DO 1040 I = 1, K
               SUM = ZERO
               DO 1020 K2 = 1, K
                  SM = ZERO
                  IF(L .LE. P) SM = PHI(I, (L - 1) * K + K2)
                  IF(L .LE. Q) SM = SM - THETA(I, (L - 1) * K + K2)
                  IF(MEAN) THEN
                     SUM = SUM + SM * (W(K2,T - 1) - MEANS(K2))
                  ELSE
                     SUM = SUM + SM * W(K2,T - 1)
                  ENDIF
 1020          CONTINUE
               TEMP(1, (L - 1) * K + I) = TEMP(1, (L -1 ) * K + I)
     *                                    + SUM
 1040       CONTINUE
 1060    CONTINUE
C
C        A(t/t-1) = TEMP + TEMPK * V(t-1)
C
 1080    DO 1160 L = 1, R
            DO 1140 L8 = 1, K
               SUM = TEMP(1, (L - 1) * K + L8)
               IF((.NOT.DELTA) .OR. (.NOT.ANOTT)) THEN
                  DO 1100 K2 = 1, K
                     SUM = SUM + TEMPK(L8, (L - 1) * K + K2) * B(K2)
 1100             CONTINUE
               ENDIF
               Z((L - 1) * K + L8) = SUM
 1140       CONTINUE
 1160    CONTINUE
C
         IF(DELTA) GOTO 2180
C
C        calculate TEMPL
C
         DO 1200 I = 1, K
            DO 1180 J = 1, K
               TEMPL(I,J) = MAT(I,J)
 1180       CONTINUE
 1200    CONTINUE
C
C        recalculate TEMP
C
         IF(.NOT.ANOTT) THEN
C
C           TEMP = T * L(t-1)
C
            DO 1300 L = 1, R
               DO 1280 I = 1, K
                  DO 1260 J = 1, K
                     SUM = ZERO
                     IF(L .LE. P) THEN
                        DO 1220 K2 = 1, K
                           SUM = SUM + PHI(I, (L - 1) * K + K2) *
     *                           MAT(K2,J)
 1220                   CONTINUE
                     ENDIF
                     IF(L .LT. R) SUM = SUM + MAT(I, L * K + J)
                     TEMP(I, (L - 1) * K + J) = SUM
C
 1260             CONTINUE
 1280          CONTINUE
 1300       CONTINUE
C
         ELSE
C
C           TEMP = A * L(t-1)
C
            DO 1420 L = 1, R
               DO 1400 I = 1, K
                  DO 1380 J = 1, K
                     SUM = ZERO
                     IF(L .LE. Q) THEN
                        DO 1340 K2 = 1, K
                           SUM = SUM + THETA(I, (L - 1) * K + K2) *
     *                           MAT(K2,J)
 1340                   CONTINUE
                     ENDIF
                     IF(L .LT. R) SUM = SUM + MAT(I, L * K + J)
                     TEMP(I, (L - 1) * K + J) = SUM
C
 1380             CONTINUE
 1400          CONTINUE
 1420       CONTINUE
         ENDIF
C
C        calculate TEMPM = M(t-1) * (h'L(t-1))'
C
C        calculate choleski factorization of MT and store as INVF
C
         CALL CHOL(MT, K, K, INVF, IFAULT)
C
C        compute h'L(t-1) * INVF
C
         DO 1500 I = 1, K
            DO 1480 J = 1, K
               SUM = ZERO
               DO 1460 K2 = J, K
                  SUM = SUM + TEMPL(I,K2) * INVF(K2,J)
 1460          CONTINUE
               WA(I,J) = SUM
 1480       CONTINUE
 1500    CONTINUE
C
C        compute INVF * WA'
C
         DO 1560 I = 1, K
            DO 1540 J = 1, K
               SUM = ZERO
               DO 1520 K2 = 1, I
                  SUM = SUM + INVF(I,K2) * WA(J,K2)
 1520          CONTINUE
               TEMPM(I,J) = SUM
 1540       CONTINUE
 1560    CONTINUE
C
C        update K(t)
C
         DO 1640 L = 1, R
            DO 1620 I = 1, K
               DO 1600 J = 1, K
                  SUM = ZERO
                  IF(PEARL .AND. L .GE. Q + 1) GAMMA(I, (L - 1) * K + J)
     *                                         = SUM
                  IF((.NOT.PEARL) .OR. L .LT. Q + 1) THEN
                     DO 1580 K2 = 1, K
                        SUM = SUM - TEMP(I, (L - 1) * K + K2) *
     *                        TEMPM(K2,J)
 1580                CONTINUE
                     GAMMA(I, (L - 1) * K + J) = GAMMA(I, (L - 1) * K
     *                                           + J) + SUM
                  ENDIF
 1600          CONTINUE
 1620       CONTINUE
 1640    CONTINUE
C
C        update F(t)
C
         DO 1700 I = 1, K
            DO 1680 J = I, K
               SUM = F(I,J)
               DO 1660 K2 = 1, K
                  SUM = SUM - WA(I,K2) * WA(J,K2)
 1660          CONTINUE
               F(I,J) = SUM
               F(J,I) = SUM
 1680       CONTINUE
 1700    CONTINUE
C
C        test for convergence of F(t)'s
C
         SUM = ZERO
         DO 1720 I = 1, K
            IF(QQ(I,I) .GT. ZERO) THEN
               SUM = MAX(SUM, ABS(F(I,I) - QQ(I,I)) / QQ(I,I))
            ELSE
               SUM = MAX(SUM, ABS(F(I,I) - QQ(I,I)))
            ENDIF
 1720    CONTINUE
         IF(SUM .LT. XTOL) DELTA = .TRUE.
C
         IF(.NOT. DELTA) GOTO 1920
C
C        set TEMPK = R * L'
C
         J7 = T
         DO 1780 I = 1, K
            DO 1760 J = 1, K
               DO 1740 L = 1, R
                  SUM = ZERO
                  IF(L .LE. P) SUM = PHI(I, (L - 1) * K + J)
                  IF(L .LE. Q) SUM = SUM - THETA(I, (L - 1) * K + J)
                  TEMPK(I, (L - 1) * K + J) = SUM
 1740          CONTINUE
 1760       CONTINUE
 1780    CONTINUE
C
         DO 1900 L = 1, R
            DO 1840 I = 1, K
               DO 1820 J = 1, K
                  SUM = ZERO
                  DO 1800 K2 = J, K
                     SUM = SUM + TEMPK(I, (L - 1) * K + K2) * A(K2,J)
 1800             CONTINUE
                  WA(I,J) = SUM
 1820          CONTINUE
 1840       CONTINUE
            DO 1880 I = 1, K
               DO 1860 J = 1, K
                  TEMPK(I, (L - 1) * K + J) = WA(I,J)
 1860          CONTINUE
 1880       CONTINUE
 1900    CONTINUE
C
C        calculate choleski decomposition of F(t) and det(F(t))
C
C        copy old A onto WA
C
 1920    DO 1960 I = 1, K
            DO 1940 J = 1, K
               WA(I,J) = A(I,J)
 1940       CONTINUE
 1960    CONTINUE
         CALL CHOL(F, K, K, A, IFAULT)
         IF(IFAULT .NE. 0) THEN
            IFAULT = 11
            RETURN
         ENDIF
         SIG = A(1,1)
         DO 1980 I = 2, K
            SIG = SIG * A(I,I)
 1980    CONTINUE
         SIG = SIG * SIG
C
         IF(DELTA) GOTO 2180
C
C        update M(t)
C
C        transpose TEMPM
C
         DO 2020 I = 1, K
            DO 2000 J = 1, I
               SUM = TEMPM(I,J)
               TEMPM(I,J) = TEMPM(J,I)
               TEMPM(J,I) = SUM
 2000       CONTINUE
 2020    CONTINUE
         CALL BKSB(A, K, K, K, .FALSE., TEMPM, IFAULT)
         DO 2080 I = 1, K
            DO 2060 J = I, K
               SUM = MT(I,J)
               DO 2040 K2 = 1, K
                  SUM = SUM + TEMPM(K2,I) * TEMPM(K2,J)
 2040          CONTINUE
               MT(I,J) = SUM
               MT(J,I) = SUM
 2060       CONTINUE
 2080    CONTINUE
C
C        update L(t)
C
         CALL BKSB(WA, K, K, K, .FALSE., TEMPL, IFAULT)
         DO 2160 L = 1, R
            DO 2140 I = 1, K
               DO 2120 J = 1, K
                  SUM = ZERO
                  IF(PEARL .AND. L .GE. Q + 1) MAT(I, (L - 1) * K + J)
     *                                         = SUM
                  IF((.NOT.PEARL) .OR. L .LT. Q + 1) THEN
                     DO 2100 K2 = 1, K
                        SUM = SUM + TEMPK(I, (L - 1) * K + K2)
     *                        * TEMPL(K2,J)
 2100                CONTINUE
                     MAT(I, (L - 1) * K + J) = TEMP(I, (L - 1) * K + J)
     *                                         - SUM
                  ENDIF
 2120          CONTINUE
 2140       CONTINUE
 2160    CONTINUE
C
 2180    DO 2200 I = 1, K
            V(I,T) = W(I,T) - Z(I)
            IF(MEAN) V(I,T) = V(I,T) - MEANS(I)
            B(I) = V(I,T)
 2200    CONTINUE
C
         CALL BKSB(A, K, K, 1, .FALSE., B, IFAULT)
         DO 2220 I = 1, K
            SSQ = SSQ + B(I) * B(I)
 2220    CONTINUE
         IF((.NOT.DELTA) .OR. T .LE. J7) DETP = DETP + LOG(SIG)
C
 2240 CONTINUE
C
      RLOGL = - CONST - (SSQ + DETP + REAL(N - J7) * LOG(TSIG)) / TWO
C
      END
C
      SUBROUTINE COVARS(K, P, Q, R, IK, PHI, THETA, QQ, GAMMA,
     *                  Z, KMAT, MAT, B, IWORK, IFAULT)
C
C     ALGORITHM AS242.3  APPL. STATIST. (1989), VOL.38, NO. 1
C
      INTEGER           P, Q, R, K, KMAT, IFAULT, M, I, J, I2, K2, J2,
     *                  KW, L, L4, IK, IWORK(KMAT)
      REAL              QQ(IK,K), PHI(IK,P*K+1), THETA(IK,Q*K+1),
     *                  GAMMA(K,K*R), Z(K*K*(R+1)), MAT(KMAT,KMAT),
     *                  B(K*K*(P+1)), SUM, ZERO, ONE
      PARAMETER         (ZERO = 0.0, ONE = 1.0)
      EXTERNAL          DECOMP, SOLVE
C
C     This auxiliary routine calculates the theoretical cross
C     covariances between the W(t)'s and between the W(t)'s and the
C     E(t)'s
C
      IFAULT = 0
C
C     generate the q gamma's
C
      DO 140 M = 1, Q
         DO 120 I = 1, K
            DO 100 J = 1, K
               SUM = ZERO
               DO 20 I2 = 1, K
                  SUM = SUM - THETA(I, (M - 1) * K + I2) * QQ(I2,J)
   20          CONTINUE
               DO 60 K2 = 1, P
                  DO 40 J2 = 1, K
                     IF(M .GT. K2) SUM = SUM + PHI(I, (K2 - 1) * K + J2)
     *                                 * GAMMA(J2, (M - K2 - 1) * K + J)
                     IF(M .EQ. K2) SUM = SUM + PHI(I, (K2 - 1) * K + J2)
     *                                   * QQ(J2,J)
   40             CONTINUE
   60          CONTINUE
               GAMMA(I, (M - 1) * K + J) = SUM
  100       CONTINUE
  120    CONTINUE
  140 CONTINUE
C
C     calculate the first p c(j)'s
C
C     first initialise MAT and Z to zero
C
      KW = K * K * (P + 1)
      DO 180 J = 1, KW
         DO 160 I = 1, KW
            MAT(I,J) = ZERO
  160    CONTINUE
         Z(J) = ZERO
  180 CONTINUE
C
      DO 360 M = 0, P
         DO 340 I = 1, K
            DO 320 J = 1, K
               L = M * K * K + (I - 1) * K + J
C
C              first calculate right-hand side vector Z
C
               IF(M .EQ. 0) Z(L) = QQ(I,J)
               IF(M .GT. 0 .AND. M .LE. Q) THEN
                  DO 200 I2 = 1, K
                     Z(L) = Z(L) - QQ(I,I2) * THETA(J, (M - 1) * K + I2)
  200             CONTINUE
               ENDIF
               DO 260 L4 = M + 1, Q
                  DO 240 I2 = 1, K
                     Z(L) = Z(L) - GAMMA(I, (L4 - M - 1) * K + I2) *
     *                     THETA(J, (L4 - 1) * K + I2)
  240             CONTINUE
  260          CONTINUE
C
C              now set up MAT array
C
               MAT(L,L) = ONE
               DO 300 I2 = 1, P
                  DO 280 K2 = 1, K
                     IF(M .GE. I2) L4 = (M - I2) * K * K + (I - 1) * K
     *                                  + K2
                     IF(M .LT. I2) L4 = (I2 - M) * K * K + (K2 - 1)
     *                                  * K + I
                     MAT(L,L4) = MAT(L,L4) - PHI(J, (I2 - 1) * K + K2)
  280             CONTINUE
  300          CONTINUE
C
  320       CONTINUE
  340    CONTINUE
  360 CONTINUE
C
      IF(P .GT. 0) THEN
         CALL DECOMP(KW, KMAT, MAT, IWORK)
         IF(IWORK(KW) .EQ. 0) THEN
            IFAULT = 1
            RETURN
         ELSE
           CALL SOLVE(KW, KMAT, MAT, Z, IWORK)
         ENDIF
      ENDIF
C
C     calculate C(j)'s for j = p + 1,...,r (if needed)
C
      DO 540 M = P + 1, R
         DO 520 I = 1, K
            DO 500 J = 1, K
               SUM = ZERO
               DO 400 L4 = 1, P
                  DO 380 I2 = 1, K
                     SUM = SUM + Z((M - L4) * K * K + (I - 1) * K + I2)
     *                     * PHI(J, (L4 - 1) * K + I2)
  380             CONTINUE
  400          CONTINUE
C
               IF(M .LE. Q) THEN
                  DO 420 I2 = 1, K
                     SUM = SUM - QQ(I,I2) * THETA(J, (M - 1) * K + I2)
  420             CONTINUE
                  DO 460 I2 = M + 1, Q
                     DO 440 L4 = 1, K
                        SUM = SUM - GAMMA(I, (I2 - M - 1) * K + L4)
     *                        * THETA(J, (I2 - 1) * K + L4)
  440                CONTINUE
  460             CONTINUE
               ENDIF
               Z(M * K * K + (I - 1) * K + J) = SUM
  500       CONTINUE
  520    CONTINUE
  540 CONTINUE
C
      END
C
      SUBROUTINE CHOL(A, IK, K, L, IFAULT)
C
C     ALGORITHM AS 242.4  APPL. STATIST. (1989), VOL.38, NO. 1
C
      INTEGER          IK, K, I, J, IFAULT, K2
      REAL             A(IK,K), L(IK,K), SUM, ZERO
      PARAMETER        (ZERO = 0.0)
      INTRINSIC        SQRT
C
      IFAULT = 0
      DO 180 J = 1, K
         SUM = A(J,J)
         DO 60 K2 = 1, J - 1
            SUM = SUM - L(J,K2) * L(J,K2)
   60    CONTINUE
         IF(SUM .LE. ZERO) THEN
            IFAULT = 1
            RETURN
         ELSE
            L(J,J) = SQRT(SUM)
         ENDIF
         DO 160 I = J + 1, K
            SUM = A(I,J)
            DO 120 K2 = 1, J - 1
               SUM = SUM - L(I,K2) * L(J,K2)
  120       CONTINUE
            L(I,J) = SUM / L(J,J)
            L(J,I) = ZERO
  160    CONTINUE
  180 CONTINUE
C
      END
C
      SUBROUTINE BKSB(L, IK, K, M, UPPER, B, IFAULT)
C
C     ALGORITHM AS 242.5  APPL. STATIST. (1989), VOL.38, NO. 1
C
      INTEGER          IK, K, IFAULT, I, J, M, J2
      REAL             L(IK,K), B(IK,M), SUM, ZERO
      LOGICAL          UPPER
      PARAMETER        (ZERO = 0.0)
C
      IFAULT = 0
      IF(UPPER) THEN
         DO 21 J2 = 1, M
            DO 40 I = K, 1, -1
               SUM = B(I,J2)
               DO 20 J = I + 1, K
                  SUM = SUM - L(I,J) * B(J,J2)
   20          CONTINUE
               IF(L(I,I) .NE. ZERO) THEN
                  B(I,J2) = SUM / L(I,I)
               ELSE
                  IFAULT = 1
                  RETURN
               ENDIF
   40       CONTINUE
   21    CONTINUE
      ELSE
         DO 22 J2 = 1, M
            DO 80 I = 1, K
               SUM = B(I,J2)
               DO 60 J = 1, I - 1
                  SUM = SUM - L(I,J) * B(J,J2)
   60          CONTINUE
               IF(L(I,I) .NE. ZERO) THEN
                  B(I,J2) = SUM / L(I,I)
               ELSE
                  IFAULT = 1
                  RETURN
               ENDIF
   80       CONTINUE
   22    CONTINUE
      ENDIF
C
      END
C
      SUBROUTINE DECOMP(N, NDIM, A, IP)
      INTEGER          N, NDIM, K, KP1, M, I, J, IP(NDIM)
      REAL             A(NDIM,NDIM), T
      INTRINSIC        ABS
C
C     matrix triangularization by Gaussian elimination.
C     INPUT..
C       N = order of matrix.
C       NDIM = declared dimension of array A.
C       A = matrix to be triangularized.
C     OUTPUT..
C       A(I,J), I .LE. J = upper triangular factor, U.
C       A(I,J), I .GT. J = multipliers = lower triangular
C                          factor, I - L.
C       IP(K),  K .LT. N = index of k-th pivot row.
C       IP(N) = (-1) ** (number of interchanges) or 0.
C     use 'SOLVE' to obtain solution of linear system.
C     determ(A) = IP(N) * A(1,1) * A(2,2) * ... * A(N,N).
C     if IP(N) = 0, A is singular, 'SOLVE' will divide by zero.
C     interchanges finished in U, only partly in L.
C
      IP(N) = 1
      DO 6 K = 1, N
         IF(K .EQ. N) GOTO 5
         KP1 = K + 1
         M = K
         DO 1 I = KP1, N
            IF(ABS(A(I,K)) .GT. ABS(A(M,K))) M = I
    1    CONTINUE
         IP(K) = M
         IF(M .NE. K) IP(N) = - IP(N)
         T = A(M,K)
         A(M,K) = A(K,K)
         A(K,K) = T
         IF(T .EQ. 0.0) GOTO 5
         DO 2 I = KP1, N
    2    A(I,K) = - A(I,K) / T
         DO 4 J = KP1, N
            T = A(M,J)
            A(M,J) = A(K,J)
            A(K,J) = T
            IF(T .EQ. 0.0) GOTO 4
            DO 3 I = KP1, N
    3       A(I,J) = A(I,J) + A(I,K) * T
    4    CONTINUE
    5    IF(A(K,K) .EQ. 0.0) IP(N) = 0
    6 CONTINUE
      RETURN
      END
C
      SUBROUTINE SOLVE(N, NDIM, A, B, IP)
      INTEGER          I, N, NDIM, IP(NDIM), NM1, K, KP1, M, KM1, KB
      REAL             A(NDIM,NDIM), B(NDIM), T
C
C     solution of linear system, A * x = b.
C     INPUT..
C       N = order of matrix.
C       NDIM = declared dimension of array A.
C       A = triangularized matrix obtained from 'DECOMP'.
C       B = right-hand side vector.
C       IP = pivot vector obtained from 'DECOMP'.
C     DO NOT USE if 'DECOMP' has set IP(N) = 0.
C     OUTPUT..
C       B = solution vector, X.
C
      IF(N .EQ. 1) GOTO 9
      NM1 = N - 1
      DO 7 K = 1, NM1
         KP1 = K + 1
         M = IP(K)
         T = B(M)
         B(M) = B(K)
         B(K) = T
         DO 7 I = KP1, N
    7 B(I) = B(I) + A(I,K) * T
      DO 8 KB = 1, NM1
         KM1 = N - KB
         K = KM1 + 1
         B(K) = B(K) / A(K,K)
         T = - B(K)
         DO 8 I = 1, KM1
    8 B(I) = B(I) + A(I,K) * T
    9 B(1) = B(1) / A(1,1)
      RETURN
      END
