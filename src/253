      SUBROUTINE RCM(TAB, FIT, FF, L, LR, LC, MU, NU, PHI, DETAIL, TOL,
     *               RWT, CWT, WRK, WRK2, ITER, NTAB, NR, NC, NM, NW,
     *               NW2, MAXIT, IFAULT)
C
C        ALGORITHM AS 253.1  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Maximum likelihood estimation of the RC(M) association model
C
      INTEGER FF(NTAB), ITER, NTAB, NR, NC, NM, NW, NW2, MAXIT, IFAULT
      REAL TAB(NTAB), FIT(NTAB), L, LR(NR), LC(NC), MU( * ), NU( * ),
     *     PHI( * ), DETAIL(MAXIT), TOL, RWT(NR), CWT(NC), WRK(NW),
     *     WRK2(NW2)
      INTEGER IPT(8), I, ICON, II, J, MNRC, JOB
      REAL SALE, TOT, ZERO
      DATA ZERO / 0.0 /
C
C        Check for errors
C
      IFAULT = 0
      IF (NM .GT. (NR - 1) .OR. NM .GT. (NC - 1)) IFAULT = 1
      IF (NM .LT. 0) IFAULT = 1
      IF (TOL .LE. ZERO) IFAULT = 2
      IF (NR .LT. 2 .OR. NC .LT. 2) IFAULT = 3
      DO 10 J = 1, NTAB
         IF (TAB(J) .LT. ZERO) IFAULT = 4
   10 CONTINUE
      IF (NTAB .LT. (NR * NC)) IFAULT = 6
      IF (NW .LT. (2 * (NR + NC) + (NR * NC))) IFAULT = 7
      IF (NW2 .LT. ((MAX(NR, NC)) + (NC + 1) + (NC ** 2) + NR +
     *    NC)) IFAULT = 7
      DO 20 I = 1, NTAB
         IF (FF(I) .NE. 1 .AND. FF(I) .NE. 0) IFAULT = 8
   20 CONTINUE
      IF (IFAULT .NE. 0) RETURN
C
C        Set pointers to locations in work space
C
      MNRC = MAX(NR, NC)
C
      IPT(1) = NR + 1
      IPT(2) = NR + NC + 1
      IPT(3) = 2 * NR + NC + 1
      IPT(4) = 2 * (NR + NC) + 1
      IPT(5) = MNRC * NC + 1
      IPT(6) = (MNRC * NC) + NC + 2
      IPT(7) = (MNRC * NC) + (NC * NC) + NC + 2
      IPT(8) = (MNRC * NC) + (NC * NC) + NR + NC + 2
C
C        Get the marginal totals for the observed table
C
      CALL MARGIN(TAB, WRK(1), WRK(IPT(1)), TOT, NR, NC)
C
C        Compute intial parameter estimates
C
      JOB = 21
      CALL START(TAB, L, LR, LC, MU, NU, PHI, RWT, CWT, WRK2(1),
     *           WRK2(IPT(5)), WRK2(IPT(6)), WRK2(IPT(7)), WRK2(IPT(8)),
     *           WRK(IPT(4)), WRK(IPT(2)), WRK(IPT(3)), NR, NC, NM,
     *           MNRC, JOB, IFAULT)
      IF (IFAULT .NE. 0) THEN
         IFAULT = 5
         RETURN
      END IF
C
C        *****  Iterations  *****
C
      ITER = MAXIT
      DO 30 II = 1, MAXIT
C
C        Update the LAMBDAs
C
         CALL LAMFIT(1, TAB, FIT, FF, L, LR, LC, MU, NU, PHI, WRK(1),
     *               WRK(IPT(1)), WRK(IPT(2)), WRK(IPT(3)), TOT, NR, NC,
     *               NR, NC, NM)
         CALL LAMFIT(2, TAB, FIT, FF, L, LC, LR, NU, MU, PHI, WRK(1),
     *               WRK(IPT(1)), WRK(IPT(2)), WRK(IPT(3)), TOT, NR, NC,
     *               NC, NR, NM)
C
         IF (NM .GT. 0) THEN
C
C         Update the MU's and NU's
C
            CALL SCRFIT(1, TAB, FIT, FF, L, LR, LC, MU, NU, PHI, NR, NC,
     *                  NM, NR, NC)
            CALL SCRFIT(2, TAB, FIT, FF, L, LC, LR, NU, MU, PHI, NR, NC,
     *                  NM, NC, NR)
C
C         Update the PHI's
C
            CALL PHIFIT(TAB, FIT, FF, L, LR, LC, MU, NU, PHI, NR, NC,
     *                  NM)
C
C         Centre the MU's and NU's about zero
C
            CALL CENTRE(L, LR, LC, MU, NU, PHI, RWT, CWT, NR, NC, NM)
C
         END IF
C
C        Check for convergence
C
         CALL CONV(TAB, FIT, FF, WRK, L, LR, LC, MU, NU, PHI, IPT, TOT,
     *             TOL, SALE, NR, NC, NM, NW, ICON)
C
         DETAIL(II) = SALE
         IF (ICON .EQ. 1 .AND. NM .GT. 0) THEN
C
C         Rescale the MU's and NU's so that they satisfy the
C         identifying restrictions
C
            CALL SCALE(MU, NU, PHI, RWT, CWT, WRK2(1), WRK2(IPT(5)),
     *                 WRK2(IPT(6)), WRK2(IPT(7)), WRK2(IPT(8)), NR, NC,
     *                 NM, MNRC, JOB, IFAULT)
            IF (IFAULT .NE. 0) IFAULT = 5
C
            ITER = II
            RETURN
         ELSE IF (ICON .EQ. 1) THEN
            ITER = II
            RETURN
         END IF
   30 CONTINUE
C
      IFAULT = 9
      RETURN
      END
      SUBROUTINE MARGIN(T, R, C, N, NR, NC)
C
C        ALGORITHM AS 253.2  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Compute total sample size, and row and column marginal totals
C
      INTEGER NR, NC
      REAL T(NR, NC), R(NR), C(NC), N
      INTEGER I, J
      REAL ZERO
      DATA ZERO / 0.0 /
C
      DO 10 I = 1, NR
         R(I) = ZERO
   10 CONTINUE
      DO 20 J = 1, NC
         C(J) = ZERO
   20 CONTINUE
      N = ZERO
C
      DO 40 J = 1, NC
         DO 30 I = 1, NR
            R(I) = R(I) + T(I, J)
            C(J) = C(J) + T(I, J)
            N = N + T(I, J)
   30    CONTINUE
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE EXPTAB(TAB, FIT, FF, LAMBDA, LR, LC, MU, NU, PHI, NR,
     *                  NC, NM)
C
C        ALGORITHM AS 253.3  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Compute expected frequencies
C
      INTEGER FF(NR, NC), NR, NC, NM
      REAL TAB(NR, NC), FIT(NR, NC), LAMBDA, LR(NR), LC(NC), MU(NR, * ),
     *     NU(NC, * ), PHI( * )
      INTEGER I, J, L
      REAL ONE, X, ZERO
      DATA ZERO / 0.0 / ONE / 1.0 /
C
      DO 30 J = 1, NC
         DO 20 I = 1, NR
            X = ZERO
            IF (NM .GT. 0) THEN
               DO 10 L = 1, NM
                  X = X + MU(I, L) * NU(J, L) * PHI(L)
   10          CONTINUE
            END IF
            FIT(I, J) = (EXP(LAMBDA + LR(I) + LC(J) + X) * FF(I, J)) +
     *                  ((ONE - FF(I, J)) * TAB(I, J))
   20    CONTINUE
   30 CONTINUE
C
      RETURN
      END
      SUBROUTINE START(TAB, LAMBDA, LR, LC, M, N, P, RWT, CWT, A, W, V,
     *                 WORK, E, LTAB, LRM, LCM, NR, NC, NM, MNRC, JOB,
     *                 IFAULT)
C
C        ALGORITHM AS 253.4  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Compute initial parameter estimates
C
      INTEGER NR, NC, NM, MNRC, JOB, IFAULT
      REAL TAB(NR, NC), LAMBDA, LR(NR), LC(NC), M(NR, * ), N(NC, * ),
     *     P( * ), RWT(NR), CWT(NC), A(MNRC, NC), W( * ), V(NC, NC),
     *     WORK(NR), E(NC), LTAB(NR, NC), LRM(NR), LCM(NC)
      INTEGER I, IA, J, L
      REAL EPS, ZERO
      DATA EPS / 0.1 / ZERO / 0.0 /
C
      DO 20 J = 1, NC
         DO 10 I = 1, NR
            IF (TAB(I, J) .LE. EPS) THEN
               LTAB(I, J) = LOG(EPS)
            ELSE
               LTAB(I, J) = LOG(TAB(I, J))
            END IF
   10    CONTINUE
   20 CONTINUE
C
      CALL MARGIN(LTAB, LRM, LCM, LAMBDA, NR, NC)
      LAMBDA = LAMBDA / (NR * NC)
C
      DO 30 I = 1, NR
         LR(I) = LRM(I) / NC - LAMBDA
   30 CONTINUE
      DO 40 J = 1, NC
         LC(J) = LCM(J) / NR - LAMBDA
   40 CONTINUE
C
      IF (NM .GE. 1) THEN
         DO 60 J = 1, NC
            DO 50 I = 1, NR
               A(I, J) = LTAB(I, J) - (LAMBDA + LR(I) + LC(J))
               A(I, J) = A(I, J) * SQRT(RWT(I)) * SQRT(CWT(J))
   50       CONTINUE
   60    CONTINUE
C
         IF (MNRC .GT. NR) THEN
            IA = NR + 1
            DO 80 J = 1, NC
               DO 70 I = IA, MNRC
                  A(I, J) = ZERO
   70          CONTINUE
   80       CONTINUE
         END IF
C
C        Call LINPACK SVD subroutine
C
         CALL SSVDC(A, MNRC, NR, NC, W, E, A, MNRC, V, NC, WORK, JOB,
     *              IFAULT)
C
         DO 110 L = 1, NM
            P(L) = W(L)
            DO 90 I = 1, NR
               M(I, L) = A(I, L) / SQRT(RWT(I))
   90       CONTINUE
C
            DO 100 J = 1, NC
               N(J, L) = V(J, L) / SQRT(CWT(J))
  100       CONTINUE
  110    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE LAMFIT(IFLAG, TAB, FIT, FF, LAMBDA, LAM1, LAM2, M1, M2,
     *                  P, RM, CM, ERM, ECM, TOT, NR, NC, N1, N2, NM)
C
C        ALGORITHM AS 253.5  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Update the estimates of the lambdas
C        IFLAG =1 => Update row parameters, else do column parameters
C
      INTEGER FF(NR, NC), IFLAG, NR, NC, NM, N1, N2
      REAL TAB(NR, NC), FIT(NR, NC), LAMBDA, LAM1(N1), LAM2(N2),
     *     M1(N1, * ), M2(N2, * ), P( * ), RM(NR), CM(NC), ERM(NR),
     *     ECM(NC), TOT
      INTEGER II, NN
      REAL ETOT
C
      IF (IFLAG .EQ. 1) THEN
         CALL EXPTAB(TAB, FIT, FF, LAMBDA, LAM1, LAM2, M1, M2, P, NR,
     *               NC, NM)
         NN = NR
      ELSE
         NN = NC
         CALL EXPTAB(TAB, FIT, FF, LAMBDA, LAM2, LAM1, M2, M1, P, NR,
     *               NC, NM)
      END IF
C
      CALL MARGIN(FIT, ERM, ECM, ETOT, NR, NC)
C
      DO 10 II = 1, NN
         IF (IFLAG .EQ. 1) THEN
            LAM1(II) = LAM1(II) + (RM(II) - ERM(II)) / ERM(II)
         ELSE
            LAM1(II) = LAM1(II) + (CM(II) - ECM(II)) / ECM(II)
         END IF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SCRFIT(IFLAG, TAB, FIT, FF, LAMBDA, LAM1, LAM2, M1, M2,
     *                  P, NR, NC, NM, N1, N2)
C
C        ALGORITHM AS 253.6  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Update the estimates of the MU's and NU's
C        IFLAG = 1 => Update MU's, else do NU's
C
      INTEGER FF(NR, NC), IFLAG, NR, NC, NM, N1, N2
      REAL TAB(NR, NC), FIT(NR, NC), LAMBDA, LAM1(N1), LAM2(N2),
     *     M1(N1, * ), M2(N2, * ), P( * )
      INTEGER I1, I2, L
      REAL C1, C2, X, XX, ZERO
      DATA ZERO / 0.0 /
C
      DO 30 L = 1, NM
         IF (IFLAG .EQ. 1) THEN
            N1 = NR
            N2 = NC
            CALL EXPTAB(TAB, FIT, FF, LAMBDA, LAM1, LAM2, M1, M2, P, NR,
     *                  NC, NM)
         ELSE
            N1 = NC
            N2 = NR
            CALL EXPTAB(TAB, FIT, FF, LAMBDA, LAM2, LAM1, M2, M1, P, NR,
     *                  NC, NM)
         END IF
C
         DO 20 I1 = 1, N1
            X = ZERO
            XX = ZERO
            DO 10 I2 = 1, N2
               IF (IFLAG .EQ. 1) THEN
                  C1 = TAB(I1, I2) - FIT(I1, I2)
                  C2 = FIT(I1, I2)
               ELSE
                  C1 = TAB(I2, I1) - FIT(I2, I1)
                  C2 = FIT(I2, I1)
               END IF
               X = X + M2(I2, L) * C1
               XX = XX + M2(I2, L) * M2(I2, L) * C2
   10       CONTINUE
            M1(I1, L) = M1(I1, L) + X / (XX * P(L))
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE PHIFIT(TAB, FIT, FF, LAMBDA, LAM1, LAM2, M1, M2, P, NR,
     *                  NC, NM)
C
C        ALGORITHM AS 253.7  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Update the estimates of the PHI's
C
      INTEGER FF(NR, NC), NR, NC, NM
      REAL TAB(NR, NC), FIT(NR, NC), LAMBDA, LAM1(NR), LAM2(NC),
     *     M1(NR, * ), M2(NC, * ), P( * )
      INTEGER I, J, L
      REAL X, XX, ZERO
      DATA ZERO / 0.0 /
C
      DO 30 L = 1, NM
         CALL EXPTAB(TAB, FIT, FF, LAMBDA, LAM1, LAM2, M1, M2, P, NR,
     *               NC, NM)
         X = ZERO
         XX = ZERO
         DO 20 J = 1, NC
            DO 10 I = 1, NR
               X = X + M1(I, L) * M2(J, L) * (TAB(I, J) - FIT(I, J))
               XX = XX + M1(I, L) * M1(I, L) * M2(J, L) * M2(J, L) *
     *              FIT(I, J)
   10       CONTINUE
   20    CONTINUE
         P(L) = P(L) + X / XX
   30 CONTINUE
C
      RETURN
      END
      SUBROUTINE CENTRE(LAMBDA, LR, LC, MU, NU, P, RWT, CWT, NR, NC, NM)
C
C        ALGORITHM AS 253.8  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Centre (with respect to the row/column weights) the MUs/NU's
C        about zero
C
      INTEGER NR, NC, NM
      REAL LAMBDA, LR(NR), LC(NC), MU(NR, * ), NU(NC, * ), P( * ),
     *     RWT(NR), CWT(NC)
      INTEGER I, J, L
      REAL C1, C2, X, XX, ZERO
      DATA ZERO / 0.0 /
C
      C1 = ZERO
      C2 = ZERO
C
      DO 10 I = 1, NR
         C1 = C1 + RWT(I)
   10 CONTINUE
      DO 20 J = 1, NC
         C2 = C2 + CWT(J)
   20 CONTINUE
C
      DO 70 L = 1, NM
         X = ZERO
         XX = ZERO
         DO 30 I = 1, NR
            X = X + MU(I, L) * RWT(I)
   30    CONTINUE
         DO 40 J = 1, NC
            XX = XX + NU(J, L) * CWT(J)
   40    CONTINUE
         DO 50 I = 1, NR
            LR(I) = LR(I) + (XX / C2 * P(L) * MU(I, L))
            MU(I, L) = MU(I, L) - X / C1
   50    CONTINUE
         DO 60 J = 1, NC
            LC(J) = LC(J) + (X / C1 * P(L) * NU(J, L))
            NU(J, L) = NU(J, L) - XX / C2
   60    CONTINUE
         LAMBDA = LAMBDA - (P(L) * X / C1 * XX / C2)
   70 CONTINUE
      RETURN
      END
      SUBROUTINE SCALE(M, N, P, RWT, CWT, A, W, V, WORK, E, NR, NC, NM,
     *                 MNRC, JOB, IFAULT)
C
C        ALGORITHM AS 253.9  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Scale the MU's/NU's so that they satisfy the
C        identifying restrictions
C
      INTEGER NR, NC, NM, MNRC, JOB, IFAULT
      REAL M(NR, * ), N(NC, * ), P( * ), RWT(NR), CWT(NC), A(MNRC, NC),
     *     W( * ), V(NC, NC), WORK(NR), E(NC)
      INTEGER I, IA, J, L
      REAL ZERO
      DATA ZERO / 0.0 /
C
      DO 30 J = 1, NC
         DO 20 I = 1, NR
            A(I, J) = ZERO
            DO 10 L = 1, NM
               A(I, J) = A(I, J) + P(L) * M(I, L) * N(J, L)
   10       CONTINUE
            A(I, J) = A(I, J) * SQRT(RWT(I)) * SQRT(CWT(J))
   20    CONTINUE
   30 CONTINUE
C
      IF (MNRC .GT. NR) THEN
         IA = NR + 1
         DO 50 J = 1, NC
            DO 40 I = IA, MNRC
               A(I, J) = ZERO
   40       CONTINUE
   50    CONTINUE
      END IF
C
C        Call LINPACK SVD subroutine
C
      CALL SSVDC(A, MNRC, NR, NC, W, E, A, MNRC, V, NC, WORK, JOB,
     *           IFAULT)
      IF (IFAULT .NE. 0) RETURN
C
      DO 80 L = 1, NM
         P(L) = W(L)
         DO 60 I = 1, NR
            M(I, L) = A(I, L) / SQRT(RWT(I))
   60    CONTINUE
C
         DO 70 J = 1, NC
            N(J, L) = V(J, L) / SQRT(CWT(J))
   70    CONTINUE
   80 CONTINUE
C
      RETURN
      END
      SUBROUTINE CONV(TAB, FIT, FF, WRK, LAMBDA, LR, LC, MU, NU, PHI,
     *                IPT, TOT, TOL, SALE, NR, NC, NM, NW, ICON)
C
C        ALGORITHM AS 253.10  APPL. STATIST. (1990) VOL.39, NO. 1
C
C        Check for convergence of the algorithm to a solution
C                of the likelihood equations
C
      INTEGER FF(NR, NC), IPT(6), NR, NC, NM, NW, ICON
      REAL TAB(NR, NC), FIT(NR, NC), WRK(NW), LAMBDA, LR(NR), LC(NC),
     *     MU(NR, * ), NU(NC, * ), PHI( * ), TOT, TOL, SALE
      INTEGER I, J, L
      REAL ETOT, X, XX, ZERO
      DATA ZERO / 0.0 /
C
      CALL EXPTAB(TAB, FIT, FF, LAMBDA, LR, LC, MU, NU, PHI, NR, NC, NM)
C
      CALL MARGIN(FIT, WRK(IPT(2)), WRK(IPT(3)), ETOT, NR, NC)
C
      X = ZERO
C
      DO 10 I = 1, NR
         X = X + ABS(WRK(I) - WRK(IPT(2) + I - 1))
   10 CONTINUE
C
      DO 20 J = 1, NC
         X = X + ABS(WRK(IPT(1) + J - 1) - WRK(IPT(3) + J - 1))
   20 CONTINUE
C
      IF (NM .GT. 0) THEN
         DO 50 L = 1, NM
            XX = ZERO
            DO 40 J = 1, NC
               DO 30 I = 1, NR
                  XX = XX + MU(I, L) * NU(J, L) *
     *                 (TAB(I, J) - FIT(I, J))
   30          CONTINUE
   40       CONTINUE
            X = X + ABS(XX)
   50    CONTINUE
C
         DO 80 L = 1, NM
            DO 70 I = 1, NR
               XX = ZERO
               DO 60 J = 1, NC
                  XX = XX + NU(J, L) * (TAB(I, J) - FIT(I, J))
   60          CONTINUE
               X = X + ABS(XX * PHI(L))
   70       CONTINUE
   80    CONTINUE
C
         DO 110 L = 1, NM
            DO 100 J = 1, NC
               XX = ZERO
               DO 90 I = 1, NR
                  XX = XX + MU(I, L) * (TAB(I, J) - FIT(I, J))
   90          CONTINUE
               X = X + ABS(XX * PHI(L))
  100       CONTINUE
  110    CONTINUE
      END IF
C
      SALE = X
      IF (SALE .LT. TOL) THEN
         ICON = 1
      ELSE
         ICON = 0
      END IF
C
      RETURN
      END
