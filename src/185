      SUBROUTINE ASRCH(NVAR, DIM, TABLE, NTAB, FIT, TFIT, INIT, CONFIG,
     *   NCON, MCON, CONFG1, CONFG2, CONFG3, MARG, NMAR, U, NU, MAXDEV,
     *   MAXIT, DEV, L1, L2, IA, MA, IB, MB, NI, NA, IFAULT)
C
C     ALGORITHM AS 185  APPL. STATIST. (1982) VOL.31, NO.3
C
C     This algorithm searches for a best log-linear model in contingency
C     tables using backward elimination procedure.
C
C     This version permits up to 7 variables.   If this limit is too
C     small, change the value of MAXVAR to its new limit, and the value
C     of the dimension of LOCMAR, MASK, & DELEF to the values of MCON,
C     MAXVAR + 1, and MAXVAR respectively.   See also changes needed in
C     subroutines LOGLIN (AS 51, Haberman 1972), LOWEFF, and DEGREF.
C
C     AS R87 has been incorporated from AS vol.40(2)
C
C     Auxiliary routines required: AS 51 and AS 66
C
      INTEGER DIM(NVAR), LOCMAR(35), CONFIG(NVAR, MCON), CONFG1(NVAR,
     *  MCON), CONFG2(NVAR, MCON), CONFG3(NVAR, MCON), DELEF(7),
     *  MASK(8), IB(NVAR, MB), SIZE
      REAL TABLE(NTAB), FIT(NTAB), TFIT(NTAB), MARG(NMAR), U(NU),
     *  MAXDEV, DEV(MAXIT), IA(4, MA), L1, L2, LR, LR1, LR3, ZERO
      DATA MAXVAR /7/, ZERO /0.0/
C
      IFAULT = 0
      NA = 0
      NB = 0
      IF (NVAR .GT. 0 .AND. NVAR .LE. MAXVAR) GO TO 10
    5 IFAULT = 4
      RETURN
C
   10 SIZE = 1
      DO 12 I = 1, NVAR
        IF (DIM(I) .LE. 0) GO TO 5
        SIZE = SIZE * DIM(I)
   12 CONTINUE
      IF (SIZE .NE. NTAB) GO TO 5
      DO 15 I = 1, NTAB
        IF (TABLE(I) .LT. ZERO .OR. FIT(I) .LT. ZERO) GO TO 5
        TFIT(I) = FIT(I)
   15 CONTINUE
C
C     Select a best initial fit.   Skip step if given initial fit by user.
C
      DO 70 NWAY = 1, NVAR
        IF (INIT .EQ. 1) GO TO 50
C
C     Construct CONFIG for model with full I-order effects
C
        NCON = 1
        DO 20 I = 1, NWAY
   20   NCON = NCON * (NVAR - I + 1) / I
        IF (MCON .LT. NCON) GO TO 400
        DO 25 J = 1, NVAR
          MASK(J) = J - 1
          DO 25 I = 1, NCON
            CONFIG(J,I) = 0
   25   CONTINUE
        MASK(NWAY+1) = NVAR
        DO 45 L = 1, NCON
          DO 35 J = 1, NWAY
   35     CONFIG(J,L) = MASK(J) + 1
          DO 40 J = 1, NWAY
            MASK(J) = MASK(J) + 1
            IF (MASK(J) .LT. MASK(J+1)) GO TO 45
            MASK(J) = J - 1
   40     CONTINUE
   45   CONTINUE
C
C     LOGLIN performs an iterative proportional fit of the marginal
C     totals (described in CONFIG) for a contingency table.
C
   50   CALL LOGLIN(NVAR, DIM, NCON, CONFIG, NTAB, TABLE, FIT, LOCMAR,
     *     NMAR, MARG, NU, U, MAXDEV, MAXIT, DEV, NLAST, IFAULT)
        IF (IFAULT .GT. 0) RETURN
        NA = NA + 1
        NI = NA
        CALL LRCHSQ(TABLE, FIT, NTAB, LR)
        CALL DEGREF(NVAR, DIM, NCON, NTAB, CONFIG, DF)
        P = XTAIL(LR, DF)
        CALL CPYLR(LR, DF, P, IA(1,NA), IA(2,NA), IA(3,NA))
        IF (NA .GT. MA .OR. (NB + NCON) .GT. MB) GO TO 400
        IA(4,NA) = NCON
        DO 55 K = 1, NCON
          NB = NB + 1
          DO 55 I = 1, NVAR
            IB(I,NB) = CONFIG(I,K)
   55   CONTINUE
        IF (INIT .EQ. 1) GO TO 90
        IF (NWAY .EQ. NVAR) RETURN
        IF (P .LE. L1) GO TO 70
        IF (NA .EQ. 1) GO TO 65
        LA = NA - 1
        IF (IA(3,LA) .LE. L1) GO TO 60
C
C     If the previous model is a good fit and the current fit does not
C     improve significantly over it, use the previous fit as the best
C     initial fit.
C
        DLR = IA(1,LA) - IA(1,NA)
        DDF = IA(2,LA) - IA(2,NA)
        DP = XTAIL(DLR,DDF)
        IF (DP .LE. L2) GO TO 60
        CALL CPYCON(CONFG1, CONFIG, NVAR, NCON1, NCON)
        CALL CPYLR(IA(1,LA), IA(2,LA), IA(3,LA), LR, DF, P)
        GO TO 90
   60   IF (NWAY .EQ. NVAR - 1) GO TO 90
   65   CALL CPYCON(CONFIG, CONFG1, NVAR, NCON, NCON1)
   70 CONTINUE
C
C     Backward elimination of effects from the best initial fit.
C
   90 PP = ZERO
      N = 0
      NA = NA + 1
      IF (NA .GT. MA .OR. NB + NCON .GT. MB) GO TO 400
      CALL CPYLR(LR, DF, P, IA(1,NA), IA(2,NA), IA(3,NA))
      IA(4,NA) = NCON
      DO 100 K = 1, NCON
        NB = NB + 1
        DO 100 I = 1, NVAR
          IB(I,NB) = CONFIG(I,K)
  100 CONTINUE
      LCON = NCON - 1
      DO 300 J = 1, NCON
C
C     Re-initialize FIT when fitting a simpler model
C
        DO 150 I = 1, NTAB
  150   FIT(I) = TFIT(I)
        CALL CPYCON(CONFIG, CONFG1, NVAR, NCON, NCON1)
C
C     Store effect to be deleted in DELEF and compress CONFG1.
C     Obtain lower order effects of DELEF not implied in CONFG1 and
C     store them in CONFG2.
C
        NCON2 = NCON1 - J + 1
        DO 200 I = 1, NVAR
          DELEF(I) = CONFG1(I,NCON2)
          CONFG1(I,NCON2) = CONFG1(I,NCON1)
  200   CONTINUE
        CALL LOWEFF(NVAR, LCON, MCON, CONFG1, CONFG2, DELEF, JJ)
        NCON1 = LCON + JJ
        IF (JJ .EQ. 0) GO TO 225
        DO 220 L = 1, JJ
          LL = LCON + L
          DO 220 I = 1, NVAR
            CONFG1(I,LL) = CONFG2(I,L)
  220   CONTINUE
  225   NA = NA + 1
        IF (NA .GT. MA .OR. NB + NCON1 .GT. MB) GO TO 400
        CALL LOGLIN(NVAR, DIM, NCON1, CONFG1, NTAB, TABLE, FIT, LOCMAR,
     *      NMAR, MARG, NU, U, MAXDEV, MAXIT, DEV, NLAST, IFAULT)
        IF (IFAULT .GT. 0) RETURN
        CALL LRCHSQ(TABLE, FIT, NTAB, LR1)
        CALL DEGREF(NVAR, DIM, NCON1, NTAB, CONFG1, DF1)
        P1 = XTAIL(LR1, DF1)
        CALL CPYLR(LR1, DF1, P1, IA(1,NA), IA(2,NA), IA(3,NA))
        IA(4,NA) = NCON1
        DO 250 K = 1, NCON1
          NB = NB + 1
          DO 250 I = 1, NVAR
            IB(I,NB) = CONFG1(I,K)
  250   CONTINUE
        DLR = LR1 - LR
        DDF = DF1 - DF
        DP = XTAIL(DLR, DDF)
C
C     Select the best fitted model that does not result in a significant
C     increase in LR chi-squared value for further simplification.
C
        IF (P1 .LT. L1 .OR. DP .LT. L2) GO TO 300
        IF (P1 .LT. PP) GO TO 300
        PP = P1
        N = N + 1
        CALL CPYCON(CONFG1, CONFG3, NVAR, NCON1, NCON3)
        CALL CPYLR(LR1, DF1, P1, LR3, DF3, P3)
  300 CONTINUE
      IF (N .EQ. 0) GO TO 500
C
C     Copy the best simplified fit of each iteration to CONFIG.
C
      CALL CPYCON(CONFG3, CONFIG, NVAR, NCON3, NCON)
      CALL CPYLR(LR3, DF3, P3, LR, DF, P)
      IF (NCON .EQ. 1 .AND. CONFIG(2,1) .EQ. 0) GO TO 500
      GO TO 90
C
  400 IFAULT = 5
      RETURN
C
C     Obtain fitted table of the best final / intermediate fit.
C
  500 CALL LOGLIN(NVAR, DIM, NCON, CONFIG, NTAB, TABLE, FIT, LOCMAR,
     *   NMAR, MARG, NU, U, MAXDEV, MAXIT, DEV, NLAST, IFAULT)
      RETURN
      END
C
      SUBROUTINE LOWEFF(NVAR, NCON, MCON, CONFG1, CONFG2, DELEF, JJ)
C
C     ALGORITHM AS 185.1  APPL. STATIST. (1982) VOL.31, NO.3
C
C     Find lower order effects of DELEF not implied in CONFG1 and sore
C     them in CONFG2.  If NVAR is to be greater than 7, the dimensions
C     of MASK, MASK1 and MASK2 must all be increased to the value NVAR.
C
      INTEGER CONFG1(NVAR,NCON), CONFG2(NVAR,MCON), DELEF(NVAR), 
     *   MASK(7), MASK1(7), MASK2(7)
C
C     Determine the order of the deleted effect.
C
      NWAY = 0
      DO 100 I = 1, NVAR
  100 IF (DELEF(I) .NE. 0) NWAY = NWAY + 1
      JJ = 0
      IF (NWAY .EQ. 1) RETURN
C
C     Find all possible lower order effects of the deleted effect.
C
      DO 300 J = 1, NWAY
        MASK(J) = J - 1
        DO 300 I = 1, NVAR
          CONFG2(I,J) = 0
  300 CONTINUE
      N = NWAY - 1
      MASK(N+1) = NWAY
  350 CONTINUE
      DO 360 I = 1, NVAR
  360 MASK1(I) = 0
      JJ = JJ + 1
      DO 400 J = 1, N
        L = MASK(J) + 1
        K = DELEF(L)
        MASK1(K) = 1
        CONFG2(J,JJ) = K
  400 CONTINUE
      IF (NCON .EQ. 0) GO TO 810
      NN = 0
C
C     Check if each lower order effect is implied in CONFG1.
C
      DO 700 J = 1, NCON
        DO 450 I = 1, NVAR
  450   MASK2(I) = 0
        DO 500 I = 1, NVAR
          K = CONFG1(I,J)
          IF (K .NE. 0) MASK2(K) = 1
  500   CONTINUE
        DO 550 K = 1, NVAR
          IF (MASK1(K) .EQ. 1 .AND. MASK2(K) .NE. 1) GO TO 650
  550   CONTINUE
        GO TO 800
  650   NN = NN + 1
  700 CONTINUE
      IF (NN .EQ. NCON) GO TO 810
  800 JJ = JJ - 1
  810 DO 850 J = 1, N
        MASK(J) = MASK(J) + 1
        IF (MASK(J) .LT. MASK(J+1)) GO TO 350
        MASK(J) = J - 1
  850 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE LRCHSQ(TABLE, FIT, NTAB, RATIO)
C
C     ALGORITHM AS 185.2  APPL. STATIST. (1982) VOL.31, NO.3
C
C     LRCHSQ calculates twice the likelihood ratio.
C
      REAL TABLE(NTAB), FIT(NTAB), ZERO
      DATA ZERO /0.0/
C
      SUM = ZERO
      DO 10 I = 1, NTAB
C
C     Note that FIT(I) = 0 implies TABLE(I) = 0.
C
   10 IF (TABLE(I) .NE. ZERO) SUM = SUM + TABLE(I) * LOG(TABLE(I) /
     *    FIT(I))
      RATIO = SUM + SUM
      RETURN
      END
C
      SUBROUTINE DEGREF(NVAR, DIM, NCON, NTAB, CONFIG, DF)
C
C     ALGORITHM AS 185.3  APPL. STATIST. (1982) VOL.31, NO.3
C
C     DEGREF computes the degrees of freedom of a model described in
C     CONFIG.   If NVAR is to be greater than 7, the dimensions of MASK
C     and HASH must be increased to the values (NVAR + 1) and
C     (2.NVAR - 1) respectively).
C
      INTEGER DIM(NVAR), CONFIG(NVAR,NCON), MASK(8), HASH(127)
C
      NP = 0
      N = 2 ** NVAR - 1
      DO 5 I = 1, N
    5 HASH(I) = 0
      DO 100 I = 1, NCON
C
C     Find the order of each effect in CONFIG.
C
        DO 10 J = 1, NVAR
        IF (CONFIG(J,I) .EQ. 0) GO TO 20
   10   CONTINUE
        J = NVAR + 1
   20   N = J - 1
C
C     Set up MASK.
C
        DO 30 J = 1, N
   30   MASK(J) = 0
C
C     Use MASK to find all possible lower order effects.
C
   40   DO 50 J = 1, N
          MASK(J) = MASK(J) + 1
          IF (MASK(J) .EQ. 1) GO TO 60
          MASK(J) = 0
   50   CONTINUE
        GO TO 100
   60   K = 0
        L = 1
        DO 80 J = 1, N
          IF (MASK(J) .NE. 1) GO TO 80
C
C     Get the lower order effect degrees of freedom.
C
          IJ = CONFIG(J,I)
          L = L * (DIM(IJ) - 1)
          K = K + 2 ** (CONFIG(J,I) - 1)
   80   CONTINUE
C
C     Check if the lower order effect already existed.
C
        IF (HASH(K) .EQ. 1) GO TO 40
        HASH(K) = 1
        NP = NP + L
        GO TO 40
  100 CONTINUE
      DF = NTAB - NP - 1
C
      RETURN
      END
C
      REAL FUNCTION XTAIL(RATIO, DF)
C
C     ALGORITHM AS 185.4  APPL. STATIST. (1982) VOL.31, NO.3
C
C     Computes probability of a chi-squared variable > RATIO using an
C     approximation from JASA (1968) vol.68, p.1420.
C
      REAL RATIO, DF, ZERO, ONE, TWO, SIX, PT04
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, SIX /6.0/, PT04 /0.04/
C
      IF (RATIO .EQ. ZERO) GO TO 50
      IF (DF .EQ. ONE) GO TO 100
      S = (DF - ONE) / TWO
      T = RATIO / TWO
      D = T - S - ONE / SIX - PT04 / DF
      Z = (D * SQRT(TWO * S * LOG(S / T) + TWO * (T - S))) / ABS(S - T)
C
C     ALNORM (AS 66, Hill 1973) evaluates the tail prob. of Z = N(0,1)
C
      XTAIL = ALNORM(Z, .TRUE.)
      RETURN
C
   50 XTAIL = ONE
      RETURN
  100 Z = SQRT(RATIO)
      XTAIL = TWO * ALNORM(Z, .TRUE.)
      RETURN
      END
C
      SUBROUTINE CPYCON(CONFG1, CONFG2, NVAR, NCON1, NCON2)
C
C     ALGORITHM AS 185.5  APPL. STATIST. (1982) VOL.31, NO.3
C
      INTEGER CONFG1(NVAR,NCON1), CONFG2(NVAR,NCON1)
C
      DO 10 J = 1, NCON1
        DO 10 I = 1, NVAR
          CONFG2(I,J) = CONFG1(I,J)
   10 CONTINUE
      NCON2 = NCON1
      RETURN
      END
C
      SUBROUTINE CPYLR(LR1, DF1, P1, LR2, DF2, P2)
C
C     ALGORITHM AS 185.6  APPL. STATIST. (1982) VOL.31, NO.3
C
      REAL LR1, DF1, P1, LR2, DF2, P2
C
      LR2 = LR1
      DF2 = DF1
      P2 = P1
      RETURN
      END
