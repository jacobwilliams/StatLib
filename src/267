      REAL FUNCTION DPROB(DIFF, NDF, NPOPS, MBEST, IFAULT)
C
C        ALGORITHM AS 267.1  APPL.STATIST. (1991), VOL.40, NO.3
C
C        Calculates the lower bound for the probability of correct
C        selection of a subset containing the best populations within
C        the set of NPOPS populations.
C
C     Auxiliary routines required: ALNORM (AS 66) and ALOGAM (CACM
C     algorithm 291 or AS 245).
C
      INTEGER IFAULT, MBEST, NDF, NPOPS
      REAL DIFF
C
      INTEGER I, J, K, JMAX, JMIN, KMAX, KMIN, L, MAXDF
      REAL DFW(2, 10), QW(2, 10), BEST, CLOG, CUTJ, CUTK, DF, DF2,
     *     EHJ, GK, GSTEP, HALF, HJ, HSTEP, ONE, PJ, PK, PROB, PZ1, PZ2,
     *     Q, R1, R2, W0, W2, ZERO
      REAL ALNORM, ALOGAM
      EXTERNAL ALNORM, ALOGAM
C
      DATA CUTK, CUTJ, GSTEP / 0.00003, 0.00001, 0.70 /
      DATA ZERO, HALF, ONE / 0.0, 0.5, 1.0 /
      DATA JMAX, JMIN, KMAX, KMIN, MAXDF / 10, 3, 15, 7, 120 /
C
C        Check input arguments
C
   10 PROB = ZERO
      IFAULT = 0
      IF (DIFF .LE. 0.0) IFAULT = 1
      IF (NDF .LE. 3) IFAULT = 20 + IFAULT
      IF (NPOPS .LT. 2 .OR. NPOPS .GT. 100) IFAULT = 300 + IFAULT
      IF (MBEST .LT. 1 .OR. MBEST .GE. NPOPS) IFAULT = 4000 + IFAULT
      IF (IFAULT .NE. 0) GO TO 100
C
C        Transfer arguments and calculate constants.
C
   20 Q = DIFF
      DF = NDF
      BEST = MBEST
      R1 = NPOPS - MBEST
      R2 = MBEST - 1
      CLOG = LOG(BEST * GSTEP) - 0.918938533
C
      IF (NDF .LE. MAXDF) THEN
         DF2 = DF * HALF
         HSTEP = GSTEP * DF ** (-HALF)
         HJ = ZERO
         CLOG = CLOG + LOG(HSTEP) - ALOGAM(DF2, IFAULT) -
     *          DF2 * (ONE - LOG(DF)) + (ONE - DF2) * 0.693147181
         DO 30 J = 1, JMAX
            HJ = HJ + HSTEP
            EHJ = EXP(HJ)
            QW(1, J) = Q * EHJ
            QW(2, J) = Q / EHJ
            DFW(1, J) = DF * (HALF + HJ - HALF * EHJ * EHJ)
            DFW(2, J) = DF * (HALF - HJ - HALF / EHJ / EHJ)
   30    CONTINUE
      END IF
C
C        Compute integral
C
   40 DO 90 K = 1, KMAX
         PK = ZERO
         DO 80 L = 1, 2
            IF (L .EQ. 1) THEN
               GK = GSTEP * REAL(K)
            ELSE
               GK = -GSTEP * REAL(K - 1)
            END IF
            W0 = CLOG - GK * GK * HALF
            W2 = ZERO
            PZ1 = ALNORM(-GK - Q, .TRUE.)
            PZ2 = ONE
            IF (MBEST .GT. 1) THEN
               PZ2 = ALNORM(GK, .TRUE.)
               IF (PZ2 .GT. ZERO) W2 = R2 * LOG(PZ2)
            END IF
            IF (PZ1 .GT. ZERO .AND. PZ2 .GT. ZERO) PK = PK +
     *          EXP(W0 + R1 * LOG(PZ1) + W2)
            IF (NDF .LE. MAXDF .AND. PZ2 .GT. ZERO) THEN
   50          DO 70 J = 1, JMAX
                  PJ = ZERO
                  DO 60 I = 1, 2
                     PZ1 = ALNORM(-GK - QW(I, J), .TRUE.)
                     IF (PZ1 .GT. ZERO) PJ = PJ +
     *                                       EXP(W0 + DFW(I, J) +
     *                                       R1 * LOG(PZ1) + W2)
   60             CONTINUE
                  PK = PK + PJ
                  IF ((J .GT. JMIN .OR. K .GT. KMIN) .AND.
     *                PJ .LT. CUTJ) GO TO 80
   70          CONTINUE
            END IF
   80    CONTINUE
         PROB = PROB + PK
         IF (PROB .GT. ONE) THEN
            PROB = ONE
            GO TO 100
         END IF
         IF (K .GT. KMIN .AND. PK .LT. CUTK) GO TO 100
   90 CONTINUE
  100 DPROB = PROB
      RETURN
      END

      REAL FUNCTION DVALU(PROB, NDF, NPOPS, MBEST, IFAULT)
C
C        ALGORITHM AS 267.2  APPL.STATIST. (1991), VOL.40, NO.3
C
C        Given a probability PROB, function DVALU calculates the
C        standardized difference required between the MBEST population
C        means and the remaining means when applying ranking and
C        selection procedures to NPOPS populations.
C
      INTEGER IFAULT, MBEST, NDF, NPOPS
      REAL PROB
      INTEGER J, JMAX
      REAL A, BEST, DF, DCUT, E1, E2, EDIF, ONE, P1, P2, PCUT1, PCUT2,
     *     POPS, Q1, Q2, TWO, ZERO
      REAL DPROB, DVALU0
      EXTERNAL DPROB, DVALU0
C
      DATA JMAX, ZERO, ONE, TWO, DCUT, PCUT1, PCUT2 / 8, 0.0, 1.0, 2.0,
     *     0.00001, 0.00015, 0.001 /
C
C        Check input arguments.
C
      DVALU = ZERO
      IFAULT = 0
      IF (PROB .LT. 0.80 .OR. PROB .GT. 0.995) IFAULT = 5
      IF (NDF .LE. 3) IFAULT = 20 + IFAULT
      IF (NPOPS .LT. 2 .OR. NPOPS .GT. 100) IFAULT = 300 + IFAULT
      IF (MBEST .LT. 1 .OR. MBEST .GE. NPOPS) IFAULT = 4000 + IFAULT
      IF (IFAULT .NE. 0) GO TO 20
C
C        Obtain initial values.
C
      DF = NDF
      POPS = NPOPS
      BEST = MBEST
      Q2 = DVALU0(PROB, DF, POPS, BEST)
      P2 = DPROB(Q2, NDF, NPOPS, MBEST, IFAULT)
      E2 = P2 - PROB
      A = (ONE - PROB) / (ONE - P2)
      IF (E2 .GT. ZERO) A = PROB / P2
      P1 = PROB - E2 * A
      Q1 = DVALU0(P1, DF, POPS, BEST)
      E1 = DPROB(Q1, NDF, NPOPS, MBEST, IFAULT) - PROB
      DVALU = Q1
C
C        Refine and test approximation.
C
      DO 10 J = 1, JMAX
         IF (ABS(E1) .LT. PCUT1) GO TO 20
         EDIF = E2 - E1
         IF (ABS(EDIF) .LT. DCUT) THEN
            DVALU = (Q1 + Q2) / TWO
            IFAULT = 6
            GO TO 20
         END IF
         DVALU = (Q1 * E2 - Q2 * E1) / EDIF
         IF (ABS(E1) .LT. PCUT2) GO TO 20
         Q2 = Q1
         E2 = E1
         Q1 = DVALU
         E1 = DPROB(Q1, NDF, NPOPS, MBEST, IFAULT) - PROB
   10 CONTINUE
      IFAULT = 6
   20 RETURN
      END
      REAL FUNCTION DVALU0(PROB, DF, POPS, BEST)
C
C        ALGORITHM AS 267.3  APPL.STATIST. (1991), VOL.40, NO.3
C
C        Calculates an initial value for the standardized difference.
C
      REAL ALBEST, ALPOPS, BEST, DF, DFMAX, PROB, POPS, T
      REAL GAUINV
      EXTERNAL GAUINV
      DATA DFMAX / 120.0 /
C
      T = GAUINV(PROB)
      IF (DF .LT. DFMAX) T = T + (T * T * T + T) / DF / 4.0
      ALBEST = LOG(AMIN1(BEST, POPS - BEST))
      ALPOPS = LOG(POPS)
      DVALU0 = (-0.5185 + 0.7927 / DF + (0.09494 - 0.1824 / DF) * T) * T
      DVALU0 = DVALU0 + 1.193 - 0.05360 * ALPOPS - 0.08017 * ALBEST
      DVALU0 = (1.0 + DVALU0 * ALPOPS * LOG(BEST + 1.7183)) * T
      RETURN
      END
