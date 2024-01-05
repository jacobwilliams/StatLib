      SUBROUTINE RMCONF(EST, CONF, M, STEPK, MODSTP, LIML, SBOUND, LIMU,
     *                  BBOUND, LSTART, SLOLIM, SHILIM, ISTART, IRSTRT,
     *                  INDIC, CLIMLO, CLIMHI, CLLOSE, CLHISE, MLO, MHI,
     *                  PTHPLT, IFAULT)
C
C        ALGORITHM AS 259  APPL. STATIST. (1990) VOL. 39, NO. 3
C
C        Calculates Monte Carlo confidence intervals using the Robbins-
C        Monro process.  The intervals are asymptotically exact as the
C        number of steps tends to infinity under general conditions.
C
      CHARACTER * 1 PTHPLT(76, 105)
      INTEGER M, ISTART, INDIC, MLO, MHI, IFAULT
      LOGICAL MODSTP, LIML, LIMU, LSTART, IRSTRT
      REAL EST, CONF, STEPK, SBOUND, BBOUND, SLOLIM, SHILIM, CLIMLO,
     *     CLIMHI, CLLOSE, CLHISE
      CHARACTER *  1 DIGIT(0:9), SYMB
      CHARACTER * 16 MIDPLT
      CHARACTER * 28 HIPLOT, LOPLOT
      CHARACTER * 30 MIDLIN
      INTEGER IX(76), MULT(2)
      REAL CVHI(20), CVLO(20), SEARCH(2, 2000)
      INTEGER I, IERR, II, IIM, IPTH, IR, ISTRT, IY, J, JJ, JX, K,
     *        KMULT, KOUNT, L, LIMIT, LL, LR, LX, LY, LYY, M0, M1, M2,
     *        MID, MM, NPT
      LOGICAL IRSTAR, ISTEP
      REAL AL2, ALPHA, BB, BEGPOS, BIG, BOOTES, CONF2, D1, D2, DIV,
     *     ENDPOS, FIFTY, HALF, HIPOS, HUN, ONE, ONEPSM, P95, POSN,
     *     PTHREE, QUART, ROOT2, SB, SEPOSN, SEVFIV, SMALL, STEP,
     *     STEP2, STEPC, STEPHI, TEN, THREP5, TRT2PI, TWO, X, YMAX,
     *     YMIN, YSCALE, Z, Z1, ZA, ZERO
      REAL PPND, SIMU
      EXTERNAL PPND, SIMU
C
      DATA DIGIT / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
      DATA ZERO, QUART, PTHREE, HALF / 0.0, 0.25, 0.3, 0.5 /
      DATA ONE, ROOT2, TWO, THREP5, TEN / 1.0, 1.414214, 2.0, 3.5,
     *     10.0 /
      DATA FIFTY, SEVFIV, P95, HUN / 50.0, 75.0, 95.0, 100.0 /
      DATA SMALL, ONEPSM, BIG / 1.0E-6, 1.00001, 1.0E20 /
      DATA HIPLOT / 'Upper limit, number of steps' /
      DATA LOPLOT / 'Lower limit, number of steps' /
      DATA MIDLIN / '------------------------------' /
      DATA MIDPLT / ' Point estimate ' /
C
C        TRT2PI is 2*SQRT(2*PI); Z1 is z st P(Z > z)=0.005,  Z is N(0,1)
C
      DATA TRT2PI / 5.013257 / , Z1 / 2.576 /
C
      IFAULT = 0
      MULT(1) = 1
      MULT(2) = 1
      IF (CONF .GE. HUN .OR. CONF .LT. ZERO) THEN
         IFAULT = 1
         RETURN
      END IF
      CONF2 = CONF
      IF (CONF .LT. SMALL) CONF2 = P95
      ALPHA = (ONE - CONF2 / HUN) / TWO
      MM = M
      IF (M .EQ. 0) MM = 1000
      M2 = MM / 2
      IF (CONF2 .GE. FIFTY) THEN
         ISTEP = .FALSE.
         AL2 = ALPHA
      ELSE
         ISTEP = .TRUE.
         AL2 = QUART
      END IF
      IF (AL2 * MM .LT. TEN / ONEPSM) THEN
         IFAULT = 5
         RETURN
      END IF
      IF (STEPK .LT. SMALL) THEN
C
C          PPND cannot fail since IFAULT=1 will have already
C          occurred if the conditions for IERR > 0 arise
C
         ZA = PPND(ONE - AL2, IERR)
         STEPC = TRT2PI * EXP(ZA ** 2 / TWO) / ZA
      ELSE
         STEPC = STEPK
      END IF
      SB = SBOUND
      BB = BBOUND
      IF (LIML) SB = -BIG
      IF (LIMU) BB = BIG
      M1 = INT(TWO / ALPHA - HALF)
      IF (M1 .LT. 39) M1 = 39
      IF (ISTART .EQ. 0) THEN
         ISTRT = INT(HALF + PTHREE * M1)
         IF (ISTRT .GT. 50) ISTRT = 50
      ELSE
         ISTRT = ISTART
      END IF
      IRSTAR = IRSTRT
      IF (CONF2 .LT. SEVFIV .AND. IRSTRT) THEN
         IFAULT = 4
         IRSTAR = .FALSE.
      END IF
      KOUNT = 0
      IF (LSTART .OR. ISTEP) THEN
C
C        Set starting points for search using the percentile method
C
         IR = INT((M1 + 1) * AL2 + HALF)
         DO 10 J = 1, IR
            CVHI(J) = -BIG
            CVLO(J) = BIG
   10    CONTINUE
C
         DO 70 I = 1, M1
C
C        Generate a bootstrap estimate, assuming the true value of
C        the parameter is EST
C
            BOOTES = SIMU(EST)
            IF (BOOTES .LT. SB / ONEPSM) IFAULT = 2
            IF (BOOTES .GT. BB * ONEPSM) IFAULT = 3
            M0 = MIN(I, IR)
C
C        Extract the IR largest bootstrap estimates
C
            DO 30 J = 1, M0
               IF (CVHI(J) .LT. BOOTES) THEN
                  IF (J .LT. M0) THEN
                     DO 20 JJ = M0, J + 1, -1
                        CVHI(JJ) = CVHI(JJ - 1)
   20                CONTINUE
                  END IF
                  CVHI(J) = BOOTES
                  GO TO 40
               END IF
   30       CONTINUE
C
C        Similarly for the lower limit
C
   40       DO 60 J = 1, M0
               IF (CVLO(J) .GT. BOOTES) THEN
                  IF (J .LT. M0) THEN
                     DO 50 JJ = M0, J + 1, -1
                        CVLO(JJ) = CVLO(JJ - 1)
   50                CONTINUE
                  END IF
                  CVLO(J) = BOOTES
                  GO TO 70
               END IF
   60       CONTINUE
   70    CONTINUE
         IF (ISTEP) THEN
            STEPHI = (CVHI(IR) - EST) * STEPC
            STEP = (CVLO(IR) - EST) * STEPC
            CLIMHI = EST + CONF2 * (CVHI(IR) - EST) / FIFTY
            CLIMLO = EST - CONF2 * (EST - CVLO(IR)) / FIFTY
         ELSE
C
C        The starting point for the upper limit will be the IRth
C        largest of the M1 bootstrap estimates, and that for the
C        lower limit the IRth smallest.  For CONF2 > 87.5%, IR=2
C
            CLIMHI = CVHI(IR)
            CLIMLO = CVLO(IR)
         END IF
      END IF
      IF ( .NOT. LSTART) THEN
C
C        Starting points supplied by the user
C
         CLIMHI = SHILIM
         CLIMLO = SLOLIM
      END IF
      IF (IFAULT .GT. 0 .AND. IFAULT .NE. 4) RETURN
C
C        Start by searching for the lower limit
C
      POSN = CLIMLO
      SEARCH(1, 1) = POSN
      BEGPOS = POSN
      IF (IRSTAR) D1 = ABS(EST - BEGPOS)
      LIMIT = 1
   80 I = ISTRT + 1
      II = 1
C
C        Generate a bootstrap estimate assuming the true value of the
C        parameter is POSN, the current position of the search
C
   90 BOOTES = SIMU(POSN)
      IF (BOOTES .LT. SB / ONEPSM) THEN
         IFAULT = 2
         RETURN
      ELSE IF (BOOTES .GT. BB * ONEPSM) THEN
         IFAULT = 3
         RETURN
      END IF
      IF (LIMIT .EQ. 1) THEN
         X = EST - BOOTES
      ELSE
         X = BOOTES - EST
      END IF
      IF (X .LT. SMALL * ABS(EST)) THEN
         X = ONE
      ELSE
         X = ZERO
      END IF
      IF ( .NOT. ISTEP) THEN
C
C        Calculate the current estimate of the steplength constant
C
         STEP = (POSN - EST) * STEPC
         IF (MODSTP) THEN
            STEP2 = BB - POSN
            IF (LIMIT .EQ. 1) STEP2 = POSN - SB
            STEP2 = (STEP2 / ALPHA) + STEP * FLOAT(ISTRT) / FLOAT(I)
            IF (STEP2 .LT. STEP) STEP = STEP2
         END IF
      END IF
C
C        Step to the next position in the search
C
      POSN = POSN + STEP * (X - ALPHA) / FLOAT(I)
C
C        Constrain the search to lie within the bounds to the
C        parameter space
C
      IF (POSN .LT. SB) POSN = SB
      IF (POSN .GT. BB) POSN = BB
      I = I + 1
      II = II + 1
      IIM = II / MULT(LIMIT)
      IF (II .EQ. MULT(LIMIT) * IIM) SEARCH(LIMIT, IIM) = POSN
C
C        If array SEARCH is full, reduce the number of points stored
C        by a factor of 10
      IF (IIM .EQ. 200) THEN
         MULT(LIMIT) = MULT(LIMIT) * 10
         DO 100 J = 1, 200
            SEARCH(LIMIT, J) = SEARCH(LIMIT, J * 10)
  100    CONTINUE
      END IF
C
C        Restart search if distance between POSN and EST has halved
C        or doubled and IRSTAR=.TRUE.
C
      IF (IRSTAR) THEN
         D2 = ABS(EST - POSN)
         IF (D1 .GT. TWO * D2 .OR. D1 .LT. HALF * D2) THEN
            I = ISTRT
            BEGPOS = POSN
            D1 = D2
         END IF
      END IF
C
C        Go to the next step if more steps are required
      IF (I .LT. MM + ISTRT .AND. INDIC .EQ. 0) GO TO 90
      IF (I .LT. M2 + ISTRT) GO TO 90
C
C        End current search
C
      IF (INDIC .GT. 0) THEN
         SEPOSN = ABS(STEP) * SQRT(ALPHA * (ONE - ALPHA) / FLOAT(I - 1))
         I = ISTRT
         BEGPOS = POSN
         IF (IRSTAR) D1 = D2
         KOUNT = KOUNT + 1
C
C        Test whether endpoint of the current search differs
C        significantly from that of previous search
C
         IF (KOUNT .GT. 1 .AND. INDIC .GT. 1) THEN
            IF (SEPOSN .LT. SMALL * POSN) THEN
               Z = ZERO
            ELSE
               Z = ABS(POSN - ENDPOS) / (ROOT2 * SEPOSN)
            END IF
            IF (Z .GT. Z1) THEN
               IF (INDIC .EQ. 2) THEN
                  KOUNT = 4
               ELSE
                  KOUNT = KOUNT - 1
               END IF
            END IF
         END IF
C
C        Is a further search required?
C
         IF (KOUNT .EQ. 1) THEN
            ENDPOS = POSN
            GO TO 90
         ELSE IF (KOUNT .NE. 4) THEN
C
C        Average the two end points selected to obtain the estimated
C        confidence limit.  If current search is for the lower limit,
C        start the search for the upper limit
C
            HIPOS = CLIMHI
            CLIMHI = (POSN + ENDPOS) / TWO
            CLHISE = SEPOSN / ROOT2
            MHI = II
            IF (LIMIT .EQ. 1) THEN
               CLIMLO = CLIMHI
               CLLOSE = CLHISE
               MLO = MHI
               LIMIT = 2
               KOUNT = 0
               POSN = HIPOS
               SEARCH(2, 1) = POSN
               BEGPOS = POSN
               IF (IRSTAR) D1 = ABS(EST - BEGPOS)
               IF (ISTEP) STEP = STEPHI
               GO TO 80
            END IF
         END IF
      END IF
C
C        Use the final endpoint as the estimate of the confidence
C        limit.  If current search is for the lower limit, start the
C        search for the upper limit
C
      IF (INDIC .EQ. 0 .OR. KOUNT .EQ. 4) THEN
         IF (INDIC .EQ. 0) SEPOSN = ABS(STEP) *
     *                              SQRT(ALPHA * (ONE - ALPHA) /
     *                              FLOAT(I - 1))
         IF (LIMIT .EQ. 1) THEN
            CLIMLO = POSN
            CLLOSE = SEPOSN
            MLO = II
            LIMIT = 2
            KOUNT = 0
            POSN = CLIMHI
            SEARCH(2, 1) = POSN
            BEGPOS = POSN
            IF (IRSTAR) D1 = ABS(EST - BEGPOS)
            IF (ISTEP) STEP = STEPHI
            GO TO 80
         ELSE
            CLIMHI = POSN
            CLHISE = SEPOSN
            MHI = II
         END IF
      END IF
C
C        Enter the path of the search for each limit into the array
C        PTHPLT so that it may be checked by the user for convergence
C
      DO 120 K = 1, 105
         DO 110 J = 1, 76
            PTHPLT(J, K) = ' '
  110    CONTINUE
  120 CONTINUE
      DO 130 J = 1, 28
         PTHPLT(J + 20, 1) = HIPLOT(J:J)
         PTHPLT(J + 20, 105) = LOPLOT(J:J)
  130 CONTINUE
      YMAX = -BIG
      YMIN = BIG
      DO 150 I = 1, 2
         IF (I .EQ. 1) THEN
            NPT = MLO / MULT(1)
         ELSE
            NPT = MHI / MULT(2)
         END IF
         DO 140 J = 1, NPT
            IF (SEARCH(I, J) .GT. YMAX) YMAX = SEARCH(I, J)
            IF (SEARCH(I, J) .LT. YMIN) YMIN = SEARCH(I, J)
  140    CONTINUE
  150 CONTINUE
      DIV = YMAX - YMIN
      IF (DIV .LE. ZERO) DIV = ONE
      MID = HUN * (YMAX - EST) / DIV + THREP5
      YSCALE = HUN / DIV
      DO 160 J = 1, 30
         PTHPLT(J, MID) = MIDLIN(J:J)
         PTHPLT(J + 46, MID) = MIDLIN(J:J)
         IF (J .LE. 16) PTHPLT(J + 30, MID) = MIDPLT(J:J)
  160 CONTINUE
      DO 210 I = 1, 2
         IF (I .EQ. 1) THEN
            NPT = MLO / MULT(1)
            LYY = YSCALE * (CLIMLO - YMIN) + HALF
            SYMB = '<'
         ELSE
            NPT = MHI / MULT(2)
            LYY = YSCALE * (CLIMHI - YMIN) + HALF
            SYMB = '>'
         END IF
         DO 170 J = 1, 76
            LX = (NPT - 1) * (J - 1) / SEVFIV + ONE + HALF
            IX(J) = LX
            LY = YSCALE * (SEARCH(I, LX) - YMIN) + HALF
            PTHPLT(J, 103 - LY) = SYMB
            IF (PTHPLT(J, 103 - LYY) .EQ. ' ')
     *          PTHPLT(J, 103 - LYY) = '.'
  170    CONTINUE
         DO 190 J = 1, 15
            JX = IX(5 * J)
            DO 180 L = 0, 4
               LL = 10 ** L
               IF (JX .GE. LL) THEN
                  LR = JX / LL - 10 * (JX / (LL * 10))
                  PTHPLT(5 * J - L, 2 + (2 - I) * 102) = DIGIT(LR)
               END IF
  180       CONTINUE
  190    CONTINUE
         IF (MULT(I) .GT. 1) THEN
            KMULT = LOG10(FLOAT(MULT(I))) + HALF
            IY = 1 + (2 - I) * 104
            PTHPLT(50, IY) = '/'
            PTHPLT(52, IY) = '1'
            IPTH = 53
            DO 200 K = 1, KMULT
               PTHPLT(IPTH, IY) = '0'
               IPTH = IPTH + 1
  200       CONTINUE
         END IF
  210 CONTINUE
      RETURN
      END
