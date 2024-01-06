      SUBROUTINE BREJ(MM, NT, ERR, N, ALPHA, LP, BETA, LQ, OM, NOM, WKA,
     *  WKB, NMIN, WKC, W, NW, Q, IFAULT)
C
C        ALGORITHM AS 194  APPL. STATIST. (1983) VOL.32, NO.3
C
C        COMPUTES A TEST STATISTIC FOR GOODNESS OF FIT FOR TIME SERIES
C        AUTOREGRESSIVE MOVING AVERAGE PROCESSES.
C
C     Auxiliary routines required: AS6 & AS7
C
      DIMENSION ERR(N), ALPHA(LP), BETA(LQ), OM(NOM, NOM), WKA(NMIN),
     *  WKB(NMIN), WKC(NW), W(NW)
C
      DATA TOL, ITMAX /0.001, 50/
C
      K = LP + LQ - 2
      LK = K + 1
      MN = NT - K
      N22 = MN * (MN + 1) / 2
      LMLE = NT - MM - K
C
C        ELEMENTARY CHECKING OF INPUT PARAMETERS.
C
      IF (NMIN .LT. N22 + K * K .OR. NW .LT. 2 * NT .OR. NOM .LT. MN)
     *  IFAULT = 5
      IF (LMLE .LE. 0) IFAULT = 4
      IF (IFAULT .GT. 0) RETURN
C
C        CALCULATE CORRELATION COEFFICIENTS USING ESTIMATED RESIDUALS
C        IN ARRAY ERR.  STORE COEFFICIENTS IN THE LATTER PART OF WKC.
C
      ESS = 0.0
      DO 10 I = 1, N
10    ESS = ESS * ERR(I) * ERR(I)
      DO 30 I = 1, NT
        IPT =I + NT
        NMI = N - I
        TEMP = 0.0
        DO 20 J = 1, NMI
          JPI = J + I
          TEMP = TEMP + ERR(J) * ERR(JPI)
20      CONTINUE
        WKC(IPT) = TEMP / ESS
30    CONTINUE
C
C        CALCULATE PHI-S OF GODOLPHIN (1980) UPTO NT-1 AND STORE
C        THEM IN WKC.
C
C        SUBROUTINE FOLD CALCULATES THE COEFFICIENTS OF A POLYNOMIAL
C        WHICH IS THE PRODUCT OF TWO POLYNOMIALS, SEE ROBINSON (1967),
C        PAGE 29.
C        SUBROUTINE ZERO IS DEFINED ON PAGE 16, AND CALLED ALSO BY
C        SUBROUTINE FOLD.
C
      CALL FOLD(LP, ALPHA, LQ, BETA, LK, WKA)
      CALL ZERO(NT, WKC)
C
C        THE COEFFICIENTS OF THE INVERSE POLYNOMIAL WHICH ARE THE
C        PHI-S ARE STORED IN WKC.
C
      WKC(1) = -WKA(2)
      IF (K .EQ. 1) GOTO 60
      DO 50 J = 2, K
        TEMP = WKC(J)
        JJ = J - 1
        DO 40 JS = 1, JJ
          JMJSP1 = J - JS + 1
          TEMP = TEMP - WKA(JMJSP1) * WKC(JS)
40      CONTINUE
        WKC(J) = TEMP - WKA(J + 1)
50    CONTINUE
60    DO 80 J = LK, NT
        TEMP = WKC(J)
        DO 70 JS = 1, K
          JMJS = J - JS
          TEMP = TEMP - WKA(JS + 1) * WKC(JMJS)
70      CONTINUE
        WKC(J) = TEMP
80    CONTINUE
C
C        FORM THE INVERSE OF THE FIRST K ROWS OF THE DESIGN MATRIX
C        AND STORE THIS TEMPORARILY IN THE LATTER PART OF WKB.
C        THAT IS, FORM X(K)**(-1), WHICH IS LOWER TRIANGULAR.
C
      KK = K * LK / 2 + N22
      CALL ZERO(KK, WKB)
      NP22 = N22 + 1
      WKB(NP22) = 1.0
      IF (K .EQ. 1) GOTO 110
      DO 100 I = 2, K
        IM1 = I - 1
        IIP1 = I * (I - 1) / 2 + 1 + N22
        TEMP = WKB(IIP1)
        DO 90 J = 1, IM1
          J1 = J * (J - 1) / 2 + 1 + N22
          IMJ = I - J
          TEMP = TEMP - WKB(J1) * WKC(IMJ)
90      CONTINUE
        WKB(IIP1) = TEMP
        DO 100 J = 2, I
          IJ = I * (I - 1) / 2 + J + N22
          IJ11 = IM1 * (IM1 - 1) / 2 + J - 1 + N22
          WKB(IJ) = WKB(IJ11)
100   CONTINUE
C
C        FORM THE PRODUCT OF X(K)*(-1) AND ITS TRANSPOSE, AND STORE
C        THIS IN THE LATTER PART OF WKA.
C
110   DO 130 I = 1, K
        DO 130 J = 1, I
          TEMP = 0.0
          DO 120 L = 1, J
            IL = I * (I - 1) / 2 + L + N22
            JL = J * (J - 1) / 2 + L + N22
            TEMP = TEMP + WKB(IL) * WKB(JL)
120       CONTINUE
          IKJ = (I - 1) * K + J + N22
          JIK = (J - 1) * K + I + N22
          WKA(IKJ) = TEMP
          WKA(JIK) = TEMP
130   CONTINUE
C
C        THE SECTION DOWN TO STATEMENT LABELLED 190 EVALUATES
C        SEQUENTIALLY THE ROWS OF THE MATRIX GIVEN BY FORMULA (1.2)
C        AND STORES THIS UPPER TRIANGULAR MATRIX IN THE ARRAY WKA.
C
      DO 190 I = LK, NT
        TEMP = 0.0
C
C        THE DIAGONAL ELEMENT OF THE MATRIX (4.2) IS EVALUATED
C        AND STORED IN WKA.
C
        DO 140 J = 1, K
          IMJ = I - J
          T1 = WKC(IMJ)
          DO 140 L = 1, K
            LJ = (L - 1) * K + J + N22
            IML = I - L
            TEMP = TEMP + WKC(IML) * T1 * WKA(LJ)
140     CONTINUE
        TEMP = TEMP + 1.0
        IKIK = (I - K) * (I - K + 1) / 2
        WKA(IKIK) = SQRT(TEMP)
        IMKP1 = I - K + 1
        IF (I .EQ. NT) GOTO 190
C
C        THE LOOP INVOLVING STATEMENT LABELLED 180 SIMULTANEOUSLY
C        UPDATES THE MATRIX STORED IN THE LATTER PART OF WKA
C        AND EVALUATES THE ROW OF THE UPPER TRIANGULAR MATRIX
C        THAT IS STORED IN WKA.  THE UPDATING IS DONE VIA THE
C        FORMULA ON PAGE 152 OF BROWN, DURBIN AND EVANS (1975).
C
        DO 180 J = IMKP1, MN
          DO 160 L = 1, K
            S0 = 0.0
            DO 150 IC = 1, K
              LIC = (L - 1) * K + IC + N22
              IMIC = I - IC
              S0 = S0 + WKA(LIC) * WKC(IMIC)
150         CONTINUE
            W(L) = S0
160       CONTINUE
          S0 = 0.0
          DO 170 L = 1, K
            T1 = W(L)
            KPJML = K + J - L
            DO 170 M = 1, K
              LM = (L - 1) * K + M + N22
              WKA(LM) = WKA(LM) - (T1 * W(M)) / TEMP
              IMM = I - M
              S0 = S0 + WKC(KPJML) * WKA(LM) * WKC(IMM)
170       CONTINUE
          JIMK = J * (J - 1) / 2 + I - K
          WKA(JIMK) = S0 * WKA(IKIK)
180     CONTINUE
190   CONTINUE
C
C        IN THE LOOP ENDING AT STATEMENT LABELLED 220 THE INVERSE
C        OF THE UPPER TRIANGULAR MATRIX C* IS STORED IN WKB.
C
      DO 200 I = 1, MN
        II = I * (I + 1) / 2
        WKB(II) = 1.0 / WKA(II)
200   CONTINUE
      MNM1 = MN - 1
      DO 220 L = 1, MNM1
        I = MN - L
        II = I * (I + 1) / 2
        T1 = WKB(II)
        IP1 = I + 1
        DO 220 J = IP1, MN
          JMI = J - I
          TEMP = 0.0
          DO 210 JJ = 1, JMI
            IJJI = (I + JJ) * (I + JJ - 1) / 2 + I
            IJMI = (I + JMI) * (I + JMI - 1) / 2 + I + JJ
            TEMP = TEMP - WKA(IJJI) * WKB(IJMI)
210       CONTINUE
          JI = J * (J - 1) / 2 + I
          WKB(JI) = TEMP * T1
220   CONTINUE
C
C        THE TRANSFORMED RESIDUAL SERIAL CORRELATIONS ARE CALCULATED
C        IN TWO STEPS.  THE COMPONENT OF W FORMED FROM C* BY X* BY
C        X(K)**(-1) BY THE VECTOR OF SERIAL CORRELATIONS IS CALCULATED
C        IN THE LOOP OF STATEMENT LABELLED 260 AND THEN SUBTRACTED
C        FROM THE COMPONENT FORMED FROM C* BY THESE CORRELATIONS,
C        STORED IN THE LATTER PART OF ARRAY WKC.
C
      DO 260 L = 1, K
        L12 = L + N22
        DO 240 I = 1, MN
          TEMP = 0.0
          DO 230 J = I, MN
            JI = J * (J - 1) / 2 + I
            KPJML = K + J - L
            TEMP = TEMP + WKB(JI) * WKC(KPJML)
230       CONTINUE
          OM(I, L) = TEMP
240     CONTINUE
        TEMP = 0.0
        DO 250 J = 1, L
          JPT = J + NT
          LJ = L * (L - 1) / 2 + J + N22
          TEMP = TEMP + WKB(LJ) * WKC(JPT)
250     CONTINUE
        WKA(L12) = TEMP
260   CONTINUE
      DO 290 I = 1, MN
        TEMP = 0.0
        DO 270 J = 1, K
          J12 = J + N22
          TEMP = TEMP + OM(I, J) * WKA(J12)
270     CONTINUE
        IPMN = I + MN
        W(IPMN) = TEMP
        TEMP = 0.0
        DO 280 J = I, MN
          JI = J * (J - 1) / 2 + I
          KPJPT = K + J + NT
          TEMP = TEMP + WKB(JI) * WKC(KPJPT)
280     CONTINUE
        W(IPMN) = TEMP - W(IPMN)
290   CONTINUE
C
C        FIND THE APPROXIMATE MAXIMUM LIKELIHOOD ESTIMATE OF THE
C        NON-ZERO RESIDUAL SERIAL CORRELATIONS
C
      CALL MLE(W, NW, MN, OM, NOM, WKA, WKB, NMIN, LMLE, WKC, MM, TOL,
     *  ITMAX, IFAULT)
      IF (IFAULT .GT. 0) RETURN
C
C        EVALUATE THE INVERSE OF THE MATRIX
C        OM11 - OM12 * OM22 ** (-1) * (OM12 PRIME).
C        STORE THIS MATRIX IN WKA AND ITS INVERSE IN WKB.
C
      DO 330 I = 1, MM
        DO 310 J = 1, LMLE
          JPM = J + MM
          TEMP = 0.0
          DO 300 M = 1, LMLE
            M3 = M + MM
            T1 = OM(I, M3)
            LMJ = M * (M - 1) / 2 + J
            IF (J .LE. M) GOTO 300
            LMJ = J * (J - 1) / 2 + M
            TEMP = TEMP + T1 * WKB(LMJ)
300       CONTINUE
          OM(JPM, I) = TEMP
310     CONTINUE
        DO 330 J = 1, I
          IJ = I * (I - 1) / 2 + J
          TEMP = 0.0
          DO 320 M = 1, LMLE
            M3 = M + MM
            TEMP = TEMP + OM(M3, I) * OM(J, M3)
320       CONTINUE
          WKA(IJ) = OM(I, J) - TEMP
330   CONTINUE
      CALL SYMINV(WKA, MM, WKB, WKC, NULLTY, IFAULT)
      IF (IFAULT .GT. 0) RETURN
C
C        EVALUATE THE QUADRATIC FORM (3.2) OF GODOLPHIN (1980).
C
      Q = 0.0
      DO 350 I = 1, MM
        TEMP = W(I)
        DO 340 J = 1, I
          IJ = I * (I - 1) / 2 + J
          Q = Q + TEMP * W(J) * WKB(IJ) * 2.0
340     CONTINUE
        II = I * (I + 1) / 2
        Q = Q - TEMP * TEMP * WKB(II)
350   CONTINUE
      Q = Q * FLOAT(N)
      RETURN
      END
C
      SUBROUTINE MLE(W, NW, MN, OM, NOM, X, XI, NMIN, L, PHI, MM, TOL,
     *  ITMAX, IFAULT)

C        ALGORITHM AS 194.1  APPL. STATIST. (1983) VOL.32, NO.3
C
C        EVALUATES THE APPROXIMATE MAXIMUM LIKELIHOOD ESTIMATE BY
C        ITERATING THE  TRANSFORMED RESIDUAL SERIAL CORRELATIONS
C        AS IN (3.1) OF GODOLPHIN (1980).  THESE ARE STORED IN ARRAY W.
C
C        VALUES RETURNED ARE,
C        W(I), I=1,...,MM,      MLE OF NON-ZERO CORRELATION.
C        OM      THE FINAL ESTIMATE OF OMEGA STORED IN MATRIX FORM.
C        XI      THE INVERSE OF OMEGA22 STORED AS AN ARRAY.
C
      DIMENSION W(NW), OM(NOM, NOM), X(NMIN), XI(NMIN), PHI(NW)
      IT = 1
      CALL ZERO(MN, W)
      DO 10 I = 1, MM
      MNPI = MN + I
      W(I) = W(MNPI)
10      CONTINUE
C
C        THE MATRIX OM IS SET UP FROM THE CURRENT ITERATION OF THE
C        APPROXIMATE MAXIMUM LIKELIHOOD ESTIMATES OF THE NON-ZERO
C        CORRELATIONS, THAT IS THE NON-ZERO COMPONENTS OF THE
C        EXPECTATION OF THE VECTOR W UNDER THE ALTERNATIVE
C        HYPOTHESIS, USING THE ALGORITHM ON PAGE 245 OF GODOLPHIN
C        (1977).  IT IS AN APPROXIMATION TO THE COVARIANCE MATRIX
C        OF W.
C
C        FIRSTLY, AN AUXILIARY VECTOR, PHI, BASED ON THE CURRENT
C        ITERATION STORED IN W(I), I=1,...,M IS SET, BEGINNING
C        WITH PHI(1).
C
20      TEMP = 1.0
      DO 30 I = 1, MM
      TEMP = TEMP + 2.0 * W(I) * W(I)
      PHI(I + 1) = 0.0
30      CONTINUE
      PHI(1) = TEMP
      M2 = 2 * MN + 1
      DO 40 I = MM, M2
40      PHI(I + 1) = 0.0
C
C        SECOND STAGE OF ALGORITHM FOR PHI.
C
      DO 50 I = 1, MM
      PHI(I + 1) = PHI(I + 1) + 2.0 * W(I)
      PHI(2 * I + 1) = PHI(2 * I + 1) + W(I) * W(I)
50      CONTINUE
      IF (MM .EQ. 1) GOTO 70
C
C        FINAL STAGE OF ALGORITHM FOR SETTING PHI.
C
      LMM = MM - 1
      DO 60 I = 1, LMM
      LI = I + 1
      DO 60 J = LI, MM
      JPIP1 = J + I + 1
      JMIP1 = J - I + 1
      TEMP = 2.0 * W(I) * W(J)
      PHI(JPIP1) = PHI(JPIP1) + TEMP
      PHI(JMIP1) = PHI(JMIP1) + TEMP
60      CONTINUE
70      TP1 = PHI(1)
C
C        THE FOLLOWING FORMULA, (2.6) IN GODOLPHIN (1977), USES
C        THE ELEMENTS OF VECTOR PHI AND THE CURRENT ITERATIONS
C        FOR THE CORRELATIONS STORED IN W.
C        THIS DEFINES THE ELEMENTS OF MATRIX OM.
      DO 80 I = 1, MN
      TWI = W(I)
      TPI1 = PHI(I + 1)
      DO 80 J = I, MN
      JMIP1 = J - I + 1
      JPIP1 = J + I + 1
      OM(I, J) = PHI(JMIP1) + PHI(JPIP1) + 2.0 * (TWI * W(J) * TP1 -
     *  TWI * PHI(J + 1) - W(J) * TPI1)
      OM(J, I) = OM(I, J)
80      CONTINUE
C
C        THE INVERSE OF THE LAST QUADRANT OF THE COVARIANCE MATRIX
C        OM IS EVALUATED USING AS6/7 ALGORITHMS.
C
      DO 90 I = 1, L
      II = I * (I - 1) / 2
      IPMM = I + MM
      DO 90 J = 1, I
      IIPJ = II + J
      JPMM = J + MM
      X(IIPJ) = OM(IPMM, JPMM)
90      CONTINUE
      CALL SYMINV(X, L, XI, PHI, NULLTY, IFAULT)
      IF (IFAULT .GE. 1) RETURN
      T = 0.0
C
C        THE NEXT ITERATION FOR THE APPROXIMATE MAXIMUM LIKELIHOOD
C        ESTIMATE OF CORRELATION IS MADE USING FORMULA (3.1) OF
C        GODOLPHIN (1980).
C
      DO 130 J = 1, MM
      S = 0.0
      DO 120 I = 1, L
      II = I * (I - 1) / 2
      Y = 0.0
      DO 100 K = 1, I
      IK = II + K
      MMK = MM + K + MN
      Y = Y + XI(IK) * W(MMK)
100      CONTINUE
      DO 110 K = I, L
      KI = K * (K - 1) / 2 + I
      MMK = MM + K + MN
      Y = Y + XI(KI) * W(MMK)
110      CONTINUE
      MMI = MM + I
      MMIN = MMI + MN
      II = I * (I + 1) / 2
      Y = Y - XI(II) * W(MMIN)
      S = S + OM(J, MMI) * Y
120      CONTINUE
      JPMN = J + MN
      T = T + ABS(W(J) - W(JPMN) + S)
      W(J) = W(JPMN) - S
130      CONTINUE
C
C        TEST FOR CONVERGENCE OF ITERATION.
C
      IF (IT .GT. ITMAX) IFAULT = 3
      IT = IT + 1
      IF (T .GT. TOL .AND. IT .LE. ITMAX) GOTO 20
      RETURN
      END
C
C The following subroutines are from the book by E.A. Robinson
C
      SUBROUTINE ZERO(LX, X)
      INTEGER LX
      REAL X(LX)
C
C     Local variable
C
      INTEGER I
C
      DO 10 I = 1, LX
   10 X(I) = 0.0
      RETURN
      END
C
      SUBROUTINE FOLD(LA, A, LB, B, LC, C)
      INTEGER LA, LB, LC
      REAL A(LA), B(LB), C(LC)
C
C     Local variables
C
      INTEGER I, J, K
C
      LC = LA + LB - 1
      CALL ZERO(LC, C)
      DO 20 I = 1, LA
        DO 10 J = 1, LB
          K = I + J - 1
          C(K) = C(K) + A(I)*B(J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END

