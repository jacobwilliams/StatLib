      SUBROUTINE STEM(K, SDF, SDN, D, T, EPS, LG, ALGAMA, KW, WK, XPT,
     *  VAR, IFAULT)
C
C        ALGORITHM AS220  APPL. STATIST. (1986) VOL. 35, NO. 2
C
C        COMPUTES EXPECTATION AND VARIANCE OF STEIN AND
C        EFRON-MORRIS LIMITED TRANSLATION ESTIMATORS.
C
      REAL SDF, SDN, D, T(K), XPT(K), VAR(K), ALGAMA(LG), WK(KW)
      REAL EPS, EPS2, EPS3, HF, ONE, TWO, P, PM, RPM, DL, BD, C1, C2,
     *  C22, PNC, PP1, PP2
      REAL SJJ, BJ, BJ2, R0, R1, R2, SLL, BV, BL, BL2, QSM, QSM2, XPG,
     *  TMP
      REAL VC, DN, QTRM, EPS4, Q(101)
      INTEGER THR
      DATA ZERO, HF, ONE, TWO /0.0, 0.5, 1.0, 2.0/
      DATA THR /3/
      DATA EPS4 /1.0E-04/
      DATA ITYPE /1/
      LENQ = 101
      IF (K .LT. THR .OR. D .LT. ZERO .OR. EPS .LE. ZERO .OR. (SDF .GT.
     *  ZERO .AND. SDN .LE. ZERO)) GOTO 21
      IFAULT = 0
      DL = D
      P = K
      PM = P - TWO
      RPM = SQRT(PM)
      IF (D .GT. RPM) DL = RPM
      EPS2 = HF * HF * EPS
      F1 = ONE
      F2 = ONE
      IF (SDF .LE. ZERO) GOTO 1
      F1 = SDF / SDN
      F2 = SDF * (TWO + SDF) / (SDN * SDN)
    1 F3 = TWO * F1 + F2
      BD = DL * DL / PM
      C1 = F1 * TWO * RPM * DL
      C22 = F2 * DL * DL
      C2 = C22 * PM
      PNC = ZERO
      DO 2 J1 = 1, K
      PNC = PNC + T(J1) * T(J1)
    2 CONTINUE
      PNC = HF * PNC
C
C        COMPUTE EXPECTATION OF THE NP-COMPONENT ESTIMATOR
C
      DO 17 NP = 1, K
      DO 3 J2 = 1, KW
      WK(J2) = ZERO
    3 CONTINUE
      PP1 = HF * T(NP) * T(NP)
      PP2 = PNC - PP1
      EPS3 = EPS2
      IF (T(NP) .LE. ONE) GOTO 4
      EPS3 = EPS2 / (T(NP) * T(NP))
    4 XPT(NP) = ZERO
      VAR(NP) = ONE + T(NP) * T(NP)
C
C        CORRECTION FOR JAMES-STEIN-EFRON-MORRIS PARTS
C        SET UP INNER POISSON SUM
C
      JJ = 1
      WK(1) = EXP(-PP1)
      SJJ = WK(1)
      JXP = 0
      R0 = ONE
      R1 = ONE
    5 JJ = JJ + 1
      IF (JJ .GT. KW) GOTO 19
      BJ = JJ - 1
      WK(JJ) = WK(JJ - 1) * PP1 / BJ
      SJJ = SJJ + WK(JJ)
      IF (SJJ .GT. ONE) SJJ = ONE
      R2 = R1
      R1 = R0
      R0 = (ONE - SJJ) + WK(JJ)
      IF (JXP .GT. 0) GOTO 6
      IF (R0 * T(NP) * (ONE + HF * C1) .LE. EPS3) JXP = JJ + 1
      IF (JXP .EQ. 0) GOTO 5
    6 IF (R0 * C22 + R1 * (C1 * PP1 + F3) + TWO * F3 * PP1 * R2 .GT.
     *  EPS2) GOTO 5
      LL = 1
C
C        SET UP OUTER POISSON SUM
C
      LP = JJ + 1
      IF (LP + 1 .GT. KW) GOTO 19
      LXP = 0
      WK(LP + 1) = EXP(-PP2)
      SLL = WK(LP + 1)
      BV = PP1 * (C1 + 2 . * F3) + C22 + F3
    7 LL = LL + 1
      IF (LP + LL .GT. KW) GOTO 19
      BL = LL - 1
      WK(LP + LL) = WK(LP + LL - 1) * PP2 / BL
      SLL = SLL + WK(LP + LL)
      IF (SLL .GT. ONE) SLL = ONE
      R0 = ONE - SLL
      IF (LXP .GT. 0) GOTO 8
      IF (R0 * T(NP) * (ONE + HF * C1) .LE. EPS3) LXP = LL + 1
      IF (LXP .EQ. 0) GOTO 7
    8 IF (R0 * BV .GT. EPS2) GOTO 7
      KSTP = K + 2 * JXP + 2 * LXP - 3
      IF (KSTP .GT. LG) GOTO 18
      IF (2 * JXP + 2 .GT. LENQ) GOTO 20
      DO 12 L2P = 1, LXP
      L2 = L2P - 1
      BL2 = L2
      JDF = 2 * JXP
      LDF = K - 1 + 2 * L2
      DO 9 JX = 1, LENQ
      Q(JX) = ZERO
    9 CONTINUE
      QSM = ZERO
      QSM2 = ZERO
      DO 10 J2P = 1, JXP
      J2 = J2P - 1
      BJ2 = J2
      JD1 = 2 * J2 + 3
      JD2 = 2 * J2 + 4
      JSB = K + 2 * J2 + 2 * L2
      XPG = EXP(ALGAMA(JSB) - ALGAMA(JSB + 1) + ALGAMA(2 * J2 + 2)
     *  - ALGAMA(2 * J2 + 3))
      Q(JD1) = -HF * DL * T(NP) * RPM * XPG * WK(J2 + 1) * WK(LP + L2 +
     *  1)
      Q(JD2) = T(NP) * WK(J2 + 1) * WK(LP + L2 + 1) * PM / (P +
     *  TWO * BL2 + TWO * BJ2)
      QSM = QSM - Q(JD1)
      QSM2 = QSM2 + Q(JD1) + Q(JD2)
   10 CONTINUE
      TMP = QSM2
      IF (ABS(DL - RPM) .LE. EPS4) GOTO 11
      TMP = BMIX(BD, JDF, LDF, Q, ITYPE, IFAULT)
      IF (IFAULT .NE. 0) RETURN
   11 XPT(NP) = XPT(NP) - QSM - TMP
   12 CONTINUE
      XPT(NP) = T(NP) + F1 * XPT(NP)
C
C        COMPUTE VARIANCE OF THE NP-COMPONENT ESTIMATOR
C
      VC = ZERO
      LLP = LL + 1
      DO 16 L2P = 1, LLP
      L2 = L2P - 1
      BL2 = L2
      JDF = 2 * JJ
      LDF = K - 1 + 2 * L2
      QSM = C2 * WK(1) * WK(LP + L2 + 1) / (PM + TWO * BL2)
      DO 13 J2 = 1, LENQ
      Q(J2) = ZERO
   13 CONTINUE
      Q(2) = -QSM
      QSM2 = -QSM
      DO 14 J2P = 1, JJ
      J2 = J2P - 1
      BJ2 = J2
      JD1 = 2 * J2 + 3
      JD2 = 2 * J2 + 4
      JSB = K + 2 * J2 + 2 * L2
      DN = JSB
      QTRM = C2 * WK(J2 + 2) * WK(LP + L2 + 1) / DN
      Q(JD2) = -QTRM + WK(J2 + 1) * WK(LP + L2 + 1) * PM * (TWO * BJ2 +
     *  ONE) / DN * (-TWO * F1 + PM * F2 / (DN - TWO))
      XPG = EXP(ALGAMA(JSB) - ALGAMA(JSB + 1) + ALGAMA(2 * J2 + 2)
     *  - ALGAMA(2 * J2 + 1))
      WK(LG + JD1) = C1 * WK(J2 + 1) * WK(LP + L2 + 1) * XPG
      QSM = QSM + QTRM - WK(LG + JD1)
      QSM2 = QSM2 + WK(LG + JD1) + WK(LG + JD2)
   14 CONTINUE
      TMP = QSM2
      IF (ABS(DL - RPM) .LE. EPS4) GOTO 15
      TMP = BMIX(BD, JDF, LDF, Q, ITYPE, IFAULT)
      IF (IFAULT .NE. 0) RETURN
   15 VC = VC + QSM + TMP
   16 CONTINUE
      VAR(NP) = VAR(NP) + VC - XPT(NP) * XPT(NP)
   17 CONTINUE
      RETURN
   18 IFAULT = 6
      RETURN
   19 IFAULT = 7
      RETURN
   20 IFAULT = 8
      RETURN
   21 IFAULT = 9
      RETURN
      END
C
      SUBROUTINE ALGAM(N, ALGAMA)
C
C        COMPUTES ALOG(GAMMA(J/2)) FOR J = 1(1)N, N EVEN.
C        PL2 IS HALF LOG(PI).
C
      DIMENSION ALGAMA(N)
      DATA ZERO, HF, PL2 /0.0, 0.5, 0.572364943/
      NHLF = N / 2
      ALGAMA(1) = PL2
      ALGAMA(2) = ZERO
      DO 1 J1 = 2, NHLF
      B1 = FLOAT(2 * J1 - 3)
      B2 = FLOAT(J1 - 1)
      ALGAMA(2 * J1) = ALGAMA(2 * J1 - 2) + ALOG(B2)
      ALGAMA(2 * J1 - 1) = ALGAMA(2 * J1 - 3) + ALOG(HF * B1)
    1 CONTINUE
      RETURN
      END
