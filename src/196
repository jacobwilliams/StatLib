C
C*********************************************************************
C This file was edited by bill bardsley on 18/08/2002 to correct four
C serious errors. To locate these, search for bill.bardsley@man.ac.uk.
C Unless the errors are corrected this code returns incorrect results.
C*********************************************************************
C
      SUBROUTINE LOGCCH(NS, NCA, NCT, NIMAX, NMAX, NMAX1, NVMAX, NVMAX1,
     *  NV, Z, IVAR, COVI, CNTR, W, WB, WDB, WD2B, U, INS, DB, D2B, DL,
     *  B, COV, CHI2, ST, IFAULT)
C
C        ALGORITHM AS 196 APPL. STATIST. (1984) VOL.33, NO.1
C
C        LOGISTIC ANALYSIS OF CASE-CONTROL STUDIES
C
c     *** WARNING  This file has been input using a scanner and may
c                  contain errors.
c
      INTEGER NCA(NS), NCT(NS), IVAR(NVMAX), INS(NS)
      REAL Z(NVMAX, NIMAX), COVI(NVMAX1), CNTR(NVMAX, NS), W(NVMAX),
     *  WB(NMAX1), WDB(NVMAX, NMAX1), WD2B(NVMAX1, NMAX1), U(NMAX),
     *  DB(NVMAX), D2B(NVMAX1), DL(NVMAX), B(NVMAX), COV(NVMAX1)
C
      REAL CHI2, EPS, ZERO, ONE, TWO, BMN, CONST, C1, FM, RLIK, RLIKP,
     *  RLIKS, ST, T, TTT
C
      REAL ABS, ALGFAC, ALOG, FLOAT
C
      DATA MAXIT /20/, EPS /0.00001/, ZERO /0.0/, ONE /1.0/, TWO /2.0/
C
C        INITIAL SETTINGS
C
      RLIKP = ONE
      IFAULT = 0
      ITS = 0
      IF (NV .GT. NVMAX) GOTO 21
      IF (NVMAX * (NVMAX + 1) / 2 .GT. NVMAX1) GOTO 21
      IF (NMAX + 1 .GT. NMAX1) GOTO 21
      INS(1) = 0
      IF (NCA(1) + NCT(1) .GT. NMAX .OR. NCA(1) + NCT(1) .GT. NIMAX)
     *  GOTO 21
      IF (NS .EQ. 1) GOTO 2
      DO 1 I = 2, NS
      IF (NCA(I) + NCT(I) .GT. NMAX) GOTO 21
      I1 = I - 1
      INS(I) = NCA(I1) + NCT(I1) + INS(I1)
    1 CONTINUE
      IF (INS(NS) + NCA(NS) + NCT(NS) .GT. NIMAX) GOTO 21
C
C        CENTRE THE INDEPENDENT VARIABLES ABOUT THE MEAN OF THE
C        COVARIATES FOR THE CASES (C.F. S.HOWARDS COMMENT TO COX (1972))
C
    2 DO 6 I = 1, NS
      IF (NCA(I) * NCT(I) .EQ. 0) GOTO 6
      M = NCA(I)
      N = M + NCT(I)
      FM = 1.0 / FLOAT(M)
      I1 = INS(I)
      DO 5 K = 1, NV
      CNTR(K, I) = ZERO
      K1 = IVAR(K)
      J1 = I1
      DO 3 J = 1, M
      J1 = J1 + 1
      CNTR(K, I) = CNTR(K, I) + Z(K1, J1)
    3 CONTINUE
      CNTR(K, I) = CNTR(K, I) * FM
      J1 = I1
      DO 4 J = 1, N
      J1 = J1 + 1
      Z(K1, J1) = Z(K1, J1) - CNTR(K, I)
    4 CONTINUE
    5 CONTINUE
    6 CONTINUE
    7 ITS = ITS + 1
      IF (ITS .GT. MAXIT) GOTO 23
      RLIK = ZERO
      K = 0
      DO 8 J = 1, NV
      DL(J) = ZERO
      DO 8 JJ = 1, J
      K = K + 1
      COVI(K) = ZERO
    8 CONTINUE
C
C        LOOP THROUGH STRATA
C
      DO 14 I = 1, NS
      IF (NCA(I) * NCT(I) .EQ. 0) GOTO 14
      M = NCA(I)
      N = M + NCT(I)
      NID = INS(I)
C
C        FIND U(J)=EXP(Z(J)*BETA) FOR EACH INDIVIDUAL J IN THE STRATA.
C        ALSO, FIHD TTT, THE TOTAL OF THE (Z(J)*BETA)S.
C
      TTT = ZERO
      DO 10 J = 1, N
      J1 = J + NID
      T = ZERO
      DO 9 K = 1, NV
      L = IVAR(K)
      T = T + B(K) * Z(L, J1)
    9 CONTINUE
      TTT = TTT + T
      U(J) = EXP(T)
   10 CONTINUE
C
C        CALC LATE THE CONSTANT CONST=
C        ((N(C)M) * EXP(M*BETA*XBAR))**(1/M)
C        WHERE XBAR IS THE MEAN OF THE COVARIATES OVER THE CASES AND
C        CONTROLS, AHD DIVIDE EACH U(J) BY THIS CONSTANT.
C        THIS KEEPS THE SUMS CALCULATED BY SUBROUTINE HOWARD FROM
C        BECOMING TOO LARGE IN ABSOLUTE VALUE.
C        NOTE - ALGFAC(X)=LN((X)FACTORIAL)
C
      C1 = ALGFAC(N) - (ALGFAC(M) + ALGFAC(N - M)) + TTT * FLOAT(M) /
     *   FLOAT(N)
      CONST = EXP(-C1 / FLOAT(M))
      DO 11 J = 1, N
   11 U(J) = U(J) * CONST
C
C        CALL TO HOWARD TO CALCULATE SUM(EXP(S(L)*BETA)) OVER ALL
C        COMBINATIONS OF N LABELS TAKEN M AT A TIME, AND ITS FIRST AND
C        SECOND DERIVATIVES WITH RESPECT TO BETA.
C        NOTE - THE VALUE FOR THE SUM AND ITS DERIVATIVES RETURNED BY
C        HOWARD ARE THE TRUE VALUES DIVIDED BY THE CONSTANT CONST**M.
C        THEREFORE RLIK IS CORRECTED FOR THIS CONSTANT SO THAT THE
C        LIKELIHOOD RATIO TEST STATISTIC WILL BE CORRECT.
C        NOTE - SC=BETA*(SUM(Z(J)) OVER THE CASES) = 0 BY DEFINITION.
C
      CALL HOWARD(M, N, U, Z, NID, IVAR, NVMAX, NIMAX, NMAX, NVMAX1, NV,
     *   NMAX1, WB, WDB, WD2B, DB, D2B, BMN)
      RLIK = RLIK - ALOG(BMN) - C1
      L = 0
      IR = 1
      IS = 0
C
C        CALCULATE THE CUMULATIVE SCORE UP TO THIS STRATUM.
C
      DO 13 K = 1, NV
      DL(K) = DL(K) - DB(K) / BMN
      DO 13 KK = 1, K
      L = L + 1
      IS = IS + 1
      IF (IS .LE. IR) GOTO 12
      IR = IR + 1
      IS = 1
C
C       CALCULATE THE CUMULATIVE INFORMATION UP TO THIS STRATUM.
C
C error 1 : bill.bardsley@man.ac.uk corrected the next line 18/08/2002
C  12 COVI(L) = COVI(L) + D2B(L) / BMN - DB(IR)
   12 COVI(L) = COVI(L) + D2B(L) / BMN - DB(IR) * DB(IS) / (BMN * MBN)
   13 CONTINUE
   14 CONTINUE
      IF (ITS .EQ. 1) RLIKS = RLIK
C
C        CALCULATE THE INVERSE OF THE INFORMATION MATRIX
C
C error 2: bill.bardsley@man.ac.uk corrected the next line 18/08/2002
C      CALL SYMINV(COVI, NV, COV, W, NULLTY, IFAULT, NVMAX1)
       CALL SYMINV(COVI, NV, NVMAX1, COV, W, NULLTY, IFAULT)
      IF (IFAULT .NE. 0) GOTO 22
C
C        CALCULATE NEW PARAMETER ESTIMATES
C
      DO 17 I = 1, NV
      W(I) = ZERO
      I2 = I * (I - 1) / 2
      DO 15 J = 1, I
      K = I2 + J
      W(I) = W(I) + DL(J) * COV(K)
   15 CONTINUE
      I1 = I + 1
      IF (I1 .GT. NV) GOTO 17
      DO 16 K = I1, NV
      J = K * (K - 1) / 2 + I
      W(I) = W(I) + DL(K) * COV(J)
   16 CONTINUE
   17 CONTINUE
      DO 18 I = 1, NV
   18 B(I) = B(I) + W(I)
      IF (ITS .NE. 1) GOTO 20
C
C        CALCULATE THE TEST SCORE
C
      ST = ZERO
      DO 19 I = 1, NV
   19 ST = ST + W(I) * DL(I)
C
C        TEST FOR CONVERGENCE
C
   20 RLIK = RLIK - RLIKS
      IF (ABS(RLIKP - RLIK) .LE. EPS) GOTO 24
      RLIKP = RLIK
      GOTO 7
   21 IFAULT = 1
      GOTO 25
   22 IFAULT = 2
      GOTO 25
   23 IFAULT = 3
      GOTO 25
   24 CHI2 = TWO * RLIK
C
C        RETURN MATRIX Z TO THE FORM IT WAS IN WHEN LOGCCH WAS CALLED
C
   25 DO 27 I = 1, NS
      IF (NCA(I) * NCT(I) .EQ. 0) GOTO 27
      N = NCA(I) + NCT(I)
      I1 = INS(I)
      DO 26 K = 1, NV
      K1 = IVAR(K)
      J1 = I1
      DO 26 J = 1, N
      J1 = J1 + 1
      Z(K1, J1) = Z(K1, J1) + CNTR(K, I)
   26 CONTINUE
   27 CONTINUE
      RETURN
      END
C
      SUBROUTINE HOWARD(M, N, U, Z, NID, IVAR, NVMAX, NIMAX, NMAX,
     *   NVMAX1, NV, NMAX1, WB, WDB, WD2B, DB, D2B, BMN)
C
C        ALGORITHM AS 122.1 APPL. STATIST. (1984) VOL.33, NO.1
C
      INTEGER IVAR(NVMAX)
      REAL U(NMAX), Z(NVMAX, NIMAX), WB(NMAX1), WDB(NVMAX, NMAX1),
     * WD2B(NVMAX1, NMAX1), DB(NVMAX), D2B(NVMAX1)
C
      REAL BMN, ZERO, ONE
C
      DATA ZERO /0.0/, ONE /1.0/
C
      NV1 = NV * (NV + 1) / 2
C
C        INITIALISE THE REQUIRED MATRICES AND VECTORS FOR THE RECURSION
C
      WB(1) = ONE
      I = 0
      DO 1 L = 1, NV
      WDB(L, 1) = ZERO
      DO 1 LL = 1, L
      I = I + 1
      WD2B(I, 1) = ZERO
    1 CONTINUE
      M1 = M + 1
      IF (M .EQ. 0) GOTO 9
      DO 3 J = 2, M1
      WB(J) = ZERO
      I = 0
      DO 2 L = 1, NV
      WDB(L, J) = ZERO
      DO 2 LL = 1, L
      I = I + 1
      WD2B(I, J) = ZERO
    2 CONTINUE
    3 CONTINUE
      NMMP1 = N - M + 1
      DO 8 IM = 1, NMMP1
      DO 7 J = 1, M
      J1 = J + 1
      I = J + IM - 1
      I1 = I + NID
      IR = 1
      IS = 0
      DO 5 L = 1, NV1
      IS = IS + 1
      IS1 = IVAR(IS)
C error 3: bill.bardsley@man.ac.uk corrected the next line 18/08/2002
C     IR1 = IVAR(IS)
      IR1 = IVAR(IR)
      IF (IS .LE. IR) GOTO 4
      IR = IR + 1
      IS = 1
      IS1 = IVAR(IS)
      IR1 = IVAR(IR)
C
C        CALCULATE SUM(EXP(S(L)*BETA)) OVER ALL COMBINATIONS OF N
C        LABELS TAKEN M AT A TIME?  AND ITS FIRST AND SECOND DERIVATIVES
C        WITH RESPECT TO BETA.
C        NOTE - THE VALUES RETURNED BT THIS ROUTINE ARE THE TRUE VALUES
C        OF THE SUMS DIVIDED BY TME CONSTANT CONST**M, WHERE CONST
C        IS AS CALCULATED IN 2UBROUTINE LOGCCH.
C
C error 4: bill.bardsley@man.ac.uk corrected the next line 18/08/2002
C   4 WD2B(L, I1) = WD2B(L, J1) + U(I) * WD2B(L, J) + Z(IR1, I1) *
    4 WD2B(L, J1) = WD2B(L, J1) + U(I) * WD2B(L, J) + Z(IR1, I1) *
     *  Z(IS1, I1) * U(I) * WB(J) + Z(IR1, I1) * U(I) * WDB(IS, J) +
     *  Z(IS1, I1) * U(I) * WDB(IR, J)
    5 CONTINUE
      DO 6 L = 1, NV
      L1 = IVAR(L)
      WDB(L, J1) = WDB(L, J1) + U(I) * WDB(L, J) + Z(L1, I1) * U(I) *
     *  WB(J)
    6 CONTINUE
      WB(J1) = WB(J1) + U(I) * WB(J)
    7 CONTINUE
    8 CONTINUE
    9 BMN = WB(M1)
      I = 0
      DO 10 L = 1, NV
      DB(L) = WDB(L, M1)
      DO 10 LL = 1, L
      I = I + 1
      D2B(I) = WD2B(I, M1)
   10 CONTINUE
      RETURN
      END

      REAL FUNCTION ALGFAC(I)
C
C        ALGORITHM AS 196.2 APPL. STATIST. (1984) VOL.33, NO.1
C
C        EVALUATES THE NATURAL LOGARITHM OF GAMMA(I+1)=LN(I-FACTORIAL)
C        (FORTRAN VERSION OF ACM291 - M.C.PIKE AND I.D.HILL,CACM,9,1966)
C
      REAL HALF, CONST, ONE, TWELVE, TREE60, AA, AK, B
C
      REAL ALOG, FLOAT
C
      DATA HALF, CONST, ONE, TWELVE, TREE60 /0.5, 0.9189385333, 1.0,
     *  12.0, 360.0/

      IF (I .GT. 7) GOTO 3
      B = ONE
      IF (I .LE. 1) GOTO 2
      DO 1 K = 2, I
    1 B = B * FLOAT(K)
    2 ALGFAC = ALOG(B)
      RETURN
    3 AA = FLOAT(I)
      AK = ONE / AA
      ALGFAC = (AA + HALF) * ALOG(AA) + AK * (ONE / TWELVE - (AK * AK)
     *   / TREE60) - AA + CONST
      RETURN
      END
