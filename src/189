C UKC NETLIB DISTRIBUTION COPYRIGHT 1990 RSS
C
      SUBROUTINE BBML(N, IX, IN, W, P, RL, MRL, ITER, CCRIT, MEW, THETA, 
     *  SEM, SETH, RNL, IFAULT)
C      
C        ALGORITHM AS189 APPL. STATIST. (1983) VOL.32, NO.2
C      
C        SUBROUTINE FOR CALCULATING THE MAXIMUM LIKELIHOOD ESTIMATES
C        OF THE PARAMETERS OF THE BETA BINOMIAL DISTRIBUTION
C      
      REAL W(N), P(N), CCRIT, MEW, THETA, SEM, SETH, RNL, INF, DUM, 
     *  FD(2), SD(3), TD(4), UB(2), DEL, EPS, A, B, C, D, E, F
      INTEGER IX(N), IN(N), RL(MRL,3), LM(3), RD1(2,2), RD2(2,3), 
     *  RD3(2,4)
      LOGICAL MC
      PARAMETER (INF = 1.0E6)
      DATA
     *  RD1(1,1), RD1(2,1), RD1(1,2), RD1(2,2)/1,-1,1,1/,
     *  RD2(1,1), RD2(2,1), RD2(1,2), RD2(2,2), 
     *  RD2(1,3), RD2(2,3)/-1,-1,-1,1,-1,-1/,
     *  RD3(1,1), RD3(2,1), RD3(1,2), RD3(2,2), RD3(1,3), 
     *  RD3(2,3), RD3(1,4), RD3(2,4)/2,-2,2,2,2,-2,2,2/
C
      I = ITER
      ITER = 0
      MC = .TRUE.
      UB(1) = 0.01
      UB(2) = 0.01
C
C        SET THE ARRAYS RL AND LM
C
      CALL SET(N, IX, IN, RL, MRL, LM, IFAULT)
      IF(IFAULT.NE.0) RETURN
      SEM = -1.0
      SETH = -1.0
      NND = 0
C
C        CALCULATION OF INITIAL ESTIMATES (BY MOMENTS)
C
      CALL BBME(N, IX, IN, W, P, INF, MEW, THETA)
      IF(THETA.EQ.INF) GOTO 50
C
C        NEWTON-RAPHSON ITERATION ON FIRST DERIVATIVES
C
    5 IF(ITER.LE.I) GOTO 10
      IFAULT = 7
      GOTO 60
C
C        CALCULATE FIRST DERIVATIVES OF LOG LIKELIHOOD
C
   10 CALL GDER(MEW, THETA, RL, MRL, LM, 2, RD1, FD)
C
C        CALCULATE SECOND DERIVATIVES OF LOG_LIKELIHOOD
C
      CALL GDER(MEW, THETA, RL, MRL, LM, 3, RD2, SD)
C
C        CALCULATE THIRD DERIVATIVES OF LOG LIKELIHOOD
C
      CALL GDER(MEW, THETA, RL, MRL, LM, 4, RD3, TD)
C
C        CALCULATE INCREMENTS
C
      DUM = SD(1)*SD(3) - SD(2)*SD(2)
      IF(SD(1).LT.0.0.AND.DUM.GT.0.0) GOTO 15
C
C        NON NEGATIVE DEFINITE MATRIX
C
      NND = NND+1
C
C        SD(1) IS ALWAYS NEGATIVE SO A GRADIENT STEP IS MADE ON MEW
C
      A = MEW - FD(1)/SD(1)
      B = THETA
      IF(FD(2).NE.0.0) B = B + SIGN(UB(2),FD(2))
      IF(A.LE.0.0) A = 0.0001
      IF(A.GE.1.0) A = 0.9999
      IF(B.LT.0.0) B = 0.0
      IF(B.GT.INF) B = INF
      CALL BBL(MEW, THETA, RL, MRL, LM, C)
      CALL BBL(A, B, RL, MRL, LM, D)
      IF(NND.GT.10.OR.C.GE.D) GOTO 40
      ITER = ITER+1
      MEW = A
      THETA = B
      GOTO 5
   15 DEL = (FD(2)*SD(2) - FD(1)*SD(3))/DUM
      EPS = (FD(1)*SD(2) - FD(2)*SD(1))/DUM
C
C        CHECK LIPSCHITZ CONDITION SATISFIED
C
      A = SD(2)*TD(2) - TD(1)*SD(3)
      B = SD(2)*TD(3) - TD(2)*SD(3)
      C = TD(1)*SD(2) - TD(2)*SD(1)
      D = SD(2)*TD(2) - SD(1)*TD(3)
      E = SD(2)*TD(4) - TD(3)*SD(3)
      F = TD(3)*SD(2) - TD(4)*SD(1)
      A = DEL*A + EPS*B
      C = DEL*C + EPS*D
      E = DEL*B + EPS*E
      F = DEL*D + EPS*F
      DUM = (A*A + C*C + E*E + F*F)/(DUM*DUM)
      IF(DUM.GE.1.0) GOTO 20
      IF(ABS(DEL).LE.CCRIT.AND.ABS(EPS).LE.CCRIT) MC = .FALSE.
      GOTO 45
C
C        FAILURE OF LIPSCHITZ CONDITION. A STEP IN THE DIRECTION OF THE
C        GRADIENT IS MADE.
C
   20 A = FD(1)*FD(1)
      B = FD(2)*FD(2)
      C = A*SD(1) + 2.0*SD(2)*FD(1)*FD(2) + B*SD(3)
      IF(C.NE.0.0) GOTO 25
      DEL = 0.0
      IF(FD(1).NE.0.0) DEL = SIGN(UB(1),FD(1))
      EPS = 0.0
      IF(FD(2).NE.0.0) EPS = SIGN(UB(2),FD(2))
      GOTO 30
   25 C = -(A+B)/C
      DEL = C*FD(1)
      EPS = C*FD(2)
      IF(ABS(DEL).GT.UB(1)) DEL = SIGN(UB(1),DEL)
      UB(1) = 2.0*ABS(DEL)
      IF(ABS(EPS).GT.UB(2)) EPS = SIGN(UB(2),EPS)
      UB(2) = 2.0*ABS(EPS)
   30 CALL BBL(MEW, THETA, RL, MRL, LM, C)
   35 A = MEW + DEL
      B = THETA + EPS
      IF(A.LE.0.0) A = 0.0001
      IF(A.GE.1.0) A = 0.9999
      DEL = A - MEW
      IF(B.LT.0.0) B = 0.0
      IF(B.GT.INF) B = INF
      EPS = B - THETA
      CALL BBL(A, B, RL, MRL, LM, D)
C
C        CHECK TO SEE IF GRADIENT STEP HAS INCREASED LOG LIKELIHOOD
C
      IF(D.GT.C) GOTO 45
      DEL = DEL/2.0
      EPS = EPS/2.0
      IF(ABS(DEL).GT.CCRIT.OR.ABS(EPS).GT.CCRIT) GOTO 35
   40 IFAULT = 8
      GOTO 60
   45 ITER = ITER + 1
      A = MEW + DEL
      B = THETA + EPS
      IF(A.GT.0.0.AND.A.LT.1.0.AND.B.GE.0.0.AND.B.LE.INF) GOTO 55
      IF(A.LE.0.0) MEW = 0.0
      IF(A.GE.1.0) MEW = 1.0
      IF(B.LT.0.0) THETA = 0.0
      IF(B.GT.INF) THETA = INF
   50 IFAULT = 6
      GOTO 60
   55 MEW = A
      THETA = B
      IF(MC) GOTO 5
C
C        CALCULATE LOG LIKELIHOOD AND S.E.S
C
      IF(SD(1).LT.0.0) SEM = SQRT(-1.0/SD(1))
      IF(SD(3).LT.0.0) SETH = SQRT(-1.0/SD(3))
   60 CALL BBL(MEW, THETA, RL, MRL, LM, RNL)
      RETURN
      END
C
      SUBROUTINE BBME(N, IX, IN, W, P, INF, MEW, THETA)
C
C        ALGORITHM AS 189.1 APPL. STATIST. (1983) VOL.32, NO.2
C
C        SUBROUTINE TO ESTIMATE MEW AND THETA OF THE BETA BINOMIAL
C        DISTRIBUTION BY THE METHOD OF MOMENTS
C
      REAL W(N), P(N), INF, MEW, THETA, D1, D2, R, S, TP, WT
      INTEGER IX(N), IN(N)
      LOGICAL J
C
      J = .FALSE.
      DO 5 I = 1,N
        W(I) = FLOAT(IN(I))
        P(I) = FLOAT(IX(I))/W(I)
    5 CONTINUE
   10 WT = 0.0
      TP = 0.0
      DO 15 I = 1,N
        WT = WT+W(I)
        TP = TP+W(I)*P(I)
   15 CONTINUE
      TP = TP/WT
      S = 0.0
      D1 = 0.0
      D2 = 0.0
      DO 20 I = 1,N
        R = P(I)-TP
        S = S+W(I)*R*R
        R = W(I)*(1.0-W(I)/WT)
        D1 = D1+R/FLOAT(IN(I))
        D2 = D2+R
   20 CONTINUE
      S = FLOAT(N-1)*S/FLOAT(N)
      R = TP*(1.0-TP)
      IF(R.EQ.0.0) GOTO 30
      R = (S-R*D1)/(R*(D2-D1))
      IF(R.LT.0.0) R = 0.0
      IF(J) GOTO 30
      DO 25 I = 1,N
   25 W(I) = W(I)/(1.0+R*(W(I)-1.0))
      J = .TRUE.
      GOTO 10
   30 MEW = TP
      IF(R.GE.1.0) GOTO 35
      THETA = R/(1.0-R)
      IF(THETA.LE.INF) RETURN
   35 THETA = INF
      RETURN
      END
C
      SUBROUTINE SET(N, IX, IN, RL, MRL, LM, IFAULT)
C
C        ALGORITHM AS 189.2 APPL. STATIST. (1983) VOL.32, NO.2
C
C        SUBROUTINE FOR SETTING UP ARRAY FOR CALCULATION OF
C        THE BETA BINOMIAL LOG LIKELIHOOD AND ITS DERIVATIVES
C
      INTEGER IX(N), IN(N), RL(MRL,3), LM(3)
C
C     TEST ADMISSIBILITY OF DATA
C
      IF(N.GT.1) GOTO 5
      IFAULT = 1
      RETURN
    5 DO 10 I = 1,N
        IF(IX(I).GT.0) GOTO 15
   10 CONTINUE
      IFAULT = 2
      RETURN
   15 DO 20 I = 1,N
        IF(IX(I).LT.IN(I)) GOTO 25
   20 CONTINUE
      IFAULT = 3
      RETURN
C
C        FORM MATRIX OF COUNTS
C
   25 IFAULT = 4
      DO 30 I = 1,3
        LM(I) = 0
        DO 30 J = 1,MRL
          RL(J,I) = 0
   30 CONTINUE
      DO 65 I = 1,N
        JJ = IX(I)
        MAR = 1
        GOTO 45
   35   JJ = IN(I)-IX(I)
        MAR = 2
        GOTO 45
   40   JJ = IN(I)
        MAR = 3
   45   IF(JJ) 50,60,55
   50   IFAULT = 5
        RETURN
   55   IF(JJ.GT.MRL) RETURN
        IF(JJ.GT.LM(MAR)) LM(MAR) = JJ
        RL(JJ,MAR) = RL(JJ,MAR)+1
   60   GOTO(35,40,65) MAR
   65 CONTINUE
      IFAULT = 0
C
C        EVALUATE NUMBER OF CALLS TO DIFFERENT TERMS OF LIKELIHOOD
C        FUNCTION
C
      DO 75 I = 1,3
        JJ = LM(I)-1
        IF(JJ.LE.0) GOTO 75
        K = JJ
        DO 70 J = 1,JJ
          RL(K,I) = RL(K,I)+RL(K+1,I)
          K = K-1
   70   CONTINUE
   75 CONTINUE
      RETURN
      END
C
      SUBROUTINE BBL(MEW, THETA, RL, MRL, LM, RNL)
C
C        ALGORITHM AS 189.3 APPL. STATIST. (1983) VOL.32, NO.2
C
C        SUBROUTINE FOR CALCULATION OF THE BETA BINOMIAL LOG
C        LIKELIHOOD
C
      REAL MEW, THETA, RNL, A
      INTEGER RL(MRL,3), LM(3)
C
      RNL = 0.0
      MLM = LM(3)
      DO 5 I = 1,MLM
        A = FLOAT(I-1)*THETA
        IF(I.LE.LM(1))RNL = RNL + FLOAT(RL(I,1))*ALOG(MEW+A)
        IF(I.LE.LM(2))RNL = RNL + FLOAT(RL(I,2))*ALOG(1.0-MEW+A)
        RNL = RNL - FLOAT(RL(I,3))*ALOG(1.0+A)
    5 CONTINUE
      RETURN
      END
C
      SUBROUTINE GDER(MEW, THETA, RL, MRL, LM, IDER, RD, PD)
C
C        ALGORITHM AS 189.4 APPL. STATIST. (1983) VOL.32, NO.2
C
C        GENERAL DERIVATIVE SUBROUTINE
C
      REAL MEW, THETA, PD(IDER), A, B, C, D
      INTEGER RL(MRL,3), LM(3), RD(2,IDER)
C
      MLM = LM(3)
      KK = IDER-1
      DO 5 I = 1,IDER
    5 PD(I) = 0.0
      DO 45 I = 1,MLM
        C = FLOAT(I-1)
        A = C*THETA
        DO 40 J = 1,3
          IF(I.GT.LM(J)) GOTO 40
          GOTO (10,15,20) J
   10     D = MEW+A
          GOTO 25
   15     D = 1.0-MEW+A
          GOTO 25
   20     D = 1.0+A
   25     B = FLOAT(RL(I,J))/D**KK
          IF(J.EQ.3) GOTO 35
          DO 30 K = 1,IDER
            PD(K) = PD(K)+FLOAT(RD(J,K))*B
            B = B*C
   30     CONTINUE
          GOTO 40
   35     D = -FLOAT(RD(1,1))*B*C**KK
          PD(IDER) = PD(IDER)+D
   40   CONTINUE
   45 CONTINUE
      RETURN
      END
