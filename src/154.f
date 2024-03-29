CSTART OF AS 154
      SUBROUTINE STARMA(IP, IQ, IR, NP, PHI, THETA, A, P, V, THETAB,
     $  XNEXT, XROW, RBAR, NRBAR, IFAULT)
C
C        ALGORITHM AS 154  APPL. STATIST. (1980) VOL.29, P.311
C
C        INVOKING THIS SUBROUTINE SETS THE VALUES OF V AND PHI, AND
C        OBTAINS THE INITIAL VALUES OF A AND P.
C        THIS ROUTINE IS NOT SUITABLE FOR USE WITH AN AR(1) PROCESS.
C        IN THIS CASE THE FOLLOWING INSTRUCTIONS SHOULD BE USED FOR
C        INITIALISATION.
C        V(1) = 1.0
C        A(1) = 0.0
C        P(1) = 1.0 / (1.0 - PHI(1) * PHI(1))
C
      REAL PHI(IR), THETA(IR), A(IR), P(NP), V(NP), THETAB(NP),
     $  XNEXT(NP), XROW(NP), RBAR(NRBAR), VJ, PHII, PHIJ, SSQERR,
     $  RECRES, YNEXT, ZERO, ONE
C
      DATA ZERO, ONE /0.0, 1.0/
C
C        CHECK FOR FAILURE INDICATION.
C
      IFAULT = 0
      IF (IP .LT. 0) IFAULT = 1
      IF (IQ .LT. 0) IFAULT = IFAULT + 2
      IF (IP .EQ. 0 .AND. IQ .EQ. 0) IFAULT = 4
      K = IQ + 1
      IF (K .LT. IP) K = IP
      IF (IR .NE. K) IFAULT = 5
      IF (NP .NE. IR * (IR + 1) / 2) IFAULT = 6
      IF (NRBAR .NE. NP * (NP - 1) / 2) IFAULT = 7
      IF (IR .EQ. 1) IFAULT = 8
      IF (IFAULT .NE. 0) RETURN
C
C        NOW SET A(0), V AND PHI.
C
      DO 10 I = 2, IR
      A(I) = ZERO
      IF (I .GT. IP) PHI(I) = ZERO
      V(I) = ZERO
      IF (I .LE. IQ + 1) V(I) = THETA(I - 1)
   10 CONTINUE
      A(1) = ZERO
      IF (IP .EQ. 0) PHI(1) = ZERO
      V(1) = ONE
      IND = IR
      DO 20 J = 2, IR
      VJ = V(J)
      DO 20 I = J, IR
      IND = IND + 1
      V(IND) = V(I) * VJ
   20 CONTINUE
C
C        NOW FIND P(0).
C
      IF (IP .EQ. 0) GOTO 300
C
C        THE SET OF EQUATIONS S * VEC(P(0)) = VEC(V)
C        IS SOLVED FOR VEC(P(0)).
C        S IS GENERATED ROW BY ROW IN THE ARRAY XNEXT.
C        THE ORDER OF ELEMENTS IN P IS CHANGED, SO AS TO
C        BRING MORE LEADING ZEROS INTO THE ROWS OF S,
C        HENCE ACHIEVING A REDUCTION OF COMPUTING TIME.
C
      IRANK = 0
      SSQERR = ZERO
      DO 40 I = 1, NRBAR
   40 RBAR(I) = ZERO
      DO 50 I = 1, NP
      P(I) = ZERO
      THETAB(I) = ZERO
      XNEXT(I) = ZERO
   50 CONTINUE
      IND = 0
      IND1 = 0
      NPR = NP - IR
      NPR1 = NPR + 1
      INDJ = NPR1
      IND2 = NPR
      DO 110 J = 1, IR
      PHIJ = PHI(J)
      XNEXT(INDJ) = ZERO
      INDJ = INDJ + 1
      INDI = NPR1 + J
      DO 110 I = J, IR
      IND = IND + 1
      YNEXT = V(IND)
      PHII = PHI(I)
      IF (J .EQ. IR) GOTO 100
      XNEXT(INDJ) = -PHII
      IF (I .EQ. IR) GOTO 100
      XNEXT(INDI) = XNEXT(INDI) - PHIJ
      IND1 = IND1 + 1
      XNEXT(IND1) = -ONE
  100 XNEXT(NPR1) = -PHII * PHIJ
      IND2 = IND2 + 1
      IF (IND2 .GT. NP) IND2 = 1
      XNEXT(IND2) = XNEXT(IND2) + ONE
      CALL INCLU2(NP, NRBAR, ONE, XNEXT, XROW, YNEXT,
     $  P, RBAR, THETAB, SSQERR, RECRES, IRANK, IFAIL)
C
C        NO NEED TO CHECK IFAIL AS WEIGHT = 1.0
C
      XNEXT(IND2) = ZERO
      IF (I .EQ. IR) GOTO 110
      XNEXT(INDI) = ZERO
      INDI = INDI + 1
      XNEXT(IND1) = ZERO
  110 CONTINUE
      CALL REGRES(NP, NRBAR, RBAR, THETAB, P)
C
C        NOW RE-ORDER P.
C
      IND = NPR
      DO 200 I = 1, IR
      IND = IND + 1
      XNEXT(I) = P(IND)
  200 CONTINUE
      IND = NP
      IND1 = NPR
      DO 210 I = 1, NPR
      P(IND) = P(IND1)
      IND = IND - 1
      IND1 = IND1 - 1
  210 CONTINUE
      DO 220 I = 1, IR
  220 P(I) = XNEXT(I)
      RETURN
C
C        P(0) IS OBTAINED BY BACKSUBSTITUTION FOR
C        A MOVING AVERAGE PROCESS.
C
  300 INDN = NP + 1
      IND = NP + 1
      DO 310 I = 1, IR
      DO 310 J = 1, I
      IND = IND - 1
      P(IND) = V(IND)
      IF (J .EQ. 1) GOTO 310
      INDN = INDN - 1
      P(IND) = P(IND) + P(INDN)
  310 CONTINUE
      RETURN
      END
C
      SUBROUTINE KARMA(IP, IQ, IR, NP, PHI, THETA, A, P,
     $  V, N, W, RESID, SUMLOG, SSQ, IUPD, DELTA, E, NIT)
C
C        ALGORITHM AS 154.1  APPL. STATIST. (1980) VOL.29, P.311
C
C        INVOKING THIS SUBROUTINE UPDATES A, P, SUMLOG AND SSQ BY
C        INCLUSION OF DATA VALUES W(1) TO W(N). THE CORRESPONDING
C        VALUES OF RESID ARE ALSO OBTAINED.
C        WHEN FT IS LESS THAN (1 + DELTA), QUICK RECURSIONS ARE USED.
C
      REAL PHI(IR), THETA(IR), A(IR), P(NP), V(NP), W(N), RESID(N),
     $  E(IR), SUMLOG, SSQ, DELTA, WNEXT, A1, DT, ET, FT, UT, G,
     $  ZERO, ZLOG, ZSQRT
C
      DATA ZERO /0.0/
C
      ZLOG(G) = ALOG(G)
      ZSQRT(G) = SQRT(G)
C
      IR1 = IR - 1
      DO 10 I = 1, IR
   10 E(I) = ZERO
      INDE = 1
C
C        FOR NON-ZERO VALUES OF NIT, PERFORM QUICK RECURSIONS.
C
      IF (NIT .NE. 0) GOTO 600
      DO 500 I = 1, N
      WNEXT = W(I)
C
C        PREDICTION.
C
      IF (IUPD .EQ. 1 .AND. I .EQ. 1) GOTO 300
C
C        HERE DT = FT - 1.0
C
      DT = ZERO
      IF (IR .NE. 1) DT = P(IR + 1)
      IF (DT .LT. DELTA) GOTO 610
      A1 = A(1)
      IF (IR .EQ. 1) GOTO 110
      DO 100 J = 1, IR1
  100 A(J) = A(J + 1)
  110 A(IR) = ZERO
      IF (IP .EQ. 0) GOTO 200
      DO 120 J = 1, IP
  120 A(J) = A(J) + PHI(J) * A1
  200 IND = 0
      INDN = IR
      DO 210 L = 1, IR
      DO 210 J = L, IR
      IND = IND + 1
      P(IND) = V(IND)
      IF (J .EQ. IR) GOTO 210
      INDN = INDN + 1
      P(IND) = P(IND) + P(INDN)
  210 CONTINUE
C
C        UPDATING.
C
  300 FT = P(1)
      UT = WNEXT - A(1)
      IF (IR .EQ. 1) GOTO 410
      IND = IR
      DO 400 J = 2, IR
      G = P(J) / FT
      A(J) = A(J) + G * UT
      DO 400 L = J, IR
      IND = IND + 1
      P(IND) = P(IND) - G * P(L)
  400 CONTINUE
  410 A(1) = WNEXT
      DO 420 L = 1, IR
  420 P(L) = ZERO
      RESID(I) = UT / ZSQRT(FT)
      E(INDE) = RESID(I)
      INDE = INDE + 1
      IF (INDE .GT. IQ) INDE = 1
      SSQ = SSQ + UT * UT / FT
      SUMLOG = SUMLOG + ZLOG(FT)
  500 CONTINUE
      NIT = N
      RETURN
C
C        QUICK RECURSIONS
C
  600 I = 1
  610 NIT = I - 1
      DO 650 II = I, N
      ET = W(II)
      INDW = II
      IF (IP .EQ. 0) GOTO 630
      DO 620 J = 1, IP
      INDW = INDW - 1
      IF (INDW .LT. 1) GOTO 630
      ET = ET - PHI(J) * W(INDW)
  620 CONTINUE
  630 IF (IQ .EQ. 0) GOTO 645
      DO 640 J = 1, IQ
      INDE = INDE - 1
      IF (INDE .EQ. 0) INDE = IQ
      ET = ET - THETA(J) * E(INDE)
  640 CONTINUE
  645 E(INDE) = ET
      RESID(II) = ET
      SSQ = SSQ + ET * ET
      INDE = INDE + 1
      IF (INDE .GT. IQ) INDE = 1
  650 CONTINUE
      RETURN
      END
C
      SUBROUTINE KALFOR(M, IP, IR, NP, PHI, A, P, V, WORK)
C
C        ALGORITHM AS 154.2  APPL. STATIST. (1980) VOL.29, P.311
C
C        INVOKING THIS SUBROUTINE OBTAINS PREDICTIONS
C        OF A AND P, M STEPS AHEAD.
C
      REAL PHI(IR), A(IR), P(NP), V(NP), WORK(IR), DT,
     $  A1, PHII, PHIJ, PHIJDT, ZERO
C
      DATA ZERO /0.0/
C
      IR1 = IR - 1
      DO 300 L = 1, M
C
C        PREDICT A.
C
      A1 = A(1)
      IF (IR .EQ. 1) GOTO 110
      DO 100 I = 1, IR1
  100 A(I) = A(I + 1)
  110 A(IR) = ZERO
      IF (IP .EQ. 0) GOTO 200
      DO 120 J = 1, IP
  120 A(J) = A(J) + PHI(J) * A1
C
C        PREDICT P.
C
  200 DO 210 I = 1, IR
  210 WORK(I) = P(I)
      IND = 0
      IND1 = IR
      DT = P(1)
      DO 220 J = 1, IR
      PHIJ = PHI(J)
      PHIJDT = PHIJ * DT
      DO 220 I = J, IR
      IND = IND + 1
      PHII = PHI(I)
      P(IND) = V(IND) + PHII * PHIJDT
      IF (J .LT. IR) P(IND) = P(IND) + WORK(J + 1) * PHII
      IF (I .EQ. IR) GOTO 220
      IND1 = IND1 + 1
      P(IND) = P(IND) + WORK(I + 1) * PHIJ + P(IND1)
  220 CONTINUE
  300 CONTINUE
      RETURN
      END
C
      SUBROUTINE INCLU2(NP, NRBAR, WEIGHT, XNEXT, XROW, YNEXT, D, RBAR,
     $  THETAB, SSQERR, RECRES, IRANK, IFAULT)
C
C        ALGORITHM AS 154.3  APPL. STATIST. (1980) VOL.29, P.311
C
C        FORTRAN VERSION OF REVISED VERSION OF ALGORITHM AS 75.1
C        APPL. STATIST. (1974) VOL.23, P.448
C        SEE REMARK AS R17 APPL. STATIST. (1976) VOL.25, P.323
C
      REAL XNEXT(NP), XROW(NP), D(NP), RBAR(NRBAR), THETAB(NP),
     $  WEIGHT, YNEXT, SSQERR, RECRES, WT, Y, DI, DPI, XI, XK,
     $  CBAR, SBAR, RBTHIS, ZERO, ZSQRT
C
      DATA ZERO /0.0/
C
      ZSQRT(Y) = SQRT(Y)
C
C        INVOKING THIS SUBROUTINE UPDATES D, RBAR, THETAB, SSQERR
C        AND IRANK BY THE INCLUSION OF XNEXT AND YNEXT WITH A
C        SPECIFIED WEIGHT. THE VALUES OF XNEXT, YNEXT AND WEIGHT WILL
C        BE CONSERVED. THE CORRESPONDING VALUE OF RECRES IS CALCULATED.
C
      Y = YNEXT
      WT = WEIGHT
      DO 10 I = 1, NP
   10 XROW(I) = XNEXT(I)
      RECRES = ZERO
      IFAULT = 1
      IF (WT .LE. ZERO) RETURN
      IFAULT = 0
C
      ITHISR = 0
      DO 50 I = 1, NP
      IF (XROW(I) .NE. ZERO) GOTO 20
      ITHISR = ITHISR + NP - I
      GOTO 50
   20 XI = XROW(I)
      DI = D(I)
      DPI = DI + WT * XI * XI
      D(I) = DPI
      CBAR = DI / DPI
      SBAR = WT * XI / DPI
      WT = CBAR * WT
      IF (I .EQ. NP) GOTO 40
      I1 = I + 1
      DO 30 K = I1, NP
      ITHISR = ITHISR + 1
      XK = XROW(K)
      RBTHIS = RBAR(ITHISR)
      XROW(K) = XK - XI * RBTHIS
      RBAR(ITHISR) = CBAR * RBTHIS + SBAR * XK
   30 CONTINUE
   40 XK = Y
      Y = XK - XI * THETAB(I)
      THETAB(I) = CBAR * THETAB(I) + SBAR * XK
      IF (DI .EQ. ZERO) GOTO 100
   50 CONTINUE
      SSQERR = SSQERR + WT * Y * Y
      RECRES = Y * ZSQRT(WT)
      RETURN
  100 IRANK = IRANK + 1
      RETURN
      END
C
      SUBROUTINE REGRES(NP, NRBAR, RBAR, THETAB, BETA)
C
C        ALGORITHM AS 154.4  APPL. STATIST. (1980) VOL.29, P.311
C
C        REVISED VERSION OF ALGORITHM AS 75.4
C        APPL. STATIST. (1974) VOL.23, P.448
C        INVOKING THIS SUBROUTINE OBTAINS BETA BY BACKSUBSTITUTION
C        IN THE TRIANGULAR SYSTEM RBAR AND THETAB.
C
      REAL RBAR(NRBAR), THETAB(NP), BETA(NP), BI
      ITHISR = NRBAR
      IM = NP
      DO 50 I = 1, NP
      BI = THETAB(IM)
      IF (IM .EQ. NP) GOTO 30
      I1 = I - 1
      JM = NP
      DO 10 J = 1, I1
      BI = BI - RBAR(ITHISR) * BETA(JM)
      ITHISR = ITHISR - 1
      JM = JM - 1
   10 CONTINUE
   30 BETA(IM) = BI
      IM = IM - 1
   50 CONTINUE
      RETURN
      END
CEND OF AS 154
