      SUBROUTINE ENTNGL(IP, L, M, NPTS, IDES, NP, NINTER, ITERMS, NTP,
     *   A, NA, X, NX, W, NW, IFAULT)
C
C        ALGORITHM AS 166  APPL. STATIST. (1981) VOL.30, NO.3
C
C        CALCULATES THE ENTANGLEMENT MATRIX FOR A SPECIFIED DESIGN
C
C     *** Warning ***
C     A scanner was used to enter this text - it may contain errors!
C
      INTEGER L(IP), M(IP), IDES(NP), ITERMS(NTP)
      REAL W(NW), A(NA), X(NX)
      LOGICAL XXONLY
      DATA ZERO, POINT4, ONE /0.0, 0.4, 1.0/
C
C        CHECK PARAMETERS FOR ERRORS
C
      IFAULT = 0
      ID = 0
      XXONLY = NINTER .EQ. 0
      IF (IP .LE. 0 .OR. NPTS .LE. 0 .OR. NINTER .LT. 0 .OR.
     *   NTP .LT. IP * NINTER) IFAULT = 5
      IF (NP .LT. NPTS * IP) IFAULT = 6
      LSUM = 0
      MSUM = 0
      LMSUM = 0
      MAXL = 0
      DO 10 JA = 1, IP
        LJA = IABS(L(JA))
        MJA = M(JA)
        IF (MJA .LE. 0 .OR. MJA .GE. LJA .OR. LJA .EQ. 0 .OR.
     *     LJA .GT. NPTS) IFAULT = 7
        IF (XXONLY) GOTO 3
        IT = JA
        DO 2 JB = 1, NINTER
          J = ITERMS(IT)
          IF (J .LT. 0 .OR. J .GT. MJA) IFAULT = 13
          IT = IT + IP
    2   CONTINUE
    3   IF (LJA .GT. MAXL) MAXL = LJA
        LSUM = LSUM + LJA
        MSUM = MSUM + MJA
        LMSUM = LMSUM + LJA * MJA
        DO 5 I = 1, NPTS
          ID = ID + 1
          LIJ = IDES(ID)
          IF (LIJ .LT. 0 .OR. LIJ .GE. LJA) IFAULT = 1
    5   CONTINUE
   10 CONTINUE
      MSIZE = MSUM + NINTER
      LASTA = MSIZE * (MSIZE + 1) / 2
      IF (NA .LT. LASTA) IFAULT = 8
      XXONLY = NX .EQ. 1
      IF (NX .LT. MSIZE * NPTS .AND. .NOT. XXONLY) IFAULT = 9
      LASTW = LSUM + MAXL + MAXL + 2 + LMSUM + MSUM + NINTER + NINTER
      IF (NW .LT. LASTW) IFAULT = 10
      IF (IFAULT .NE. 0) RETURN
C
C        INITIALIZE LAST PART OF W TO ZERO
C
      IWG = LASTW - NINTER + 1
      DO 15 IW = IWG, LASTW
   15 W(IW) = ZERO
C
C        CALCULATE ORTHOGONAL POLYNOMIAL COEFFICIENTS FOR EACH FACTOR
C
      LIMWD = LSUM + MAXL + MAXL + 2
      IWD = LIMWD
      ID = 0
      LIMWB = LSUM + 1
      LIMWC = LIMWB + MAXL
      DO 100 JA = 1, IP
        KR = 0
        LJA = IABS(L(JA))
        MJA = M(JA)
        IF (LJA .EQ. 2) GOTO 52
C
C        (1)  FACTORS WITH THREE OR MORE LEVELS
C
        DO 16 J = 1, LJA
          IW = LSUM + J
          W(IW) = ZERO
   16   CONTINUE
        DO 18 I = 1, NPTS
          ID = ID + 1
          IW = IDES(ID) + LIMWB
          W(IW) = W(IW) + ONE
   18   CONTINUE
        CONST = FLOAT(NPTS / LJA)
        IFAULT = 12
        DO 19 J = 1, LJA
          IW = LSUM + J
          IF (ABS(W(IW) - CONST) .GT. POINT4) RETURN
   19   CONTINUE
        CONST = ONE / SQRT(CONST)
        IFAULT = 0
        IF (L(JA) .LT. 0) GOTO 40
      DO 20 J = 1, LJA
        IW = IWA + J
        W(IW) = FLOAT(J)
   20 CONTINUE
   40 IW = IWA + 1
      CALL ORTHON(0, W(IW), W(LIMWB), W(LIMWC), LJA, LJA+2, IFAULT)
      IF (IFAULT .NE. 0 .AND. IFAULT .NE. 2) RETURN
      IF (IFAULT .EQ. 2) IFAIL = 2
   50 KR = KR + 1
      IF (KR .GT. MJA) GOTO 80
      CALL ORTHON(KR, W(IW), W(LIMWB), W(LIMWC), LJA, LJA+2, IFAULT)
      IF (IFAULT .NE. 0 .AND. IFAULT . NE . 2) RETURN
      GOTO 58
C
C        (2) FACTORS WITH TWO LEVELS ONLY (POSSIBLY OCCURRING UNEQUALLY)
C
   52 ISUM = 0
      DO 54 I = 1, NPTS
        ID = ID + 1
        ISUM = ISUM + IDES(ID)
   54 CONTINUE
      CONST = FLOAT(ISUM) / FLOAT(NPTS)
      W(LIMWB) = -CONST
      W(LIMWB + 1) = ONE - CONST
      CONST = 1.0 / SQRT(CONST * FLOAT(NPTS - ISUM))
C
C        COPY COEFFICIENTS FROM WORKSPACE AREA TO STORAGE AREA
C
   58 DO 60 J = 1, LJA
        IWD = IWD + 1
        IWB = LSUM + J
        W(IWD) = W(IWB) + CONST
   60 CONTINUE
      IF (LJA .NE. 2) GOTO 50
   80 IWA = IWA + LJA
  100 CONTINUE
C
C        PREPARE TO SET UP A
C
      LIMWE = IWD
      DO 120 IA = 1, LASTA
  120 A(IA) = ZERO
C
C       ADD EACH DESIGN POINT INTO A
C
      LIMWD = LIMWD + 1
      DO 240 I = 1, NPTS
      ID = 0
C
C        (1) MAIN EFFECTS
C
      IWD = LIMWD
      IWE = LIMWE
      DO 140 JA = 1, IP
        LJA = IABS(L(JA))
        MJA = M(JA)
        LIJ = ID + I
        LIJ = IDES(LIJ)
        DO 130 J = 1, MJA
          IWE = IWE + 1
          IW = IWD + LIJ
           W(IWE) = W(IW)
           IWD = IWD + LJA
  130   CONTINUE
        ID = ID + NPTS
  140 CONTINUE
C
C       (2)  INTERACTIONS
C
      IF  (NINTER .EQ. 0)  GOTO  200
      IT = 0
      IWF = IWE
      IWG = IWF + NINTER
      DO 180 JB = 1, NINTER
        IWF = IWF + 1
        IWG = IWG + 1
        WW = ONE
        IWE = LIMWE
        DO 160 JA = 1, IP
          IT = IT + 1
          IF (ITERMS(IT) .EQ. 0) GOTO 150
          IW = IWE + ITERMS(IT)
          WW = WW * W(IW)
  150     IWE = IWE + M(JA)
  160   CONTINUE
        W(IWF) = WW
        W(IWG) = W(IWG) + WW
  180 CONTINUE
C
C        CALCULATE THE TERMS IN THE SSP MATRIX
C
  200 IA = 1
      DO 220 JA = 1, MSIZE
        IW = LIMWE + JA
        DO 210 JB = 1, JA
          IW2 = LIMWE + JB
          A(IA) = A(IA) + W(IW) * W(IW2)
          IA = IA + 1
  210   CONTINUE
  220 CONTINUE
C                                                                                                 C
C        STORE CONTRASTS IN X
C
      IF (XXONLY) GOTO 240
      I X  =  I
      DO 230 JA = 1, MSIZE
        IW = LIMWE + JA
        X(IX) = W(IW)
        IX = IX + NPTS
  230 CONTINUE
  240 CONTINUE
C
C        ADJUST INT X INT TERMS IN A ABOUT THEIR MEANS
C        (TERMS INVOLVING MAIN EFFECTS ARE ALREADY CENTRED)
C
      IF (NINTER .EQ. 0)  RETURN
      J = MSUM + 1
      IA = J * MSUM / 2
      IA2 = IA
      IWG = IWG - MSIZE
      CONST = ONE / FLOAT(NPTS)
      DO 260 JA = J, MSIZE
        IA = IA + J - 1
        IW = IWG + JA
        DO 250 JB = J, JA
          IA = IA + 1
          IW2 = IWG + JB
          A(IA) = A(IA) - W(IW) * W(IW2) * CONST
  250   CONTINUE
  260 CONTINUE
C
C        CALCULATE 1 / SQRT(DIAG(A)) AND STORE
C
      IWE = LIMWE
      IA = 0
      DO 320 JA = 1, MSIZE
        IA = IA + JA
        IWE = IWE + 1
        W(IWE) = 1.0 / SQRT(LA(IA))
  320 CONTINUE
C
C        CONVERT A FROM A SSP-TYPE MATRIX TO A CO
C
      DO 360 JA = J, MSIZE
        IW = LIMWE + JA
        WW = W(IW)
        DO 350 JB = 1, JA
          IW2 = LIMWE + JB
          IA2 = IA2 + 1
          A(IA2) = A(IA2) * WW * W(IW2)
  350   CONTINUE
  360 CONTINUE
C
C        STANDARDIZE INTERACTION CONTRASTS
C
      IF (XXONLY) RETURN
      IWG = IWE
      IWE = IWE - NINTER
      IX = (J - 1) * NPTS
      DO 400 JA = 1, NINTER
        IWE = IWE +  1
        IWG = IWG + 1
        XBAR = W(IWG) + 1
        WW = W(IWE)
        DO 380 I = 1, NPTS
          IX = IX + 1
          X(IX) = WW * (X(IX) - XBAR)
  380   CONTINUE
  400 CONTINUE
      RETURN
      END

