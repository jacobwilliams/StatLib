      SUBROUTINE EXTR(BP,EP,M,N,F,YM,XM,NE,E1,E2,N1,NP,X1,Y1,X2,Y2,
     +                IFAULT)
C
C        ALGORITHM AS 133 APPL. STATIST. (1978) VOL.27, NO.3
C
C        OPTIMIZATION OF ONE DIMENSIONAL MULTIMODAL FUNCTIONS
C
C     Auxiliary routines required: The user must provide a real function
C     F(X) to return the value of the function to be minimized.
C     Brent's LOCALM has been added to this routine.
C
      DIMENSION X1(NP),Y1(NP),X2(NP),Y2(NP)
      EXTERNAL YL,F
      REAL LOCALM
      REAL ONE,TEN,TWENT5,ZERO,SEVEN,SMALL,NINET8,PT99,TWO,THREE,FIVE
      DOUBLE PRECISION EE1, EE2, EE3, DM17
C
      COMMON /EXTR1/A, EB, YH, XJ1, XJ2, YJ1, YJ2
C
C     K is a machine-dependant constant.   10^(-K) should be set to the
C     smallest number which can be represented in the machine.
C
      DATA K/37/, EPSL/0.0/
      DATA ONE/1.0/, TEN/10.0/, TWENT5/25.0/, ZERO/0.0/, SEVEN/7.0/,
     +    SMALL/3.0E-4/, NINET8/98.0/, PT99/0.99/, TWO/2.0/, THREE/3.0/,
     +    FIVE/5.0/
C
      K1 = K/3
      CUND = -2*K
      CP5 = TEN** (K1-1)
      CP6 = TEN*CP5
      CM6 = ONE/CP6
      CM5 = ONE/CP5
      DM17 = 10.0D0** (-K+2)
      CM17 = DM17
      IFAULT = 0
      IF (NE.GE.6) GO TO 5
      IFAULT = 3
      RETURN

    5 N1 = NE
      N2 = N1 + 2
      CM = M
      A = MIN(BP,EP)
      EB = MAX(BP,EP) - A
      EB2 = EB/TWENT5
      IF (EB.GE.CM5 .AND. EB.LE.CP5) GO TO 10
      IFAULT = 1
      RETURN

   10 DA = EB/REAL(N1-1)
      N6 = N1 + 3
      DO 20 J = 4,N6
        X1(J) = A
        Y1(J) = F(A)*CM
        A = A + DA
   20 CONTINUE
      A = ZERO
      YM = Y1(4)
      XM = X1(4)
C
C        A IS ESTIMATE OF PARAMETER OF WIENER PROCESS
C        (PARAMETER OF MODEL OF FUNCTION F)
C
      DO 30 J = 5,N6
        A = A + (Y1(J)-Y1(J-1))**2
        IF (YM.LE.Y1(J)) GO TO 30
        YM = Y1(J)
        XM = X1(J)
   30 CONTINUE
      YH = YM
      A = SQRT(A/EB)*SEVEN
      IF (A.GE.CM6 .AND. A.LE.CP6) GO TO 40
      IFAULT = 2
      RETURN
C
C        A IS EVALUATED.  N5 IS THE NUMBER OF LOCAL OPTIMA
C        COMPUTED WITH ACCURACY E1, E2
C
   40 N5 = 0
      DO 50 J = 4,N2
   50 Y2(J) = ZERO
C
C        BAYESIAN STEP - FIND INTERVAL OF BIGGEST IMPROVEMENT (J)
C
   60 DO 70 J = 4,N2
        IF (Y2(J).GT.ZERO) GO TO 70
        XJ1 = X1(J)
        XJ2 = X1(J+1)
        YJ1 = Y1(J)
        YJ2 = Y1(J+1)
        TL = (XJ2-XJ1)*SMALL
        Y2(J) = LOCALM(XJ1,XJ2,EPSL,TL,YL,X2(J))
   70 CONTINUE
   80 P1 = Y2(4)
      J = 4
      DO 90 N3 = 5,N2
        IF (Y2(N3).GE.P1) GO TO 90
        J = N3
        P1 = Y2(N3)
   90 CONTINUE
C
C        SET UP X1, Y1, X2, Y2 FOR NEXT STEP
C
      N4 = N1 - J + 3
      N7 = N1 + 4
      DO 100 N3 = 1,N4
        N7 = N7 - 1
        X1(N7+1) = X1(N7)
        Y1(N7+1) = Y1(N7)
        X2(N7) = X2(N7-1)
        Y2(N7) = Y2(N7-1)
  100 CONTINUE
      X1(J+1) = X2(J)
      Y1(J+1) = F(X2(J))*CM
      N4 = 1
      IF (Y1(J+1).GT.YM) GO TO 110
      YM = Y1(J+1)
      XM = X1(J+1)
C
C        IF NEW EVALUATION IS BETTER THAN YM,
C        X2 AND Y2 WILL BE UPDATED (N4 = 2)
C
      N4 = 2
  110 N1 = N1 + 1
      N2 = N1 + 2
      YH = YM
C
C        IF NUMBER OF EVALUATIONS OF F EXCEEDS (N-1) GO TO END
C
      IF (N1.LT.N) GO TO 115
      IFAULT = 4
      GO TO 310
C
C        P4 IS PROBABILITY OF EVALUATING
C        XM, YM WITH ACCURACY E1, E2
C
  115 P4 = ONE
      DO 120 N3 = 4,N2
        IF (Y2(N3).GT.ZERO) GO TO 120
        P1 = -NINET8* ((Y1(N3)-YM+E1)/A)* ((Y1(N3+1)-YM+E1)/A)/
     +         (X1(N3+1)-X1(N3))
        IF (P1.GT.CUND) P4 = P4* (ONE-EXP(P1))
  120 CONTINUE
      IF (P4.GE.PT99) GO TO 310
C
C        IF N5 = 0 A IS UPDATED IN EACH FIFTH STEP
C
      DO 130 N3 = 1,N1,5
        IF (N3.EQ.N1) GO TO 140
  130 CONTINUE
      GO TO 160

  140 IF (N5.GT.0) GO TO 160
      A = ZERO
      N4 = 2
      DO 150 N3 = 4,N2
  150 A = A + ((Y1(N3+1)-Y1(N3))**2)/ (X1(N3+1)-X1(N3))
      A = SQRT(A/REAL(N1-1))*SEVEN
C
C        PREPARATION OF BOUNDS OF X1, Y1 FOR TESTING -
C        IS INTERVAL OF LOCAL OPTIMUM FOUND
C
  160 X1(3) = (THREE*X1(4)-X1(5))/TWO
      X1(1) = X1(3)
      X1(2) = X1(3)
      X1(N1+4) = (THREE*X1(N1+3)-X1(N1+2))/TWO
      X1(N1+5) = X1(N1+4)
      X1(N1+6) = X1(N1+4)
      NR = N1 + 3
      DO 170 N3 = 1,3
        Y1(N3) = Y1(5)
        NR = NR + 1
        Y1(NR) = Y1(N1+2)
  170 CONTINUE
      N7 = MAX(1,J-5)
      N8 = MIN(N1,J+1)
C
C        TEST IF INTERVAL OF LOCAL OPTIMUM FOUND
C
      DO 180 N3 = N7,N8
        IF (Y1(N3).GE.Y1(N3+1) .AND. Y1(N3+1).GE.Y1(N3+2) .AND.
     +        Y1(N3+2).GT.Y1(N3+3) .AND. Y1(N3+3).LT.Y1(N3+4) .AND.
     +        Y1(N3+4).LE.Y1(N3+5) .AND. Y1(N3+5).LE.Y1(N3+6) .AND.
     +        (X1(N3+6)-X1(N3)).LE. (EB/FIVE)) GO TO 200
  180 CONTINUE
      IF (N4.EQ.2) GO TO 60
C
C        IF INTERVAL OF LOCAL OPTIMUM IS FOUND GO TO
C        LOCAL OPTIMIZATION, ELSE BAYESIAN STEP
C
      N7 = J + 1
      DO 190 N3 = J,N7
        XJ1 = X1(N3)
        XJ2 = X1(N3+1)
        YJ1 = Y1(N3)
        YJ2 = Y1(N3+1)
        IF ((XJ2-XJ1).LT.EB2 .AND. Y2(N3-1).GT.ZERO .AND.
     +        YJ1.LT.YJ2) GO TO 185
        TL = (XJ2-XJ1)*SMALL
        Y2(N3) = LOCALM(XJ1,XJ2,EPSL,TL,YL,X2(N3))
        GO TO 190

  185   Y2(N3) = ONE
  190 CONTINUE
      IF (Y2(N7+1).GT.ZERO .AND. Y1(N7).GT.Y1(N7+1) .AND.
     +    X1(N7+1)-X1(N7).LT.EB2) Y2(N7) = ONE
      IF (Y2(N7).GT.ZERO .AND. Y1(N7-1).GT.Y1(N7) .AND.
     +    X1(N7)-X1(N7-1).LT.EB2) Y2(N7-1) = ONE
      GO TO 80
C
C        START OF LOCAL OPTIMIZATION
C
  200 N6 = N3 + 5
      DO 210 N7 = N3,N6
  210 Y2(N7) = ONE
      N3 = N3 + 3
      N5 = N5 + 1
      N7 = 1
  220 EE1 = DBLE(Y1(N3-1)-Y1(N3))
      EE2 = DBLE(X1(N3-1)-X1(N3))
      IF (ABS(EE1).LE.DM17 .OR. ABS(EE2).LE.DM17) GO TO 60
      EE3 = DBLE(TWO* (Y1(N3-1)-Y1(N3+1)))/EE1 -
     +      DBLE(TWO* (X1(N3-1)-X1(N3+1)))/EE2
      IF (ABS(EE3).LE.DM17) GO TO 60
      XA = (DBLE(X1(N3-1)+X1(N3))*DBLE(Y1(N3-1)-Y1(N3+1))/EE1-
     +     DBLE(X1(N3-1)-X1(N3+1))*DBLE(X1(N3-1)+X1(N3+1))/EE2)/EE3
      EE1 = DBLE(X1(N3)-XA)
      EE2 = DBLE(X1(N3-1)-XA)
      EE3 = DBLE(Y1(N3-1)-Y1(N3))
      IF (ABS(EE2).GE.DM17) GO TO 230
      P1 = Y1(N3-1)
      GO TO 250

  230 EE1 = ONE - (EE1/EE2)**2
      IF (ABS(EE1).GE.DM17) GO TO 240
      IF (EE3.LE.ABS(EE1)*1.0D6) GO TO 240
      P1 = Y1(N3-1) - 1.0E5
      GO TO 250

  240 EE3 = EE3/EE1
      P1 = Y1(N3-1) - EE3
  250 YA = F(XA)*CM
      P3 = ABS(XA-X1(N7))
      IF (YA.GE.YM) GO TO 260
      YM = YA
      XM = XA
  260 N1 = N1 + 1
      N2 = N1 + 2
      YH = YM
      N7 = N3
      IF (XA.GT.X1(N3)) N7 = N7 + 1
      N6 = N1 - N7 + 3
      N4 = N1 + 3
      DO 290 N8 = 1,N6
        N4 = N4 - 1
        X1(N4+1) = X1(N4)
        Y1(N4+1) = Y1(N4)
        X2(N4) = X2(N4-1)
        Y2(N4) = Y2(N4-1)
  290 CONTINUE
      X1(N7) = XA
      Y1(N7) = YA
      Y2(N7) = ONE
      Y2(N7-1) = ONE
      X1(3) = (THREE*X1(4)-X1(5))/TWO
      X1(1) = X1(3)
      X1(2) = X1(3)
      X1(N1+4) = (THREE*X1(N1+3)-X1(N1+2))/TWO
      X1(N1+5) = X1(N1+4)
      X1(N1+6) = X1(N1+4)
      NR = N1 + 3
      DO 300 N8 = 1,3
        Y1(N8) = Y1(5)
        NR = NR + 1
        Y1(NR) = Y1(N1+2)
  300 CONTINUE
      IF (N1.LT.N) GO TO 305
      IFAULT = 4
      GO TO 310

  305 P2 = ABS(P1-YA)
      IF (P2.LE.E1 .AND. P3.LE.E2 .OR. P2.LE.CM17 .OR.
     +    P3.LE.CM17) GO TO 60
      IF (N3.EQ.N7 .AND. Y1(N7).GT.Y1(N7+1) .OR.
     +    N3+1.EQ.N7 .AND. Y1(N7).LT.Y1(N7-1)) N3 = N3 + 1
      GO TO 220

  310 YM = YM*CM
      RETURN

      END
C
C
      FUNCTION YL(C)
C
C        ALGORITHM AS 133.1 APPL. STATIST. (1978) VOL.27, NO.3
C
C        -YL IS THE MEAN OF THE IMPROVEMENT
C
      REAL VSMALL,ZERO
      COMMON /EXTR1/A,EB,YH,XJ1,XJ2,YJ1,YJ2
      DATA VSMALL/1.0E-17/,ZERO/0.0/
C
      P1 = EB*VSMALL
      IF (C-XJ1.LE.P1 .OR. XJ2-C.LE.P1) GO TO 1
      P1 = XJ2 - XJ1
      P2 = A*SQRT((C-XJ1)* (XJ2-C)/P1)
      P1 = (YJ1* (XJ2-C)+YJ2* (C-XJ1))/P1 - YH
      P3 = P1/P2
      YL = P1*0.65*EXP(-0.443* (0.75+P3)**2) -
     +     P2*0.3989*EXP(- (P3**2)/2.0)
      RETURN

    1 YL = ZERO
      RETURN

      END
C
C
      REAL FUNCTION LOCALM(A, B, EPS, T, F, X)
C
C     Entered from pages 188-190 of 'Algorithms for minimization without
C     derivatives' by Richard P. Brent, Prentice-Hall, 1973
C     Comments added from the Algol version on pages 79-80.
C
      REAL A, B, EPS, T, F, X
      EXTERNAL F
C
C     Local variables
C
      REAL CONST
      REAL SA, SB, D, E, M, P, Q, R, TOL, T2, U, V, W, FU, FV, FW, FX
C
C     CONST = (3 - sqrt(5))/2
C
      DATA CONST/0.381966/
C
      SA = A
      SB = B
      X = SA + CONST*(SB - SA)
      W = X
      V = W
      E = 0.0
      FX = F(X)
      FW = FX
      FV = FW
C
C     Main loop
C
   10 M = 0.5*(SA + SB)
      TOL = EPS*ABS(X) + T
      T2 = 2.0*TOL
C
C     Check stopping criterion
C
      IF (ABS(X-M) .LE. T2-0.5*(SB-SA)) GO TO 190
      R = 0.0
      Q = R
      P = Q
      IF (ABS(E) .LE. TOL) GO TO 40
C
C     Fit parabola
C
      R = (X - W)*(FX - FV)
      Q = (X - V)*(FX - FW)
      P = (X - V)*Q - (X - W)*R
      Q = 2.0*(Q - R)
      IF (Q .LE. 0.0) GO TO 20
      P = -P
      GO TO 30
   20 Q = -Q
   30 R = E
      E = D
   40 IF (ABS(P) .GE. ABS(0.5*Q*R)) GO TO 60
      IF ((P .LE. Q*(SA-X)) .OR. (P .GE. Q*(SB-X))) GO TO 60
C
C     A parabolic interpolation step
C
      D = P/Q
      U = X + D
C
C     F must not be evaluated too close to A or B
C
      IF ((U-SA .GE. T2) .AND. (SB-U .GE. T2)) GO TO 90
      IF (X .GE. M) GO TO 50
      D = TOL
      GO TO 90
   50 D = -TOL
      GO TO 90
C
C     A golden section step
C
   60 IF (X .GE. M) GO TO 70
      E = SB - X
      GO TO 80
   70 E = SA - X
   80 D = CONST*E
C
C     F must not be evaluated too close to X
C
   90 IF (ABS(D) .LT. TOL) GO TO 100
      U = X + D
      GO TO 120
  100 IF (D .LE. 0.0) GO TO 110
      U = X + TOL
      GO TO 120
  110 U = X - TOL
  120 FU = F(U)
C
C     Update A, B, V, W and X
C
      IF (FU .GT. FX) GO TO 150
      IF (U .GE. X) GO TO 130
      SB = X
      GO TO 140
  130 SA = X
  140 V = W
      FV = FW
      W = X
      FW = FX
      X = U
      FX = FU
      GO TO 10
  150 IF (U .GE. X) GO TO 160
      SA = U
      GO TO 170
  160 SB = U
  170 IF ((FU .GT. FW) .AND. (W .NE. X)) GO TO 180
      V = W
      FV = FW
      W = U
      FW = FU
      GO TO 10
  180 IF ((FU .GT. FV) .AND. (V .NE. X) .AND. (V .NE. W)) GO TO 10
      V = U
      FV = FU
      GO TO 10
C
  190 LOCALM = FX
      RETURN
      END
