      REAL FUNCTION QF(ALB,ANC,N,IRR,SIGMA,CC,LIM1,ACC,ITH,TRACE,IFAULT)
C
C     ALGORITHM AS 155  APPL. STATIST. (1980) VOL.29, NO.3
C
C     Distribution of a linear combination of non-central chi-squared
C     random variables.
C
      INTEGER IRR,LIM1,IFAULT
      REAL SIGMA,CC,ACC
      REAL TRACE(7),ALB(IRR),ANC(IRR)
      INTEGER N(IRR),ITH(IRR)
      INTEGER J,NJ,NT,NTM
      REAL ACC1,ALMX,XLIM,XNT,XNTM
      REAL UTX,TAUSQ,SD,AINTV,AINTV1,X,UP,UN,D1,D2,ALJ,ANCJ
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,ZERO,HALF,ONE,TWO,FOUR,
     1  SIXTN,FOURP5,PT07,PT2,QUART,TEN,PT33,PT67,PT75,ONEPT5,THREE,
     2  ONEPT1
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA ZERO /0.0/,HALF/0.5/,ONE/1.0/,TWO/2.0/,FOUR/4.0/,SIXTN/16.0/,
     1   FOURP5/4.5/,PT07/0.07/,PT2/0.2/,QUART/0.25/,TEN/10.0/,
     2   PT33/0.33/,PT67/0.67/,PT75/0.75/,ONEPT5/1.5/,THREE/3.0/,
     3   ONEPT1/1.1/
C
      IROUND(X) = INT(X + SIGN(HALF,X))
C
C     Setting constants in COMMON.   ALN28 = ln(2) / 8.
C
      PI = 3.14159265358979
      ALN28 = 0.0866
C
      C = CC
      IR = IRR
      LIM = LIM1
      DO 10 J = 1,7
	TRACE(J) = ZERO
 10   CONTINUE
      IFAULT = 0
      ICOUNT = 0
      AINTL = ZERO
      ERSM = ZERO
      QF = -ONE
      ACC1 = ACC
      NDTSRT = .TRUE.
      FAIL = .FALSE.
C
C     Find mean, sd, max & min of ALB.
C     Check that parameter values are valid.
C
      XLIM = LIM
      SIGSQ = SIGMA**2
      SD = SIGSQ
      ALMAX = ZERO
      ALMIN = ZERO
      AMEAN = ZERO
      J = 1
 20   IF (.NOT.(J.LE.IR)) GO TO 60
      NJ = N(J)
      ALJ = ALB(J)
      ANCJ = ANC(J)
      IF (.NOT.(NJ.LT.0.OR.ANCJ.LT.ZERO)) GO TO 30
      IFAULT = 3
      GO TO 260
 30   SD = SD + ALJ**2*(2*NJ + FOUR*ANCJ)
      AMEAN = AMEAN + ALJ*(NJ + ANCJ)
      IF (.NOT.(ALMAX.LT.ALJ)) GO TO 40
      ALMAX = ALJ
      GO TO 50
 40   IF (.NOT.(ALMIN.GT.ALJ)) GO TO 50
      ALMIN = ALJ
 50   J = J + 1
      GO TO 20
 60   IF (.NOT.(SD.EQ.ZERO)) GO TO 80
      IF (.NOT.(C.GT.ZERO)) GO TO 70
      QF = ONE
      GO TO 260
 70   QF = ZERO
      GO TO 260
 80   IF (.NOT.(ALMIN.EQ.ZERO.AND.ALMAX.EQ.ZERO.AND.SIGMA.EQ.ZERO))
     1  GO TO 90
      IFAULT = 3
      GO TO 260
 90   SD = SQRT(SD)
      IF (.NOT.(ALMAX.LT.-ALMIN)) GO TO 100
      ALMX = -ALMIN
      GO TO 110
 100  ALMX = ALMAX
C
C     Starting values for FINDU * CTFF.
C
 110  UTX = SIXTN/SD
      UP = FOURP5/SD
      UN = -UP
C
C     Truncation point with no convergence factor.
C
      CALL FINDU (N,ALB,ANC,UTX,HALF*ACC1)
C
C     Does convergence factor help ?
C
      IF (.NOT.(C.NE.ZERO.AND.ALMX.GT.PT07*SD)) GO TO 130
      TAUSQ = QUART*ACC1/CFE(N,ALB,ANC,ITH,C)
      IF (.NOT.(FAIL)) GO TO 120
      FAIL = .FALSE.
      GO TO 130
 120  IF (.NOT.(TRUNCN(N,ALB,ANC,UTX,TAUSQ).LT.PT2*ACC1)) GO TO 130
      SIGSQ = SIGSQ + TAUSQ
      CALL FINDU (N,ALB,ANC,UTX,QUART*ACC1)
      TRACE(6) = SQRT(TAUSQ)
 130  TRACE(5) = UTX
      ACC1 = HALF*ACC1
C
C     Find 'range' of distribution, quit if outside of this.
C
 140  D1 = CTFF(N,ALB,ANC,ACC1,UP) - C
      IF (.NOT.(D1.LT.ZERO)) GO TO 150
      QF = ONE
      GO TO 260
 150  D2 = C - CTFF(N,ALB,ANC,ACC1,UN)
      IF (.NOT.(D2.LT.ZERO)) GO TO 160
      QF = ZERO
      GO TO 260
C
C     Find integration interval.
C
 160  IF (.NOT.(D1.GT.D2)) GO TO 170
      AINTV = D1
      GO TO 180
 170  AINTV = D2
 180  AINTV = TWO*PI/AINTV
C
C     Calculate number of terms required for main & auxiliary
C     integrations.
C
      XNT = UTX/AINTV
      XNTM = THREE/SQRT(ACC1)
      IF (.NOT.(XNT.GT.XNTM*ONEPT5)) GO TO 220
      IF (.NOT.(XNTM.GT.XLIM)) GO TO 190
      IFAULT = 1
      GO TO 260
C
C     Parameters for auxiliary integration.
C
 190  NTM = IROUND(XNTM)
      AINTV1 = UTX/XNTM
      X = TWO*PI/AINTV1
      IF (.NOT.(X.LE.ABS(C))) GO TO 200
      GO TO 220
C
C     Calculate convergence factor.
C
 200  TAUSQ = CFE(N,ALB,ANC,ITH,C - X) + CFE(N,ALB,ANC,ITH,C + X)
      TAUSQ = PT33*ACC1/(ONEPT1*TAUSQ)
      IF (.NOT.(FAIL)) GO TO 210
      GO TO 220
 210  ACC1 = PT67*ACC1
C
C     Auxiliary integration.
C
      CALL INTEGR (N,ALB,ANC,NTM,AINTV1,TAUSQ,.FALSE.)
      XLIM = XLIM - XNTM
      SIGSQ = SIGSQ + TAUSQ
      TRACE(3) = TRACE(3) + 1
      TRACE(2) = TRACE(2) + NTM + 1
C
C     Find truncation point with new convergence factor.
C
      CALL FINDU (N,ALB,ANC,UTX,QUART*ACC1)
      ACC1 = PT75*ACC1
      GO TO 140
C
C     Main integration.
C
 220  TRACE(4) = AINTV
      IF (.NOT.(XNT.GT.XLIM)) GO TO 230
      IFAULT = 1
      GO TO 260
 230  NT = IROUND(XNT)
      CALL INTEGR (N,ALB,ANC,NT,AINTV,ZERO,.TRUE.)
      TRACE(3) = TRACE(3) + 1
      TRACE(2) = TRACE(2) + NT + 1
      QF = HALF - AINTL
      TRACE(1) = ERSM
      UP = ERSM
C
C     Test whether round-off error could be significant.
C     Allow for radix 8 or 16 machines.
C
      X = UP + ACC/TEN
      J = 1
 240  IF (.NOT.(J.LE.8)) GO TO 260
      IF (.NOT.(J*X.EQ.J*UP)) GO TO 250
      IFAULT = 2
 250  J = J*2
      GO TO 240
 260  TRACE(7) = ICOUNT
      RETURN
      END
C
      SUBROUTINE COUNTR
C
C     Count number of calls to ERRBD, TRUNCN & CFE.
C
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
C 
      ICOUNT = ICOUNT + 1
      IF (.NOT.(ICOUNT.GT.LIM)) GO TO 20
      WRITE (6,10)
 10   FORMAT (' qf: cannot locate integration parameters'/)
      STOP
 20   RETURN
      END
C
      REAL FUNCTION ALOG1(X,FIRST)
C
C     If FIRST then return ln(1 + x) else ln(1 + x) - x.
C
      REAL X
      LOGICAL FIRST
      REAL S,S1,TERM,Y,AK,PT1,ONE,TWO,THREE
      DATA PT1/0.1/,ONE/1.0/,TWO/2.0/,THREE/3.0/
C 
      F1(I) = S + TERM/AK
C 
      IF (.NOT.(ABS(X).GT.PT1)) GO TO 20
      IF (.NOT.(FIRST)) GO TO 10
      ALOG1 = LOG(ONE + X)
      GO TO 70
 10   ALOG1 = LOG(ONE + X) - X
      GO TO 70
 20   Y = X/(TWO + X)
      TERM = TWO*Y**3
      AK = THREE
      IF (.NOT.(FIRST)) GO TO 30
      S = TWO
      GO TO 40
 30   S = -X
 40   S = S*Y
      Y = Y**2
      S1 = F1(0)
 50   IF (.NOT.(S1.NE.S)) GO TO 60
      AK = AK + TWO
      TERM = TERM*Y
      S = S1
      S1 = F1(0)
      GO TO 50
 60   ALOG1 = S
 70   RETURN
      END
C
      REAL FUNCTION EXP1(X)
      REAL X
      REAL ZERO,NEG50
      DATA ZERO/0.0/,NEG50/-50.0/
C 
      IF (.NOT.(X.LT.NEG50)) GO TO 10
      EXP1 = ZERO
      GO TO 20
 10   EXP1 = EXP(X)
 20   RETURN
      END
C
      SUBROUTINE ORDER (ALB,ITH)
C
C     Find order of absolute values of ALB.
C
      REAL ALB(*)
      INTEGER ITH(*)
      INTEGER J,K,K1,ITHK
      REAL ALJ
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
C 
      J = 1
 10   IF (.NOT.(J.LE.IR)) GO TO 60
      ALJ = ABS(ALB(J))
      K = J - 1
 20   IF (.NOT.(K.GT.0)) GO TO 40
      ITHK = ITH(K)
      K1 = K + 1
      IF (.NOT.(ALJ.GT.ABS(ALB(ITHK)))) GO TO 50
      ITH(K1) = ITHK
      GO TO 30
 30   K = K - 1
      GO TO 20
 40   K = 0
      K1 = 1
 50   ITH(K1) = J
      J = J + 1
      GO TO 10
 60   NDTSRT = .FALSE.
      RETURN
      END
C
      REAL FUNCTION ERRBD(N,ALB,ANC,UU,CX)
C
C     Find bound on tail probability using mgf.
C     Cut-off point returned to CX.
C
      REAL U,UU,CX
      INTEGER N(*)
      REAL ALB(*),ANC(*)
      REAL SUM1,ALJ,ANCJ,X,Y,CONST
      INTEGER J,NJ
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,HALF,ONE,TWO
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA HALF/0.5/,ONE/1.0/,TWO/2.0/
C 
      CALL COUNTR
      U = UU
      CONST = U*SIGSQ
      SUM1 = U*CONST
      U = TWO*U
      J = IR
 10   IF (.NOT.(J.GT.0)) GO TO 20
      NJ = N(J)
      ALJ = ALB(J)
      ANCJ = ANC(J)
      X = U*ALJ
      Y = ONE - X
      CONST = CONST + ALJ*(ANCJ/Y + NJ)/Y
      SUM1 = SUM1 + ANCJ*(X/Y)**2
      SUM1 = SUM1 + NJ*(X**2/Y + ALOG1(-X,.FALSE.))
      J = J - 1
      GO TO 10
 20   ERRBD = EXP1(-HALF*SUM1)
      CX = CONST
      RETURN
      END
C
      REAL FUNCTION CTFF(N,ALB,ANC,ACCX,UPN)
C
C     Find CTFF so that P(QF > CTFF) < ACCX if UPN > 0;
C     P(QF < CTFF) < ACCX otherwise.
C
      REAL ACCX,UPN
      INTEGER N(*)
      REAL ALB(*),ANC(*)
      REAL U1,U2,U,RB,CONST,C1,C2
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,ZERO,ONE,TWO
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA ZERO/0.0/,ONE/1.0/,TWO/2.0/
C 
      F1(I) = U2/(ONE + U2*RB)
      F2(I) = (C1 - AMEAN)/(C2 - AMEAN)
C 
      U2 = UPN
      U1 = ZERO
      C1 = AMEAN
      IF (.NOT.(U2.GT.ZERO)) GO TO 10
      RB = ALMAX
      GO TO 20
 10   RB = ALMIN
 20   RB = TWO*RB
      U = F1(0)
 30   IF (.NOT.(ERRBD(N,ALB,ANC,U,C2).GT.ACCX)) GO TO 40
      U1 = U2
      C1 = C2
      U2 = TWO*U2
      U = F1(0)
      GO TO 30
 40   U = F2(0)
 50   IF (.NOT.(U.LT.0.9)) GO TO 80
      U = (U1 + U2)/TWO
      IF (.NOT.(ERRBD(N,ALB,ANC,U/(ONE + U*RB),CONST).GT.ACCX)) GO TO 60
      U1 = U
      C1 = CONST
      GO TO 70
 60   U2 = U
      C2 = CONST
 70   U = F2(0)
      GO TO 50
 80   CTFF = C2
      UPN = U2
      RETURN 
      END
C
      REAL FUNCTION TRUNCN(N,ALB,ANC,UU,TAUSQ)
C
C     Bound integration error due to truncation at U.
C
      REAL U,UU,TAUSQ
      INTEGER N(*)
      REAL ALB(*),ANC(*)
      REAL SUM1,SUM2,PROD1,PROD2,PROD3,ALJ,ANCJ,X,Y,ERR1,ERR2
      INTEGER J,NJ,NS
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,ZERO,QUART,HALF,ONE,TWO
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA ZERO/0.0/,QUART/0.25/,HALF/0.5/,ONE/1.0/,TWO/2.0/
C 
      CALL COUNTR
      U = UU
      SUM1 = ZERO
      PROD2 = ZERO
      PROD3 = ZERO
      NS = 0
      SUM2 = (SIGSQ + TAUSQ)*U**2
      PROD1 = TWO*SUM2
      U = TWO*U
      J = 1
 10   IF (.NOT.(J.LE.IR)) GO TO 40
      ALJ = ALB(J)
      ANCJ = ANC(J)
      NJ = N(J)
      X = (U*ALJ)**2
      SUM1 = SUM1 + ANCJ*X/(ONE + X)
      IF (.NOT.(X.GT.ONE)) GO TO 20
      PROD2 = PROD2 + NJ*LOG(X)
      PROD3 = PROD3 + NJ*ALOG1(X,.TRUE.)
      NS = NS + NJ
      GO TO 30
 20   PROD1 = PROD1 + NJ*ALOG1(X,.TRUE.)
 30   J = J + 1
      GO TO 10
 40   SUM1 = HALF*SUM1
      PROD2 = PROD1 + PROD2
      PROD3 = PROD1 + PROD3
      X = EXP1(-SUM1 - QUART*PROD2)/PI
      Y = EXP1(-SUM1 - QUART*PROD3)/PI
      IF (.NOT.(NS.EQ.0)) GO TO 50
      ERR1 = ONE
      GO TO 60
 50   ERR1 = X*TWO/NS
 60   IF (.NOT.(PROD3.GT.ONE)) GO TO 70
      ERR2 = 2.5*Y
      GO TO 80
 70   ERR2 = ONE
 80   IF (.NOT.(ERR2.LT.ERR1)) GO TO 90
      ERR1 = ERR2
 90   X = HALF*SUM2
      IF (.NOT.(X.LE.Y)) GO TO 100
      ERR2 = ONE
      GO TO 110
 100  ERR2 = Y/X
 110  IF (.NOT.(ERR1.LT.ERR2)) GO TO 120
      TRUNCN = ERR1
      GO TO 130
 120  TRUNCN = ERR2
 130  RETURN
      END
C
      SUBROUTINE FINDU (N,ALB,ANC,UTX,ACCX)
C
C     Find U such that TRUNCN(U) < ACCX & TRUNCN(U / 1.2) > ACCX.
C
      REAL UTX,ACCX
      INTEGER N(*)
      REAL ALB(*),ANC(*)
      REAL U,UT
      REAL DIVIS(4)
      INTEGER I
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,FOUR
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA DIVIS/2.0,1.4,1.2,1.1/, FOUR/4.0/
C 
      UT = UTX
      U = UT/FOUR
      IF (.NOT.(TRUNCN(N,ALB,ANC,U,ZERO).GT.ACCX)) GO TO 20
      U = UT
 10   IF (.NOT.(TRUNCN(N,ALB,ANC,U,ZERO).GT.ACCX)) GO TO 40
      UT = UT*FOUR
      U = UT
      GO TO 10
 20   UT = U
      U = U/FOUR
 30   IF (.NOT.(TRUNCN(N,ALB,ANC,U,ZERO).LE.ACCX)) GO TO 40
      UT = U
      U = U/FOUR
      GO TO 30
 40   DO 50 I = 1,4
	U = UT/DIVIS(I)
	IF (.NOT.(TRUNCN(N,ALB,ANC,U,ZERO).LE.ACCX)) GO TO 50
	UT = U
 50   CONTINUE
      UTX = UT
      RETURN 
      END
C
      SUBROUTINE INTEGR (N,ALB,ANC,NTERM,AINTRV,TAUSQ,MAIN)
C
C     Carry out integration with NTERM terms, at interval AINTRV.
C     If not MAIN then multiply integrand by  1 - exp(-0.5 * TAUSQ *
C     U**2).
C
      INTEGER NTERM
      REAL AINTRV,TAUSQ
      LOGICAL MAIN
      INTEGER N(*)
      REAL ALB(*),ANC(*)
      REAL AINPI,U,SUM1,SUM2,SUM3,X,Y,Z
      INTEGER K,J,NJ
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,QUART,HALF,ONE,TWO
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA QUART/0.25/,HALF/0.5/,ONE/1.0/,TWO/2.0/
C 
      AINPI = AINTRV/PI
      K = NTERM
 10   IF (.NOT.(K.GE.0)) GO TO 50
      U = (K + HALF)*AINTRV
      SUM1 = -TWO*U*C
      SUM2 = ABS(SUM1)
      SUM3 = -HALF*SIGSQ*U**2
      J = IR
 20   IF (.NOT.(J.GT.0)) GO TO 30
      NJ = N(J)
      X = TWO*ALB(J)*U
      Y = X**2
      SUM3 = SUM3 - QUART*NJ*ALOG1(Y,.TRUE.)
      Y = ANC(J)*X/(ONE + Y)
      Z = NJ*ATAN(X) + Y
      SUM1 = SUM1 + Z
      SUM2 = SUM2 + ABS(Z)
      SUM3 = SUM3 - HALF*X*Y
      J = J - 1
      GO TO 20
 30   X = AINPI*EXP1(SUM3)/U
      IF (.NOT.(.NOT.MAIN)) GO TO 40
      X = X*(ONE - EXP1(-HALF*TAUSQ*U**2))
 40   SUM1 = SIN(HALF*SUM1)*X
      SUM2 = HALF*SUM2*X
      AINTL = AINTL + SUM1
      ERSM = ERSM + SUM2
      K = K - 1
      GO TO 10
 50   RETURN
      END
C
      REAL FUNCTION CFE(N,ALB,ANC,ITH,X)
C
C     Coefficient of TAUSQ in error when convergence factor of
C     exp(-0.5 * TAUSQ * U**2) is used when df is evaluated at X.
C
      REAL X
      INTEGER N(*),ITH(*)
      REAL ALB(*),ANC(*)
      REAL AXL,AXL1,AXL2,SXL,SUM1,ALJ
      INTEGER J,K,IT,ITK
      DOUBLE PRECISION AINTL,ERSM
      REAL PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,ZERO,ONE,FOUR,HUNDRD
      INTEGER ICOUNT,IR,LIM
      LOGICAL NDTSRT,FAIL
      COMMON /QFCOM/ AINTL,ERSM,PI,ALN28,SIGSQ,ALMAX,ALMIN,AMEAN,C,
     1  ICOUNT,IR,NDTSRT,FAIL,LIM
      DATA ZERO/0.0/,ONE/1.0/,FOUR/4.0/,HUNDRD/100.0/
C 
      CALL COUNTR
      IF (.NOT.(NDTSRT)) GO TO 10
      CALL ORDER (ALB,ITH)
 10   AXL = ABS(X)
      SXL = SIGN(ONE,X)
      SUM1 = ZERO
      J = IR
 20   IF (.NOT.(J.GT.0)) GO TO 70
      IT = ITH(J)
      IF (.NOT.(ALB(IT)*SXL.GT.ZERO)) GO TO 60
      ALJ = ABS(ALB(IT))
      AXL1 = AXL - ALJ*(N(IT) + ANC(IT))
      AXL2 = ALJ/ALN28
      IF (.NOT.(AXL1.GT.AXL2)) GO TO 30
      AXL = AXL1
      GO TO 60
 30   IF (.NOT.(AXL.GT.AXL2)) GO TO 40
      AXL = AXL2
 40   SUM1 = (AXL - AXL1)/ALJ
      K = J - 1
 50   IF (.NOT.(K.GT.0)) GO TO 70
      ITK = ITH(K)
      SUM1 = SUM1 + (N(ITK) + ANC(ITK))
      K = K - 1
      GO TO 50
 60   J = J - 1
      GO TO 20
 70   IF (.NOT.(SUM1.GT.HUNDRD)) GO TO 80
      CFE = ONE
      FAIL = .TRUE.
      GO TO 90
 80   CFE = 2**(SUM1/FOUR)/(PI*AXL**2)
 90   RETURN
      END
