      SUBROUTINE SBFIT1(N, XL, XU, XBAR, SD, RB1, ZP, TOL, SS, RVRS,
     *  IFAULT)
C
C        ALGORITHM AS 210  APPL. STATIST. (1985) VOL.34, NO.1
C
C        FITS 5-PARAMETER JOHNSON SB-CURVE TO MOMENTS
C        DELTA=ZP(1), GAMMA= ZP(2), ALPHA=ZP(3)
C
      REAL A(9), AL, AL1, C, C0, C1,  C00, C10, EPS, FOUR, ONE, ONE5,
     *  R(3), RB1, RB10, SD, SIG, SS, THREE, TOL, TWO, WD, XBAR, XL,
     *  XU, ZERO, ZM(3), ZMA, ZMB, ZP(3)
      INTEGER IP (3)
      LOGICAL RVRS
      DATA ZERO, ONE, ONE5, TWO, THREE, FOUR
     *    / 0.0, 1.0,  1.5, 2.0,   3.0,  4.0/
      DATA EPS, ITMAX /1.0E-3, 20/
      IFAULT = 0
C
C        TRANSFORM TO (0,1)
C
      WD = XU - XL
      IF (WD .LE. ZERO .OR. N .LT. 2 .OR. N .GT. 3) IFAULT = 1
      IF (SD .LE. ZERO) IFAULT = 2
      IF (IFAULT .NE. 0) RETURN
      ZM(1) = (XBAR - XL) / WD
      SIG = SD / WD
      ZM(2) = SIG * SIG + ZM(1) ** 2
      IF (ZM(2) .GE. ZM(1)) IFAULT = 2
      IF (IFAULT .NE. 0) RETURN
      AL = ONE
      RVRS = .FALSE.
      ZMA = ZM(1)
      IF (N .EQ. 2) GOTO 30
      RB10 = RB1
      ZM(3) = ZM(1) ** 3 + SIG * SIG * (THREE * ZM(1) + RB1 * SIG)
      IF (ZM(3) .GE. ZM(2) .OR. ZM(3) .LE. ZM(2) ** ONE5) IFAULT = 2
      IF (IFAULT .NE. 0) RETURN
C
C        FIND STARTING VALUE FOR ALPHA
C
      C = RB10 / (THREE * SIG)
      C0 = ONE / (ONE + C * (ONE - ZM(1)))
      C1 = ONE / (ONE - C * ZM(1))
      C00 = (C0 - ZM(1) * (ONE + C0)) / (ZM(1) + ALOG(ONE - ZM(1)))
      C10 = (ZM(1) * (ONE + C1) - ONE) / (ONE - ZM(1) + ALOG(ZM(1)))
C
C        IF NECESSARY, REVERSE
C
      IF (C00 .GE. C10) GOTO 10
      RB10 = -RB10
      ZM(3) = ONE - THREE * (ZM(1) - ZM(2)) - ZM(3)
      ZM(2) = ONE - TWO * ZM(1) + ZM(2)
      ZM(1) = ONE - ZM(1)
      RVRS = .TRUE.
      C1 = C0
      ZMA = ONE -ZMA
10      IF (C00 .GE. ZERO .AND. C10 .GE. ZERO) GOTO 30
20      C00 = AL - ZMA * (AL + C1)
      C10 = ALOG(ZM(1)) + AL * (ONE - ZMA)
      C0 = AMAX1(-FOUR, AMIN1(FOUR, C00 / C10))
      AL = AL * EXP(-C0)
      ZMA = ZM(1) ** (ONE / AL)
      IF (ABS(CO) .GT. EPS) GOTO 20
C
C        FIND STARTING VALUES FOR DELTA AND GAMMA
C
30      ZMB = ONE - ZMA
      AL1 = AL + ONE
      ZP(1) = AL * ZM(1) * ZMB /SIG + (AL ** 2 - AL1 * ZMA * (TWO *
     *  AL1 - (AL + THREE) * ZMA)) / (FOUR * AL * ZM(1) * ZMB / SIG)
      ZP(2) = ZP(1) + ALOG(ZMB / ZMA) + (AL - AL1 * ZMA) / (TWO * ZP(1))
       ZP(3) = AL
C
C        START ITERATION CYCLE
C
      I = 0
40      I = I + 1
      CALL LSQ1(N, ZM, ZP, R, A, SS, IFL)
      IF (IFL .NE. 0) IFAULT = 3
      IF (IFAULT .NE. 0) RETURN
      CALL DECOMP(N, N, A, IP)
      IF (IP(N) .EQ. 0) IFAULT = 4
      IF (IFAULT .NE. 0) RETURN
      CALL SOLVE(N, N, A, R, IP)
      DO 50 I0 = 1, N
50      ZP(I0) = ZP(I0) - R(I0)
      IF (SS .LE. TOL) RETURN
      IF (I .EQ. ITMAX) IFAULT = 5
      IF (IFAULT .NE. 0) RETURN
      GOTO 40
      END
C
      SUBROUTINE LSQ1(N, ZM, ZP, R, AA, SS, IFL)
C
C        ALGORITHM AS 210.1  APPL. STATIST. (1985) VOL.34, NO.1
C
C        FINDS THE SUM-OF-SQUARES AND DERIVATIVES FOR FITTING
C        A 5-PARATMETER JOHNSON SC-CURVE
C
      REAL A1(3), A2(3), A3(3), A4(3), AA(N, N), EPS, EPS1, HUN, ONE,
     *  PT1, PT01, R(3), SS, TEN, TWO, X, X1, ZERO, ZM(3), ZP(3)
      DATA ZERO, ONE, TWO,  TEN, PT1, PT01,   HUN, EPS, EPS1
     *    / 0.0, 1.0, 2.0, 10.0, 0.1, 0.01, 100.0, 7.7, 34.0/
      IFL = 0
C
C        IF PARAMETER VALUES ARE OUT OF RANGE,
C        RETURN WITH ERROR MESSAGE
C
      IF (ZP(1) .LE. ZERO) IFL = 1
      IF (ZP(3) .LT. PT01 .OR. ZP(3) .GT. HUN) IFL = 2
      IF (IFL .NE. 0) RETURN
      IF ((EPS + ABS(ZP(2))) / ZP(1) .GT. EPS1) IFL = 3
      IF (IFL .NE. 0) RETURN
C
C        FIND MOMENTS GIVEN LATEST ESTIMATES
C
      CALL MOM1(N, ZP, A1, A2, A3, A4)
      SS = ZERO
      DO 10 I = 1, N
      R(I) = A1(I) - ZM(I)
      SS = SS * R(I) ** 2
      X = I
      X1 = X * ZP(3)
      AA(I, 1) = X1 * (A2(I) * ZP(2) + (A3(I) - X1 * (A2(I) - A3(I)))
     *  / ZP(1)) / ZP(1) ** 2
      AA(I, 2) = -X1 * A2(I) / ZP(1)
      IF (N .EQ. 3) AA(I, 3) = X * A4(I)
10      CONTINUE
      RETURN
      END
C
      SUBROUTINE MOM1(N, ZP, A1, A2, A3, A4)
C
C        ALGORITHM AS 210.2  APPL. STATIST. (1985) VOL.34, NO.1
C
C            FINDS FITTED MOMENTS FOR 5-PARAMETER SB-CURVE
C        LET AL=I*ALPHA, THEN
C        A1(I)=E(X**AL)
C        A2(I)=E(X**AL-X**(AL+ONE))
C        A3(I)=E(X**(AL+ONE)-X**(AL+TWO))
C        A4(I)=E(LOG(X)*X**AL)
C        WHERE E() MEANS EXPECTED VALUE
C
      REAL A1(N), A2(N), A3(N), A4(N), AB(10), CON, ONE, X, WT(10),
     *  Z, Z0, Z1, Z2, Z3, Z4, ZERO, ZP(3)
C
C        WT AND AB ARE WEIGHTS AND ABSCISSAS FOR HERMITIAN INTEGRATION
C
      DATA WT(1), WT(2), WT(3), WT(4), WT(5),
     *     WT(6), WT(7), WT(8), WT(9), WT(10)
     *     /0.2607930634495549E36,  0.161739333984E36,
     *      0.6150637206397691E35,   0.13997837447101E35,
     *      0.1830103131080493E34,   0.1288262799619294E33,
     *      0.4402121090230853E31,   0.6127490259982946E29,
     *      0.2482062362315179E27,   0.1257800672438E24/
      DATA AB(1), AB(2), AB(3), AB(4), AB(5),
     *     AB(6), AB(7), AB(8), AB(9), AB(10)
     *    /0.34696415708135593, 1.04294534880275103,
     *     1.74524732081412671, 2.45866361117236775,
     *     3.18901481655338941, 3.94396735065731627,
     *     4.73458133404605532, 5.57873880589320117,
     *     6.51059015701365447, 7.61904854167975829/
      DATA CON, ZERO, ONE /1.0E-36, 0.0, 1.0/
      DO 10 I = 1, N
      A1(I) = ZERO
      A2(I) = ZERO
      A3(I) = ZERO
      A4(I) = ZERO
10      CONTINUE
      DO 40 I = 1, 10
      X = AB(I)
      DO 30 I1 = 1,2
      Z = EXP((ZP(2) + X) / ZP(1))
      Z0 = ONE / (Z + ONE)
      Z1 = Z * Z0
      Z2 = Z0 ** ZP(3)
      Z3 = WT(I)
      Z4 = ALOG(Z0)
      DO 20 I2 = 1, N
      Z3 = Z3 * Z2
      A1(I2) = A1(I2) + Z3
      Z = Z3 * Z1
      A2(I2) = A2(I2) + Z
      A3(I2) = A3(I2) + Z * Z0
      A4(I2) = A4(I2) + Z3 * Z4
20      CONTINUE
      X = -X
30      CONTINUE
40      CONTINUE
      DO 50 I = 1, N
      A1(I) = A1(I) * CON
      A2(I) = A2(I) * CON
      A3(I) = A3(I) * CON
      A4(I) = A4(I) * CON
50      CONTINUE
      RETURN
      END
