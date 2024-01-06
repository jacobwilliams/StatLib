      SUBROUTINE LGSTN(ZM, ZP, TOL, IFAULT)
C
C        ALGORITHM AS 208  APPL. STATIST. (1985) VOL.34, NO.1
C
C        FITS LOGISTIC NORMAL FOR 3 DIMENSIONAL CASE, USING
C        MOMENTS OF FIRST AND SECOND ORDER
C
      REAL DJ(5, 5), EPS, ONE, R(5), SS, TOL, TWO, WW, WW0, WW1, WW2,
     *  WW3, ZERO, ZM(5), ZM1, ZM2, ZP(5)
      INTEGER IP(5)
      DATA ZERO, ONE, TWO, EPS, ITMAX / 0.0, 1.0, 2.0, 0.01, 20/
C
C        FIND STARTING VALUES
C
      IFAULT = 0
      N = 5
      R(1) = ZM(3) - ZM(1) ** 2
      R(2) = ZM(5) - ZM(2) ** 2
      R(3) = ZM(4) - ZM(1) * ZM(2)
      ZM1 = ONE - ZM(1)
      ZM2 = ONE - ZM(2)
      WW = ONE - ZM(1) - ZM(2)
      WW0 = ZM1 * ZM2 + ZM(1) * ZM(2)
      WW1 = R(1) * ZM2 ** 2 + R(2) * ZM(1) ** 2 + TWO * R(3) * ZM2 *
     *  ZM(1)
      WW2 = R(1) * ZM(2) * ZM2 + R(2) * ZM1 * ZM(1) + R(3) * WW0
      WW3 = R(1) * ZM(2) ** 2 + R(2) * ZM1 ** 2 + TWO * R(3) * ZM1 *
     *  ZM(2)
C
C        CHECK THAT MOMENTS ARE IN RANGE
C
      IF (AMIN1(R(1), R(2), ZM(1), ZM(2), WW, WW1, WW3) .LE. ZERO)
     *  IFAULT = 1
      IF (IFAULT .NE. 0) RETURN
      WW1 = SQRT(ONE / WW1)
      WW3 = SQRT(ONE / WW3)
      ZP(1) = ALOG(ZM(1) / WW)
      ZP(2) = ALOG(ZM(2) / WW)
      ZP(3) = ZM(1) * WW * WW1
      ZP(4) = ZM(2) * WW * WW3
C
C        FIND INITIAL VALUE FOR ZP(5) THAT IS INSIDE ALLOWED RANGE
      ZP(5) = AMIN1(AMAX1(WW1 * WW2 * WW3, -ONE + EPS), ONE - EPS)
C
C        INITIAL VALUES IN ZP SO START ITERATION CYCLE
C
      IT = 0
100      IT = IT + 1
C
C        FIND MOMENTS AND DERIVATIVES GIVEN LATEST PARAMETER VALUES
C
      CALL LSQ2(ZM, ZP, R, DJ, SS, IER)
      IF (IER .NE. 0) IFAULT = 4
      IF (IFAULT .NE. 0) RETURN
C
C        DECOMPOSE 2ND DERIVATIVE MATRIX
C
      CALL DECOMP(N, N, DJ, IP)
      IF (IP(N) .EQ. 0) IFAULT = 3
      IF (IFAULT .NE. 0) RETURN
C
C        SOLVE LINEAR EQUATION
C
      CALL SOLVE(N, N, DJ, R, IP)
C
C        FIND UPDATED PARAMETER ESTIMATES
C
      ZP(1) = ZP(1) - R(1)
      ZP(2) = ZP(2) - R(2)
      ZP(3) = ZP(3) * EXP(-R(3))
      ZP(4) = ZP(4) * EXP(-R(4))
      WW = EXP(-R(5)) * (ONE + ZP(5)) / (ONE - ZP(5))
      ZP(5) = (WW - ONE) / (WW + ONE)
      IF (SS .LE. TOL) RETURN
      IF (IT .EQ. ITMAX) IFAULT = 2
      IF (IFAULT .NE. 0) RETURN
      GOTO 100
      END
C
      SUBROUTINE LSQ2(ZM, ZP, R, DJ, SS, IER)
C
C        ALGORITHM AS 208.1  APPL. STATIST. (1985) VOL.34, NO.1
C
C        LSQ2 FINDS SS AND DERIVATIVES FOR FITTING LOGISTIC NORMAL
C        FOR S2 BY MOMENTS
C
      REAL D(5, 5, 6), DC, DJ(5, 5), DR, DT, D34, EPS, EPS1, HALF, ONE,
     *  R(5), SS, ZERO, ZI, ZJ, ZM(5), ZMOM(5, 5), ZP(5)
      INTEGER IMAP(6), KMAP(6)
      DATA KMAP(1), KMAP(2), KMAP(3), KMAP(4), KMAP(5), KMAP(6)
     *  /0, 1, 2, 1, 4, 2/
      DATA IMAP(1), IMAP(2), IMAP(3), IMAP(4), IMAP(5), IMAP(6)
     *  /0, 4, 3, 4, 3, 3/
      DATA ZERO, HALF, ONE, EPS, EPS1 /0.0, 0.5, 1.0, 7.7, 18.0/
C
      IER = 0
C
C      CHECK THAT PARAMETERS ARE IN RANGE
C
      IF (AMIN1(ZP(3), ZP(4)) .LE. ZERO) IER = 1
      IF (IER .NE. 0) RETURN
      IF (ABS(ZP(5)) .GT. ONE) IER = 2
      IF (AMAX1(EPS / ZP(3) + ABS(ZP(1)), EPS / ZP(4) + ABS(ZP(2)))
     *  .GT. EPS1) IER = 3
      IF (IER .NE. 0) RETURN
C
C        GET MOMENTS
C
      CALL MOM2(ZP, ZMOM)
C
C        INITIALISE D WITH MOMENTS OR ZEROS
C
      DO 100 K = 2, 4, 2
      DO 100 I = 1, 5
      J = 6 - I
100      D(I, J, K) = ZERO
      DO 120 I = 1, 5
      J0 = 6 - I
      DO 120 J = 1, J0
120      D(I, J, 1) = ZMOM(I, J)
C
C        DO DIFFERENCING ON D TO GET DERIVATIVES
C        KMPA AND IMAP TELL WHICH DIFFERENCES TO USE
C
      DO 200 K = 2, 6
      K0 = KMAP(K)
      I0 = IMAP(K)
      DO 200 I = 1, I0
      ZI = I - 1
      J0 = I0 - I + 1
      DO 200 J = 1, J0
      ZJ = J - 1
      IF (K .LE. 3) D(I, J, K) = ZI * D(I, J, K0) - (ZI + ZJ) *
     *  D(I + 1, J, K0)
      IF (K .GT. 3) D(I, J, K) = ZJ * D(I, J, K0) - (ZI + ZJ) *
     *  D(I, J + 1, K0)
200       CONTINUE
      L = 0
      SS = ZERO
      D34 = ONE / (ZP(3) * ZP(4))
      DC = HALF * (ONE - ZP(5) ** 2) * D34
      DR = ZP(3) / ZP(4)
      DO 400 I = 2, 3
      DO 400 J = 1, I
      K = I - J + 1
      L = L + 1
C
C        FIND SS, RESIDUALS AND DERIVATIVES W.R.T. PARAMETERS
C        D(K,J,1) CONTAINS FITTED MOMENT CORRESPONDING TO ZM(L)
C
      R(L) = D(K, J, 1) - ZM(L)
      SS = SS + R(L) ** 2
      DT = ZP(5) * D(K, J, 6)
      DJ(L, 1) = D(K, J, 2)
      DJ(L, 2) = D(K, J, 4)
      DJ(L, 3) = -D34 * (D(K, J, 3) / DR + DT)
      DJ(L, 4) = -D34 * (D(K, J, 5) * DR + DT)
      DJ(L, 5) = DC * D(K, J, 6)
400      CONTINUE
      RETURN
      END
C
      SUBROUTINE MOM2(ZP, A)
C
C        ALGORITHM AS 208.2  APPL. STATIST. (1985) VOL.34, NO.1
C
C        MOM2 CALCULATES E((X**(I-1))(Y**(J-1)))  I,J=1 TO 5
C          =A(I,J)  FOR LOGISTIC NORMAL
C           I+J LESS THAN 6
C           USING (20,20) POINT HERMITIAN INTEGRATION
C
      REAL A(5, 5), AB(10), CON, ONE, W, WT(10), W0, X, Y, Z, ZERO,
     *  ZP(5), Z1, Z10, Z2, Z20, Z3
C
C        WT AND AB CONTAIN THE WEIGHTS AND ABCISSAS FOR THE
C        HERMITIAN INTEGRATION
C
      DATA WT(1), WT(2), WT(3), WT(4), WT(5), WT(6), WT(7), WT(8),
     *  WT(9), WT(10)
     *  /0.2607930634495549E18, 0.161739333984E18,
     *  0.6150637206397691E17,0.13997837447101E17,
     *  0.1830103131080493E16, 0.1288262799619294E15,
     *  0.4402121090230853E13, 0.6127490259982946E11,
     *  0.2482062362315179E9, 0.1257800672437927E6/
      DATA AB(1), AB(2), AB(3), AB(4), AB(5), AB(6), AB(7), AB(8),
     *  AB(9), AB(10)
     *  /0.34696415708135593, 1.04294534880275103, 1.74524732081412671,
     *  2.45866361117236775, 3.18901481655338941, 3.94396735065731627,
     *  4.73458133404605532, 5.57873880589320117, 6.51059015701365447,
     *  7.61904854167975829/
      DATA ZERO, ONE, CON /0.0, 1.0, 1.0E-36/
C
      Z = SQRT(ONE - ZP(5) ** 2)
C
C        INITIALIZE A() TO ZERO
C
      DO 100 J = 1, 5
      J0 = 6 - J
      DO 100 J1 = 1, J0
100      A(J, J1) = ZERO
      DO 500 I = 1, 10
      X = AB(I)
      W0 = WT(I)
      DO 500 J = 1, 10
      Y = AB(J)
      W = W0 * WT(J)
C
C        GIVEN (X,Y)= GRID POINT AND W= WEIGHT DO FOR (+-X,+-Y)
C
      DO 500 I1 = 1, 2
      DO 400 J1 = 1, 2
      Z10 = W
      Z1 = EXP(X / ZP(3) + ZP(1))
      Z2 = EXP((Y * Z + ZP(5) * X) / ZP(4) + ZP(2))
      Z3 = ONE / (ONE + Z1 + Z2)
      Z1 = Z1 * Z3
      Z2 = Z2 * Z3
      DO 150 I0 = 1, 5
      Z20 = Z10
      J10 = 6 - I0
      DO 140 J0 = 1, J10
      A(I0, J0) = A(I0, J0) + Z20
      Z20 = Z20 * Z2
140      CONTINUE
      Z10 = Z10 * Z1
150      CONTINUE
      Y = -Y
400      CONTINUE
      X = -X
500      CONTINUE
C
C        RESCALE MOMENTS TO AVOID UNDERFLOW
C
      DO 600 J = 1, 5
      J0 = 6 - J
      DO 600 J1 = 1, J0
600      A(J, J1) = A(J, J1) * CON
      RETURN
      END

