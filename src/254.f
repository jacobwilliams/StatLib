C This file starts with a test driver program.
C It includes all auxiliary routines required by the algorithm
C
      PROGRAM TESTMGT
      DIMENSION X(19),N(19),XBAR(10),PI(10),VAR(10),PSE(10),XSE(10),
     *          VSE(10),WK(500)
      DATA N/10,21,56,79,114,122,110,85,85,
     *       61,47,49,47,44,31,20,11,4,4/
C
C        TEST DATA FOR MGT OBTAINED FROM:
C
C        EVERITT B.S. AND HAND D.J. (1981) FINITE MIXTURE DISTRIBUTIONS.
C        CHAPMAN AND HALL. P.46.  
C      
C      PUBLISHED PARAMETER ESTIMATES
C     
C      MEANS:       19.96  ,26.16
C      VARIANCES:    4.6225, 7.6176
C      PROPORTIONS:  0.65  , 0.35
C
C
      LOUT = 6
      NG = 2
      NR = 19
      X0 = 14.5
      DO 10 I = 1,19
	FI = FLOAT(I)
	X(I) = X0+FI
   10 CONTINUE
      KW = 500
      PI(1) = .5
      PI(2) = .5
      XBAR(1) = 20.0
      XBAR(2) = 26.0
      VAR(1) = 4.0
      VAR(2) = 8.0
      TOL = 1.0E-8
      ITRUNC = 1
      ITER = 200
      CALL MGT(NG,NR,ITRUNC,X0,X,NL,NU,N,WK,KW,
     *  PI,XBAR,VAR,PSE,XSE,VSE,TOL,XLOGL,ITER,IFAULT)
      WRITE(LOUT,35)
      WRITE(LOUT,20)IFAULT
      WRITE(LOUT,20)ITER
   20 FORMAT(I10,//)
   30 FORMAT(2F14.7)
   40 FORMAT(E14.7)
   35 FORMAT(//)
      WRITE(LOUT,35)
      WRITE(LOUT,30)(XBAR(I),XSE(I),I=1,2)
      WRITE(LOUT,35)
      WRITE(LOUT,30)(VAR(I),VSE(I),I=1,2)
      WRITE(LOUT,35)
      WRITE(LOUT,30)(PI(I),PSE(I),I=1,2)
      WRITE(LOUT,35)
      WRITE(LOUT,40)XLOGL
      WRITE(LOUT,35)
      STOP
      END
C
      SUBROUTINE MGT(NG, NR, ITRUNC, X0, X, NL, NU, N, WK, KW, PI, XBAR,
     *               VAR, PSE, XSE, VSE, TOL, XLOGL, ITER, IFAULT)
C
C     ALGORITHM AS 254 APPL. STATIST. (1990), VOL. 39, NO. 2
C
C     Subroutine for fitting a mixture of normal distributions to
C     grouped and truncated data.
C   
C     Auxiliary routines required: ALNORM (= AS66), CHOL (= AS6) and
C     SYMINV (= AS7)
C
      INTEGER NR, N(NR), NG, ITRUNC, NL, NU, KW, ITER, IFAULT
      REAL X(NR), WK(KW), PI(NG), XBAR(NG), VAR(NG), PSE(NG), XSE(NG),
     *     VSE(NG), X0, TOL, XLOGL
      INTEGER I1, I2, I3, IFA, IG, IG1, II1, IL, IP, IQ, IR, IR1,
     *        IT, IU, J1, J2, J3, JG, JG1, K1, K2, KG, KG1, LOC1, LOC2,
     *        LOC3, LOC4, LOC5, LOC6, MAXITR, NG1, NG3, NN, NPAR, NR1,
     *        NR2, NROW, NTOT, NULL
      REAL FNT(2), AA, AB, ARG, BA, BB, CA, CB, CONST, FN, ONE, R1, R2,
     *     RATIO, SMALL, SMALLM, SUMPI, T, TOT, TOTN, TRUNC, TWO, TXBAR,
     *     TXLOGL, UFLO, WW,XS,ZERO
      EXTERNAL ALNORM
      REAL ALNORM
      DATA ZERO, ONE, TWO / 0.0, 1.0, 2.0 /
      DATA CONST / 0.398942280401432 /
      DATA UFLO / -87.0 /
      DATA SMALL / -1.0E30 /
      DATA SMALLM / 1.0E-37 /
C
C         Check on initial parameter estimates
C
      IFAULT = 0
      IF(NG.GT.NR)IFAULT = 1
      IF(IFAULT.NE.0)RETURN
      SUMPI = ZERO
      DO 1 IG = 1,NG
	IF(PI(IG).LT.ZERO.OR.PI(IG).GT.ONE)IFAULT = 2
        IF(IFAULT.NE.0)RETURN
	IF(VAR(IG).LE.ZERO)IFAULT = 4
        IF(IFAULT.NE.0)RETURN
	SUMPI = SUMPI+PI(IG)
    1 CONTINUE
      IF(SUMPI.NE.ONE)IFAULT = 3
      IF(IFAULT.NE.0)RETURN
C
C       Set storage location pointers
C
      NG3 = 3*NG
      NPAR = NG3-1
      NPAR = 3*NG-1
      NG1 = NG-1
      NR1 = NR+1
      NR2 = NR+2
      LOC1 = NG*NR2
      LOC2 = 2*LOC1
      LOC3 = 3*LOC1
      LOC4 = LOC3+3*NG
      LOC5 = LOC4+NR2
      LOC6 = LOC5+NPAR+1
      IF(ITRUNC.EQ.1)GOTO 20
      FNT(1) = NL
      FNT(2) = NU
   20 CONTINUE
      NTOT = 0
      DO 30 IR = 1,NR
   30 NTOT = NTOT+N(IR)
      FN = NTOT
C
C        Begin iteration loop
C
      XLOGL = SMALL
      MAXITR = ITER
      ITER = 1
C
C        Calculate the A, B and C's
C
   40 DO 60 IG = 1,NG
	I1 = LOC3+IG
	I2 = I1+NG
	I3 = I2+NG
	WK(I1) = ZERO
	WK(I2) = ZERO
	WK(I3) = ZERO
	AA = ZERO
	BA = ZERO
	CA = ZERO
	DO 50 IR = 1,NR1
          IF(IR.EQ.1)THEN
	    XS = X0
          ELSE
	    IR1 = IR-1
	    XS = X(IR1)
          ENDIF
	  WW = SQRT(VAR(IG))
	  ARG = (XS-XBAR(IG))/WW
          IF(ARG.LE.ZERO)THEN
	    AB = ALNORM(ARG,.FALSE.)
          ELSE
	    AB = ONE-ALNORM(ARG,.TRUE.)
          ENDIF
	  ARG = -ARG*ARG/TWO
          IF(ARG.GT.UFLO)THEN
	    ARG = EXP(ARG)
          ELSE
	    ARG = ZERO
          ENDIF
	  BB = ARG*CONST/WW
	  CB = XS*BB
	  J1 = (IG-1)*NR2+IR
	  J2 = LOC1+J1
	  J3 = LOC2+J1
	  WK(J1) = AB-AA
	  WK(J2) = BB-BA
	  WK(J3) = CB-CA
	  AA = AB
	  BA = BB
	  CA = CB
   50   CONTINUE
	J1 = (IG-1)*NR2+NR2
	J2 = LOC1+J1
	J3 = LOC2+J1
	WK(J1) = ONE-AA
	WK(J2) = -BA
	WK(J3) = -CA
   60 CONTINUE
C
C         Calculate the P's
C
      TRUNC = ZERO
      DO 80 IR = 1,NR2
	K1 = LOC4+IR
	WK(K1) = ZERO
	DO 70 IG = 1,NG
	  J1 = (IG-1)*NR2+IR
	  WK(K1) = WK(K1)+PI(IG)*WK(J1)
   70   CONTINUE
	IF(WK(K1).LE.ZERO)IFAULT = 6
        IF(IFAULT.NE.0)RETURN
	TRUNC = TRUNC+WK(K1)
   80 CONTINUE
C
C        Calculate the sums of the A, B and C's
C
      K1 = LOC4+1
      K2 = LOC4+NR2
      TRUNC = TRUNC-WK(K1)-WK(K2)
      TOT = ZERO
      TXLOGL = ZERO
      DO 110 IR = 1,NR2
	K1 = LOC4+IR
        IF(IR.EQ.1.OR.IR.EQ.NR2)THEN
	  K2 = IR/NR2+1
          IF(ITRUNC.EQ.0)THEN
	      T = FNT(K2)
          ELSE
	    T = FN*WK(K1)/TRUNC
	    FNT(K2) = T
          ENDIF
        ELSE
	  IR1 = IR-1
	  T = FLOAT(N(IR1))
        ENDIF
        IF(ITRUNC.EQ.1.AND.(IR.EQ.1.OR.IR.EQ.NR2))GOTO 90
	TXLOGL = TXLOGL+T*LOG(WK(K1))
   90   CONTINUE
	TOT = TOT+T
	T = T/WK(K1)
	DO 100 IG = 1,NG
	  I1 = LOC3+IG
	  I2 = I1+NG
	  I3 = I2+NG
	  J1 = (IG-1)*NR2+IR
	  J2 = LOC1+J1
	  J3 = LOC2+J1
	  WK(I1) = WK(I1)+T*WK(J1)
	  WK(I2) = WK(I2)+T*WK(J2)
	  WK(I3) = WK(I3)+T*WK(J3)
  100   CONTINUE
  110 CONTINUE
      IF(ITRUNC.NE.0) TXLOGL = TXLOGL-FN*LOG(TRUNC)
C
C        Calculate the new estimates
C
      DO 120 IG = 1,NG
	I1 = LOC3+IG
	I2 = I1+NG
	I3 = I2+NG
	WW = WK(I1)
	R1 = WK(I2)/WW
	R2 = WK(I3)/WW
	ARG = VAR(IG)
	TXBAR = XBAR(IG)-R1*ARG
	VAR(IG) = ARG*(ONE-R2+
     *          (TWO*TXBAR-XBAR(IG))*R1)+
     *          (TXBAR-XBAR(IG))**2
	XBAR(IG) = TXBAR
	PI(IG) = PI(IG)*WW/TOT
  120 CONTINUE
C
C        If convergence criterion is not satisfied,
C        perform another iteration
C
      ARG = ABS((TXLOGL-XLOGL)/XLOGL)
      XLOGL = TXLOGL
      IF(ARG.LE.TOL)GOTO 130
      ITER = ITER+1
      IF(ITER.LE.MAXITR)GOTO 40
      IFAULT = 5
      RETURN
C
C        Calculate the standard errors of the estimates
C
  130 IT = 0
      DO 150 IP = 1,NG3
	I1 = LOC5+IP
	WK(I1) = ZERO
	DO 140 IQ = 1,IP
	  IT = IT+1
	  I2 = LOC6+IT
	  WK(I2) = ZERO
  140   CONTINUE
  150 CONTINUE
      IF(ITRUNC.EQ.0)THEN
	IL = 1
	IU = NR2
	TOTN = FN+FNT(1)+FNT(2)
      ELSE
	IL = 2
	IU = NR1
	TOTN = FN
      ENDIF
C
      DO 190 IR = IL,IU
	IR1 = IR-1
	IF(IR.EQ.1)NN = NL
	IF(IR.EQ.NR2)NN = NU
	IF(IR.NE.1.AND.IR.NE.NR2)NN = N(IR1)
	K1 = LOC4+IR
	K2 = NG1*NR2+IR
	WW = WK(K2)
C
	DO 160 IG = 1,NG
	  J1 = (IG-1)*NR2+IR
	  J2 = LOC1+J1
	  J3 = LOC2+J1
	  RATIO = PI(IG)/WK(K1)
	  WK(J1) = (WK(J1)-WW)/WK(K1)
	  WK(J2) = -RATIO*WK(J2)
	  WK(J3) = -(WK(J2)*XBAR(IG)+RATIO*WK(J3))/(TWO*VAR(IG))
  160   CONTINUE
C
	IT = 0
	DO 180 IG = 1,NG3
	  IG1 = IG-1
	  I1 = IG1*NR2+IR
	  II1 = LOC5+IG
	  WK(II1) = WK(II1)+NN*WK(I1)
	  DO 170 JG = 1,IG
	    IT = IT+1
	    I2 = (JG-1)*NR2+IR
	    I3 = LOC6+IT
            IF(LOG(ABS(WK(I1))+SMALLM)+
     *         LOG(ABS(WK(I2))+SMALLM).LT.
     *         LOG(SMALLM))GOTO 170
	    WK(I3) = WK(I3)+NN*WK(I1)*WK(I2)
  170     CONTINUE
  180   CONTINUE
  190 CONTINUE
C
      IF(ITRUNC.EQ.0)GOTO 220
C
      NL = FNT(1)
      NU = FNT(2)
      IT = 0
      DO 210 IG = 1,NG3
	I1 = LOC5+IG
	DO 200 JG = 1,IG
	  IT = IT+1
	  I2 = LOC5+JG
	  I3 = LOC6+IT
	  WK(I3) = WK(I3)-WK(I1)*WK(I2)/TOTN
  200   CONTINUE
  210 CONTINUE
C
  220 IT = 0
      IP = 0
      DO 240 IG = 1,NG3
	DO 230 JG = 1,IG
	  IT = IT+1
          IF(IG.EQ.NG.OR.JG.EQ.NG)GO TO 230
	  IP = IP+1
	  I1 = LOC6+IT
	  WK(IP) = WK(I1)
  230   CONTINUE
  240 CONTINUE
C
      NROW = NPAR
      NN = NPAR*(NPAR+1)/2
      CALL SYMINV(WK,NROW,NN,WK,WK(LOC3),NULL,IFA)
      IF(NG.EQ.1)THEN
	PSE(1) = ZERO
      ELSE
	IT = 0
	PSE(NG) = ZERO
	DO 260 IG = 1,NG1
	  KG = IG*(IG+1)/2
	  PSE(IG) = WK(KG)
	  PSE(NG) = PSE(NG)-PSE(IG)
	  PSE(IG) = SQRT(PSE(IG))
	  DO 250 JG = 1,IG
	    IT = IT+1
	    PSE(NG) = PSE(NG)+TWO*WK(IT)
  250     CONTINUE
  260   CONTINUE
	PSE(NG) = SQRT(PSE(NG))
      ENDIF
C
      DO 270 IG = 1,NG
	JG = IG+NG-1
	JG1 = (JG+1)*JG/2
	KG = JG+NG
	KG1 = (KG+1)*KG/2
	XSE(IG) = SQRT(WK(JG1))
	VSE(IG) = SQRT(WK(KG1))
  270 CONTINUE
      RETURN
      END
C
      REAL FUNCTION ALNORM(X,UPPER)
C
C         ALGORITHM AS66 APPLIED STATISTICS (1973) VOL22 NO.3
C
C       EVALUATES THE TAIL AREA OF THE STANDARDISED NORMAL CURVE
C       FROM X TO INFINITY IF UPPER IS .TRUE. OR
C       FROM MINUS INFINITY TO X IF UPPER IS .FALSE.
C
      REAL LTONE, UTZERO, ZERO, HALF, ONE, CON,
     $   A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6,
     $   B7, B8, B9, B10, B11, B12, X, Y, Z, ZEXP
      LOGICAL UPPER,UP
C
C         LTONE AND UTZERO MUST BE SET TO SUIT THE PARTICULAR COMPUTER
C         (SEE INTRODUCTORY TEXT)
C
      DATA LTONE, UTZERO /7.0, 12.5/
      DATA ZERO, HALF, ONE, CON /0.0, 0.5, 1.0, 1.28/
      DATA           A1,               A2,                A3,
     $               A4,               A5,                A6,
     $               A7
     $   /0.398942280444, 0.399903438504, 5.75885480458,
     $     29.8213557808,  2.62433121679, 48.6959930692,
     $     5.92885724438/
      DATA           B1,               B2,                B3,
     $               B4,               B5,                B6,
     $               B7,               B8,                B9,
     $              B10,              B11,               B12
     $   /0.398942280385,       3.8052E-8,  1.00000615302,
     $       3.98064794E-4, 1.98615381364, 0.151679116635,
     $     5.29330324926,  4.8385912808,  15.1508972451,
     $    0.742380924027,  30.789933034,  3.99019417011/
C
      ZEXP(Z) = EXP(Z)
C
      UP=UPPER
      Z=X
      IF(Z.GE.ZERO)GOTO 10
      UP=.NOT.UP
      Z=-Z
   10 IF(Z.LE.LTONE.OR.UP.AND.Z.LE.UTZERO)GOTO 20
      ALNORM=ZERO
      GOTO 40
   20 Y=HALF*Z*Z
      IF(Z.GT.CON) GOTO 30
C
      ALNORM=HALF-Z*(A1-A2*Y/(Y+A3-A4/(Y+A5+A6/(Y+A7))))
      GOTO 40
   30 ALNORM=B1*ZEXP(-Y)/(Z-B2+B3/(Z+B4+B5/(Z-B6+B7/(Z+B8-B9/
     $   (Z+B10+B11/(Z+B12))))))
   40 IF(.NOT.UP)ALNORM=ONE-ALNORM
      RETURN
      END
C
        SUBROUTINE SYMINV(A, N, NN, C, W, NULLTY, IFAULT)
C
C       ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968, P.198.
C
C       FORMS IN C( ) AS LOWER TRIANGLE, A GENERALISED INVERSE
C       OF THE POSITIVE SEMI-DEFINITE SYMMETRIC MATRIX A( )
C       ORDER N, STORED AS LOWER TRIANGLE.
C
        REAL A(NN), C(NN), W(N), X, ZERO, ONE
C
        DATA ZERO, ONE /0.0, 1.0/
C
C       CHOLESKY FACTORIZATION OF A, RESULT IN C
C
        CALL CHOL(A, N, NN, C, NULLTY, IFAULT)
        IF(IFAULT.NE.0) RETURN
C
C       INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
C       OF C, ROW BY ROW STARTING WITH THE LAST ROW.
C       IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.
C
        IROW=N
        NDIAG=NN
   10   L=NDIAG
        IF (C(NDIAG) .EQ. ZERO) GOTO 60
        DO 20 I=IROW,N
        W(I)=C(L)
        L=L+I
   20   CONTINUE
        ICOL=N
        JCOL=NN
        MDIAG=NN
   30   L=JCOL
        X=ZERO
        IF(ICOL.EQ.IROW) X=ONE/W(IROW)
        K=N
   40   IF(K.EQ.IROW) GO TO 50
        X=X-W(K)*C(L)
        K=K-1
        L=L-1
        IF(L.GT.MDIAG) L=L-K+1
        GO TO 40
   50   C(L)=X/W(IROW)
        IF(ICOL.EQ.IROW) GO TO 80
        MDIAG=MDIAG-ICOL
        ICOL=ICOL-1
        JCOL=JCOL-1
        GO TO 30
   60   DO 70 J=IROW,N
        C(L)=ZERO
        L=L+J
   70   CONTINUE
   80   NDIAG=NDIAG-IROW
        IROW=IROW-1
        IF(IROW.NE.0) GO TO 10
        RETURN
        END
C
        SUBROUTINE CHOL(A, N, NN, U, NULLTY, IFAULT)
C
C       ALGORITHM AS6, APPLIED STATISTICS, VOL.17, (1968)
C
C       GIVEN A SYMMETRIC MATRIX ORDER N AS LOWER TRIANGLE IN A( )
C       CALCULATES AN UPPER TRIANGLE, U( ), SUCH THAT UPRIME * U = A.
C       A MUST BE POSITIVE SEMI-DEFINITE.  ETA IS SET TO MULTIPLYING
C       FACTOR DETERMINING EFFECTIVE ZERO FOR PIVOT.
C
        REAL A(NN), U(NN), ETA, ETA2, X, W, ZERO, ZABS, ZSQRT
C
C       THE VALUE OF ETA WILL DEPEND ON THE WORD-LENGTH OF THE 
C       COMPUTER BEING USED.  SEE INTRODUCTORY TEXT.
C
        DATA ETA, ZERO /1.E-5, 0.0/
C
        ZABS(X) = ABS(X)
        ZSQRT(X) = SQRT(X)
C
        IFAULT=1
        IF(N.LE.0) GO TO 100
        IFAULT=3
        IF (NN .NE. N*(N+1)/2) RETURN
        IFAULT=2
        NULLTY=0
        J=1
        K=0
        ETA2 = ETA*ETA
        II=0
C
C       FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.
C
        DO 80 ICOL=1,N
        II=II+ICOL
        X=ETA2*A(II)
        L=0
        KK=0
C
C       IROW = ROW NUMBER WITHIN COLUMN ICOL
C
        DO 40 IROW=1,ICOL
        KK=KK+IROW
        K=K+1
        W=A(K)
        M=J
        DO 10 I=1,IROW
        L=L+1
        IF(I.EQ.IROW) GO TO 20
        W=W-U(L)*U(M)
        M=M+1
   10   CONTINUE
   20   IF(IROW.EQ.ICOL) GO TO 50
        IF(U(L).EQ.ZERO) GO TO 30
        U(K)=W/U(L)
        GO TO 40
   30   IF (W*W .GT. ZABS(X*A(KK))) RETURN
        U(K)=ZERO
   40   CONTINUE
   50   IF (ZABS(W) .LE. ZABS(ETA*A(K))) GOTO 60
        IF(W.LT.ZERO) RETURN
        U(K)=ZSQRT(W)
        GO TO 70
   60   U(K)=ZERO
        NULLTY=NULLTY+1
   70   J=J+ICOL
   80   CONTINUE
        IFAULT=0
  100   RETURN
        END
