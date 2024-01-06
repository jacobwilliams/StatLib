      program t_as139
c
c     Test AS139
c
      IMPLICIT NONE
      INTEGER n, p(100), m, mplone, rowx, colx, lenw, lenwrk, maxits,
     *        ifault, i, i1, j
      PARAMETER (n=100, m=5, mplone=m+1, rowx=n, colx=m, lenw=mplone+n,
     *        lenwrk=m*n, maxits=25)
      DOUBLE PRECISION y, y1(n), y2(n), x(rowx, colx), w(lenw),
     *        vcov(lenwrk), work(lenwrk), alpha(mplone), tol(mplone),
     *        temp, se(m)
      REAL rand, randn
      INTEGER IX, IY, IZ
      COMMON /RANDC/ IX, IY, IZ
c
c     Generate artificial data
c
      WRITE(*, *)'Enter 3 integers for random number seeds: '
      READ(*, *) ix, iy, iz

      DO i = 1, n
        x(i, 1) = 1.0
        x(i, 2) = i/10
        x(i, 3) = MOD(i, 10)
        x(i, 4) = x(i, 2) + rand() - 0.5
        x(i, 5) = x(i, 3) + rand() - 0.5
        y = x(i, 1) + x(i, 2) + x(i, 3) + x(i, 4) + x(i, 5) + randn()

        temp = rand()
        IF (temp .LT. 0.1) THEN
          p(i) = 1
          y1(i) = INT(y)
        ELSE IF (temp .LT. 0.2) THEN
          p(i) = -1
          y1(i) = INT(y) + 1
        ELSE IF (temp .LT. 0.3) THEN
          p(i) = 2
          y1(i) = INT(y)
          y2(i) = INT(y) + 1
        ELSE
          p(i) = 0
          y1(i) = y
        END IF
      END DO

      DO i = 1, mplone
        tol(i) = 1.d-06
      END DO

      alpha(mplone) = -1.0             ! No initial estimate for sigma

      CALL REGRES(N, Y1, Y2, P, MPLONE, X, ROWX, COLX, W, LENW,
     *    VCOV, WORK, LENWRK, ALPHA, TOL, MAXITS, IFAULT)

      IF (ifault .LT. 0) THEN
        WRITE(*, *)'IFAULT = ', ifault
        STOP
      ELSE
        WRITE(*, *)'No. of iterations = ', ifault
      END IF

      WRITE(*, 900) (alpha(i), i=1,m)
  900 FORMAT(1x, 'Regression coefficients:'/ 1x, 5f10.4)
      DO i = 1, m
        se(i) = SQRT( vcov((i-1)*mplone + i) )
      END DO
      WRITE(*, 910) (se(i), i=1,m)
  910 FORMAT(1x, 'Estimated standard errors:'/ 1x, 5f10.4)
      WRITE(*, 920) alpha(mplone)
  920 FORMAT(1x, 'Residual variance estimate = ', f10.4)
      WRITE(*, *) 'Variance-covariance matrix:-'
      i1 = 1
      DO i = 1, m
        WRITE(*, '(1x, 7f10.4)') (vcov(j),j=i1,i1+m-1)
        i1 = i1 + mplone
      END DO
      end


      REAL FUNCTION RAND()
C
C     The Wichmann & Hill random number generator
C     Algorithm AS183, Appl. Statist., 31, 188-190, 1982.
C     The cycle length is 6.95E+12.
C     This random number generator is very slow compared with most
C     others, but it is dependable, and the results are reproducible.
C
      INTEGER IX, IY, IZ
      COMMON /RANDC/ IX, IY, IZ
C
C     Initialize IX, IY & IZ if necessary
C
      IF(IX .LE. 0) IX = 777
      IF(IY .LE. 0) IY = 777
      IF(IZ .LE. 0) IZ = 777
C
   10 IX = MOD(171*IX, 30269)
      IY = MOD(172*IY, 30307)
      IZ = MOD(170*IZ, 30323)
      RAND = MOD(FLOAT(IX)/30269. + FLOAT(IY)/30307. +
     1       FLOAT(IZ)/30323. , 1.0)
      RETURN
      END


      REAL FUNCTION RANDN()
C
C      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
C
C  The function RANDN() returns a normally distributed pseudo-random
C  number with zero mean and unit variance.  Calls are made to a
C  function subprogram RAND() which must return independent random
C  numbers uniform in the interval (0,1).
C
C  The algorithm uses the ratio of uniforms method of A.J. Kinderman
C  and J.F. Monahan augmented with quadratic bounding curves.
C
      DATA S,T,A,B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA R1,R2/ 0.27597, 0.27846/
C
C  Generate P = (u,v) uniform in rectangle enclosing acceptance region
 50   U = RAND()
      V = RAND()
      V = 1.7156 * (V - 0.5)
C  Evaluate the quadratic form
      X = U - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C  Accept P if inside inner ellipse
      IF (Q .LT. R1) GO TO 100
C  Reject P if outside outer ellipse
      IF (Q .GT. R2) GO TO 50
C  Reject P if outside acceptance region
      IF (V**2 .GT. -4.0*ALOG(U)*U**2) GO TO 50
C  Return ratio of P's coordinates as the normal deviate
 100  RANDN = V/U
      RETURN
      END



      SUBROUTINE REGRES(N, Y1, Y2, P, MPLONE, X, ROWX, COLX, W, LENW,
     *    VCOV, WORK, LENWRK, ALPHA, TOL, MAXITS, IFAULT)
C***** NOTE: this routine uses the auxiliary routine
C***** AS7  with the modified argument lists suggested by Freeman 
C***** Vol 31 No 3 p336--339.
C
C  as r31 Vol 29 No 2 1980 p228 -- incorporated
C  as r32 Vol 30 No 1 1981 p105 -- incorporated
C  as r91 Vol 42 No 3 1993 p583 -- incorporated
C
C        ALGORITHM AS139 APPL. STATIST. (1979) VOL.28, NO.2
C
C        COMPUTE MAXIMUM LIKELIHOOD ESTIMATES
C        FROM A LINEAR MODEL WITH NORMAL HETEROGENOUS VARIANCE.
C        THE DESIGN MATRIX MUST BE NON-SINGULAR.  THE DEPENDENT
C        VARIABLE MAY INCLUDE OBSERVATIONS CENSORED IN EITHER TAIL
C        AND/OR OBSERVATIONS CONFINED BETWEEN FINITE LIMITS.
C
C     Auxiliary routine required: RMILLS which is the last routine in
C     AS138.
C
      IMPLICIT NONE

      INTEGER N, MPLONE, ROWX, COLX, P(N), LENW, LENWRK, MAXITS, IFAULT
      DOUBLE PRECISION X(ROWX,COLX), TOL(MPLONE), Y1(N), Y2(N),
     *    ALPHA(MPLONE)
      DOUBLE PRECISION VCOV(LENWRK), WORK(LENWRK), W(LENW)
C
C     Local variables
C
      INTEGER M, I, J, K, II, NDIMC, NUL, JJ, IIK, IJ, IPT, NJ, IIJ,
     *      JJI, KK
      DOUBLE PRECISION QLIMIT, RLIMIT, C, ZERO, HALF, ONE
      PARAMETER (QLIMIT = 0.00001D0, RLIMIT=0.00001D0, C=0.39894228D0)
      PARAMETER (zero=0.d0, half=0.5D0, one=1.D0)
      DOUBLE PRECISION TEMP, SUM2, DEMP, R, R2, XSIG, TD, YMEAN, F, A,
     *      FU, YN, YQ, TMPU, YNU, YQU, TINT, TEMP2, YS, YSU, B
C
C        INITIALIZATION
C
      M = MPLONE-1
C
C        CHECK ARRAY SIZES, ETC
C
      IFAULT = -7
      IF(ROWX.LT.N) RETURN
      IFAULT = -8
      IF(COLX.LT.M) RETURN
      IFAULT = -9
      IF(LENW.LT.(MPLONE+N)) RETURN
      IFAULT = -10
      IF(LENWRK.LT.(M*N)) RETURN
C
C        COMPUTE X'X IN LOWER TRIANGULAR FORM
C
      II = 0
      DO 53 I = 1, M
        DO 50 J = 1, I
          TEMP = ZERO
          DO 40 K = 1, N
   40     TEMP = TEMP + X(K,I)*X(K,J)
          II = II + 1
          VCOV(II) = TEMP
   50   CONTINUE
   53 CONTINUE
      NDIMC = M*(M+1)/2
      CALL SYMINV(VCOV, M, NDIMC, WORK, W, NUL, IFAULT)
      IF(IFAULT.NE.0) GOTO 60
      IF(NUL.EQ.0) GOTO 70
   60 VCOV(2) = IFAULT
      VCOV(1) = NUL
      IFAULT = -5
      RETURN
C
C        MATRIX NON-SINGULAR AND INVERSE OBTAINED
C        COMPUTE (X'X)INVERSE*X
C        FOLLOWING SCHEME USED TO REDUCE NUMBER OF STORAGE ARRAYS
C        NEEDED.   EXPAND FROM TRIANGULAR TO SQUARE MATRIX
C
   70 CALL UNPACK(WORK, M, LENWRK)
C
C        DO MULTIPLICATION - ONE ROW AT A TIME - STARTING
C        WITH THE LAST ONE
C
      JJ = N*M
      II = M*M
      DO 220 I = 1, M
        II = II - M
        DO 200 J = 1, N
          TEMP = ZERO
          DO 170 K = 1, M
            IIK = II + K
            TEMP = TEMP + WORK(IIK)*X(J,K)
  170     CONTINUE
          W(J) = TEMP
  200   CONTINUE
        DO 210 J = 1, N
          IJ = N+1-J
          WORK(JJ) = W(IJ)
          JJ = JJ-1
  210   CONTINUE
  220 CONTINUE
C
      XSIG = ALPHA(MPLONE)
      IF(XSIG.GT.ZERO) GOTO 500
C
C        NO ACCEPTABLE INITIAL VALUE FOR SIGMA HAS BEEN INPUT,
C        OBTAIN INITIAL ESTIMATES FROM EXACTLY SPECIFIED
C        OBSERVATIONS ONLY (ALTHOUGH MATRIX BASED ON ALL
C        OBSERVATIONS) AND CONFINED OBSERVATIONS
C
      II = -N
      DO 300 I = 1, M
        II = II + N
        TEMP = ZERO
        DO 280 J = 1, N
          IIJ = II + J
          IPT = P(J)
          IF(IPT.EQ.0) GOTO 270
          IF(IPT.EQ.2) TEMP = TEMP + WORK(IIJ)*(Y1(J) + Y2(J))*HALF
          GOTO 280
  270     TEMP = TEMP + WORK(IIJ)*Y1(J)
  280   CONTINUE
        ALPHA(I) = TEMP
  300 CONTINUE
C
C        CALCULATE INITIAL ESTIMATE OF SIGMA
C
      SUM2 = ZERO
      TEMP = ZERO
      DO 350 I = 1, N
        IPT = P(I)
        IF(IABS(IPT).EQ.1) GOTO 350
        DEMP = Y1(I)
        IF(IPT.EQ.2) DEMP = (DEMP + Y2(I))*HALF
        DO 320 J = 1, M
  320   DEMP = DEMP - ALPHA(J)*X(I, J)
        SUM2 = SUM2 + DEMP**2
        TEMP = TEMP + ONE
  350 CONTINUE
      XSIG = SQRT(SUM2/TEMP)
C
C        COMPUTE SOME CONSTANTS NEEDED THROUGHOUT
C
  500 R = ZERO
      R2 = ZERO
      IFAULT = -2
      DO 600 I = 1, N
        IPT = P(I)
        IF(IPT.EQ.0) GOTO 550
        IF(IPT.EQ.2 .AND. ABS(Y1(I)-Y2(I)).LE. QLIMIT*ABS(Y1(I)))
     *                     GOTO 540
        IF(IPT.NE.2) GOTO 600
        R2 = R2 + ONE
        IF(Y1(I).LT.Y2(I)) GOTO 600
        RETURN
  540   P(I) = 0
  550   R = R + ONE
        W(I) = Y1(I)
  600 CONTINUE
      I = R + R2 + 0.01
      IFAULT = -4
      IF(I.LT.MPLONE) RETURN
      IFAULT = 0
C
C        START OF ITERATION PROCEDURE
C
  620 TD = R
      SUM2 = ZERO
C
C        COMPLETE W-VECTOR
C
      DO 1000 I = 1, N
        IPT = P(I)
        YMEAN = ZERO
        DO 650 J = 1, M
  650   YMEAN = YMEAN + ALPHA(J)*X(I,J)
        IF(IPT.EQ.0) GOTO 990
C
C        OBSERVATION NOT EXACTLY SPECIFIED
C
        TEMP = (Y1(I) - YMEAN)/XSIG
        IF(IPT-1)750, 700, 800
C
C        OBSERVATION CENSORED FROM ABOVE - LOWER BOUND KNOWN
C
  700   CALL RMILLS(TEMP, F, RLIMIT)
        W(I) = YMEAN + XSIG*F
        TD = TD + F*(F-TEMP)
        GOTO 990
C
C        OBSERVATION CENSORED FROM BELOW - UPPER BOUND KNOWN
C
  750   CALL RMILLS(-TEMP, F, RLIMIT)
        W(I) = YMEAN - XSIG*F
        TD = TD + F*(F+TEMP)
        GOTO 990
C
C        OBSERVATION CONFINED TO LIE BETWEEN TWO FINITE LIMITS
C
  800   YN = EXP(-HALF*TEMP**2)*C
        CALL RMILLS(TEMP, F, RLIMIT)
        YQ = YN/F
        TMPU = (Y2(I) - YMEAN)/XSIG
        YNU = EXP(-HALF*TMPU**2)*C
        CALL RMILLS(TMPU, FU, RLIMIT)
        YQU = YNU/FU
        TINT = YQ - YQU
        IF(TINT.GE.QLIMIT) GOTO 820
C
C        AFTER STANDARDIZING, UPPER AND LOWER LIMITS RESULT IN
C        SAME PROBABILITY INTEGRAL
C
        IFAULT = -3
        RETURN
  820   A = (YN - YNU)/TINT
        W(I) = YMEAN + XSIG*A
        TD = TD + (A**2 + (TMPU*YNU - TEMP*YN)/TINT)
C
C        CALCULATE RESIDUAL SUM OF SQUARES
C
  990   SUM2 = SUM2 + (W(I) - YMEAN)**2
 1000 CONTINUE
C
C        UPDATE PARAMETER ESTIMATES - STORE IN END OF W-VECTOR
C
      JJ = -N
      DO 1200 J = 1, M
        JJ = JJ + N
        TEMP = ZERO
        DO 1100 I = 1, N
          JJI = JJ + I
          TEMP = TEMP + WORK(JJI)*W(I)
 1100   CONTINUE
        NJ = N + J
        W(NJ) = TEMP
 1200 CONTINUE
      NJ = N + MPLONE
      W(NJ) = SQRT(SUM2/TD)
C
C        TEST FOR CONVERGENCE
C
      DO 1300 J = 1, MPLONE
        NJ = N + J
        IF(ABS(ALPHA(J)-W(NJ)).GE.TOL(J)) GOTO 1400
 1300 CONTINUE
C
C        IF WE REACH HERE, CONVERGENCE OBTAINED
C
      IJ = IFAULT
      IFAULT = -1
C
C        UPDATE VALUES
C
 1400 DO 1450 J = 1, MPLONE
        NJ = N + J
        ALPHA(J) = W(NJ)
 1450 CONTINUE
      XSIG = ALPHA(MPLONE)
      IFAULT = IFAULT + 1
      IF(IFAULT.EQ.0) GOTO 1600
      IF(IFAULT.LE.MAXITS) GOTO 620
      IFAULT = -1
      RETURN
C
C        CONVERGENCE OBTAINED - COMPUTE VARIANCE-COVARIANCE
C        MATRIX. INITIALIZE WORK ARRAY
C
 1600 II = MPLONE*(MPLONE + 1)/2
      DO 1650 I = 1, II
 1650 WORK(I) = ZERO
      DO 2500 I = 1, N
        IPT = P(I)
        YS = Y1(I)
        DO 1680 J = 1, M
 1680   YS = YS - ALPHA(J)*X(I,J)
        YS = YS/XSIG
        JJ = 0
        IF(IPT.NE.0) GOTO 1900
C
C        EXACTLY SPECIFIED OBSERVATION
C
        DO 1750 K = 1, M
          DO 1720 J = 1, K
            JJ = JJ + 1
            WORK(JJ) = WORK(JJ) + X(I,K)*X(I,J)
 1720     CONTINUE
          KK = II - MPLONE + K
          WORK(KK) = WORK(KK) + YS*X(I,K)
 1750   CONTINUE
        WORK(II) = WORK(II) + ONE + YS**2
        GOTO 2500
 1900   IF(IPT-1) 2100, 2000, 2300
C
C        OBSERVATION CONSORED FROM ABOVE - LOWER BOUND KNOWN
C
 2000   CALL RMILLS(YS, F, RLIMIT)
        TEMP = F*(F - YS)
        GOTO 2150
C
C        OBSERVATION CONSORED FROM BELOW - UPPER BOUND KNOWN
C
 2100   CALL RMILLS(-YS, F, RLIMIT)
        TEMP = F*(F + YS)
C
C        ROUTINE FOR CENSORED OBSERVATIONS
C
 2150   DO 2190 K = 1, M
          DO 2170 J = 1, K
            JJ = JJ + 1
            WORK(JJ) = WORK(JJ) + X(I,J)*X(I,K)*TEMP
 2170     CONTINUE
          KK = II - MPLONE + K
          WORK(KK) = WORK(KK) + YS*X(I,K)*TEMP
 2190   CONTINUE
        WORK(II) = WORK(II) + YS**2*TEMP
        GOTO 2500
C
C        OBSERVATION CONFINED BETWEEN TWO FINITE LIMITS
C
 2300   YN = EXP(-HALF*YS**2)*C
        CALL RMILLS(YS, F, RLIMIT)
        YQ = YN/F
        YSU = YS + (Y2(I) - Y1(I))/XSIG
        CALL RMILLS(YSU, FU, RLIMIT)
        YNU = EXP(-HALF*YSU**2)*C
        YQU = YNU/FU
        TINT = YQ - YQU
        A = (YN - YNU)/TINT
        B = (YNU*YSU - YN*YS)/TINT
        TEMP = A**2 + B
        TEMP2 = A*B + (YS**2 * YN - YSU**2 * YNU) / TINT
        DO 2350 K = 1, M
          DO 2330 J = 1, K
            JJ = JJ + 1
            WORK(JJ) = WORK(JJ) + X(I,J)*X(I,K)*TEMP
 2330     CONTINUE
          KK = II - MPLONE + K
          WORK(KK) = WORK(KK) - X(I,K) * TEMP2
 2350   CONTINUE
        TEMP = (YS**3*YN - YSU**3*YNU)/TINT
        WORK(II) = WORK(II) - TEMP + B**2
 2500 CONTINUE
C
C        INVERT THE MATRIX
C
      CALL SYMINV(WORK, MPLONE, II, VCOV, W, NUL, IFAULT)
      IF(IFAULT.EQ.0 .AND. NUL.EQ.0) GOTO 2550
      VCOV(2) = IFAULT
      VCOV(1) = NUL
      IFAULT = -6
      RETURN
C
C        RESTORE ITERATION COUNTER
C
 2550 IFAULT = IJ
C
C        MULTIPLY BY SIGMA-SQUARED
C
      TEMP = XSIG**2
      DO 2580 I = 1, II
 2580 VCOV(I) = VCOV(I)*TEMP
C
C        UNPACK THE MATRIX
C
      CALL UNPACK(VCOV, MPLONE, LENWRK)
      RETURN
      END
*****************************************************************
      SUBROUTINE UNPACK(X, N, LENX)
C
C        ALGORITHM AS139.1 APPL. STATIST. (1979) VOL.28 NO.2
C
C        THIS SUBROUTINE EXPANDS A SYMMETRIC MATRIX STORED IN LOWER
C        TRIANGULAR FORM IN THE FIRST N*(N+1)/2 POSITIONS OF X
C        INTO A MATRIX USING THE FIRST N*N POSITIONS
C
C        LENX = LENGTH OF VECTOR X - MUST BE LESS THAN N*N
C
      INTEGER N, LENX
      DOUBLE PRECISION X(LENX)
C
C     Local variables
C
      INTEGER NSQ, II, JJ, I, IJ, KK, J

      NSQ = N*N
      II = NSQ
      JJ = N*(N+1)/2
C
C        STORE LAST ROW
C
      DO 10 I = 1, N
        X(II) = X(JJ)
        II = II-1
        JJ = JJ-1
   10 CONTINUE
      DO 80 I = 2, N
C
C        OBTAIN UPPER PART OF MATRIX FROM PART ALREADY SHIFTED
C
        IJ = I - 1
        KK = NSQ+1-I
        DO 50 J = 1, IJ
          X(II) = X(KK)
          II = II - 1
          KK = KK - N
   50   CONTINUE
C
C        OBTAIN LOWER PART OF MATRIX FROM
C        ORIGINAL TRIANGULAR STORAGE
C
        IJ = N - IJ
        DO 70 J = 1, IJ
          X(II) = X(JJ)
          II = II - 1
          JJ = JJ - 1
   70   CONTINUE
   80 CONTINUE
      RETURN
      END


      SUBROUTINE RMILLS(X, FUNC, TOL)
C
C       ALGORITHM AS 138.1 APPL. STATIST. (1979) VOL.28 NO.2
C
C       COMPUTE THE RECIPROCAL OF MILLS RATIO
C
      DOUBLE PRECISION X, FUNC, TOL
C
C     Local variables
C
      DOUBLE PRECISION Y, SGN, S, A, T, R, B, B1, A1, A2, B2, A0, B0
      DOUBLE PRECISION FPI, FPII, ZERO, HALF, ONE, TWO, TEN, SMALL
      DATA FPI /1.2533141D0/, FPII /0.7978846D0/, ZERO /0.D0/,
     *     HALF /0.5D0/, ONE /1.D0/, TWO /2.D0/, TEN /10.D0/,
     *     SMALL /0.000001D0/
C
      FUNC = ZERO
      IF (X .LT. -TEN) RETURN
      FUNC = FPII
      Y = ABS(X)
      IF (Y .LT. SMALL) RETURN
      SGN = ONE
      IF (X .LT. ZERO) SGN = -ONE
      IF (Y .GT. TWO) GOTO 100
      S = ZERO
      A = ONE
      T = Y
      R = Y
      B = Y ** 2
40    A = A + TWO
      S = T
      R = R * B / A
      T = T + R
      IF (R .GT. TOL) GOTO 40
      FUNC = ONE / (FPI * EXP(HALF * B) - SGN * T)
      RETURN
C
100   A = TWO
      B1 = Y
      S = Y
      A1 = Y ** 2 + ONE
      A2 = Y * (A1 + TWO)
      B2 = A1 + ONE
      T = A2 / B2
140   A = A + ONE
      A0 = A1
      A1 = A2
      A2 = Y * A1 + A * A0
      B0 = B1
      B1 = B2
      B2 = Y * B1 + A * B0
      R = S
      S = T
      T = A2 / B2
      IF (T - R .GT. TOL .OR. T - S .GT. TOL) GOTO 140
      FUNC = T
      IF (SGN .LT. ZERO) FUNC = T / (TWO*FPI*EXP(HALF * Y**2) * T - ONE)
      RETURN
      END



        subroutine syminv(a, n, nn, c, w, nullty, ifault)
c
c       Algorithm AS7, Applied Statistics, vol.17, 1968, p.198.
c
c       Forms in c( ) as lower triangle, a generalised inverse
c       of the positive semi-definite symmetric matrix a( )
c       order n, stored as lower triangle.
c
c       arguments:-
c       a()     = input, the symmetric matrix to be inverted, stored in
c                 lower triangular form
c       n       = input, order of the matrix
c       nn      = input, the size of the a and c arrays     n*(n+1)/2
c       c()     = output, the inverse of a (a generalized inverse if c is
c                 singular), also stored in lower triangular.
c                 c and a may occupy the same locations.
c       w()     = workspace, dimension at least n.
c       nullty  = output, the rank deficiency of a.
c       ifault  = output, error indicator
c                       = 1 if n < 1
c                       = 2 if a is not +ve semi-definite
c                       = 3 if nn < n*(n+1)/2
c                       = 0 otherwise
c
c***************************************************************************
c
        double precision a(nn), c(nn), w(n), x, zero, one
c
        data zero, one /0.0d0, 1.0d0/
c
c       cholesky factorization of a, result in c
c
        call chol(a, n, nn, c, nullty, ifault)
        if(ifault.ne.0) return
c
c       invert c & form the product (cinv)'*cinv, where cinv is the inverse
c       of c, row by row starting with the last row.
c       irow = the row number, ndiag = location of last element in the row.
c
        irow=n
        ndiag=nn
   10   l=ndiag
        if (c(ndiag) .eq. zero) goto 60
        do 20 i=irow,n
          w(i)=c(l)
          l=l+i
   20   continue
        icol=n
        jcol=nn
        mdiag=nn
   30   l=jcol
        x=zero
        if(icol.eq.irow) x=one/w(irow)
        k=n
   40   if(k.eq.irow) go to 50
        x=x-w(k)*c(l)
        k=k-1
        l=l-1
        if(l.gt.mdiag) l=l-k+1
        go to 40
   50   c(l)=x/w(irow)
        if(icol.eq.irow) go to 80
        mdiag=mdiag-icol
        icol=icol-1
        jcol=jcol-1
        go to 30
   60   do 70 j=irow,n
          c(l)=zero
          l=l+j
   70   continue
   80   ndiag=ndiag-irow
        irow=irow-1
        if(irow.ne.0) go to 10
        return
        end



      SUBROUTINE CHOL (A, N, NN, U, NULLTY, IFAULT)
C 
C       Algorithm AS6, Applied Statistics, vol.17, (1968)
C 
C       Given a symmetric matrix order n as lower triangle in a( )
C       calculates an upper triangle, u( ), such that uprime * u = a.
C       a must be positive semi-definite.  eta is set to multiplying
C       factor determining effective zero for pivot.
C 
C       arguments:-
C       a()     = input, a +ve definite matrix stored in lower-triangula
C                 form.
C       n       = input, the order of a
C       nn      = input, the size of the a and u arrays      n*(n+1)/2
C       u()     = output, a lower triangular matrix such that u*u' = a.
C                 a & u may occupy the same locations.
C       nullty  = output, the rank deficiency of a.
C       ifault  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 3 if nn < n*(n+1)/2
C                       = 0 otherwise
C 
C***********************************************************************
C 
      DOUBLE PRECISION A(NN), U(NN), ETA, ETA2, X, W, ZERO
C 
C       The value of eta will depend on the word-length of the
C       computer being used.  See introductory text.
C 
      DATA ETA, ZERO/1.D-9, 0.0D0/
C 
      IFAULT = 1
      IF (N.LE.0) RETURN
      IFAULT = 3
      IF (NN.LT.N*(N+1)/2) RETURN
      IFAULT = 2
      NULLTY = 0
      J = 1
      K = 0
      ETA2 = ETA*ETA
      II = 0
C 
C       Factorize column by column, icol = column no.
C 
      DO 80 ICOL = 1,N
        II = II+ICOL
        X = ETA2*A(II)
        L = 0
        KK = 0
C 
C       IROW = row number within column ICOL
C 
        DO 40 IROW = 1,ICOL
          KK = KK+IROW
          K = K+1
          W = A(K)
          M = J
          DO 10 I = 1,IROW
            L = L+1
            IF (I.EQ.IROW) GO TO 20
            W = W-U(L)*U(M)
            M = M+1
 10       CONTINUE
 20       IF (IROW.EQ.ICOL) GO TO 50
          IF (U(L).EQ.ZERO) GO TO 30
          U(K) = W/U(L)
          GO TO 40
 30       IF (W*W.GT.ABS(X*A(KK))) RETURN
          U(K) = ZERO
 40     CONTINUE
 50     IF (ABS(W).LE.ABS(ETA*A(K))) GO TO 60
        IF (W.LT.ZERO) RETURN
        U(K) = SQRT(W)
        GO TO 70
 60     U(K) = ZERO
        NULLTY = NULLTY+1
 70     J = J+ICOL
 80   CONTINUE
      IFAULT = 0
      END
