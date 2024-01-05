      SUBROUTINE ISRE(K, NH, XO, WO, X, Y, Z, WX, WY, WZ, XA, IFAULT)
C
C        ALGORITHM AS 257.1  APPL.STATIST. (1990), VOL.39, NO.3
C
C        Subroutine to compute isotonic regression for umbrella
C        ordering with known peak
C        Algorithm AMALGM (AS 149, Cran, 1980) is called
C
      INTEGER K, NH, IFAULT
      REAL XO(K), WO(K), X(K), Y(K), Z(K), WX(K), WY(K), WZ(K), XA(K)
      INTEGER I, I1, IH, IH1, IHIGH, ILOW, IXMAX, J, L, N
      REAL S, TOL, XMAX
      DATA TOL / 1.0E-6 /
C
C        Check input data K, NH and WO
C
      IFAULT = 1
      IF (K .LT. 2) RETURN
      IFAULT = 3
      IF (NH .LT. 1 .OR. NH .GT. K) RETURN
      IFAULT = 2
      DO 10 I = 1, K
         IF (WO(I) .LE. TOL) RETURN
         X(I) = XO(I)
         WX(I) = WO(I)
   10 CONTINUE
      IFAULT = 0
      N = K
      IH = NH
C
C        Pool the active block with the next lower block
C        from X(1) to X(IH-1)
C
      IF (IH .LE. 2) GO TO 40
      DO 20 I = 1, IH - 1
         Y(I) = X(I)
         WY(I) = WX(I)
   20 CONTINUE
      CALL AMALGM(IH - 1, Y, WY, Z, WZ, XA, IFAULT)
      DO 30 I = 1, IH - 1
         X(I) = XA(I)
   30 CONTINUE
C
C        Pool the active block with the next higher block
C        from X(IH+1) to X(N)
C
   40 IF (N - IH .LE. 1) GO TO 70
      L = N - IH
      DO 50 I = IH + 1, N
         Y(L) = X(I)
         WY(L) = WX(I)
         L = L - 1
   50 CONTINUE
      CALL AMALGM(N - IH, Y, WY, Z, WZ, XA, IFAULT)
      L = N - IH
      DO 60 I = IH + 1, N
         X(I) = XA(L)
         L = L - 1
   60 CONTINUE
C
C        Pool adjacent violators such that X(IH) is largest
C
   70 ILOW = MAX(1, IH - 1)
      IHIGH = MIN(N, IH + 1)
      XMAX = X(IH)
      IXMAX = IH
      DO 80 I = ILOW, IHIGH
         IF (I .EQ. IH .OR. XMAX .GE. X(I)) GO TO 80
         XMAX = X(I)
         IXMAX = I
   80 CONTINUE
      IF (IXMAX .EQ. IH) GO TO 100
      IF (IXMAX .LT. IH) IH = IH - 1
      IH1 = IH + 1
      X(IH) = (WX(IH) * X(IH) + WX(IH1) * X(IH1)) / (WX(IH) + WX(IH1))
      WX(IH) = WX(IH) + WX(IH1)
      DO 90 I = IH1, N - 1
         X(I) = X(I + 1)
         WX(I) = WX(I + 1)
   90 CONTINUE
      N = N - 1
      GO TO 70
C
C        Obtain XA(1), ..., XA(K) from X(1), ..., X(N)
C
  100 I1 = 1
      DO 130 I = 1, N
         S = 0.0
         DO 110 J = I1, K
            S = S + WO(J)
            XA(J) = X(I)
            IF (ABS(S - WX(I)) .LT. TOL) GO TO 120
  110    CONTINUE
  120    I1 = J + 1
  130 CONTINUE
      RETURN
      END
      SUBROUTINE ISREU(MAXK, K, XO, WO, X, Y, Z, WX, WY, WZ, MINVAL,
     *                 XMIN, XB, RES, IFAULT)
C
C        ALGORITHM AS 257.2  APPL.STATIST. (1990), VOL.39, NO.3
C
C        Subroutine to compute isotonic regression for umbrella
C        ordering with unknown peak
C        A subroutine ISRE given in this paper is called
C
      INTEGER MAXK, K, MINVAL, IFAULT
      REAL XO(K), WO(K), X(K), Y(K), Z(K), WX(K), WY(K), WZ(K), XMIN(K),
     *     XB(MAXK, K), RES(K)
      INTEGER I, NH
      REAL RESMIN
C
C        Check that K .GE. 2
C
      IF (K .GE. 2) GO TO 10
      IFAULT = 1
      RETURN
C
C        For all cases with the peak NH = 1, ..., K
C
   10 DO 30 NH = 1, K
         CALL ISRE(K, NH, XO, WO, X, Y, Z, WX, WY, WZ, XB(1, NH),
     *             IFAULT)
         IF (IFAULT .NE. 0) RETURN
C
C        Calculate square error to obtain min
C
         RES(NH) = 0.0
         DO 20 I = 1, K
            RES(NH) = RES(NH) + WO(I) * (XB(I, NH) - XO(I)) ** 2
   20    CONTINUE
         IF (NH .NE. 1 .AND. RES(NH) .GT. RESMIN) GO TO 30
         MINVAL = NH
         RESMIN = RES(NH)
   30 CONTINUE
      DO 40 I = 1, K
         XMIN(I) = XB(I, MINVAL)
   40 CONTINUE
      RETURN
      END
