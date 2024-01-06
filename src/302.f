      SUBROUTINE MISRE(K, IP, IPP, NH, INV, COV, X, EPS, ITMAX,
     *                 IB, COVINV, XA, Y, Z, AII, A, W, WK, IT, SITA, 
     *                 STAT, IFAULT)
C
C        ALGORITHM AS 302.1 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Subroutine to compute a multiple isotonic regression
C        for umbrella ordering with known peak
C
      INTEGER I, IFAULT, II, II1, INC, IP, IP1, IPP, IT, ITMAX, J,
     *        J1, K, L, NULLTY
      REAL DET, EPS, RABS, RMAX, STAT
      REAL A(IP), AII(K), COV(IPP, K), COVINV(IPP, K), SITA(IP, K),
     *     W(IP), WK(K, 6), X(IP, K), XA(K), Y(K), Z(K)
      INTEGER IB(IP), NH(IP)
C
      EXTERNAL ISRE, SUBINV
C
      IFAULT = 0
      IF (IP .LE. 1) IFAULT = 1
      IF (IPP .NE. IP * (IP + 1) / 2) IFAULT = 6
      IF (IFAULT .NE. 0) RETURN
C
C        Find inversions of covariances
C
      IF (INV .EQ. 1)THEN 
        DO 10 I = 1, IP
          IB(I) = I
   10   CONTINUE
        DO 20 I = 1, K
          CALL SUBINV(COV(1, I), IP, IB, IP, COVINV(1, I), W, NULLTY,
     *     IFAULT, IPP, DET)
          IF (IFAULT .NE. 0) RETURN
   20   CONTINUE
        INV = 0
      ENDIF
C
      IT = 1
      IP1 = IP - 1
      DO 40 I = 2, IP
        DO 30 J = 1, K
          SITA(I, J) = X(I, J)
   30   CONTINUE
   40 CONTINUE
C
C        Find multivariate isotonic regression iteratively
C
   50 RMAX = 0.0
      DO 120 I = 1, IP
        II = I * (1 + I) / 2
C
C        For the i-th univariate isotonic regression
C
        DO 100 L = 1, K
          AII(L) = COVINV(II, L)
C
C        Set A(1),..., A(IP - 1) to the elements of A21 
C
          IF (I .GT. 1) THEN
            II1 = II - I
            DO 60 J = 1, I - 1
              II1 = II1 + 1
              A(J) =  COVINV(II1, L)
   60       CONTINUE
          ENDIF
          IF (I .LT. IP) THEN
            II1 = II + I
            INC = I + 1
            DO 70 J = I, IP1
              A(J) = COVINV(II1, L)
              II1 = II1 + INC
              INC = INC + 1
   70       CONTINUE
          ENDIF
C
C        Prepare the input data Z(1),...,Z(IP - 1) for
C        univariate isotonic regression
C
          J1 = 1
          DO 80 J = 1, IP
            IF (J .NE. I) THEN
              Y(J1) = X(J, L) - SITA(J, L)
              J1 = J1 + 1
            ENDIF
   80     CONTINUE
          Z(L) = X(I, L)
          DO 90 J = 1, IP1
            Z(L) = Z(L) + Y(J) * A(J) / AII(L)
   90     CONTINUE
  100   CONTINUE
C
C        Find a univariate isotonic regression
C
        CALL ISRE(K, NH(I), Z, AII, WK(1, 1), WK(1, 2), WK(1, 3),
     *            WK(1, 4), WK(1, 5), WK(1, 6), XA, IFAULT)
        IF (IFAULT .NE. 0) THEN
          IFAULT = IFAULT + 2
          RETURN
        ENDIF
        DO 110 J = 1, K
          IF (IT .NE. 1) THEN
            RABS = ABS(XA(J) - SITA(I,J))
            IF (RABS .GT. RMAX) RMAX = RABS
          ENDIF
          SITA(I, J) = XA(J)
  110   CONTINUE
  120 CONTINUE
      IT = IT + 1
      IF (IT .EQ. 2 .OR. RMAX .GT. EPS .AND. IT .LE. ITMAX) GOTO 50
      IF (RMAX .GT. EPS .AND. IT .GT. ITMAX) IFAULT = 7
C
C        Calculate likelihood ratio statistics
C
      STAT = 0.0
      DO 150 L = 1, K
        DO 140 I = 1, IP
          II = I * (1 + I) / 2
          STAT = STAT + COVINV(II, L) * (X(I, L) - SITA(I, L)) ** 2 
          IF (I .GT. 1) THEN
            II1 = II - I
            DO 130 J = 1, I - 1
              STAT = STAT + 2.0 * COVINV(II1+J,L) * (X(I,L) - SITA(I,L))
     *                          * (X(J, L) - SITA(J, L))
  130       CONTINUE
          ENDIF
  140   CONTINUE
  150 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MISREU(K, IP, IPP, KP, NH, COV, X, EPS, ITMAX, IB,
     *                  COVINV, XA, Y, Z, AII, A, W, WK, ITM, SITA,
     *                  SITAMIN, STATMIN, NHMIN, STATMAT, IFAULT)
C
C        ALGORITHM AS 302.2 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Subroutine to compute a multiple isotonic regression
C        for umbrella ordering with unknown peak
C
      INTEGER I, IFAULT, II, II1, INC, INV, IP, IP1, IPP, IT, ITM, 
     *        ITMAX, J, J1, K, KP, L, NULLTY
      REAL DET, EPS, RABS, RMAX, STAT, STATMIN
      REAL A(IP), AII(K), COV(IPP, K), COVINV(IPP, K), SITA(IP, K),
     *     SITAMIN(IP,K), STATMAT(KP), W(IP), WK(K, 6), X(IP, K), 
     *     XA(K), Y(K), Z(K)
      INTEGER IB(IP), NH(IP), NHMIN(IP)
C
      EXTERNAL MISRE
C
      IF (KP .LT. K ** IP) THEN
        IFAULT = 8
        RETURN
      ENDIF
C
      INV = 1
      ITM = 0
      STATMIN = 0.99E10
      NUM = 1
      DO 10 I = 1, IP
        NH(I) = 1
   10 CONTINUE
C
C        Calculate multiple isotonic regression for all possible 
C        peaks
C
   20 CALL MISRE(K, IP, IPP, NH, INV, COV, X, EPS, ITMAX,
     *           IB, COVINV, XA, Y, Z, AII, A, W, WK, IT, SITA, STAT, 
     *           IFAULT)
      IF (IFAULT .NE. 0) RETURN
      STATMAT(NUM) = STAT
      NUM = NUM + 1
      IF (IT .GT. ITM) ITM = IT
C
C        Save the best regression and peak vector
C
      IF (STAT .LE. STATMIN) THEN
        STATMIN = STAT
        DO 40 I = 1, IP
          DO 30 J = 1, K
            SITAMIN(I, J) = SITA(I, J)
   30     CONTINUE
          NHMIN(I) = NH(I)
   40   CONTINUE
      ENDIF
C
C        Set the next peak
C
      DO 50 I = 1, IP
        NH(I) = NH(I) + 1
        IF (NH(I) .LE. K) GOTO 20
        NH(I) = 1
   50 CONTINUE
C
      RETURN
      END
