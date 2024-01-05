      SUBROUTINE SSPINV(X, XM, XSSPI, E, WT, SUMWT, NP, NPP, IFAULT)
C
C        ALGORITHM AS240  APPL. STATIST. (1988) VOL. 37, NO. 3
C
C        THIS SUBROUTINE UPDATES THE MEAN VECTOR XM, AND THE INVERSE
C        MATRIX OF CORRECTED SUMS OF SQUARES AND PRODUCTS XSSPI.
C        A DATA VECTOR X WITH WEIGHT WT .GT. 0 IS INCLUDED,
C        OR EXCLUDED WHEN WT .LT. 0.
C        IFAULT = 2, IF UPDATED SUMWT .LE. 0,
C               = 1, IF MANIPULATED MATRIX XSSPI DOESN'T EXIST,
C               = 0, OTHERWISE.
C
      REAL X(NP), XM(NP), E(NP), XSSPI(NPP), WT, SUMWT, ZERO, ONE, EPS,
     *  A, B, C, ZABS
C
      DATA ZERO, ONE, EPS /0.0, 1.0, 1.0E-5/
C
      ZABS(A) = ABS(A)
C
      IFAULT = 2
      SUMWT = SUMWT + WT
      IF (SUMWT .LE. ZERO) RETURN
      IFAULT = 1
C
C        UPDATE MEANS AND INVERSE MATRIX OF SUMS OF SQUARES AND
C        PRODUCTS.
C
      B = WT / SUMWT
      C = WT - B * WT
      DO 10 I = 1, NP
      X(I) = X(I) - XM(I)
   10 XM(I) = XM(I) + B * X(I)
C
      DO 30 I = 1, NP - 1
      E(I) = ZERO
      K = I * (I - 1) / 2
      DO 20 J = 1, I
      K = K + 1
   20 E(I) = E(I) + XSSPI(K) * X(J)
      DO 25 J = I + 1, NP
      K = K + J - 1
   25 E(I) = E(I) + XSSPI(K) * X(J)
   30 CONTINUE
      E(NP) = ZERO
      K = NP * (NP - 1) / 2
      DO 35 J = 1, NP
      K = K + 1
   35 E(NP) = E(NP) + XSSPI(K) * X(J)
C
      A = ZERO
      DO 40 I = 1, NP
   40 A = A + X(I) * E(I)
      A = A + ONE / C
      IF (ZABS(A) .LE. EPS) RETURN
      IFAULT = 0
C
      K = 0
      DO 70 I = 1, NP
      DO 70 J = 1, I
      K = K + 1
      XSSPI(K) = XSSPI(K) - E(I) * E(J) / A
   70 CONTINUE
      RETURN
      END
