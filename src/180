      REAL FUNCTION TRIMSD(X, N, K, IW, IFAULT)
C
C     ALGORITHM AS 180  APPL. STATIST. (1982) VOL.31, NO.2
C
C     The value returned is a linear estimate of the standard deviation
C     from a sample of N values of X symmetrically trimmed to retain K
C     values.
C
      INTEGER IW(N)
      REAL X(N), QTR, HALF, P74, P98, AMID, B, P, SUM, C(6)
      DATA C/1.06683, 1.72534, 0.403049, 1.32458, -3.91154, 17.5265/,
     *  QTR/0.25/, HALF/0.5/, P74/0.74/, P98/0.98/
C
      IFAULT = 1
      TRIMSD = -1.0
      IF (N .LT. 4) RETURN
      IFAULT = 2
      P = K / FLOAT(N)
      IF (P .LT. HALF .OR. P .GT. P98) RETURN
      IFAULT = 3
      IF (MOD(N-K, 2) .NE. 0) RETURN
      IFAULT = 0
C
      ILO = (N-K)/2 + 1
      IHI = N - ILO + 1
      AMID = N + 1
C
C     Get ranks of X in IW
C
      IW(1) = 1
      DO 10 I = 2, N
	IW(I) = 1
	DO 10 J = 1, I-1
	  IF (X(I) .GE. X(J)) GO TO 9
	  IW(J) = IW(J) + 1
	  GO TO 10
    9     IW(I) = IW(I) + 1
   10 CONTINUE
C
      SUM = 0.0
      DO 12 I = 1, N
   12 IF (IW(I) .GE. ILO .AND. IW(I) .LE. IHI) SUM = 
     *            SUM + (2.0 * IW(I) - AMID) * X(I)
C
C     Calculate the unbiassing factor
C
      P = P74 - P
      B = ((((C(6)*P + C(5))*P + C(4))*P + C(3))*P + C(2))*P + C(1)
C
      TRIMSD = EXP(B) * SUM / (QTR * (2*K*(2*K - 1)))
      RETURN
      END
