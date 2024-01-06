      REAL FUNCTION BETANC(X, A, B, LAMBDA, IFAULT)
C
C     ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
C     Incorporates modification AS R84 from AS vol. 39, pp311-2, 1990
C
C     Returns the cumulative probability of X for the non-central beta
C     distribution with parameters A, B and non-centrality LAMBDA
C
C     Auxiliary routines required: ALOGAM - log-gamma function (ACM
C     291 or AS 245), and BETAIN - incomplete-beta function (AS 63)
C
      REAL A, AX, B, BETA, C, ERRBD, ERRMAX, GX, HALF, LAMBDA, ONE, Q,
     *     SUMQ, TEMP, X, XJ, ZERO
      REAL A0, X0, UALPHA
C
C     Change ERRMAX and ITRMAX if desired ...
C
      DATA ERRMAX, ITRMAX /1.0E-6, 100/, UALPHA /5.0/
C
      DATA ZERO, HALF, ONE /0.0, 0.5, 1.0/
C
      BETANC = X
C
      IFAULT = 2
      IF (LAMBDA .LT. ZERO .OR. A .LE. ZERO .OR. B .LE. ZERO) RETURN
      IFAULT = 3
      IF (X .LT. ZERO .OR. X .GT. ONE) RETURN
      IFAULT = 0
      IF (X .EQ. ZERO .OR. X .EQ. ONE) RETURN
C
      C = LAMBDA * HALF
C
C     Initialize the series ...
C
      X0 = INT(MAX(C - UALPHA*SQRT(C), ZERO)
      A0 = A + X0
      BETA = ALOGAM(A0, IFAULT) + ALOGAM(B, IFAULT) - 
     *       ALOGAM(A0+B, IFAULT)
      TEMP = BETAIN(X, A0, B, BETA, IFAULT)
      GX = EXP(A0 * LOG(X) + B * LOG(ONE - X) - BETA - LOG(A0))
      IF (A0 .GT. A) THEN
	Q = EXP(-C + X0*LOG(C) - ALOGAM(X0 + ONE, IFAULT)
      ELSE
        Q = EXP(-C)
      END IF
      XJ = ZERO
      AX = Q * TEMP
      SUMQ = ONE - Q
      BETANC = AX
C
C     Recur over subsequent terms until convergence is achieved...
C
   10 XJ = XJ + ONE
      TEMP = TEMP - GX
      GX = X * (A + B + XJ - ONE) * GX / (A + XJ)
      Q = Q * C / XJ
      SUMQ = SUMQ - Q
      AX = TEMP * Q
      BETANC = BETANC + AX
C
C     Check for convergence and act accordingly...
C
      ERRBD = (TEMP - GX) * SUMQ
      IF ((INT(XJ) .LT. ITRMAX) .AND. (ERRBD .GT. ERRMAX)) GO TO 10
      IF (ERRBD .GT. ERRMAX) IFAULT = 1
C
      RETURN
      END
