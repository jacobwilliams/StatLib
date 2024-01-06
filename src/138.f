	SUBROUTINE EM(N, X1, X2, P, XMEAN, XSIGMA, E1, E2, MAXITS,
     *  COV, NOBS, K, IFAULT)
C
C	  ALGORITHM AS 138 APPL. STATIST. (1979) VOL.28, NO.2
C
C	  COMPUTE THE MAXIMUM LIKELIHOOD ESTIMATES OF THE MEAN
C	  AND STANDARD DEVIATION OF A SINGLE NORMAL POPULATION,
C	  THE DATA MAY CONTAIN CENSORED OR CONFINED OBSERVATIONS.
C
	DIMENSION X1(N), X2(N), COV(2, 2)
	INTEGER P(N), NOBS(4)
C
	DATA C /0.39894228/, ONEPLS /1.0001/, TOL /0.00001/,
     *  TOLL /0.00001/
C
	IFAULT = -2
	IF (N .LT. 2) RETURN
C
C	  INITIALIZE COUNTERS
C
	K = 0
	SUM = 0.0
	SUM2 = 0.0
	SUMG = 0.0
	SUMG2 = 0.0
	IP = 0
	IQ = 0
	IR = 0
	IS = 0
C
C	  EXACTLY SPECIFIED OBSERVATIONS ARE REMOVED,
C	  THE REMAINING DATA PACKED INTO FIRST PART OF ARRAY X
C
	IFAULT = -4
	DO 200 I = 1, N
	IPT = P(I)
	IF (IPT .EQ. 0) GOTO 100
	IF (IPT .EQ. 2 .AND. ABS(X1(I) - X2(I)) .LE. ABS(X1(I) * TOLL))
     *  GOTO 100
C
C	  OBSERVATION NOT EXACTLY SPECIFIED
C
	IS = IS + 1
	P(IS) = IPT
	X1(IS) = X1(I)
C
C	  HANDLE GROUPED DATA
C
	IF (IPT .NE. 2) GOTO 50
	IQ = IQ + 1
	IF (X1(I) .GT. X2(I)) RETURN
	X2(IS) = X2(I)
	XTEMP = 0.5 * (X1(I) + X2(I))
	SUMG = SUMG + XTEMP
	SUMG2 = SUMG2 + XTEMP ** 2
	GOTO 200
C	  ACCUMULATE NUMBER OF OBSERVATIONS CENSORED ON THE RIGHT
C
50	IF (IPT .EQ. 1) IR = IR + 1
	GOTO 200
C
C	  HANDLE EXACTLY-SPECIFIED OBSERVATIONS
C
100	IP = IP + 1
	XTEMP = X1(I)
	SUM = SUM + XTEMP
	SUM2 = SUM2 + XTEMP ** 2
200	CONTINUE
C
C	  INITIAL PASS THROUGH DATA COMPLETED
C
	NOBS(1) = IP
	NOBS(2) = IR
	NOBS(3) = N - IP - IR - IQ
	NOBS(4) = IQ
	RIM = IP + IQ
	IF (IP .EQ. N) GOTO 230
	IF (XSIGMA .GT. 0.0) GOTO 350
	IF (RIM .GT. ONEPLS) GOTO 250
C
C	  AT MOST ONE OBSERVATION HAS BEEN EXACTLY
C	  SPECIFIED OR CONFINED
C
	XMEAN = 1.0
	XSIGMA = 1.0
	GOTO 350
C
C	  ALL OBSERVATIONS EXACTLY SPECIFIED
C
230	XMEAN = SUM / RIM
	XSIGMA = SQRT((SUM2 - RIM * XMEAN ** 2) / RIM)
	COV(1, 1) = XSIGMA ** 2 / RIM
	COV(2, 2) = COV(1, 1) * 0.5
	COV(1, 2) = 0.0
	COV(2, 1) = 0.0
C
C	  NORMAL RETURN
C
240	IFAULT = 0
	RETURN
C
C	  OBTAIN INITIAL ESTIMATES
C
250	XMEAN = (SUM + SUMG) / RIM
	XSIGMA = SQRT((SUM2 + SUMG2 - RIM * XMEAN ** 2) / RIM)
C
C	  INITIALIZE BEFORE STARTING FIRST ITERATION
C
350	RP = IP
	RN = N
C
C	  START OF ITERATION CYCLE,
C	  ESTIMATE CONDITIONAL EXPECTATION OF CONFINED AND
C	  CENSORED OBSERVATIONS
C
	IFAULT = -3
400	TS = SUM
	SUMG2 = SUM2
	TD = RP
	DO 610 I = 1, IS
	YS = (X1(I) - XMEAN) / XSIGMA
	IF (P(I) - 1) 500, 450, 550
C
C	  OBSERVATION CENSORED ON THE RIGHT
C
450	CALL RMILLS(YS, F, TOL)
	W = XMEAN + XSIGMA * F
	TD = TD + F * (F - YS)
	GOTO 600
C
C	  OBSERVATION CENSORED ON THE LEFT
C
500	CALL RMILLS(-YS, F, TOL)
	W = XMEAN - XSIGMA * F
	TD = TD + F * (F + YS)
	GOTO 600
C
C	  CONFINED OBSERVATION.
C	  USE MILLS RATIO RECIPROCAL TO COMPUTE PROBABILITY
C	  INTEGRALS THAT ARE REQUIRED,
C	  AS IN ORIGINAL ALGORITHM ASSUMING X1(I) IS
C	  NEVER GREATER THAN X2(I) FOR CONFINED OBSERVATIONS
C
550	YN = EXP(-0.5 * YS ** 2) * C
	CALL RMILLS(YS, F, TOL)
	YQ = YN / F
	YSU = (X2(I) - XMEAN) / XSIGMA
	YNU = EXP(-0.5 * YSU ** 2) * C
	CALL RMILLS(YSU, FU, TOL)
	YQU = YNU / FU
	YD = YQ - YQU
C
C	  IF INTEGRAL NOT EQUAL TO ZERO, CARRY ON
C
	IF (YD .LT. TOLL) RETURN
	A = (YN - YNU) / YD
	W = XMEAN + XSIGMA * A
	TD = TD + (A ** 2 + (YSU * YNU - YS * YN) / YD)
C
600	TS = TS + W
	SUMG2 = SUMG2 + W ** 2
610	CONTINUE
C
C	  CALCULATE NEW ESTIMATES
C
	XNEW = TS / RN
	YNEW = SQRT((SUMG2 + RN * XMEAN ** 2 - 2.0 * TS * XMEAN) / TD)
	K = K + 1
	IF (ABS(XNEW - XMEAN) .LT. E1 .AND. ABS(YNEW - XSIGMA) .LT. E2)
     *  GOTO 700
	IF (K .GE. MAXITS) GOTO 650
	XMEAN = XNEW
	XSIGMA = YNEW
	GOTO 400
C
C	  MAXIMUM NUMBER OF ITERATIONS EXCEEDED
C
650	IFAULT = -1
	COV(1, 1) = 0.0
	COV(2, 2) = 0.0
	COV(1, 2) = XNEW - XMEAN
	COV(2, 1) = YNEW - XSIGMA
	RETURN
C
C	  CONVERGENCE OBTAINED
C
700	XMEAN = XNEW
	XSIGMA = YNEW
	XSIG2 = XSIGMA ** 2
C
C	  CALCULATE VARIANCE-COVARIANCE MATRIX
C
	X11 = RP
	X12 = (SUM - RP * XMEAN) / XSIGMA
	X22 = RP + (SUM2 + RP * XMEAN ** 2 - 2.0 * SUM * XMEAN) / XSIG2
	DO 800 I = 1, IS
	YS = (X1(I) - XMEAN) / XSIGMA
	IF (P(I) - 1) 740, 710, 770
710	CALL RMILLS(YS, F, TOL)
C
C	  OBSERVATION CENSORED ON THE RIGHT
C
	FL = F * (F - YS)
730	X11 = X11 + FL
	FL = FL * YS
	X12 = X12 + FL
	FL = FL * YS
	X22 = X22 + FL
	GOTO 800
C
740	CALL RMILLS(-YS, F, TOL)
C
C	  OBSERVATION CENSORED ON THE LEFT
C
	FL = F * (F + YS)
	GOTO 730
C
770	CALL RMILLS(YS, F, TOL)
C
C	  OBSERVATION CONFINED BETWEEN 2 FINITE LIMITS
C
	YN = EXP(-0.5 * YS ** 2) * C
	YQ = YN / F
	YSU = (X2(I) - XMEAN) / XSIGMA
	CALL RMILLS(YSU, FU, TOL)
	YNU = EXP(-0.5 * YSU ** 2) * C
	YQU = YNU / FU
	YD = YQ - YQU
	A = (YN - YNU) / YD
	B = (YNU * YSU - YN * YS) / YD
	X11 = X11 + A ** 2 + B
	B1 = (YS ** 2 * YN - YSU ** 2 * YNU) / YD
	X12 = X12 - A * B - B1
	B1 = (YS ** 3 * YN - YSU ** 3 * YNU) / YD
	X22 = X22 - B1 + B ** 2
800	CONTINUE
	CONST = XSIG2 / (X11 * X22 - X12 * X12)
	COV(1, 1) = CONST * X22
	COV(2, 2) = CONST * X11
	COV(1, 2) = -CONST * X12
	COV(2, 1) = COV(1, 2)
	GOTO 240
	END
C
	SUBROUTINE RMILLS(X, FUNC, TOL)
C
C	  ALGORITHM AS 138.1 APPL. STATIST. (1979) VOL.28 NO.2
C
C	  COMPUTE THE RECIPROCAL OF MILLS RATIO
C
	DATA FPI /1.2533141/, FPII /0.7978846/
C
	FUNC = 0.0
	IF (X .LT. -10.0) RETURN
	FUNC = FPII
	Y = ABS(X)
	IF (Y .LT. 0.000001) RETURN
	SGN = 1.0
	IF (X .LT. 0.0) SGN = -1.0
	IF (Y .GT. 2.0) GOTO 100
	S = 0.0
	A = 1.0
	T = Y
	R = Y
	B = Y ** 2
40	A = A + 2.0
	S = T
	R = R * B / A
	T = T + R
	IF (R .GT. TOL) GOTO 40
	FUNC = 1.0 / (FPI * EXP(0.5 * B) - SGN * T)
	RETURN
100	A = 2.0
	B1 = Y
	S = Y
	A1 = Y ** 2 + 1.0
	A2 = Y * (A1 + 2.0)
	B2 = A1 + 1.0
	T = A2 / B2
140	A = A + 1.0
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
	IF (SGN .LT. 0.0) FUNC =
     *  T / (2.0 * FPI * EXP(0.5 * Y ** 2) * T - 1.0)
	RETURN
	END
