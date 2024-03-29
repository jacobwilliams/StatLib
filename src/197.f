      SUBROUTINE FLIKAM(P,MP,Q,MQ,W,E,N,SUMSQ,FACT,VW,VL,
     *                  MRP1,VK,MR,TOLER,IFAULT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     ALGORITHM AS 197  APPLIED STATISTICS (1984) VOL.33, NO.1
C
C     COMPUTES THE LIKELIHOOD FUNCTION OF AN AUTOREGRESSIVE-MOVING
C     AVERAGE PROCESS, EXPRESSED AS FACT * SUMSQ.
C
C     3/96 CAC  FROM STATLIB; CONVERTED TO DOUBLE PRECISION, NO TABS,
C     UPPERCASE, MIN --> MIN0, MAX --> MAX0 TO AGREE WITH JOURNAL
C     LISTING.  FIXED 4 BUGS (SEE COMMENTS WITH ORIGINAL BELOW).
C
C     TESTED SUCCESSFULLY ON THE FOLLOWING MODELS AND BENCHMARK RESULTS:
C      MA(1)  BOX-JENKINS(1976) SERIES A  -  ARIMA(0,1,1).
C      MA(2)  OSBORN(1976) USING 28-OBS SECTIONS OF B-J SERIES C
C      ARMA(1,1)  BOX-JENKINS(1976) SERIES A  (MEAN REMOVED FIRST).
C      ARMA(2,1)  DAVID CUSHMAN'S SERIES.  COMPARED WITH TSP PROC CODE
C                 THAT IMPLEMENTED NEWBOLD(1974) METHOD.
C      ARMA(3,1), ARMA(4,1), ARMA(5,1)   COMPARED WITH GAUSS RESULTS USING
C                 ANSLEY(1979) METHOD, WHICH IS NOT COMPLETELY EXACT ML
C                 (IS MISSING THE QUADRATIC FORM IN THE INITIAL RESIDUALS
C                  MENTIONED IN JUDGE, ET AL (1980) THEORY AND PRACTICE).
C
C     ORIGINAL CODE FROM STATLIB USED P(MP),Q(MQ) BELOW, BUT THIS
C     FAILS ON COMPUTERS LIKE VAX/VMS WHEN MP OR MQ ARE ZERO.
      DIMENSION P(1),Q(1),W(N),E(N),VW(MRP1),VL(MRP1),VK(MR)
C
      DATA EPSIL1 /1.D-10/
      DATA ZERO, P0625, ONE, TWO, FOUR, SIXTEN /
     *   0.D0, 0.0625D0, 1.D0, 2.D0, 4.D0, 16.D0/
C
      FACT = ZERO
      DETMAN = ONE
      DETCAR = ZERO
      SUMSQ = ZERO
      MXPQ = MAX0(MP, MQ)
      MXPQP1 = MXPQ + 1
      MQP1 = MQ + 1
      MPP1 = MP + 1
C
C      CALCULATION OF THE AUTOCOVARIANCE FUNCTION OF A PROCESS WITH
C      UNIT INNOVATION VARIANCE (VW) AND THE COVARIANCE BETWEEN THE
C      VARIABLE AND THE LAGGED INNOVATIONS (VL).
C
      CALL TWACF(P,MP,Q,MQ,VW,MXPQP1,VL,MXPQP1,VK,MXPQ,IFAULT)
      IF (MR .NE. MAX0(MP, MQP1)) IFAULT = 6
      IF (MRP1 .NE. MR + 1) IFAULT = 7
      IF (IFAULT .GT. 0) RETURN
C
C      COMPUTATION OF THE FIRST COLUMN OF MATRIX P (VK)
C
      VK(1) = VW(1)
      IF (MR .EQ. 1) GO TO 150
      DO 140 K = 2, MR
        VK(K) = ZERO
        IF (K .GT. MP) GO TO 120
        DO 110 J = K, MP
          JP2MK = J + 2 - K
C          ORIGINAL CODE FROM STATLIB
C          VK(K) = VK(K) + P(J) + VW(JP2MK)
C          FIXED TO AGREE WITH JOURNAL LISTING AND EQN (10.2)
          VK(K) = VK(K) + P(J) * VW(JP2MK)
  110   CONTINUE
  120   IF (K .GT. MQP1) GO TO 140
        DO 130 J = K, MQP1
          JP1MK = J + 1 - K
          VK(K) = VK(K) - Q(J-1) * VL(JP1MK)
  130   CONTINUE
  140 CONTINUE
C
C      COMPUTATION OF THE INITIAL VECTORS L AND K (VL, VK).
C
  150 R = VK(1)
      VL(MR) = ZERO
      DO 160 J = 1, MR
        VW(J) = ZERO
        IF (J .NE. MR) VL(J) = VK(J+1)
        IF (J .LE. MP) VL(J) = VL(J) + P(J) * R
        VK(J) = VL(J)
  160 CONTINUE
C
C      INITIALIZATION
C
      LAST = MPP1 - MQ
      LOOP = MP
      JFROM = MPP1
      VW(MPP1) = ZERO
      VL(MXPQP1) = ZERO
C
C      EXIT IF NO OBSERVATION, ELSE LOOP ON TIME.
C
      IF (N .LE. 0) GO TO 500
      DO 290 I = 1, N
C
C      TEST FOR SKIPPED UPDATING
C
        IF (I .NE. LAST) GO TO 170
        LOOP = MIN0(MP, MQ)
        JFROM = LOOP + 1
C
C      TEST FOR SWITCHING
C
        IF (MQ .LE. 0) GO TO 300
  170   IF (R .LE. EPSIL1) GO TO 400
        IF (DABS(R - ONE) .LT. TOLER .AND. I .GT. MXPQ) GO TO 300
C
C      UPDATING SCALARS
C
        DETMAN = DETMAN * R
  190   IF (DABS(DETMAN) .LT. ONE) GO TO 200
        DETMAN = DETMAN * P0625
        DETCAR = DETCAR + FOUR
        GO TO 190
  200   IF (DABS(DETMAN) .GE. P0625) GO TO 210
        DETMAN = DETMAN * SIXTEN
        DETCAR = DETCAR - FOUR
        GO TO 200
  210   VW1 = VW(1)
        A = W(I) - VW1
        E(I) = A / DSQRT(R)
        AOR = A / R
        SUMSQ = SUMSQ + A * AOR
        VL1 = VL(1)
        ALF = VL1 / R
        R = R - ALF * VL1
        IF (LOOP .EQ. 0) GO TO 230
C
C      UPDATING VECTORS
C
        DO 220 J = 1, LOOP
          FLJ = VL(J+1) + P(J) * VL1
          VW(J) = VW(J+1) + P(J) * VW1 + AOR * VK(J)
          VL(J) = FLJ - ALF * VK(J)
          VK(J) = VK(J) - ALF * FLJ
  220   CONTINUE
  230   IF (JFROM .GT. MQ) GO TO 250
        DO 240 J = JFROM, MQ
          VW(J) = VW(J+1) + AOR * VK(J)
          VL(J) = VL(J+1) - ALF * VK(J)
          VK(J) = VK(J) - ALF * VL(J+1)
  240   CONTINUE
  250   IF (JFROM .GT. MP) GO TO 270
        DO 260 J = JFROM, MP
  260   VW(J) = VW(J+1) + P(J) * W(I)
  270   CONTINUE
  290 CONTINUE
      GO TO 390
C
C      QUICK RECURSIONS
C
  300 NEXTI = I
      IFAULT = -NEXTI
      DO 310 I = NEXTI, N
  310 E(I) = W(I)
      IF (MP .EQ. 0) GO TO 340
      DO 330 I = NEXTI, N
        DO 320 J = 1, MP
          IMJ = I - J
          E(I) = E(I) - P(J) + W(IMJ)
  320   CONTINUE
  330 CONTINUE
  340 IF (MQ .EQ. 0) GO TO 370
      DO 360 I = NEXTI, N
        DO 350 J = 1, MQ
          IMJ = I - J
          E(I) = E(I) + Q(J) * E(IMJ)
  350   CONTINUE
  360 CONTINUE
C
C      RETURN SUM OF SQUARES AND DETERMINANT
C
  370 DO 380 I = NEXTI, N
  380 SUMSQ = SUMSQ + E(I) * E(I)
C
C      CODE FOR CONDITIONAL SUM OF SQUARES.
C      REPLACES ALL EXECUTABLE UP TO AND INCLUDING THAT LABELLED 380.
C
C            FACT = ZERO
C            DETMAN = ONE
C            DETCAR = ZERO
C            SUMSQ = ZERO
C            MXPQ = MAX0(MP, MQ)
C            DO 380 I = MXPQ, N
C              E(I) = W(I)
C              IF (MP .LE. 0) GO TO 340
C              DO 320 J = 1, MP
C                IMJ = I - J
C                E(I) = E(I) - P(J) * W(IMJ)
C        320        CONTINUE
C        340        IF (MQ .LE. 0) GO TO 380
C              DO 350 J = 1, MQ
C                IMJ = I - J
C                E(I) = E(I) + Q(J) * E(IMJ)
C        350        CONTINUE
C        380      SUMSQ = SUMSQ + E(I) * E(I)
C
  390 FN = N
      FACT = DETMAN ** (ONE / FN) * TWO ** (DETCAR / FN)
      RETURN
C
C      EXECUTION ERRORS
C
  400 IFAULT = 8
      RETURN
  500 IFAULT = 9
      RETURN
      END
C-----------------------
      SUBROUTINE TWACF(P,MP,Q,MQ,ACF,MA,CVLI,MXPQP1,ALPHA,MXPQ,IFAULT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C      ALGORITHM AS 197.1 APPLIED STATISTICS (1984) VOL.33, NO. 1
C
C      IMPLEMENTATION OF THE ALGORITHM OF G. TUNNICLIFFE WILSON
C      (J. STATIST. COMPUT. SIMUL. 8, 1979, 301-309) FOR THE COMPUTATION
C      OF THE AUTOCOVARIANCE FUNCTION (ACF) OF AN ARMA PROCESS OF ORDER
C      (MP, MQ) AND UNIT INNOVATION VARIANCE.   THE AUTOREGRESSIVE AND
C      MOVING AVERAGE COEFFICIENTS ARE STORED IN VECTORS P AND Q, USING
C      BOX AND JENKINS NOTATION.   ON OUTPUT, VECTOR CVLI CONTAINS THE
C      COVARIANCES BETWEEN THE VARIABLE AND THE (K-1)-LAGGED INNOVATION
C      FOR K = 1, ..., MQ+1.
C
C     ORIGINAL CODE FROM STATLIB USED P(MP),Q(MQ) BELOW, BUT THIS
C     FAILS ON COMPUTERS LIKE VAX/VMS WHEN MP OR MQ ARE ZERO.
      DIMENSION P(1), Q(1), ACF(MA), CVLI(MXPQP1), ALPHA(MXPQ)
C
      DATA EPSIL2 /1.D-10/
      DATA ZERO, HALF, ONE, TWO /0.D0, 0.5D0, 1.D0, 2.D0/
C
      IFAULT = 0
      IF (MP .LT. 0 .OR. MQ .LT. 0) IFAULT = 1
      IF (MXPQ .NE. MAX0(MP, MQ)) IFAULT = 2
      IF (MXPQP1 .NE. MXPQ + 1) IFAULT = 3
      IF (MA .LT. MXPQP1) IFAULT = 4
      IF (IFAULT .GT. 0) RETURN
C
C      INITIALIZATION, AND RETURN IF MP = MQ = 0
C
      ACF(1) = ONE
      CVLI(1) = ONE
      IF (MA .EQ. 1) RETURN
      DO 10 I = 2, MA
   10 ACF(I) = ZERO
      IF (MXPQP1 .EQ. 1) RETURN
      DO 20 I = 2, MXPQP1
   20 CVLI(I) = ZERO
      DO 90 K = 1, MXPQ
   90 ALPHA(K) = ZERO
C
C      COMPUTATION OF THE A.C.F. OF THE MOVING AVERAGE PART,
C      STORED IN ACF.
C
      IF (MQ .EQ. 0) GO TO 180
      DO 130 K = 1, MQ
        CVLI(K+1) = -Q(K)
        ACF(K+1) = -Q(K)
        KC = MQ - K
        IF (KC .EQ. 0) GO TO 120
        DO 110 J = 1, KC
          JPK = J + K
          ACF(K+1) = ACF(K+1) + Q(J) * Q(JPK)
  110   CONTINUE
  120   ACF(1) = ACF(1) + Q(K) * Q(K)
  130 CONTINUE
C
C      INITIALIZATION OF CVLI = T.W.-S PHI -- RETURN IF MP = 0.
C
  180 IF (MP .EQ. 0) RETURN
      DO 190 K = 1, MP
        ALPHA(K) = P(K)
        CVLI(K) = P(K)
  190 CONTINUE
C
C      COMPUTATION OF T.W.-S ALPHA AND DELTA (DELTA STORED IN ACF
C      WHICH IS GRADUALLY OVERWRITTEN).
C
      DO 290 K = 1, MXPQ
        KC = MXPQ - K
        IF (KC .GE. MP) GO TO 240
        DIV = ONE - ALPHA(KC+1) * ALPHA(KC+1)
        IF (DIV .LE. EPSIL2) GO TO 700
        IF (KC .EQ. 0) GO TO 290
        DO 230 J = 1, KC
          KCP1MJ = KC + 1 - J
          ALPHA(J) = (CVLI(J) + ALPHA(KC+1) * CVLI(KCP1MJ)) / DIV
  230   CONTINUE
  240   IF (KC .GE. MQ) GO TO 260
        J1 = MAX0(KC + 1 - MP, 1)
        DO 250 J = J1, KC
          KCP1MJ = KC + 1 - J
C          ORIGINAL CODE FROM STATLIB
C          ACF(J+1) = ACF(J+1) + ACF(K+2) * ALPHA(KCP1MJ)
C          FIXED TO AGREE WITH JOURNAL LISTING
          ACF(J+1) = ACF(J+1) + ACF(KC+2) * ALPHA(KCP1MJ)
  250   CONTINUE
  260   IF (KC .GE. MP) GO TO 290
        DO 270 J = 1, KC
  270   CVLI(J) = ALPHA(J)
  290 CONTINUE
C
C      COMPUTATION OF T.W.-S NU (NU IS STORED IN CVLI, COPIED INTO
C      ACF).
C
      ACF(1) = HALF * ACF(1)
      DO 330 K = 1, MXPQ
        IF (K .GT. MP) GO TO 330
        KP1 = K + 1
        DIV = ONE - ALPHA(K) * ALPHA(K)
        DO 310 J = 1, KP1
          KP2MJ = K + 2 - J
          CVLI(J) = (ACF(J) + ALPHA(K) * ACF(KP2MJ)) / DIV
  310   CONTINUE
        DO 320 J = 1, KP1
  320   ACF(J) = CVLI(J)
  330 CONTINUE
C
C      COMPUTATION OF ACF (ACF IS GRADUALLY OVERWRITTEN).
C
      DO 430 I = 1, MA
        MIIM1P = MIN0(I-1, MP)
        IF (MIIM1P .EQ. 0) GO TO 430
        DO 420 J = 1, MIIM1P
          IMJ = I - J
          ACF(I) = ACF(I) + P(J) * ACF(IMJ)
  420   CONTINUE
  430 CONTINUE
      ACF(1) = ACF(1) * TWO
C
C      COMPUTATION OF CVLI - RETURN WHEN MQ = 0.
C
      CVLI(1) = ONE
      IF (MQ .LE. 0) GO TO 600
      DO 530 K = 1, MQ
        CVLI(K+1) = -Q(K)
        IF (MP .EQ. 0) GO TO 530
        MIKP = MIN0(K, MP)
        DO 520 J = 1, MIKP
          KP1MJ = K + 1 - J
          CVLI(K+1) = CVLI(K+1) + P(J) * CVLI(KP1MJ)
  520   CONTINUE
  530 CONTINUE

  600 RETURN
C
C      EXECUTION ERROR DUE TO (NEAR) NON-STATIONARITY.
C
  700 IFAULT = 5
      RETURN
      END
