      SUBROUTINE NPMLG(NI, NJ, ITER, ICYC, JN, DEPS, Y, S, R, P, DFL, X,
     *  W, M, V, AFI, LP, IFAULT)
C
C        ALGORITHM AS221  APPL. STATIST. (1986) VOL. 35, NO. 3
C
C        SUBROUTINE FOR CALCULATING THE NONPARAMETRIC MAXIMUM
C        LIKELIHOOD ESTIMATE OF THE MIXING DISTRIBUTION IN
C        A MIXTURE OF NORMALS
C
C     N.B. The dimension of array R below has been changed in accordance
C          with the correction on page 176 of vol. 39 (1990) of Applied
C          Statistics.
C
      REAL R(0:NI), P(NI), Y(NI), S(NI), M(NI), V(NI), AFI(NI), X, W,
     *  DEPS, DFL, RMIN, RMAX, RCONS, A, Q, D, DENOM, PR, FL, DIFF,
     *  POLD, ROLD, PQ, RR, AI, ANB, TEMP, DFOLD, AN
      INTEGER LP(ICYC), NI, NJ, ITER, ICYC, JN, IFAULT, JJ, IT, IC, IJ,
     *  JNN, NEXLAS, MORE
      DATA DCONS /1.0E-6/, RCONS /1.0E-2/, DINF /1.0E-6/
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, PI /6.2831853/
C
C        CHECK INPUT VALUES FOR SAMPLING VARIANCES
C
      IFAULT = 0
      DO 1 I = 1, NI
      IF (S(I) .LE. ZERO) IFAULT = 2
      IF (IFAULT .EQ. 2) GOTO 999
    1 CONTINUE
C
C        SET THE INITIAL ESTIMATE : UNIFORM DISTRIBUTION OVER
C        THE RANGE OF OBSERVED VALUES
C
      FNJ = FLOAT(NJ)
      FNI = FLOAT(NI)
      POLD = ONE / FNJ
      RMIN = Y(1)
      RMAX = Y(1)
      DO 5 I = 2, NI
      IF (RMIN .GT. Y(I)) RMIN = Y(I)
      IF (RMAX .LT. Y(I)) RMAX = Y(I)
    5 CONTINUE
      PQ = (RMAX - RMIN) / (FNJ - ONE)
      DO 10 J = 1, NJ
      FJ = FLOAT(J) - ONE
      R(J) = RMIN + PQ * FJ
      P(J) = POLD
   10 CONTINUE
C
C        CALCULATE GRID SIZE FOR EVALUATING CONVERGENCE
C        TO GLOBAL MAXIMUM
C
      CALL GRIDN(NI, Y, JN, RMIN, RMAX, RCONS, DEPS)
C
C        START THE ITERATIONS
C
      DFOLD = DINF
C
C        DFOLD REPRESENTS THE VALUE OF THE LOG LIKELIHOOD
C        FUNCTION FROM THE PREVIOUS ITERATION
C
      IT = 0
      IC = ICYC + 1
   30 IT = IT + 1
C
C        'IT' COUNTS THE NUMBER OF EM CYCLES REQUIRED
C        TO REACH CONVERGENCE AT GLOBAL MAXIMUM
C
      DO 75 LPOUT = 1, ITER
      DFL = ZERO
      LP(IT) = LPOUT
      DO 40 I = 1, NI
      FL = ZERO
      DO 35 J = 1, NJ
      PR = (Y(I) - R(J)) * (Y(I) - R(J)) / (TWO * S(I))
      A = EXP(-PR) / SQRT(PI * S(I))
      FL = FL + P(J) * A
   35 CONTINUE
C
C        A REPRESENTS THE SAMPLING DENSITY
C        FL REPRESENTS THE MARGINAL DENSITY
C        DFL REPERESENTS THE LOG LIKELIHOOD FUNCTION
C
      DFL = DFL + LOG(FL)
      AFI(I) = FL
   40 CONTINUE
      DIFF = ABS((DFL - DFOLD) / DFOLD)
C
C        EVALUATE WHETHER THE RELATIVE CHANGE IN LOG
C        LIKELIHOOD IS LESS THAN DCONS
C
      IF (DIFF .LE. DCONS) GOTO 80
      DFOLD = DFL
C
C        CALCULATE THE NEW ESTIMATE
C
      DO 50 J = 1, NJ
      POLD = P(J)
      ROLD = R(J)
      R(J) = ZERO
      P(J) = ZERO
      DENOM = ZERO
      DO 45 I = 1, NI
      PR = (Y(I) - ROLD) * (Y(I) - ROLD) / (TWO * S(I))
      A = EXP(-PR) / SQRT(PI * S(I))
      Q = A / (S(I) * AFI(I))
      R(J) = R(J) + Y(I) * Q
      P(J) = P(J) + A / AFI(I)
      DENOM = DENOM + Q
   45 CONTINUE
      R(J) = R(J) / DENOM
      P(J) = POLD * P(J) / FNI
   50 CONTINUE
C
C        COLLAPSE THE NUMBER OF SUPPORT POINTS
C
      JJ = 1
      DO 70 J = 2, NJ
      DIFF = ABS(R(J) - R(JJ))
      IF (DIFF .LE. DEPS) GOTO 60
      JJ = JJ + 1
      R(JJ) = R(J)
      P(JJ) = P(J)
      GOTO 65
   60 CONTINUE
      P(JJ) = P(JJ) + P(J)
   65 CONTINUE
   70 CONTINUE
      NJ = JJ
   75 CONTINUE
C
C        END OF ITERATIONS
C
   80 CONTINUE
C
C        CHECK IF THE ESTIMATE IS THE GLOBAL MAXIMUM
C
      IF (IT .EQ. IC) GOTO 169
      R(0) = RMIN
      IJ = 0
      JNN = 1
      FJN = FLOAT(JN)
      PQ = (RMAX - RMIN) / (FJN - ONE)
      DO 130 J = 1, JN
      FJ = FLOAT(J) - ONE
      RR = RMIN + PQ * FJ
      IF (RR .LT. R(IJ)) GOTO 130
      AI = ZERO
      DO 90 I = 1, NI
      PR = (Y(I) - RR) ** 2
      PR = PR / (TWO * S(I))
      A = EXP(-PR) / SQRT(PI * S(I))
      AI = AI + (A / AFI(I))
   90 CONTINUE
      DIFF = AI - FNI
C
C        CHECK IF THE CRITERION OF GLOBAL MAXIMUM
C        HOLDS FOR ALL VALUES IN THE RANGE OF THE
C        OBSERVATIONS
C
      IF (DIFF .LE. DINF) GOTO 130
      DO 100 IJ = 1, JJ
      ANB = RR - R(IJ)
      AN = ABS(ANB)
      IF (AN .LE. DEPS) GOTO 130
      IF (ANB .LE. ZERO) GOTO 105
  100 CONTINUE
  105 CONTINUE
C
C        IF THE CRITERION OF GLOBAL MAXIMUM IS NOT MET,
C        NEW SUPPORT IS ADDED AT THE POINT OF FAILURE
C
      JNN = 0
      IF (NJ .GE. NI) GOTO 152
      NJ = NJ + 1
      FNJ = FLOAT(NJ)
      DO 110 I = 1, NJ
      P(I) = ONE / FNJ
  110 CONTINUE
      R(NJ) = (R(IJ) + RR) / TWO
      IF (RR .GE. R(JJ)) R(NJ) = (RMAX + RR) / TWO
      IF (RR .GE. R(JJ)) GOTO 140
  130 CONTINUE
  140 CONTINUE
C
C        CHECK IF THE CRITERION OF GLOBAL MAXIMUM
C        HOLDS AT THE ESTIMATED POINTS OF SUPPORT
C
      DO 150 J = 1, JJ
      AI = ZERO
      DO 145 I = 1, NI
      PR = (Y(I) - R(J)) * (Y(I) - R(J)) / (TWO * S(I))
      A = EXP(-PR) / SQRT(PI * S(I))
      A = A / AFI(I)
      AI = AI + A
  145 CONTINUE
      DIFF = ABS(AI - FNI)
      IF (DIFF .LE. DINF) GOTO 150
      IF (P(J) .LE. DINF) GOTO 150
      JNN = 0
  150 CONTINUE
  152 JJ = NJ
C
C        IF NO NEW SUPPORT IS ADDED, THEN CONVERGENCE TO
C        GLOBAL MAXIMUM IS ASSUMED
C
      IF (JNN .EQ. 1) GOTO 170
C
C        SORT THE POINTS OF SUPPORT IN AN ASCENDING ORDER
C
      NEXLAS = NJ - 1
      MORE = 0
  155 IF (MORE .EQ. 1) GOTO 165
      MORE = 1
      DO 160 I = 1, NEXLAS
      J = I + 1
      IF (R(I) .LE. R(J)) GOTO 160
      TEMP = R(I)
      R(I) = R(J)
      R(J) = TEMP
      MORE = 0
  160 CONTINUE
      GOTO 155
  165 CONTINUE
C
C        IF THE GLOBAL MAXIMUM IS NOT REACHED AND IF THE
C        MAXIMUM ALLOWED NUMBER OF CYCLES IS NOT EXHAUSTED,
C        THEN A NEW CYCLE IS RESTARTED WITH THE NEW SUPPORT
C
      IF (IT .LE. ICYC) GOTO 30
C
C        IF THE GLOBAL MAXIMUM IS NOT REACHED AND IF THE
C        MAXIMUM ALLOWED NUMBER OF CYCLES IS EXHAUSTED,
C        THEN IFAULT IS SET TO 1
C
  169 IFAULT = 1
  170 CONTINUE
      ICYC = IT
      IF (IFAULT .EQ. 1) ICYC = ICYC - 1
      X = ZERO
      W = ZERO
C
C        CALCULATE THE PRIOR MEAN AND VARIANCE
C
      DO 175 J = 1, NJ
      X = X + P(J) * R(J)
  175 CONTINUE
      DO 180 J = 1, NJ
      W = W + P(J) * ((R(J) - X) ** 2)
  180 CONTINUE
C
C        CALCULATE THE POSTERIOR MEANS AND VARIANCES
C
      DO 195 I = 1, NI
      Q = ZERO
      D = ZERO
      DO 185 J = 1, NJ
      PR = (Y(I) - R(J)) * (Y(I) - R(J)) / (TWO * S(I))
      A = EXP(-PR) / SQRT(PI * S(I))
      A = A / AFI(I)
      Q = Q + P(J) * R(J) * A
  185 CONTINUE
      DO 190 J = 1, NJ
      PR = (Y(I) - R(J)) * (Y(I) - R(J)) / (TWO * S(I))
      A = EXP(-PR) / SQRT(PI * S(I))
      A = A / AFI(I)
      D = D + P(J) * A * ((R(J) - Q) ** 2)
  190 CONTINUE
      M(I) = Q
      V(I) = D
  195 CONTINUE
  999 CONTINUE
      RETURN
      END
C
      SUBROUTINE GRIDN(NI, Y, JN, RMIN, RMAX, RCONS, DEPS)
C
C        ALGORITHM AS221  APPL. STATIST. (1986) VOL. 35, NO. 3
C
C        AUXILIARY SUBROUTINE FOR CALCULATING THE GRID SIZE TO
C        BE USED IN EVALUATING CONVERGENCE TO THE GLOBAL MAXIMUM
C
      REAL Y(NI), RMIN, RMAX, RCONS, DEPS, SUM, AVE, VAR, SD
      INTEGER NI, JN
      DATA ZERO /0.0/
      SUM = ZERO
      VAR = ZERO
      DO 5 I = 1, NI
      SUM = SUM + Y(I)
    5 CONTINUE
      FNI = FLOAT(NI)
      AVE = SUM / FNI
      DO 10 I = 1, NI
      VAR = VAR + (Y(I) - AVE) ** 2
   10 CONTINUE
      VAR = VAR / (FNI - 1.0)
      SD = SQRT(VAR)
      DEPS = RCONS * SD
      JN = (RMAX - RMIN) / DEPS
      RETURN
      END
