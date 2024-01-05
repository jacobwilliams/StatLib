      SUBROUTINE OCCALC(SSEQ, L, MAXL, OC)
C
C        ALGORITHM AS 291.1 APPL.STATIST. (1994), VOL.43, NO.2
C
C        Calculates overlap capability of subsequence SSEQ of length L
C        OC(I) = 1 if SSEQ can overlap its own last I letters
C                  (I = 1, ... , L)
C              = 0 otherwise
C
      INTEGER L, MAXL, OC(MAXL), SSEQ(MAXL)
C
      INTEGER I, INDEX, LM1, LMI
C
      OC(L) = 1
      LM1 = L - 1
      DO 20 I = 1, LM1
         LMI = L - I
         OC(I) = 0
         DO 10 INDEX = 1, I
            IF (SSEQ(INDEX) .EQ. SSEQ(INDEX + LMI)) GO TO 10
            GO TO 20
   10    CONTINUE
         OC(I) = 1
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE EXPVAR(M, L, MAXL, OC, P, EXPECT, VARIAN, IFAULT)
C
C        ALGORITHM AS 291.2 APPL.STATIST. (1994), VOL.43, NO.2
C
C        Compute expected value and variance of frequency of
C        occurrence
C
      INTEGER M, L, MAXL, OC(MAXL), IFAULT
      REAL P(MAXL), EXPECT, VARIAN
C
      INTEGER I, INDEX, LM1, NN
      REAL TERM
C
      IFAULT = 1
      IF (M .LT. L .OR. L .LT. 2 .OR. L .GT. MAXL) RETURN
      IFAULT = 0
      NN = M - L + 1
      EXPECT = NN * P(L)
      VARIAN = EXPECT * (1.0 - EXPECT)
      IF ((NN - L) .GT. 0) VARIAN = VARIAN +
     *   ((P(L) ** 2) * (NN - L) * (NN - L + 1))
      TERM = 0.0
      LM1 = L - 1
      DO 10 I = 1, LM1
         IF ((NN - I) .LE. 0) GO TO 10
         INDEX = L - I
         IF (OC(INDEX) .EQ. 0) GO TO 10
         TERM = TERM + OC(INDEX) * (NN - I) * P(I)
   10 CONTINUE
      TERM = 2.0 * P(L) * TERM
      VARIAN = VARIAN + TERM
      RETURN
      END
C
      SUBROUTINE PROW(M, N, MAXNP1, L, MAXL, OC, INDPP, C, POWER, PP,
     *                ROWNEW, IFAULT)
C
C        ALGORITHM AS 291.3 APPL.STATIST. (1994), VOL.43, NO.2
C
C        Calculate prob(J occurrences of a subsequence of length L
C                     within a sequence of length M) for J = 0, ... , N.
C        Prob(J occurrences) is stored in ROWNEW(J + 1)
C
      INTEGER M, N, MAXNP1, L, MAXL, OC(MAXL), INDPP(MAXL), IFAULT
      REAL C(MAXL, 2), POWER(MAXL), PP(MAXL, MAXNP1),
     *     ROWNEW(MAXNP1)
C
      INTEGER I, INDOLD, IP1, JP1, LM1, NP1, SEQLEN
C
C        Calculate IFAULT
C
      IFAULT = 0
      NP1 = N + 1
      IF (M .GE. L .AND. N .GE. 0 .AND. N .LE. M - L + 1
     *  .AND. NP1 .LE. MAXNP1 .AND. L .GE. 2 .AND. L .LE. MAXL) GO TO 10
      IFAULT = 1
      GO TO 80
C
C        Initialize PP using boundary conditions
C        (PP = previous L rows of probabilities needed to calculate
C              ROWNEW recursively)
C
   10 DO 30 IP1 = 1, L
         PP(IP1, 1) = 1.0
         IF (N .LT. 1) GO TO 30
         DO 20 JP1 = 2, NP1
            PP(IP1, JP1) = 0.0
   20    CONTINUE
   30 CONTINUE
C
C        Initialize INDPP
C        (INDPP = address in PP of rows 1 to L of P)
C
      DO 40 I = 1, L
         INDPP(I) = I
   40 CONTINUE
C
C        Loop to generate new rows of PP recursively
C
      DO 70 SEQLEN = L, M
C
C        Generate next row of probabilities
C
         CALL ROWGEN(N, L, MAXL, OC, MAXNP1, PP, INDPP, C,
     *               POWER, ROWNEW)
         IF (SEQLEN .EQ. M) GO TO 80
C
C        Shift rows of PP up one, overwriting first row
C        (effectively, using indirect addresses)
C
         LM1 = L - 1
         INDOLD = INDPP(1)
         DO 50 I = 1, LM1
            INDPP(I) = INDPP(I + 1)
   50    CONTINUE
         INDPP(L) = INDOLD
C
C        Move ROWNEW into bottom row of PP
C        (effectively, using indirect addresses)
C
         DO 60 JP1 = 1, NP1
            PP(INDOLD, JP1) = ROWNEW(JP1)
   60    CONTINUE
   70 CONTINUE
   80 RETURN
      END
C
      SUBROUTINE ROWGEN (N, L, MAXL, OC, MAXNP1, PP,
     *                  INDPP, C, POWER, ROWNEW)
C
C        ALGORITHM AS 291.4 APPL.STATIST. (1994), VOL.43, NO.2
C
C        Used by SUBROUTINE PROW to calculate one row of values:
C           prob(J occurrences of subsequence of length L)
C                for J = 0, ... , N.
C        Prob(J occurrences) is stored in ROWNEW(J + 1)
C
      INTEGER N, L, MAXL, OC(MAXL), MAXNP1, INDPP(MAXL)
      REAL PP(MAXL, MAXNP1), C(MAXL, 2), POWER(MAXL),
     *     ROWNEW(MAXNP1)
C
      INTEGER I, INDEX, J, JNDEX, JP1, K, LIMJ, LM1, LM2, NP1
      REAL Q, Q1
C
      NP1 = N + 1
      LM1 = L - 1
      LM2 = L - 2
      DO 10 I = 1, NP1
         ROWNEW(I) = 0.0
   10 CONTINUE
C
C        Define C and POWER:
C           - C is coefficients in recursive formula;
C           - POWER is powers in recursive formula
C
      DO 20 I = 1, L
         POWER(I) = 4 ** I
   20 CONTINUE
      Q = OC(L - 1)
      C(1, 1) = (4.0 - Q) / POWER(1)
      C(1, 2) = Q / POWER(1)
      Q = OC(1)
      C(L, 1) = (4.0 * Q - 1.0) / POWER(L)
      C(L, 2) = -C(L, 1)
      IF (L .EQ. 2) GO TO 40
      DO 30 K = 1, LM2
         Q = OC(K)
         Q1 = OC(K + 1)
         C(L - K, 1) = ( - Q + 4.0 * Q1) / POWER(L - K)
         C(L - K, 2) = -C(L - K, 1)
   30 CONTINUE
C
C        Compute new probabilities from old ones
C
   40 DO 70 JP1 = 1, NP1
         ROWNEW(JP1) = 0.0
         LIMJ = 2
         IF (JP1 .EQ. 1) LIMJ = 1
         DO 60 J = 1, LIMJ
            DO 50 I = 1, L
               INDEX = INDPP(L - I + 1)
               JNDEX = JP1 - J + 1
               ROWNEW(JP1) = ROWNEW(JP1) + C(I, J) * PP(INDEX, JNDEX)
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
      RETURN
      END
