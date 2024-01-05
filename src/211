c  This file includes ASR 71 and ASR 74 which provide replacements for
c  subroutines FALG and GALG respectively.   The replacements have been
c  placed at the end.
c
      SUBROUTINE FGALG(A, RN, B, IP, K, IFAULT, NF)
C
C        ALGORITHM AS211  APPL. STATIST. (1985) VOL. 34, NO. 2
C        F-G DIAGONALIZATION ALGORITHM
C
      REAL A(5, 10, 10), H(5, 10, 10), RN(5), B(10, 10)
      REAL ZERO, ONE, P8, P6, G1, G2
      DATA ZERO /0.0/, ONE /1.0/, P8 /0.8/, P6 /0.6/
C
C        CHECK INPUT PARAMETERS
C
      IFAULT = 1
      IF (K .LT. 1 .OR. K .GT. 5) RETURN
      IFAULT = 2
      IF (IP .LT. 2 .OR. IP .GT. 10) RETURN
      IFAULT = 3
      DO 3 I = 1, K
      IF (RN(I) .LE. ZERO) RETURN
    3 CONTINUE
      DO 4 I = 1, K
      IFAULT = -I
      DO 4 L = 2, IP
      JEND = L - 1
      DO 4 J = 1, JEND
      IF (A(I, L, J) .NE. A(I, J, L)) RETURN
    4 CONTINUE
      IFAULT = 0
C
      DO 5 L = 1, IP
      DO 5 J = 1, IP
      B(L, J) = ZERO
      IF (L .EQ. J) B(L, J) = ONE
    5 CONTINUE
C
      DO 7 I = 1, K
      DO 7 L = 1, IP
      DO 7 J = 1, IP
    7 H(I, L, J) = A(I, L, J)
      CALL FALG(A, RN, B, IP, K, IFAULT, NF)
      IF (IFAULT .NE. 0 .OR. NF .GT. 1) RETURN
C
C        IF F-ALGORITHM STOPPED AFTER ONLY 1 ITERATION, TRY
C        ANOTHER INITIAL APPROXIMATION TO THE ORTHOGONAL MATRIX B.
C
      DO 8 I = 1, K
      DO 8 L = 1, IP
      DO 8 J = 1, IP
    8 A(I, L, J) = H(I, L, J)
      IP1 = IP - 1
      DO 6 L = 1, IP1
      DO 6 J = 1, IP
      G1 = P8 * B(L, J) + P6 * B(L + 1, J)
      G2 = P8 * B(L + 1, J) - P6 * B(L, J)
      B(L, J) = G1
      B(L + 1, J) = G2
    6 CONTINUE
      CALL FALG(A, RN, B, IP, K, IFAULT, NF)
      RETURN
      END
C
      SUBROUTINE FALG(A, RN, B, IP, K, IFAULT, NF)
      REAL A(5, 10, 10), RN(5), B(10, 10), BOLD(10, 10)
      REAL AUX(10, 10), T(5, 2, 2), Q(2, 2), G(10, 10)
      REAL EPS, EPSF, DIFF, ZERO
      DATA ZERO /0.0/, EPSF /0.0001/, MAXF /15/
C
      NF = 0
      IFAULT = 0
      IP1 = IP - 1
C
C        START ITERATION STEP OF F-ALGORITHM
C
    5 NF = NF + 1
      DO 6 L = 1, IP
      DO 6 J = 1, IP
    6 BOLD(L, J) = B(L, J)
      DO 7 L = 1, IP1
      JSTART = L + 1
      DO 7 J = JSTART, IP
      DO 8 M = 1, IP
      G(M, 1) = B(M, L)
    8 G(M, 2) = B(M, J)
      DO 9 I = 1, K
      DO 10 M1 = 1, IP
      DO 10 M2 = 1, IP
   10 AUX(M1, M2) = A(I, M1, M2)
      CALL MULT(AUX, G, IP, 2)
      T(I, 1, 1) = AUX(1, 1)
      T(I, 2, 1) = AUX(2, 1)
      T(I, 1, 2) = AUX(1, 2)
      T(I, 2, 2) = AUX(2, 2)
    9 CONTINUE
C
      CALL GALG(T, Q, RN, K, IFAULT)
      IF (IFAULT .NE. 0) RETURN
C
      DO 7 M = 1, IP
      B(M, L) = G(M, 1) * Q(1, 1) + G(M, 2) * Q(2, 1)
      B(M, J) = G(M, 1) * Q(1, 2) + G(M, 2) * Q(2, 2)
    7 CONTINUE
      EPS = ZERO
      DO 11 L = 1, IP
      DO 11 J = 1, IP
      DIFF = ABS(B(L, J) - BOLD(L, J))
      IF (DIFF .GT. EPS) EPS = DIFF
   11 CONTINUE
      IF (EPS .GE. EPSF .AND. NF .LT. MAXF) GOTO 5
      IF (EPS .GE. EPSF) IFAULT = 4
      DO 12 I = 1, K
      DO 13 L = 1, IP
      DO 13 J = 1, IP
   13 AUX(L, J) = A(I, L, J)
      CALL MULT(AUX, B, IP, IP)
      DO 12 L = 1, IP
      DO 12 J = 1, IP
      A(I, L, J) = AUX(L, J)
   12 CONTINUE
      RETURN
      END
C
      SUBROUTINE GALG(T, Q, RN, K, IFAULT)
      REAL T(5, 2, 2), Q(2, 2), RN(5), DELTA(5, 2), COEF(5), U(2, 2)
      REAL EPSG, ZERO, ONE, TWO, CO, SI, CO2, SI2, COSI, SINOLD, CP
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, EPSG /0.0001/
      MAXG = 5
C
      IF (K .EQ. 1) MAXG = 1
      CO = ONE
      SI = ZERO
C
C        START ITERATION STEP OF G-ALGORITHM
C
      DO 1 NG = 1, MAXG
      SINOLD = SI
      CO2 = CO * CO
      SI2 = SI * SI
      COSI = TWO * CO * SI
C
      DO 2 I = 1, K
      DELTA(I, 1) = CO2 * T(I, 1, 1) + SI2 * T(I, 2, 2) + COSI *
     *  T(I, 1, 2)
      DELTA(I, 2) = CO2 * T(I, 2, 2) + SI2 * T(I, 1, 1) - COSI *
     *  T(I, 2, 1)
      IF (DELTA(I, 1) .GT. ZERO .AND. DELTA(I, 2) .GT. ZERO) GOTO 4
      IFAULT = -I
      RETURN
    4 COEF(I) = RN(I) * (DELTA(I, 1) - DELTA(I, 2)) / (DELTA(I, 1)
     *  * DELTA(I, 2))
    2 CONTINUE
      DO 3 L = 1, 2
      DO 3 J = 1, 2
      U(L, J) = ZERO
      DO 3 I = 1, K
      U(L, J) = U(L, J) + COEF(I) * T(I, L, J)
    3 CONTINUE
C
      CALL EIGVEC(U, Q)
C
C        REORDER MATRIX Q AND/OR MULTIPLY COLUMN(S) BY -1 IF
C        NECESSARY, SUCH THAT Q IS A ROTATION, AND Q(1,1) IS
C        THE LARGEST ELEMENT
C
      INDEX = 1
      IF(ABS(Q(1, 1)) .LT. ABS(Q(1, 2))) INDEX = 2
      CP = Q(1, INDEX)
      CO = ABS(CP)
      SI = Q(2, INDEX)
      IF (CP .LT. ZERO) SI = -SI
      Q(1, 1) = CO
      Q(2, 1) = SI
      Q(1, 2) = -SI
      Q(2, 2) = CO
      IF (ABS(SI - SINOLD) .LT. EPSG) RETURN
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE EIGVEC(U, Q)
      REAL U(2, 2), Q(2, 2)
      REAL ZERO, ONE, TWO, ROOT, DENOM, EP
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, EP /1 .E - 10/
C
      ROOT = (U(1, 1) + U(2, 2)) / TWO
      ROOT = ROOT + SQRT(((U(1, 1) - U(2, 2)) / TWO) ** 2 + U(1, 2)
     *  ** 2)
      DENOM = SQRT((ROOT - U(1, 1)) ** 2 + U(1, 2) ** 2)
      Q(1, 1) = ONE
      Q(2, 1) = ZERO
      IF (DENOM .LT. EP) GOTO 1
      Q(1, 1) = U(1, 2) / DENOM
      Q(2, 1) = (ROOT - U(1, 1)) / DENOM
    1 Q(1, 2) = -Q(2, 1)
      Q(2, 2) = Q(1, 1)
      RETURN
      END
C
      SUBROUTINE MULT(A, B, IP, IQ)
      REAL A(10, 10), B(10, 10), H(10, 10)
      REAL ZERO
      DATA ZERO /0.0/
C
      DO 1 I = 1, IQ
      DO 1 J = 1, IP
      H(I, J) = ZERO
      DO 1 K = 1, IP
      H(I, J) = H(I, J) + B(K, I) * A(K, J)
    1 CONTINUE
      DO 2 I = 1, IQ
      DO 2 J = 1, I
      A(I, J) = ZERO
      DO 3 K = 1, IP
    3 A(I, J) = A(I, J) + H(I, K) * B(K, J)
      A(J, I) = A(I, J)
    2 CONTINUE
      RETURN
      END
c
c---------------------------------------------------------------------
c
      SUBROUTINE FALG(A, RN, B, IP, K, IFAULT, NF)
C
C        ASR 71 (REMARK ON AS 211) APPL. STATIST. (1988) VOL. 37, NO. 1
C
      INTEGER I, IFAULT, IP, IP1, J, JSTART, K, L, M, NF, MAXF
      REAL A(5, 10, 10), AUX(10, 10), B(10, 10), B1, B2, BOLD(10, 10)
      REAL C, C2, DIFF, EPS, EPSF, Q(2, 2), R1, RN(5), S, T(5, 2, 2)
      REAL T1, T2, ZERO, ONE, TWO
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, EPSF /0.0001/, MAXF /15/
C
      NF = 0
      IFAULT = 0
      IP1 = IP - 1
C
C        INITIAL MULTIPLICATION
C
      DO 2 I = 1, K
      DO 1 L = 1, IP
      DO 1 J = 1, IP
    1 AUX(L, J) = A(I, L, J)
      CALL MULT(AUX, B, IP, IP)
      DO 2 L = 1, IP
      DO 2 J = 1, IP
    2 A(I, L, J) = AUX(L, J)
C
C        START ITERATION STEP OF F-ALGORITHM
C
    4 NF = NF + 1
      DO 5 L = 1, IP
      DO 5 J = 1, IP
    5 BOLD(L, J) = B(L, J)
      DO 8 L = 1, IP1
      JSTART = L + 1
      DO 8 J = JSTART, IP
      DO 6 I = 1, K
      T(I, 1, 1) = A(I, L, L)
      T(I, 1, 2) = A(I, L, J)
      T(I, 2, 1) = A(I, J, L)
      T(I, 2, 2) = A(I, J, J)
    6 CONTINUE
C
C        GET ROTATION SINE AND COSINE
C
      CALL GALG(T, Q, RN, K, IFAULT)
      IF (IFAULT .NE. 0) RETURN
C
C        ROTATE B
C
      C = Q(1, 1)
      S = Q(2, 1)
      R1 = S / (ONE + C)
      DO 7 M = 1, IP
      B1 = B(M, L)
      B2 = B(M, J)
      B(M, L) = B1 + S * (B2 - R1 * B1)
      B(M, J) = B2 - S * (B1 + R1 * B2)
    7 CONTINUE
C
C        UPDATE THE A MATRICES
C
      T1 = S / C
      T2 = T1 * T1
      C2 = C * C
      DO 8 I = 1, K
      A(I, L, L) = C2 * (T(I, 1, 1) + TWO * T1 * T(I, 1, 2) + T2 *
     *  T(I, 2, 2))
      A(I, J, J) = C2 * (T2 * T(I, 1, 1) - TWO * T1 * T(I, 1, 2)
     *  + T(I, 2, 2))
      A(I, L, J) = C2 * (T1 * (T(I, 2, 2) - T(I, 1, 1)) + (ONE - T2)
     *  * T(I, 1, 2))
      A(I, J, L) = A(I, L, J)
      DO 8 M = 1, IP
      IF (M .EQ. J .OR. M .EQ. L) GOTO 8
      B1 = A(I, M, L)
      B2 = A(I, M, J)
      A(I, M, L) = B1 + S * (B2 - R1 * B1)
      A(I, M, J) = B2 - S * (B1 + R1 * B2)
      A(I, L, M) = A(I, M, L)
      A(I, J, M) = A(I, M, J)
    8 CONTINUE
C
C        CHECK FOR CONVERGENCE
C
      EPS = ZERO
      DO 9 L = 1, IP
      DO 9 J = 1, IP
      DIFF = ABS(B(L, J) - BOLD(L, J))
      IF (DIFF .GT. EPS) EPS = DIFF
    9 CONTINUE
      IF (EPS .GE. EPSF .AND. NF .LT. MAXF) GOTO 4
      IF (EPS .GE. EPSF) IFAULT = 4
      RETURN
      END
c
c---------------------------------------------------------------------
c
c  This is ASR 74
c
      SUBROUTINE GALG (T, Q, RN, K, IFAULT)
      REAL C, CO, EIGHT, EPSG, FMG, FOUR, H, ONE, PI, Q(2,2), RN(5), SI, 
     &     T(5,2,2), THETA, TMP, TWO, U, WT, ZERO
C
      DATA EPSG/1.0E-8/, ZERO/0.0/, ONE/1.0/, TWO/2.0/, FOUR/4.0/, 
     &     EIGHT/8.0/, PI/3.141592653589793238/
C
      FMG = ZERO
      H = ZERO
C                                  COMPUTE ITERATION CONSTANTS
      DO 1  I=1, K
      U  = (T(I,2,2)-T(I,1,1))/TWO
      C  = T(I,1,2)
      WT = RN(I)
      FMG = FMG + WT*(U*U-C*C)
    1 H = H + WT*U*C
C                                  CHECK FOR DEGENERATE CASES
      IF (ABS(H) .LT. EPSG .OR. ABS(FMG) .LT. EPSG) GOTO 2
C
C                                  NODEGENERATE CASES (NO. 1 TO 4)
C
      THETA = ATAN(-TWO*H/FMG)/FOUR
C                                  CASES 2 AND 4
      IF (H .GT. ZERO .AND. FMG .LT. ZERO) THETA = THETA - PI/FOUR
      IF (H .LT. ZERO .AND. FMG .LT. ZERO) THETA = THETA + PI/FOUR
      GO TO 3
C
C                                  DEGENERATE CASES (NO. 5 TO 9)
C
    2 THETA = ZERO
      IF (FMG .LE. -EPSG) THETA = -PI/FOUR
      IF (H .GE. EPSG) THETA = -PI/EIGHT
      IF (H .LE. -EPSG) THETA = PI/EIGHT
C                                  GET COSINES AND SINES FROM ANGLE 
    3 SI     = SIN(THETA)
      CO     = COS(THETA)
      Q(1,1) = CO
      Q(1,2) = -SI
      Q(2,1) = SI
      Q(2,2) = CO
      RETURN
      END
