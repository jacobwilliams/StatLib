      SUBROUTINE GG1WT(WAIT, DUMM, N, H, CTRL, IFAULT)
C
C        ALGORITHM AS 265.1  APPL.STATIST. (1991), VOL.40, NO.2
C
C        Obtains the distribution of the stationary waiting time for
C        a G/G/1 queue on using the fast Fourier transform algorithm
C
C     Auxiliary routines required: FORRT and REVRT from AS 97.
C
      INTEGER N, IFAULT
      REAL WAIT(N), DUMM(N), H, CTRL(3)
      INTEGER I, J, M, MM
      REAL AN, CA, FOUR, HALF, ONE, PI, Q, SA, TINY, TSI, TSR, TWI,
     *     TWO, TWR, XI, XR, XOLD, XNEW, YR, YI, Z, ZERO
      REAL ARRIVE, SERVE
      EXTERNAL ARRIVE, SERVE
      DATA ZERO, HALF, ONE, TWO, FOUR / 0.0, 0.5, 1.0, 2.0, 4.0 /
      DATA TINY / 1.0E-6 /
C
      PI = FOUR * ATAN(ONE)
      M = N / 2
      MM = M + 1
C
C        Check admissibility of H and N
C
      IFAULT = 1
      IF (H .LE. ZERO) RETURN
      IFAULT = 2
      J = 4
      DO 10 I = 2, 20
         IF (J .EQ. N) GO TO 20
         J = J * 2
   10 CONTINUE
      RETURN
C
C        Find tail transforms of interarrival
C        and service time distributions
C
   20 IFAULT = 3
      XOLD = ONE
      DO 30 I = 1, M
         DUMM(I) = ZERO
         XNEW = ARRIVE((REAL(I) - HALF) * H)
         IF ((XNEW .LT. ZERO) .OR. (XNEW .GT. XOLD)) RETURN
         XOLD = XNEW
         WAIT(I) = XNEW
   30 CONTINUE
      IFAULT = 4
      XOLD = ZERO
      DO 40 I = MM, N
         WAIT(I) = ZERO
         XNEW = SERVE((REAL(N - I) + HALF) * H)
         IF ((XNEW .GT. ONE) .OR. (XNEW .LT. XOLD)) RETURN
         XOLD = XNEW
         DUMM(I) = -XNEW
   40 CONTINUE
      CTRL(1) = WAIT(M)
      CTRL(2) = -DUMM(M + 1)
      CALL TRANSF(WAIT, N)
      CALL TRANSF(DUMM, N)
C
C        Find transform of normalized harmonic
C        renewal sequence and invert
C
      IFAULT = 5
      Z = WAIT(1) + DUMM(1)
      IF (Z .LE. ZERO) RETURN
      WAIT(1) = -LOG(Z)
      Z = WAIT(MM) + DUMM(MM) - TWO * WAIT(MM) * DUMM(MM)
      IF (Z .LE. ZERO) RETURN
      WAIT(MM) = -LOG(Z)
      DO 50 I = 2, M
         AN = PI * (REAL(I) - ONE) / REAL(M)
         SA = SIN(AN)
         CA = COS(AN)
         TSR = DUMM(I)
         TSI = DUMM(I + M)
         TWR = WAIT(I)
         TWI = WAIT(I + M)
         XR = ONE + (CA - ONE) * TSR - SA * TSI
         XI = (CA - ONE) * TSI + SA * TSR
         YR = TSR + XR * TWR - XI * TWI
         YI = TSI + XR * TWI + XI * TWR
         Z = YR * YR + YI * YI
         IF (Z .LE. TINY) RETURN
         Z = HALF * LOG(Z)
         AN = ATAN2(YI, YR)
         WAIT(I) = -Z
         WAIT(I + M) = -AN
   50 CONTINUE
      CALL REVERT(WAIT, N)
C
C        Change appropriate elements to zero and transform
C
      DO 60 I = 2, M
         WAIT(I) = ZERO
   60 CONTINUE
      CTRL(3) = ABS(WAIT(MM)) + ABS(WAIT(MM + 1)) + ABS(WAIT(MM + 5))
      CALL TRANSF(WAIT, N)
C
C        Find transform of first descending ladder height
C
      WAIT(1) = ONE - EXP(-WAIT(1))
      WAIT(MM) = ONE - EXP(-WAIT(MM))
      DO 70 I = 2, M
         TWR = WAIT(I)
         TWI = WAIT(I + M)
         WAIT(I) = ONE - EXP(-TWR) * COS(TWI)
         WAIT(I + M) = EXP(-TWR) * SIN(TWI)
   70 CONTINUE
C
C        Find transform of associated random walk minimum
C
      Q = ONE - WAIT(1)
      WAIT(1) = ONE
      WAIT(MM) = Q / (ONE - WAIT(MM))
      DO 80 I = 2, M
         TWR = ONE - WAIT(I)
         TWI = WAIT(I + M)
         Z = Q / (TWR * TWR + TWI * TWI)
         WAIT(I) = Z * TWR
         WAIT(I + M) = Z * TWI
   80 CONTINUE
C
C        Invert and shift
C
      CALL REVERT(WAIT, N)
      DO 90 I = 2, M
         WAIT(I) = WAIT(N - I + 2)
   90 CONTINUE
      DO 100 I = MM + 1, N
         WAIT(I) = ZERO
  100 CONTINUE
      IFAULT = 0
      RETURN
      END

      SUBROUTINE TRANSF(X, N)
C
C        ALGORITHM AS 265.2  APPL.STATIST. (1991), VOL.40, NO.2
C
C        Adapts transform to definition used in algorithm AS 97
C
      INTEGER N
      REAL X(N)
      INTEGER I, MM
      REAL Z
C
      Z = REAL(N)
      MM = N / 2 + 1
      DO 10 I = 1, N
         X(I) = Z * X(I)
   10 CONTINUE
      CALL FORRT(X, N)
      DO 20 I = MM + 1, N
         X(I) = -X(I)
   20 CONTINUE
      RETURN
      END

      SUBROUTINE REVERT(X, N)
C
C        ALGORITHM AS 265.3  APPL.STATIST. (1991), VOL.40, NO.2
C
C        Adapts transform to definition used in algorithm AS 97
C
      INTEGER N
      REAL X(N)
      INTEGER I, MM
      REAL ONE, Z
      DATA ONE / 1.0 /
C
      Z = ONE / REAL(N)
      MM = N / 2 + 1
      DO 10 I = MM + 1, N
         X(I) = -X(I)
   10 CONTINUE
      CALL REVRT(X, N)
      DO 20 I = 1, N
         X(I) = Z * X(I)
   20 CONTINUE
      RETURN
      END

