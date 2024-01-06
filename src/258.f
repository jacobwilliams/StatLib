      SUBROUTINE ARL2(DELTA, K, H, S0, ARL, ARLFIR, IFAULT)
C
C        ALGORITHM AS 258.1  APPL.STATIST. (1990), VOL.39, NO.3
C
C        Computes the average run length for a cumulative
C        sum control scheme
C
      REAL DELTA, K, H, S0, ARL, ARLFIR
      INTEGER IFAULT
      REAL ARLH, ARLHF, ARLL, ARLLF, BIGARL, BIGDEL
      INTEGER JFAULT
      DATA BIGARL / 1.E30 / , BIGDEL / 5.0 /
C
      IFAULT = 0
      IF (DELTA .LT. 0.0) THEN
         IFAULT = 1
      ELSE
C
C        Compute ARL's for upper tail.
C
         CALL ARL1(DELTA, K, H, S0, ARLH, ARLHF, IFAULT)
         IF (IFAULT .EQ. 0) THEN
C
C        If DELTA=0, then ARL's for lower tail are the same as for
C        the upper.
C
            IF (DELTA .EQ. 0.0) THEN
               ARLLF = ARLHF
               ARLL = ARLH
C
C        If DELTA is too large, skip the low-side ARL calculation.
C
            ELSE IF (DELTA .GT. BIGDEL) THEN
               ARLL = BIGARL
               ARLLF = BIGARL
            ELSE
C
C        Otherwise compute ARL's for lower tail.
C
               CALL ARL1(-DELTA, K, H, S0, ARLL, ARLLF, JFAULT)
C
C        Set lower ARL's large if negative JFAULT .GT. 0
C
               IF (ARLL .LE. ARLH .OR. ARLLF .LE. ARLHF .OR.
     *             ARLL .LT. ARLLF .OR. JFAULT .GT. 0) THEN
                  ARLL = BIGARL
                  ARLLF = BIGARL
               END IF
            END IF
C
C        Compute two-sided ARL for S0=0.0
C
            ARL = ARLH / (1.0 + ARLH / ARLL)
C
C        Compute two-sided ARL for specified value of S0.
C
            ARLFIR = ARLHF / (1.0 + ARLH / ARLL) +
     *               ARLH / (ARLH / ARLLF + ARLL / ARLLF) - ARL
C
C        Set IFAULT=3 if two-sided ARL's are lower bounds.
C
            IF (IFAULT .EQ. 0 .AND. S0 .GT. H / 2.0 + K) IFAULT = 3
         END IF
      END IF
      RETURN
      END
      SUBROUTINE ARL1(DELTA, K, H, S0, ARL, ARLFIR, IFAULT)
C
C        ALGORITHM AS 258.2  APPL.STATIST. (1990), VOL.39, NO.3
C
      REAL DELTA, K, H, S0, ARL, ARLFIR
      INTEGER IFAULT
      INTEGER N, N1, N2, I, J
      REAL XN
      DOUBLE PRECISION XCOND
      PARAMETER (N=12, N1=N + 1, N2=N + 2, XN=N, XCOND=100.D0)
      INTEGER IPVT(N1)
      REAL P1, P2, ALNORM
      DOUBLE PRECISION A(N1, N1), B(N1), R(N1), W(N2),
     *                 C, E1, E2, RCOND, S, T
      EXTERNAL ALNORM
C
C        N is the degree of the polynomial approximation.
C        XCOND defines the criterion for singularity:
C              XCOND+RCOND .LE. XCOND,
C        where RCOND is the reciprocal of the condition number.
C
      IFAULT = 0
      IF (K .LT. 0.0 .OR. H .LT. 0.0 .OR. S0 .LT. 0.0 .OR.
     *    S0 .GT. H) THEN
         IFAULT = 1
      ELSE IF (H .EQ. 0.0) THEN
         ARL = 1.0 / ALNORM(DELTA - K, .FALSE.)
         ARLFIR = ARL
      ELSE
C
C        Set C.
C
         C = MAX(0.0, K - DELTA)
C
C        For each point S at which the polynomial approximation is to be
C        evaluated...
C
         DO 40 I = 0, N
C
C        Compute S
C
            S = H * I / XN
C
C        Calculate necessary exponentials in S.
C
            E1 = EXP(C * S)
            E2 = EXP((S + DELTA - K) * C + C * C / 2.0)
C
C        Apply left-hand-side of integral equation.
C
            T = E1
            DO 10 J = 1, N + 1
               A(I + 1, J) = T
               T = T * S
   10       CONTINUE
C
C        Apply lower integration limit.
C
            CALL MOMENT(-S - DELTA - C + K, -S - DELTA - C + K, N, R, W)
            DO 20 J = 1, N + 1
               A(I + 1, J) = A(I + 1, J) - R(J) * E2
   20       CONTINUE
C
C        Apply upper integration limit.
C
            CALL MOMENT(H - S - DELTA - C + K, -S - DELTA - C + K, N, R,
     *                  W)
            DO 30 J = 1, N + 1
               A(I + 1, J) = A(I + 1, J) + R(J) * E2
   30       CONTINUE
C
C        Apply term '1 + L(0) F(-S-DELTA+K)'.
C
            A(I + 1, 1) = A(I + 1, 1) - ALNORM(REAL(-S - DELTA + K),
     *                    .FALSE.)
            B(I + 1) = 1.0
   40    CONTINUE
C
C        Normalize the simultaneous equations
C
         DO 70 I = 1, N + 1
            S = 0.0
            DO 50 J = 1, N + 1
               S = MAX(S, ABS(A(I, J)))
   50       CONTINUE
            B(I) = B(I) / S
            DO 60 J = 1, N + 1
               A(I, J) = A(I, J) / S
   60       CONTINUE
   70    CONTINUE
         DO 100 J = 1, N + 1
            W(J) = 0.0
            DO 80 I = 1, N + 1
               W(J) = MAX(W(J), ABS(A(I, J)))
   80       CONTINUE
            DO 90 I = 1, N + 1
               A(I, J) = A(I, J) / W(J)
   90       CONTINUE
  100    CONTINUE
C
C        Factor matrix A.  If equations are singular to working
C        precision, IFAULT=2.
C
C        ***************************************
C        SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
C        on entry:
C        A:     the matrix to be factored.
C        LDA:   the leading dimension of array A.
C        N:     the order of the matrix A.
C        on return:
C        A:     the lu factorization of A.
C        IPVT:  pivot indices.
C        RCOND: an estimate of the reciprocal condition of A.
C        Z:     a working vector.
C        ***************************************
C
         CALL DGECO(A, N + 1, N + 1, IPVT, RCOND, R)
         IF (XCOND + RCOND .EQ. XCOND) THEN
            IFAULT = 2
         ELSE
C
C        Solve for the polynomial coefficients
C
C        ***************************************
C        SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
C        on entry:
C        A:     the output from dgeco.
C        LDA:   the leading dimension of array A.
C        N:     the order of the matrix A.
C        IPVT:  the pivot vector from dgeco.
C        B:     the right hand side vector.
C        JOB:   = 0       to solve A*X=B.
C               = nonzero to solve trans(A)*X=B.
C        on return:
C        B:     the solution vector X.
C        ***************************************
C
            CALL DGESL(A, N + 1, N + 1, IPVT, B, 0)
C
C        Get ARL and ARLFIR.
C
            ARL = B(1) / W(1)
            ARLFIR = 0.0
            DO 110 I = 0, N
               ARLFIR = S0 * ARLFIR + B(N - I + 1) / W(N - I + 1)
  110       CONTINUE
            ARLFIR = ARLFIR * EXP(C * S0)
         END IF
      END IF
      RETURN
      END
      SUBROUTINE MOMENT(X, Y, N, R, W)
C
C        ALGORITHM AS 258.3  APPL.STATIST. (1990), VOL.39, NO.3
C
C        For k=0,...,n,  computes the integral from x to
C        infinity of the quantity
C
C                 R(k+1) = ( t - y )**k z(t) dt,
C            where
C                 z(t) = 1/sqrt(2 pi) exp( -t**2 / 2 ) .
C
      INTEGER N
      DOUBLE PRECISION X, Y, R( * ), W( * )
      INTEGER I, K
      REAL ALNORM
      DOUBLE PRECISION XMY, SQR2PI, FACT(19)
      EXTERNAL ALNORM
      DATA SQR2PI / 2.506628274631000502415765D0 /
      DATA FACT / 2 * 1.D0, 2.D0, 6.D0, 24.D0, 120.D0, 720.D0, 5040.D0,
     *     40320.D0, 362880.D0, 3628800.D0, 39916800.D0, 479001600.D0,
     *     6227020800.D0, 87178291200.D0, 1307674368000.D0,
     *     20922789888000.D0, 355687428096000.D0, 6402373705728000.D0 /
C
C        Compute first term of R.
C
      W(1) = EXP(-X * X / 2.0) / SQR2PI
      W(2) = ALNORM(REAL(-X), .FALSE.)
      R(1) = W(2)
      IF (N .GT. 0) THEN
         DO 10 I = 1, N
            W(I + 2) = (W(I) - X * W(I + 1)) / I
            R(I + 1) = W(I + 2) * FACT(I + 1)
   10    CONTINUE
C
C        If X=Y, then R is already computed.
C
         IF (X .NE. Y) THEN
C
C        Compute R.
C
            DO 30 K = 0, N
               R(K + 1) = W(2) / FACT(K + 1)
               XMY = X - Y
               DO 20 I = 1, K
                  R(K + 1) = R(K + 1) * XMY + W(I + 2) / FACT(K - I + 1)
   20          CONTINUE
               R(K + 1) = R(K + 1) * FACT(K + 1)
   30       CONTINUE
         END IF
      END IF
      RETURN
      END
