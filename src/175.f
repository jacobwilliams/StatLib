      SUBROUTINE FACSYM(A, B, W, N1, ITER0, ITER1, ITER2, EPS, IFAULT)
C
C     ALGORITHM AS 175  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Cramer-Wold factorization; 3-stage algorithm.
C     Factorizes: A(x) = a0 + a1(x + 1/x) + a2(x^2 + 1/x^2) + ...
C                      = B(x).B(1/x)
C     where       B(x) = b0 + b1.x + b2.x^2 + ... + bn.x^n
C
C     Auxiliary routine required: SDSDOT from IMSL; an equivalent function
C     has been added to this file.
C
      REAL A(N1), B(N1), W(N1,2), HALF
      DATA HALF /0.5E0/
C
C     Generate initial values for B using Bauer algorithm
C
      IFAULT = 0
      DO 10 I = 1, N1
   10 W(I,1) = A(I)
      CALL BAUER(W, B, N1, ITER0, IER)
      IF (IER .NE. 0) GO TO 130
C
C     Iterate to convergence using Wilson algorithm
C
      IF (ITER1 .LE. 0) GO TO 50
      DO 40 ITER = 1, ITER1
	DO 20 I = 1, N1
	  W(I,1) = A(I)
	  W(I,2) = B(I)
   20   CONTINUE
	CALL RECURS(B, W, N1, IER)
	IF (IER .NE. 0) GO TO 120
	DO 30 I = 1, N1
   30   B(I) = B(I) + HALF * W(I,2)
	IF (B(1) .LE. ABS(B(N1))) GO TO 120
	IF (W(1,2) - B(1) .LE. EOS * W(1,2) .AND. ITER .GT. 1) GO TO 50
   40 CONTINUE
      GO TO 110
C
C     Iterative improvement using incremental form
C
   50 IF (ITER2 .EQ. 0) GO TO 140
      DO 80 ITER = 1, ITER2
	DO 60 I = 1, N1
	  W(I,1) = - A(I)
	  W(I,2) = B(I)
   60   CONTINUE
	CALL RESIDU(W, B, N1)
	CALL RECURS(B, W, N1, IER)
	IF (IER .NE. 0) GO TO 120
	DO 70 I = 1, N1
   70   B(I) = W(I,2) - B(I)
	IF (B(1) .LE. ABS(B(N1))) GO TO 120
	IF (B(1) .GE. W(1,2) .AND. ITER .GT. 1) GO TO 140
   80 CONTINUE
      GO TO 140
C
C     Exits from subroutine
C
  110 IFAULT = 1
      GO TO 140
  120 IFAULT = 2
      GO TO 140
  130 IFAULT = 3
  140 RETURN
      END
C
      SUBROUTINE BAUER(A, B, N1, MAXIT, IFAULT)
C
C     ALGORITHM AS 175.1  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Bauer-Rissanen algorithm.
C
      REAL A(N1), B(N1), C, ONE, ZERO
      DATA ONE /1.0E0/, ZERO /0.0E0/
C
      N = N1 - 1
      IFAULT = 0
      DO 10 I = 1, N
   10 B(I) = A(I+1)
      B(N1) = ZERO
      IF (MAXIT .EQ. 0) GO TO 30
      DO 20 ITER = 1, MAXIT
	IF (A(1) .LE. ABS(A(N1))) GO TO 100
	C = - B(1) / A(1)
	IF (ABS(C) .GT. ONE) GO TO 100
	DO 20 I = 1, N
	  A(I) = A(I) + C * B(I)
	  B(I) = B(I+1) + C * A(I+1)
   20 CONTINUE
   30 IF (A(1) .LE. ZERO) GO TO 100
      C = SQRT(A(1))
      DO 40 I = 1, N1
   40 B(I) = A(I) / C
      RETURN
C
  100 IFAULT = 1
      RETURN
      END
C
      SUBROUTINE RECURS(B, W, N1, IFAULT)
C
C     ALGORITHM AS 175.2  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Recursive solution for polynomial coefficients
C
      REAL B(N1), W(N1), S, HALF, ZERO
      DATA HALF /0.5E0/, ZERO /0.0E0/
C
      N = N1 - 1
      M = N1
      W(1) = HALF * W(1)
      DO 20 I = 1, N
	IF (B(1) .LE. ZERO) GO TO 40
	S = W(M) / B(1)
	K = M
	DO 10 J = 1, M
	  W(J) = W(J) - S * B(K)
	  K = K - 1
   10   CONTINUE
	W(M) = S
	S = B(M) / B(1)
	CALL LINCOM(B, S, M)
	B(M) = S
	M = M - 1
   20 CONTINUE
C
      IF (B(1) .LE. ZERO) GO TO 40
      B(1) = W(1) / B(1)
      DO 30 I = 2, N1
	S = B(I)
	B(I) = ZERO
	CALL LINCOM(B, S, I)
	B(I) = B(I) + W(I)
   30 CONTINUE
      IFAULT = 0
      RETURN
C
   40 IFAULT = 1
      RETURN
      END
C
      SUBROUTINE LINCOM(X, P, M)
C
C     ALGORITHM AS 175.3  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Replaces X with X - P * XSTAR
C
      REAL X(M), P, S
C
      I = 1
      J = M
   10 S = X(I)
      X(I) = X(I) - P * X(J)
      IF (I .NE. J) X(J) = X(J) - P * S
      I = I + 1
      J = J - 1
      IF (I .LE. J) GO TO 10
C
      RETURN
      END
C
      SUBROUTINE RESIDU(A, B, N1)
C
C     ALGORITHM AS 175.4  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Replaces A with A + B * BSTAR
C
      REAL A(N1), B(N1), SDSDOT
C
      DO 10 I = 1, N1
   10 A(I) = SDSDOT(N1-I+1, A(I), B(1), 1, B(I), 1)
C
      RETURN
      END
C
C
      REAL FUNCTION SDSDOT(N, SB, SX, INCX, SY, INCY)
C
C     Evaluate the single-precision dot-product B + X'Y using double
C     precision internally.
C     Modified from the BLAS function SDOT by Alan Miller, CSIRO
C     Division of Mathematics & Statistics, Melbourne, Australia
C
      INTEGER N, INCX, INCY
      REAL SB, SX(*), SY(*)
C
C     Local variables
C
      DOUBLE PRECISION DDOT
      INTEGER I, IX, IY, M, MP1, NS
C***FIRST EXECUTABLE STATEMENT  SDSDOT
      DDOT = SB
      IF(N.LE.0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5, 20, 60
    5 CONTINUE
C
C        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0) IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0) IY = (-N+1)*INCY + 1 
      DO 10 I = 1, N 
        DDOT = DDOT + DBLE(SX(IX)) * DBLE(SY(IY))
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      SDSDOT = DDOT
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. 
C
   20 M = MOD(N, 5) 
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1, M 
        DDOT = DDOT + DBLE(SX(I)) * DBLE(SY(I))
   30 CONTINUE
      IF( N .LT. 5 ) THEN
	SDSDOT = DDOT
        RETURN
      END IF
   40 MP1 = M + 1
      DO 50 I = MP1, N, 5
        DDOT = DDOT + DBLE(SX(I)) * SY(I) + DBLE(SX(I+1)) * SY(I+1) +
     1   DBLE(SX(I+2)) * SY(I+2) + DBLE(SX(I+3)) * SY(I+3) + 
     2   DBLE(SX(I+4)) * SY(I+4)
   50 CONTINUE
      SDSDOT = DDOT
      RETURN
C
C        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
      DO 70 I = 1, NS, INCX
        DDOT = DDOT + DBLE(SX(I)) * DBLE(SY(I))
   70   CONTINUE
      SDSDOT = DDOT
      RETURN
      END 
