ALGORITHM AS 285
----------------
This file contains in order:
1. Some notes supplied by the author of the algorithm.
2. A MAKEFILE.
3. A driver program.
4. The algorithm.
5. User functions for various regions of integration - a cube, an ellipsoid,
   a rectangle, a sphere, a convex polynomial and a star shape.
6. The input data file FORT.1 needed by the driver program.

The sections are separated by the following:

C--------------------------------------------------------------------------

1. Some notes supplied by the author of the algorithm.

The  subroutines  in  this file compute the probability that a multivariate
normal random vector falls within a prescribed region A.  A is input to the
program via a user-defined function subprogram.

The function should begin as follows:

      real function f(n,v)
      integer n
      real v(n)

      lines calculating f, which should take the value:
            = 0 if v is on the boundary of A
            < 0 if v is inside A
            > 0 if v is outside A

The following functions in this directory  may  be  used  to  describe  the
boundary.   Each  must be altered to provide the specific boundary (or, the
covariance  matrix  input  into  the  driving   program   may   be   scaled
appropriately).

sphere      A is sphere of radius sqrt(radsq)

cube        A is a cube with infinity-norm cubsz

rect        x is in A if a(i) < x(i) < b(i) for all i

ellipsoid   The boundary of A is x'ax + bx - c = 0

convpoly    A is a convex polytope with faces described by
            the inequalities
            a(1,i)*x(1) + a(2,i)*x(2) +...+ a(n,i)*x(n) <= 1
            for i=1,...,p.

pentagram   A is a five-pointed star

Once  an  appropriate boundary function is selected or written, the program
should be  reloaded.   This  may  be  done  quickly  by  editing  the  file
"Makefile" so that the filename listed after "FNCN =" is the object file of
the boundary function, e.g.

      FNCN = sphere.o

To load the file, type

      make FFLAGS=-O

and type

      mulnor

to execute the program.

All input to the driving program should be in file fort.1, of the form

n          order of covariance matrix
stdev      desired standard deviation of estimate
rhigh      radius of sphere circumscribing region A
init       initial number of simulations to be run
ix,iy,iz   three integers between 1 and 50000 to initialize
           random normal generator
covar      lower triangle of the symmetric covariance matrix,
           stored row-wise as a one-dimensional array.

C--------------------------------------------------------------------------

2. A MAKEFILE

FNCN = cube.o 

MULNOR =  mulnor.o chol.o zerovc.o chiprb.o orthp.o 

RANDOM = rnortm.o normal.o alnorm.o random.o ppnd.o 

AUX =  drive.o tdiag.o lrvt.o 

OBJECTS = $(MULNOR) $(RANDOM) $(AUX) $(FNCN)

mulnor: $(OBJECTS)
      f77  $(OBJECTS) -O0  -o mulnor

clean:
      rm *.o

C--------------------------------------------------------------------------

3. A driver program.

C  DRIVING PROGRAM
C
      REAL COVAR(1275),WORK(50,52),T(1275)      
      INTEGER ISEED,N,IFAULT,J,IX,IY,IZ,ITER,INIT
      REAL F,SDEV,PROB,RHIGH,ETA
      COMMON /RAND/ IX,IY,IZ
      EXTERNAL F
      DATA ETA /1.0E-5/
C
C READ IN COVARIANCE AND INITIALIZE RANDOM NORMAL GENERATOR.
C
C With many Fortran compilers, you will need to insert a line such as:
C     OPEN(1, file='FORT.1', status='old')
C before reading from unit 1.
C
      READ (1,*) N,SDEV,RHIGH,INIT,IX,IY,IZ
      READ (1,*) (COVAR(J),J=1,N*(N+1)/2)
C
C  FIND SMALLEST NONZERO EIGENVALUE OF COVARIANCE MATRIX
C
      NN = N*(N+1)/2
      CALL CHOL(COVAR,N,NN,T,NULLTY,IFAULT)
      CALL TDIAG(N,ETA,COVAR,WORK(1,1),WORK(1,2),WORK(1,3),IFAULT)
      CALL LRVT(N,ETA,WORK(1,1),WORK(1,2),WORK(1,3),IFAULT)
      RHIGH = RHIGH/SQRT(WORK(NULLTY+1,1))
C
      CALL MULNOR(COVAR,N,F,SDEV,RHIGH,INIT,ISEED,T,WORK,PROB,
     *    ITER,IFAULT)
      PRINT *, 'ESTIMATED PROBABILITY IS', PROB,'WITH IFAULT = ',IFAULT
      PRINT *, 'ITER = ',ITER
      STOP
      END

      SUBROUTINE TDIAG(N,TOL,A,D,E,Z,IFAULT)
C
      REAL A(*),D(N),E(N),F,G,H,HH,TOL,Z(N,N),ZERO,
     *            ONE,ZSQRT
      INTEGER I,J,K,L,N,IFAULT,I1,J1
C
      DATA ZERO, ONE /0.0, 1.0/
C
      ZSQRT(H) = SQRT(H)
C
      IFAULT = 1
      IF (N .LE. 1) RETURN
      IFAULT = 0
      K = 0
      DO 10 I=1,N
      DO 10 J=1,I
      K = K+1
10    Z(I,J) = A(K)
      I=N
      DO 70 I1=2,N
      L = I-2
      F = Z(I,I-1)
      G=ZERO
      IF (L .LT. 1) GOTO 25
      DO 20 K=1,L
20    G=G+Z(I,K)**2
25    H=G+F*F
C
      IF (G .GT. TOL) GOTO 30
      E(I)=F
      D(I)=ZERO
      GOTO 65
30    L=L+1
      G=ZSQRT(H)
      IF (F .GE. ZERO) G=-G
      E(I) = G
      H = H-F*G
      Z(I,I-1)=F-G
      F=ZERO
      DO 50 J=1,L
      Z(J,I) = Z(I,J)/H
      G = ZERO
C
      DO 40 K=1,J
40    G = G+Z(J,K) *Z(I,K)
      IF (J .GE. L)GOTO 47
      J1  = J+1
      DO 45 K=J1,L
45    G=G + Z(K,J)*Z(I,K)
C
47    E(J)=G/H
      F = F+G*Z(J,I)
50    CONTINUE
C
      HH = F/(H+H)
C
      DO 60 J=1,L
      F=Z(I,J)
      G=E(J)-HH*F
      E(J)=G
      DO 60 K=1,J
      Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
60    CONTINUE
      D(I)=H
65    I=I-1
70    CONTINUE
      D(1) = ZERO
      E(1) = ZERO
C
      DO 110 I=1,N
      L=I-1
      IF (D(I) .EQ. ZERO .OR. L .EQ. 0) GO TO 100
      DO 90 J=1,L
      G=ZERO
      DO 80 K=1,L
80    G=G + Z(I,K)*Z(K,I)
      DO 90 K=1,L
      Z(K,J) = Z(K,J) - G*Z(K,I)
90    CONTINUE
100   D(I) = Z(I,I)
      Z(I,I) = ONE
      IF (L .EQ. 0) GO TO 110
      DO 105 J=1,L
      Z(I,J) = ZERO
      Z(J,I) = ZERO
105   CONTINUE
110   CONTINUE
      RETURN
      END


      SUBROUTINE LRVT(N,PRECIS,D,E,Z,IFAULT)
C
      REAL B,C,D(N),E(N),F,G,H,P,PR,PRECIS,R,S,
     *          Z(N,N),ZERO,ONE,TWO,ZABS,ZSQRT
      INTEGER I,J,K,L,N,IFAULT,MITS,N1,JJ,M1,M,I1
C
      DATA MITS,ZERO,ONE,TWO /30,0.0,1.0,2.0/
      ZABS(B) = ABS(B)
      ZSQRT(B) = SQRT(B)
C
      IFAULT = 2
      IF (N .LE. 1) RETURN
      IFAULT = 1
      N1 = N-1
      DO 10 I=2,N
10    E(I-1) = E(I)
      E(N) = ZERO
      B = ZERO
      F = ZERO
      DO 90 L=1,N
      JJ=0
      H = PRECIS * (ZABS(D(L)) + ZABS(E(L)))
      IF (B .LT. H) B = H
C
      DO 20 M1=L,N
      M=M1
      IF (ZABS(E(M)) .LE. B) GOTO 30
20    CONTINUE
30    IF (M .EQ. L) GOTO 85
40    IF (JJ .EQ. MITS) RETURN
      JJ = JJ+1
C
      P = (D(L+1) - D(L))/(TWO*E(L))
      R = ZSQRT(P*P + ONE)
      PR = P + R
      IF (P .LT. ZERO) PR = P - R
      H = D(L) - E(L)/PR
      DO 50 I=L,N
50    D(I) = D(I)-H
      F=F+H
C
      P=D(M)
      C=ONE
      S=ZERO
      M1=M-1
      I=M
      DO 80 I1=L,M1
      J=I
      I=I-1
      G=C*E(I)
      H=C*P
      IF(ZABS(P) .GE. ZABS(E(I))) GOTO 60
      C = P/E(I)
      R = ZSQRT(C*C+ONE)
      E(J) = S*E(I)*R
      S = ONE/R
      C = C/R
      GOTO 70
60    C=E(I)/P
      R=ZSQRT(C*C+ONE)
      E(J) = S*P*R
      S=C/R
      C=ONE/R
70    P=C*D(I)-S*G
      D(J)=H+S*(C*G+S*D(I))
C
      DO 80 K=1,N
      H=Z(K,J)
      Z(K,J)=S*Z(K,I) + C*H
      Z(K,I) = C*Z(K,I) - S*H
80    CONTINUE
      E(L)=S*P
      D(L)=C*P
      IF (ZABS(E(L)) .GT. B)GOTO 40
85    D(L) = D(L) + F
90    CONTINUE
C
      DO 120 I=1,N1
      K=I
      P=D(I)
      I1=I+1
      DO 100 J=I1,N
      IF (D(J) .GE. P) GOTO 100
      K=J
      P=D(J)
100   CONTINUE
      IF (K .EQ. I) GOTO 120
      D(K) = D(I)
      D(I) = P
      DO 110 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
      Z(J,K)=P
110   CONTINUE
120   CONTINUE
      IFAULT=0
      RETURN
      END


C--------------------------------------------------------------------------

4. The algorithm.

      SUBROUTINE MULNOR(COVAR,N,F,SDEV,RHIGH,INIT,ISEED,T,WORK,PROB,
     *   ITER,IFAULT)
C
C     ALGORITHM AS285 APPL. STATIST. (1993) VOL.42, NO.3
C
C  FINDS THE PROBABILITY THAT A NORMALLY DISTRIBUTED RANDOM
C  N-VECTOR WITH MEAN 0 AND COVARIANCE COVAR FALLS
C  IN AREA ENCLOSED BY THE EXTERNAL USER-DEFINED FUNCTION F.
C
      INTEGER N,ISEED,IFAULT,INIT,INIT1,
     *          K,II,KK,Q,ITER,NN,NULLTY,I,J,MAXITR
      REAL F,SDEV,RHIGH,PROB,P,VAR,ZERO,COX,
     *       COVAR(*), WORK(N,*),T(*)
C
      EXTERNAL F
      DATA ZERO,MAXITR /0.0,10000/
C
      IFAULT = 0
      INIT1 = INIT
      IF (INIT .LT. 25) INIT1 = 25
      IF (INIT .GT. 1000) INIT1 = 1000
      COX = 1.0 + 2.0/INIT1
C
C  FIND CHOLESKY DECOMPOSITION OF COVARIANCE MATRIX
C
      NN = N*(N+1)/2
      CALL CHOL(COVAR,N,NN,T,NULLTY,IFAULT)
      Q = N - NULLTY
      IFAULT = IFAULT*IFAULT
      IF (IFAULT .NE. 0) RETURN
C
C  TRANSPOSE CHOLESKY FACTOR, OMITTING COLUMNS OF ZEROES.
C  RESULTING VECTOR IS PACKED FORM OF T'J,
C  WRITTEN ROWWISE BEGINNING WITH LAST ROW.
C
      DO 10 I = 1,N
      DO 10 J = 1,I
         WORK(J,I) = ZERO
10    CONTINUE
      II = 0
      DO 20 I = 1,N
      DO 20 J = 1,I
         II = II + 1
         WORK(I,J) = T(II)
20    CONTINUE
      J = 0
      DO 40 I = N,1,-1
         IF (J .EQ. NULLTY) GO TO 50
         IF (WORK(I,I) .EQ. ZERO) THEN
            J = J+1
            DO 30 K = I+1,N
            DO 30 KK = I, N-J
              WORK(K,KK) = WORK(K,KK+1)
30          CONTINUE
         END IF
40    CONTINUE
50    CONTINUE
      K = 0
      DO 60 J = N,1,-1
      DO 60 I = 1,MIN0(J,Q)
         K = K + 1
         T(K) = WORK(J,I)
60    CONTINUE
C
C  TAKE PILOT SAMPLE; ESTIMATE VARIANCE OF PROB FROM ORTHP.
C
      PROB = ZERO
      VAR = ZERO
      DO 70 K=1,INIT1
         CALL ORTHP(F,RHIGH,N,Q,T,WORK,P,ISEED,IFAULT)
         IF (IFAULT .NE. 0) RETURN
         PROB = PROB + P 
         VAR = VAR + P*P
70    CONTINUE
      PROB = PROB/INIT1
      VAR = VAR/INIT1 - PROB*PROB
C
      ITER = INT(COX*VAR/(SDEV*SDEV)) + 1
      IF (ITER .GT. MAXITR) IFAULT = 5
C
      DO 80 K=INIT1+1,ITER
         CALL ORTHP(F,RHIGH,N,Q,T,WORK,P,ISEED,IFAULT)
         IF (IFAULT .NE. 0) RETURN
         PROB = PROB + (P - PROB)/K
80    CONTINUE
C
      RETURN
      END

      SUBROUTINE CHOL(A,N,NN,U,NULLTY,IFAULT)
      REAL A(NN),U(NN),ETA,ETA2,X,W,ZERO,ZABS,ZSQRT
      INTEGER N,NN,NULLTY,IFAULT
      INTEGER J,K,II,ICOL,L,KK,M,IROW,I
C
      DATA ETA,ZERO /1.0E-05,0.0/
      ZABS(X) = ABS(X)
      ZSQRT(X) = SQRT(X)
C
      IFAULT = 1
      IF (N .LE. 0) RETURN
      IFAULT = 3
      IF (NN .NE. N*(N+1)/2) RETURN
      IFAULT = 2
      NULLTY = 0
      J = 1
      K = 0
      ETA2 = ETA*ETA
      II = 0
      DO 80 ICOL = 1,N
      II = II + ICOL
      X = ETA2 * A(II)
      L = 0
      KK = 0
      DO 40 IROW = 1, ICOL
      KK = KK + IROW
      K = K+1
      W = A(K)
      M = J
      DO 10 I = 1, IROW
      L = L + 1
      IF (I .EQ. IROW) GO TO 20
      W = W - U(L) * U(M)
      M = M + 1
10    CONTINUE
20    IF (IROW .EQ. ICOL) GO TO 50
      IF (U(L) .EQ. ZERO) GO TO 30
      U(K) = W / U(L)
      GO TO 40
30    IF (W * W .GT. ZABS(X * A(KK))) RETURN
      U(K) = ZERO
40    CONTINUE
50    IF (ZABS(W) .LE. ZABS(ETA * A(K))) GO TO 60
      IF (W .LT. ZERO) RETURN
      U(K) = ZSQRT(W)
      GO TO 70
60    U(K) = ZERO
      NULLTY = NULLTY + 1
70    J = J + ICOL
80    CONTINUE
      IFAULT = 0
      RETURN
      END

      SUBROUTINE ORTHP(F,RHIGH,N,Q,T,Z,PROB,ISEED,IFAULT)
C
C  FINDS THE PSEUDO-MONTE CARLO ESTIMATE OF PROB USING AN
C  ORTHONORMAL SYSTEM OF VECTORS.
C
      INTEGER N,ISEED,IFAULT,Q,I,J,K,JJ,II,QP1,QP2
      REAL T(*),RHIGH,Z(N,*),PROB,S,A,B,R,RT2RCP,ZERO,
     *           F,CHIPRB,FAB
C
      EXTERNAL F
      DATA ZERO,RT2RCP /0.0,0.707106781/
C
      PROB = ZERO
      QP1 = Q+1
      QP2 = Q+2
C
C  PUT RANDOMLY ORIENTED ORTHONORMAL SYSTEM OF Q-VECTORS IN Z.
C
      CALL RNORTM(N,Q, QP1, 1, 1, Z(1,QP1), 1, Z(1,QP1), ISEED, Z,
     *   Z(1,QP1), FAB, Z(1,QP2), Z(1,QP2), IFAULT)
C
C  REPLACE EACH COLUMN OF Z BY TRANS(CHOLESKY FACTOR)*J*COLUMN OF Z.
C
      DO 30 I=1,Q
        JJ = 0
        DO 20 J = N,1,-1
           S = ZERO
           DO 10 K=1,MIN0(J,Q)
              JJ = JJ + 1
              S = S + T(JJ)*Z(K,I)
10         CONTINUE
           Z(J,I) = S
20      CONTINUE
30    CONTINUE
C
C  ESTIMATE THE PROBABILITY THAT N(0,COVAR) IS WITHIN THE
C  BOUNDARY GIVEN BY F.
C
      DO 40 I = 1,Q
         A = ZERO 
         B = RHIGH
         CALL ZEROVC (F,R,A,B,N,Z(1,I),Z(1,QP2),IFAULT)
         PROB = PROB - CHIPRB(R*R,Q)
40    CONTINUE
C
      DO 60 I = 1,Q
         A = ZERO 
         B = RHIGH
         DO 50 K = 1,N
           Z(K,QP1) = -Z(K,I)
50       CONTINUE
         CALL ZEROVC (F,R,A,B,N,Z(1,QP1),Z(1,QP2),IFAULT)
         PROB = PROB - CHIPRB(R*R,Q)
60    CONTINUE
C
      DO 80 I = 1,Q-1
      DO 80 J = I+1,Q
      DO 80 JJ= -1,1,2
      DO 80 II = -1,1,2
          DO 70 K=1,N
            Z(K,QP1) = RT2RCP*(Z(K,I)*II + Z(K,J)*JJ)
70        CONTINUE
          A = ZERO
          B = RHIGH
          CALL ZEROVC (F,R,A,B,N,Z(1,QP1),Z(1,QP2),IFAULT)
          PROB = PROB - CHIPRB(R*R,Q)
80    CONTINUE
      PROB = PROB/(2.0*Q*Q) + 1.0
C
      RETURN
      END

      SUBROUTINE RNORTM(LDA,N, NP1, IBCONF, NB, B, NDBI, DBI, ISEED, A,
     *   CHISQ, FAB, U, W, IFAULT)
C
C      ALGORITHM AS 127  APPL. STATIST.  (1978) VOL.27, NO.2
C
C      RNORTM GENERATES ORTHOGONAL MATRICES A FROM A DISTRIBUTION
C      WITH DENSITY FUNCTION DEFINED ON THE GROUP OF ORTHOGONAL
C      MATRICES.  WHEN THE INPUT PARAMETER IBCONF = 1 THE 
C      MATRICES ARE GENERATED FROM THE INVARIANT HAAR MEASURE.
C
      REAL B(NB), DBI(NDBI), A(LDA,N), CHISQ(N), U(NP1), W(N), FAB, UL,
     *   WL, WIL, WWIL, T2, HV, HW, ZERO, ONE
      INTEGER N, NP1, IBCONF, NB, NDBI, ISEED, IFAULT, I, J,IBSUB,L,NP1L
      INTEGER LP1,LDA
C
      DATA ZERO /0.0/, ONE /1.0/
C
C      STATEMENT FUNCTIONS FOR DOUBLE PRECISION
C      DOUBLE PRECISION SQRT, SIGN, ABS
C      SQRT(UL) = DSQRT(UL)
C      SIGN(UL, WL) = DSIGN(UL, WL)
C      ABS(UL) = DABS(UL)
C
C      SECTION 1 INITIALIZATION
C
      IFAULT = 0
      GO TO 4
C
C      ERROR RETURNS
C
2     IFAULT = 2
      RETURN
C
C      INITIALIZE DENSITY, H(STORED IN A), AND IBSUB
C
4     FAB = ONE
      DO 100 I=1,N
      DO 50 J=1,N
50    A(I,J) = ZERO
      A(I,I) = ONE
100   CONTINUE
      IBSUB = 0
C
C      CONSTRUCT  A  ONE COLUMN AT A TIME
C
      DO 315 L=1,N
      NP1L = NP1 - L
C
C      SECTION 2
C      CONSTRUCT U(L...N) UNIFORMLY DISTRIBUTED ON THE SURFACE OF THE
C      (N+1-L)-SPHERE OF RADIUS UL = SQRT(CHISQ(L)).
C
      CALL NORMAL(U(L),NP1L+1,ISEED, CHISQ(L))
      UL = SQRT(CHISQ(L))
      IF (UL .EQ. ZERO) GO TO 2
      IF (L .EQ. N) GO TO 301
C
C      SECTION 3
C      CALCULATE  W = B*U  AND STORE ITS LENGTH IN WL.
C      SIGN OF WL PREVENTS LOSS OF SIGNIFICANCE IN WWIL.
C
C      B IS IDENTITY
C
      WL = SIGN(UL,-U(L))
      DO 219 J=L,N
219   W(J) = U(J)
C
C      SECTION 4
C      PROJECT  W  ONTO ORTHOGOANL COMPLEMENT OF FIRST L-1
C      COLUMNS OF A, NORMALIZE, AND STORE IN L COLUMN OF A.
C      A(L) = H(L) * W * WIL
C      CALCULATE PROJECTION MATRIX H(L) FOR ORTHOGOANL
C      COMPLEMENT OF FIRST L COLUMNS OF A AND STORE IN LAST
C      N-L COLUMNS OF A.
C
      LP1 = L+1
      WIL = ONE/WL
      WWIL = ONE/(WL - W(L))
      DO 250 I=1,N
      T2 = ZERO
      DO 230 J=L,N
230   T2 = T2 + A(I,J) * W(J)
      HV = T2 * WIL
      HW = (HV - A(I,L)) * WWIL
      A(I,L) = HV
      DO 240 J = LP1, N
240   A(I,J) = A(I,J) - HW*W(J)
250   CONTINUE
C
C      END OF L COLUMN OF MATRIX A
C
C      SECTION 5
C      WHEN L = N ONLY SIGN NEEDS TO BE CHOSEN
C      REMEMBER MATRIX H(N) IS STORED IN N COLUMN OF A.
C
301   IF (U(L) .GT. ZERO) GOTO 315
      DO 310 I=1,N
310   A(I,L) = -A(I,L)
C
C      FAB = ABS(B(LAST) * DBI(N)*FAB)
C
315   CONTINUE
320   FAB = ABS(FAB)
C
      RETURN
      END

      SUBROUTINE ZEROVC (F,R,R0,R1,N,V,Y,IFAULT)
C
C  CALCULATES THE SOLUTION R TO F(R*V) = 0.
C
      INTEGER N,IFAULT,I,J,NITER
      REAL F,R,R0,R1,ETA,V(N),Y(N),F0,F1,FR,HALF,ZABS,X
      EXTERNAL F
      DATA HALF, ETA /0.5,1.0E-05/ NITER /500/
      ZABS(X) = ABS(X)
C
      DO 10 I=1,N
         Y(I) = R0*V(I)
10    CONTINUE
      F0 = F(N,Y)
      DO 20 I=1,N
         Y(I) = R1*V(I)
20    CONTINUE
      F1 = F(N,Y)
      IF (SIGN(HALF,F1) .EQ. SIGN(HALF,F0)) THEN
       IFAULT = 3
       RETURN
      END IF
       DO 40 J = 1,NITER
         R = R1 - F1*(R1-R0)/(F1-F0)
         DO 30 I=1,N
               Y(I) = R*V(I)
30       CONTINUE
         FR = F(N,Y)
         IF (ZABS(R-R1) .LT. ETA) GO TO 50
         IF (SIGN(HALF,FR) .NE. SIGN(HALF,F1)) THEN
            R0=R1
            F0=F1
         ELSE
            F0 = F0 * HALF
         END IF
         R1 = R
         F1 = FR
40     CONTINUE
       IFAULT = 6
50    CONTINUE
      RETURN
      END

      REAL FUNCTION CHIPRB(X,DF)
C
C  FINDS THE PROBABILITY THAT A CHI-SQUARED RANDOM VARIABLE
C  WITH DF DEGREES OF FREEDOM EXCEEDS X.
C  ADAPTED FROM HILL AND PIKE (1967) ALG. 299, CHI-SQUARED INTEGRAL,
C  COMM. ACM, VOL.10, 243-244.
C
      REAL X,ZERFP1,ZSQRT,SQRT2,
     *     A,SQRTA,Y,C,E,Z,S,BIG,SMALL,RTPIRP,ZERO,ONE,HALF,TWO
      INTEGER DF,ODD,N2,I
      REAL ALNORM
C     REAL ERF
C
      DATA RTPIRP /0.564189583547756/, SQRT2 /1.41421356237309/
      DATA BIG/88.03/,SMALL/-85.2/,ZERO/0.0/,ONE /1.0/,HALF/0.5/,
     *     TWO /2.0/
      ZSQRT(A) = SQRT(A)
      ZERFP1(A) = TWO*ALNORM(SQRT2*A,.FALSE.)
C
C  NOTE: IF THE FUNCTION "ERF" IS AVAILABLE, THE ABOVE STATEMENT
C  FUNCTION MAY BE REPLACED BY
C     ZERFP1(A) = ERF(A) + ONE
C
      S = ZERO
      IF (X .GT. BIG) GO TO 50
      S = ONE
      IF (X .LT. SMALL) GO TO 50
C
      A = HALF*X
      Y = EXP(-A)
      ODD = MOD(DF,2)
      IF (ODD .EQ. 0) THEN
         S = Y
         E = ONE
         Z = ZERO
      ELSE
         SQRTA = ZSQRT(A)
         S = ZERFP1(-SQRTA)
         E = RTPIRP/SQRTA
         Z = -HALF
      END IF
      IF (DF .LT. 3) GO TO 50
C
      N2 = DF/2 - 1 + ODD
      C = ZERO
      DO 40 I = 1,N2
      Z = Z + 1
      E = E*(A/Z)
      C = C + E
40    CONTINUE
      S = S + C*Y
C
50    CONTINUE
      CHIPRB = S
      END

      SUBROUTINE NORMAL(U,NP1,ISEED, CHISQ)
C
C  GENERATES AN (NP1-1)-VECTOR OF INDEPENDENT NORMAL VARIATES IN U.
C  ON OUTPUT, CHISQ CONTAINS THE SQUARED LENGTH OF THE VECTOR.
C
      REAL U(NP1)
      INTEGER ISEED,NP1,I,IFAULT,IX,IY,IZ
      REAL RANDOM,PPND,CHISQ,ZERO
C
      COMMON /RAND/ IX,IY,IZ
      DATA ZERO /0.0/
C
      CHISQ = ZERO
      DO 10 I = 1, NP1-1
      U(I) = PPND(RANDOM(ISEED),IFAULT)
      CHISQ = CHISQ + U(I)*U(I)
10    CONTINUE
      RETURN
      END

      REAL FUNCTION PPND(P, IFAULT)
C
C  ALGORITHM AS 111  APPL. STATIST. (1977), VOL.26, NO.1
C
C  PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA OF P
C  REAL VERSION FOR EPS = 2 **(-31)
C  THE HASH SUMS ARE THE SUMS OF THE MODULI OF THE COEFFICIENTS
C  THEY HAVE NO INHERENT MEANINGS BUT ARE INCLUDED FOR USE IN
C  CHECKING TRANSCRIPTIONS
C  STANDARD FUNCTIONS ABS, ALOG AND SQRT ARE USED
C
      REAL ZERO, SPLIT, HALF, ONE
      REAL A0, A1, A2, A3, B1, B2, B3, B4, C0, C1, C2, C3, D1, D2
      REAL P, Q, R
      INTEGER IFAULT
      DATA ZERO /0.0E0/, HALF/0.5E0/, ONE/1.0E0/
      DATA SPLIT /0.42E0/
      DATA A0 / 2.50662823884E0/
      DATA A1 / -18.61500062529E0/
      DATA A2 / 41.39119773534E0/
      DATA A3 / -25.44106049637E0/
      DATA B1 / -8.47351093090E0/
      DATA B2 / 23.08336743743E0/
      DATA B3 / -21.06224101826E0/
      DATA B4 / 3.13082909833E0/
      DATA C0 / -2.78718931138E0/
      DATA C1 / -2.29796479134E0/
      DATA C2 / 4.85014127135E0/
      DATA C3 / 2.32121276858E0/
      DATA D1 / 3.54388924762E0/
      DATA D2 / 1.63706781897E0/
C
      IFAULT = 0
      Q = P - HALF
      IF (ABS(Q) .GT. SPLIT) GOTO 1
      R = Q*Q
      PPND = Q * (((A3*R + A2)*R + A1) * R + A0) /
     *  ((((B4*R + B3)*R + B2) * R + B1) * R + ONE)
      RETURN
1     R = P
      IF (Q .GT. ZERO)R = ONE - P
      IF (R .LE. ZERO) GOTO 2
      R = SQRT(-ALOG(R))
      PPND = (((C3 * R + C2) * R + C1) * R + C0)/
     *  ((D2*R + D1) * R + ONE)
      IF (Q .LT. ZERO) PPND = -PPND
      RETURN
2     IFAULT = 1
      PPND = ZERO
      RETURN
      END

      REAL FUNCTION RANDOM(L)
C
C  ALGORITHM AS 183  APPL. STATIST. (1982) VOL. 31, NO. 2
C
C  RETURNS A PSEUDO-RANDOM NUMBER RECTANGULARLY DISTRIBUTED
C  BETWEEN 0 AND 1.
C
C  IX, IY, AND IZ SHOULD BE SET TO INTEGER VALUES BETWEEN
C  1 AND 30000 BEFORE FIRST ENTRY
C
C  INTEGER ARITHMETIC UP TO 30323 IS REQUIRED
C
      INTEGER IX, IY, IZ, L
      COMMON /RAND/ IX,IY,IZ
      IX = 171 * MOD(IX, 177) - 2 * (IX / 177)
      IY = 172 * MOD(IY, 176) - 35 * (IY / 176)
      IZ = 170 * MOD(IZ,178) - 63 * (IZ / 178)
C
      IF (IX .LT. 0) IX = IX + 30269
      IF (IY .LT. 0) IY = IY + 30307
      IF (IZ .LT. 0) IZ = IZ + 30323
C
      RANDOM = MOD(FLOAT(IX) / 30269.0 + FLOAT(IY) / 30307.0 +
     *              FLOAT(IZ) / 30323.0, 1.0)
      RETURN
      END

      REAL FUNCTION ALNORM(X,UPPER)
C
C  ALGORITHM AS 66  APPL. STATIST. (1973) VOL. 22, NO. 3
C
C  EVALUATES THE TAIL AREA OF THE STANDARDIZED NORMAL CURVE
C  FROM X TO INFINITY IF UPPER IS .TRUE. OR
C  FROM MINUS INFINITY TO X IF UPPER IS .FALSE.
C
      REAL LTONE, UTZERO, ZERO, HALF, ONE, CON, Z, Y, X
      LOGICAL UPPER, UP
C
C  LTONE AND UTZERO MUST BE SET TO SUIT THE PARTICULAR COMPUTER
C  (SEE INTRODUCTORY TEXT)
C
      DATA LTONE, UTZERO /7.0, 18.66/
      DATA ZERO, HALF, ONE, CON /0.0,0.5,1.0,1.28/
      UP = UPPER
      Z = X
      IF (Z .GE. ZERO) GOTO 10
      UP = .NOT. UP
      Z = -Z
10    IF (Z .LE. LTONE .OR. UP .AND. Z .LE. UTZERO) GOTO 20
      ALNORM = ZERO
      GOTO 40
20    Y = HALF *Z * Z
      IF (Z .GT. CON) GOTO 30
C
      ALNORM = HALF - Z*(0.398942280444 - 0.399903438504 * Y /
     1  (Y + 5.75885480458 - 29.8213557808 /
     2  (Y + 2.62433121679 + 48.6959930692 /
     3  (Y + 5.92885724438))))
      GO TO 40
C
30    ALNORM = 0.398942280385 * EXP(-Y) /
     1  (Z - 3.8052E-08 + 1.00000615302 /
     2  (Z + 3.98064794E-4 + 1.98615381364 /
     3  (Z - 0.151679116635 + 5.29330324926 /
     4  (Z + 4.8385912808 - 15.1508972451 /
     5  (Z + 0.742380924027 + 30.789933034 / (Z + 3.99019417011))))))
C
40    IF (.NOT. UP) ALNORM = ONE - ALNORM
      RETURN
      END

C--------------------------------------------------------------------------

5. Various user functions for different shapes of region of integration.

      REAL FUNCTION F(N,V)
      REAL V(N),CUBSZ
      INTEGER N,I
      DATA CUBSZ /3.00/
C  
C  BOUNDARY OF REGION IS AN N-DIMENSIONAL CUBE.
C  CUBSZ IS THE DISTANCE FROM THE ORIGIN TO THE SIDE OF THE CUBE
C
      F = 0.0
      DO 20 I = 1,N
      F = AMAX1(F,ABS(V(I)))
20    CONTINUE
      F = F - CUBSZ
      RETURN
      END


      REAL FUNCTION F(N,V)
      REAL V(N)
      INTEGER N
      REAL A(200),B(10),C
      INTEGER I,J,K
C
C  THE BOUNDARY IS AN ELLIPSOID OF THE FORM X'AX + BX + C = 0
C  THE UPPER TRIANGLE OF A IS STORED COLUMN-WISE IN PACKED FORM.
C  IN THIS EXAMPLE, A = INVERSE (1  .3), B = 0, C = -1.
C                               (.3  1)
C
      DATA A(1),A(2),A(3) /1.0989011,-.32967033,1.0989011/
     *    B(1),B(2) /0.0,0.0/ C /-1.0/
C
      F = C
      K = 0
      DO 20 I = 1,N
      K = K + 1
      F = F + (V(I)**2)*A(K)
      DO 10 J = I+1,N
      K = K + 1
      F = F + 2.0*V(I)*V(J)*A(K)
10    CONTINUE
20    CONTINUE
      DO 30 I = 1,N
      F = F + B(I)*V(I)
30    CONTINUE
      RETURN
      END


      REAL FUNCTION F(N,V)
      REAL V(N)
      INTEGER N,I
      REAL TEMP,A(30),B(30)
C
C  RECTANGLE IS GIVEN BY A < V < B FOR EACH COORDINATE.
C  THIS SUBROUTINE MAY ALSO BE USED TO CALCULATE THE
C  CDF BY MAKING THE APPROPRIATE SIDES OF THE RECTANGLE
C  SUFFICIENTLY LARGE.
C  THIS EXAMPLE IS THE FOUR-DIMENSIONAL RECTANGLE DEFINED
C  BY  -INFINITY (ALMOST) < V(I) < 2 FOR EACH I.
C
      DATA A(1),A(2),A(3),A(4)
     *  /-20.0,-20.0,-20.0,-20.0/,
     *       B(1),B(2),B(3),B(4)
     *  /2.0,2.0,2.0,2.0/
C
      F = -1.0
      DO 10 I=1,N
         IF (V(I) .GT. 0.0) THEN
            TEMP = V(I) - B(I)
         ELSE
            TEMP = A(I) - V(I)
         END IF
         F  = AMAX1(F,TEMP)
10    CONTINUE
      RETURN
      END


      REAL FUNCTION F(N,V)
      REAL V(N),RADSQ
      INTEGER N,I
      DATA RADSQ /4.0/
C  
C  BOUNDARY IS SPHERE OF RADIUS SQRT(RADSQ).
C
      F = 0.0
      DO 20 I = 1,N
      F = F + V(I)*V(I)
20    CONTINUE
      F = F - RADSQ
      RETURN
      END


      REAL FUNCTION F(N,V)
      REAL V(N)
      INTEGER N
      INTEGER I,J,P
      REAL TEMP, A(2,6)
C
C  A CONVEX POLYTOPE WITH P FACES IS DEFINED BY THE INEQUALITIES
C  A(1,I)*V(1) + A(2,I)*V(2) +...+ A(N,I)*V(N) <= 1
C  FOR I=1,...,P.
C
C  THE FIGURE IN THIS DEMONSTRATION PROGRAM IS THE 
C  REGULAR PENTAGON WITH DISTANCE FROM CENTER TO POINT = 1.
C
      DATA P/5/ 
     *       A(1,1),A(2,1) /0.726542528, -1.0/
     *       A(1,2),A(2,2) /1.175570504, 0.381966011/
     *       A(1,3),A(2,3) /0.0, 1.236067978/
     *       A(1,4),A(2,4) /-0.726542528, -1.0/
     *       A(1,5),A(2,5) /-1.175570504, 0.381966011/
C
      F = -1.0
      DO 10 J=1,N
         F = F + A(J,1)*V(J)
10    CONTINUE
      DO 30 I=2,P
        TEMP = -1.0
        DO 20 J=1,N
           TEMP = TEMP + A(J,I)*V(J)
20      CONTINUE
        F = AMAX1(F,TEMP)
30    CONTINUE
      RETURN
      END


      REAL FUNCTION F(N,V)
      REAL V(N)
      INTEGER I,J,P,INDEX,N
      REAL A(2,6),TANG(5),X(2),TANGNT
C
C  THE FIGURE IN THIS DEMONSTRATION PROGRAM IS THE 
C  REGULAR PENTAGRAM WITH DISTANCE FROM CENTER TO POINT = 2.618.
C
      DATA P/5/ 
     *       A(1,3),A(2,3) /0.726542528, -1.0/
     *       A(1,2),A(2,2) /1.175570504, 0.381966011/
     *       A(1,4),A(2,4) /0.0, 1.236067978/
     *       A(1,1),A(2,1) /-0.726542528, -1.0/
     *       A(1,5),A(2,5) /1.175570504, 0.381966011/
     *           TANG(1) / -1.376381921/
     *           TANG(2) / -0.3249196962329/
     *           TANG(3) / 0.3249196962329/
     *           TANG(4) / 1.376381921/
C
      X(1) = ABS(V(1))
      X(2) = V(2)
      INDEX = P
      IF (X(1) .EQ. 0.0) THEN
         INDEX = 3
      ELSE
         TANGNT = X(2)/X(1)
         DO 10 I=1,P
            IF (TANGNT .LE. TANG(I)) THEN
            INDEX = I
            GO TO 20
            END IF
10       CONTINUE
20       CONTINUE
      END IF
      F = -1.0
      DO 30 J=1,2
         F = F + A(J,INDEX)*X(J)
30    CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------

6. The data input file FORT.1 needed by the driver program.


7   .0005   107.0   100   81501   19837   39470
  1
  0
  1
  0
  0
  1
  0
  0
  0
  1
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
