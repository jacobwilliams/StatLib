C This file contains AS6 and the enhanced version ASR44.   See AS7 also.
C 
C 
      SUBROUTINE CHOL (A,N,NN,U,NULLTY,IFAULT)
C 
C       Algorithm AS6, Applied Statistics, vol.17, (1968)
C 
C       Given a symmetric matrix order n as lower triangle in a( )
C       calculates an upper triangle, u( ), such that uprime * u = a.
C       a must be positive semi-definite.  eta is set to multiplying
C       factor determining effective zero for pivot.
C 
C       arguments:-
C       a()     = input, a +ve definite matrix stored in lower-triangula
C                 form.
C       n       = input, the order of a
C       nn      = input, the size of the a and u arrays      n*(n+1)/2
C       u()     = output, a lower triangular matrix such that u*u' = a.
C                 a & u may occupy the same locations.
C       nullty  = output, the rank deficiency of a.
C       ifault  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 3 if nn < n*(n+1)/2
C                       = 0 otherwise
C 
C***********************************************************************
C 
      DOUBLE PRECISION A(NN),U(NN),ETA,ETA2,X,W,ZERO
C 
C       The value of eta will depend on the word-length of the
C       computer being used.  See introductory text.
C 
      DATA ETA,ZERO/1.D-9,0.0D0/
C 
      IFAULT=1
      IF (N.LE.0) RETURN
      IFAULT=3
      IF (NN.LT.N*(N+1)/2) RETURN
      IFAULT=2
      NULLTY=0
      J=1
      K=0
      ETA2=ETA*ETA
      II=0
C 
C       Factorize column by column, icol = column no.
C 
      DO 80 ICOL=1,N
        II=II+ICOL
        X=ETA2*A(II)
        L=0
        KK=0
C 
C       IROW = row number within column ICOL
C 
        DO 40 IROW=1,ICOL
          KK=KK+IROW
          K=K+1
          W=A(K)
          M=J
          DO 10 I=1,IROW
            L=L+1
            IF (I.EQ.IROW) GO TO 20
            W=W-U(L)*U(M)
            M=M+1
 10       CONTINUE
 20       IF (IROW.EQ.ICOL) GO TO 50
          IF (U(L).EQ.ZERO) GO TO 30
          U(K)=W/U(L)
          GO TO 40
 30       IF (W*W.GT.ABS(X*A(KK))) RETURN
          U(K)=ZERO
 40     CONTINUE
 50     IF (ABS(W).LE.ABS(ETA*A(K))) GO TO 60
        IF (W.LT.ZERO) RETURN
        U(K)=SQRT(W)
        GO TO 70
 60     U(K)=ZERO
        NULLTY=NULLTY+1
 70     J=J+ICOL
 80   CONTINUE
      IFAULT=0
      END
C
C
C
C
      SUBROUTINE SUBCHL (A,B,N,U,NULLTY,IFAULT,NDIM,DET)
C 
C     REMARK ASR 44  APPL. STATIST. (1982) VOL. 31, NO. 3
C 
C     A revised and enhanced version of
C     ALGORITHM AS 6  APPL. STATIST. (1968) VOL. 17, NO. 2
C 
C     Given a symmetric matrix of order N as lower triangle in A(),
C     calculates an upper triangle, U(), such that U'U = the sub-matrix
C     of A whose rows and columns are specified in the integer array
C     B().
C     U() may coincide with A().   A() must be +ve semi-definite.
C     ETA is set to multiplying factor determining effective zero for
C     a pivot.
C     NULLTY is returned as number of effective zero pivots.
C     IFAULT is returned as 1 if N <= 0, 2 if A() is not +ve semi-
C     definite, otherwise 0 is returned.
C 
      DOUBLE PRECISION A(NDIM),U(NDIM),DET
      INTEGER B(N)
C 
C     Local variables
C 
      DOUBLE PRECISION ETA,ONE,ZERO,W,ETA2
C 
C     The value of ETA below will depend upon the word length of the
C     computer being used.
C 
      DATA ETA/1.D-14/,ONE/1.D0/,ZERO/0.D0/
C 
      IFAULT=1
      IF (N.LE.0) GO TO 90
      IFAULT=2
      NULLTY=0
      DET=ONE
      J=1
      K=0
      ETA2=ETA*ETA
      DO 80 ICOL=1,N
        IJ=B(ICOL)*(B(ICOL)-1)/2
        II=IJ+B(ICOL)
        X=ETA2*A(II)
        L=0
        DO 40 IROW=1,ICOL
          KK=B(IROW)*(B(IROW)+1)/2
          K=K+1
          JJ=IJ+B(IROW)
          W=A(JJ)
          M=J
          DO 10 I=1,IROW
            L=L+1
            IF (I.EQ.IROW) GO TO 20
            W=W-U(L)*U(M)
            M=M+1
 10       CONTINUE
 20       IF (IROW.EQ.ICOL) GO TO 50
          IF (U(L).EQ.ZERO) GO TO 30
          U(K)=W/U(L)
          GO TO 40
 30       IF (W*W.GT.ABS(X*A(KK))) GO TO 90
          U(K)=ZERO
 40     CONTINUE
 50     IF (ABS(W).LE.ABS(ETA*A(KK))) GO TO 60
        IF (W.LT.ZERO) GO TO 90
        U(K)=SQRT(W)
        GO TO 70
 60     U(K)=ZERO
        NULLTY=NULLTY+1
 70     J=J+ICOL
        DET=DET*U(K)*U(K)
 80   CONTINUE
C 
      IFAULT=0
 90   RETURN
      END
