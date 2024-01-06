      SUBROUTINE TDIAG (N,TOL,A,D,E,Z,MAXDIM,IFAULT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C 
      DIMENSION A(MAXDIM,MAXDIM), D(MAXDIM), E(MAXDIM), Z(MAXDIM,MAXDIM)
      DATA ZERO/0.0/,ONE/1.0/
      DATA ETA/1.0D-37/,EPS/1.0D-14/
C 
C       Algorithm as 60.1 appl.statist. (1973) vol.22 no.2
C 
C        reduces real symmetric matrix to tridiagonal form
C 
C       tol is a machine dependent constant , tol = eta/eps , where
C       eta is the smallest positive number representable in the
C       computer and eps is the smallest positive number for which
C       1+eps.ne.1.
C 
C         eta=eps*tol
C         eps=0.7105427358e-14
C         tol=0.3131513063e-293
C         precis=1.0e-14
C 
C         nb
C           real constants must be le 15 decimal digits
C           the range of a real constant is from 1.0e-293 to 1.0e+322
C 
      TOL=ETA/EPS
      IFAULT=1
      IF (N.LE.1) RETURN
      IFAULT=0
      DO 10 I=1,N
        DO 10 J=1,I
 10     Z(I,J)=A(I,J)
      I=N
      DO 110 I1=2,N
        L=I-2
        F=Z(I,I-1)
        G=ZERO
        IF (L.LT.1) GO TO 30
        DO 20 K=1,L
 20       G=G+Z(I,K)**2
 30     H=G+F*F
C 
C       if g is too small for orthogonality to be guaranteed, the
C       transformation is skipped
C 
        IF (G.GT.TOL) GO TO 40
        E(I)=F
        D(I)=ZERO
        GO TO 100
 40     L=L+1
	G=SQRT(H)
        IF (F.GE.ZERO) G=-G
        E(I)=G
        H=H-F*G
        Z(I,I-1)=F-G
        F=ZERO
        DO 80 J=1,L
          Z(J,I)=Z(I,J)/H
          G=ZERO
C 
C       form element of a * u
C 
          DO 50 K=1,J
 50         G=G+Z(J,K)*Z(I,K)
          IF (J.GE.L) GO TO 70
          J1=J+1
          DO 60 K=J1,L
 60         G=G+Z(K,J)*Z(I,K)
C 
C       form element of p
C 
 70       E(J)=G/H
          F=F+G*Z(J,I)
 80     CONTINUE
C 
C       form k
C 
        HH=F/(H+H)
C 
C       form reduced a
C 
        DO 90 J=1,L
          F=Z(I,J)
          G=E(J)-HH*F
          E(J)=G
          DO 90 K=1,J
          Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
 90     CONTINUE
        D(I)=H
 100    I=I-1
 110  CONTINUE
      D(1)=ZERO
      E(1)=ZERO
C 
C       accumulation of transformation matrices
C 
      DO 160 I=1,N
        L=I-1
        IF (D(I).EQ.ZERO.OR.L.EQ.0) GO TO 140
        DO 130 J=1,L
          G=ZERO
          DO 120 K=1,L
 120        G=G+Z(I,K)*Z(K,J)
          DO 130 K=1,L
          Z(K,J)=Z(K,J)-G*Z(K,I)
 130    CONTINUE
 140    D(I)=Z(I,I)
        Z(I,I)=ONE
        IF (L.EQ.0) GO TO 160
        DO 150 J=1,L
          Z(I,J)=ZERO
          Z(J,I)=ZERO
 150    CONTINUE
 160  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE LRVT (N,PRECIS,D,E,Z,IFAULT,MAXDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C 
C       algorithm as 60.2 appl.statist. (1973) vol.22, no.2
C 
C        finds latent roots and vectors of tridiagonal matrix
C 
      DIMENSION D(MAXDIM), E(MAXDIM), Z(MAXDIM,MAXDIM)
      DATA MITS/30/,ZERO/0.0/,ONE/1.0/,TWO/2.0/
C 
      PRECIS=1.0D-14
      IFAULT=2
      IF (N.LE.1) RETURN
      IFAULT=1
      N1=N-1
      DO 10 I=2,N
 10     E(I-1)=E(I)
      E(N)=ZERO
      B=ZERO
      F=ZERO
      DO 100 L=1,N
        JJ=0
	H=PRECIS*(ABS(D(L))+ABS(E(L)))
        IF (B.LT.H) B=H
C 
C       look for small sub-diagonal element
C 
        DO 20 M1=L,N
          M=M1
	  IF (ABS(E(M)).LE.B) GO TO 30
 20     CONTINUE
 30     IF (M.EQ.L) GO TO 90
 40     IF (JJ.EQ.MITS) RETURN
        JJ=JJ+1
C 
C       form shift
C 
        P=(D(L+1)-D(L))/(TWO*E(L))
	R=SQRT(P*P+ONE)
        PR=P+R
        IF (P.LT.ZERO) PR=P-R
        H=D(L)-E(L)/PR
        DO 50 I=L,N
 50       D(I)=D(I)-H
        F=F+H
C 
C       ql transformation
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
	  IF (ABS(P).GE.ABS(E(I))) GO TO 60
          C=P/E(I)
	  R=SQRT(C*C+ONE)
          E(J)=S*E(I)*R
          S=ONE/R
          C=C/R
          GO TO 70
 60       C=E(I)/P
          R=SQRT(C*C+ONE)
          E(J)=S*P*R
          S=C/R
          C=ONE/R
 70       P=C*D(I)-S*G
          D(J)=H+S*(C*G+S*D(I))
C 
C       form vector
C 
          DO 80 K=1,N
          H=Z(K,J)
          Z(K,J)=S*Z(K,I)+C*H
          Z(K,I)=C*Z(K,I)-S*H
 80     CONTINUE
        E(L)=S*P
        D(L)=C*P
	IF (ABS(E(L)).GT.B) GO TO 40
 90     D(L)=D(L)+F
 100  CONTINUE
C 
C       order latent roots and vectors
C 
      DO 130 I=1,N1
        K=I
        P=D(I)
        I1=I+1
        DO 110 J=I1,N
          IF (D(J).LE.P) GO TO 110
          K=J
          P=D(J)
 110    CONTINUE
        IF (K.EQ.I) GO TO 130
        D(K)=D(I)
        D(I)=P
        DO 120 J=1,N
          P=Z(J,I)
          Z(J,I)=Z(J,K)
          Z(J,K)=P
 120    CONTINUE
 130  CONTINUE
      IFAULT=0
      RETURN
      END
