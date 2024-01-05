C---- SUBROUTINE ARAVSS -----------------------------------------------
C
C  Computes sufficient statistics for analysis of variance table when
C  the errors follow AR(1) normal process, given the design matrix X(N,P),
C  observation matrix Y(N,T) and assuming that XXINV(P,P) contains
C  the inverse matrix of Xtranspose*X.
C
C----------------------------------------------------------------------
      SUBROUTINE ARAVSS(Y,N,T,P,QS,XD1,XD2,YBAR,X,XXINV,IFAULT)
C
C     ALGORITHM AS 246.1  APPL. STATIST. (1989), VOL.38, NO. 2
C
      INTEGER T,P
      REAL Y(N,T),X(N,P),DLIN(2),DQUAD(3),YBAR(T)
      REAL XD1(P),XD2(P),XXINV(P,P),QS(14),ZERO,PROD
      DATA ZERO/0.0/
C
C      Computes admissibility of parameters
C
      IFAULT=0
      IF(P.LT.1) IFAULT=1
      IF(T.LT.3) IFAULT=IFAULT+2
      IF(N.LT.P) IFAULT=IFAULT+4
      J=0
      DO 10 I=1,N
  10    IF(X(I,1).NE.1.0) J=1
      IF(J.NE.0) IFAULT=IFAULT+8
      IF(IFAULT.GT.0) RETURN
C
C      Sets initial values to zero
C
      DO 20 I=1,14
  20    QS(I)=ZERO
      DO 30 I=1,P
        XD1(I)=ZERO
        XD2(I)=ZERO
  30    CONTINUE
C
C      For each individual time series, sufficient
C      statistics are computed and accumulated
C
      DO 60 I=1,N
        DO 40 J=1,T
  40      YBAR(J)=Y(I,J)
        CALL ARSUFF(YBAR,T,1,DLIN,DQUAD,IFAULT)
        QS(1)=QS(1)+DQUAD(1)
        QS(2)=QS(2)+DQUAD(2)
        QS(3)=QS(3)+DQUAD(3)
        QS(4)=QS(4)+DLIN(1)**2
        QS(5)=QS(5)+DLIN(1)*DLIN(2)
        QS(6)=QS(6)+DLIN(2)**2
        DO 50 J=1,P
          XD1(J)=XD1(J)+X(I,J)*DLIN(1)
          XD2(J)=XD2(J)+X(I,J)*DLIN(2)
  50      CONTINUE
  60    CONTINUE
C
C      Computes d(i)'X(X'X)**(-1)X'd(j) for i,j=1,2
C
      QS(7)=PROD(XD1,XXINV,XD1,P)
      QS(8)=PROD(XD1,XXINV,XD2,P)
      QS(9)=PROD(XD2,XXINV,XD2,P)
C
C      Computing column means and their sufficient statistics
C
      DO 80 J=1,T
        YBAR(J)=ZERO
        DO 70 I=1,N
  70      YBAR(J)=YBAR(J)+Y(I,J)
        YBAR(J)=YBAR(J)/N
  80    CONTINUE
      CALL ARSUFF(YBAR,T,1,DLIN,DQUAD,IFAULT)
      QS(10)=DQUAD(1)
      QS(11)=DQUAD(2)
      QS(12)=DQUAD(3)
      QS(13)=DLIN(1)
      QS(14)=DLIN(2)
      RETURN
      END


C---- SUBROUTINE ARAVLK ---------------------------------------------
C
C  Computes the log-likelihood and the components of the analysis of
C  variance table supposing that the errors follow an AR(1) stationary
C  normal process, and that the autoregressive parameter is equal to
C  alpha. The sufficient statistics must have been previously computed
C  using SUBROUTINE ARAVSS.
C
C----------------------------------------------------------------------
      SUBROUTINE ARAVLK(ALPHA,N,T,P,QS,XD1,XD2,XXINV,LOGLIK,GM,
     +                  BETA,SE,SIGMA2,PSI,SS,DF,RANDOM,TIME,IFAULT)
C
C     ALGORITHM AS 246.2  APPL. STATIST. (1989), VOL.38, NO. 2
C
      INTEGER T,P,DF(6),DFR
      REAL LOGLIK,SS(6),BETA(P),SE(P)
      REAL XD1(P),XD2(P),XXINV(P,P),QS(14),DELTA,RSS,DA,HA,BA
      LOGICAL RANDOM,TIME
      DELTA(B11,B12,B22)=B11-2.0*ALPHA*B12+ALPHA*ALPHA*B22
C
C      Checks whether alpha is between -1 and 1
C
      IF(ALPHA**2.LT.1.0) THEN
        IFAULT=0
      ELSE
        IFAULT=1
        RETURN
        ENDIF
C
      FN=FLOAT(N)
      TAU=T-(T-2)*ALPHA
      CTAU=(1.0-ALPHA)/TAU
C
C      Computes D(alpha), H(alpha), B(alpha)
C
      DA=DELTA(QS(1),QS(2),QS(3))
      HA=DELTA(QS(4),QS(5),QS(6))
      BA=DELTA(QS(7),QS(8),QS(9))
      RSS=DA
      DFR=N*T
C
C      Subtracts required sum of squares from RSS
C      and computes profile log-likelihood
C
      GM=(QS(13)-ALPHA*QS(14))/TAU
      SS(1)=CTAU*FN*(GM*TAU)**2
      DF(1)=1
      SS(6)=DA
      DF(6)=N*T
      DO 10 I=2,5
        SS(I)=0.0
        DF(I)=0
  10    CONTINUE
      IF(P.GT.1) THEN
        SS(2)=CTAU*BA-SS(1)
        DF(2)=P-1
        ENDIF
C
C      Branch for time effect
C
      IF(TIME) THEN
        SS(4)=FN*DELTA(QS(10),QS(11),QS(12))-SS(1)
        DF(4)=T-1
        RSS=RSS-SS(4)
        DFR=DFR-DF(4)
        ENDIF

      IF(RANDOM) THEN
C
C         Random effect model
C
        SS(3)=(HA-BA)*CTAU
        DF(3)=N-P
        PSI=SS(3)/FN
        RSS=RSS-CTAU*HA
        DFR=DFR-N
        SIGMA2=RSS/((T-1)*FN)
        IF(PSI.LT.SIGMA2) THEN
          SIGMA2=(RSS+CTAU*(HA-BA))/(FN*T)
          PSI=SIGMA2
          ENDIF
        LOGLIK=(T-1)*ALOG(SIGMA2)+ALOG(PSI)-ALOG(1.0-ALPHA**2)
      ELSE
C
C         Fixed effect model
C
        RSS=RSS-BA*CTAU
        DFR=DFR-P
        SIGMA2=RSS/(N*T)
        PSI=SIGMA2
        LOGLIK=T*ALOG(SIGMA2)-ALOG(1.0-ALPHA**2)
        ENDIF

      SS(5)=RSS
      DF(5)=DFR
      LOGLIK=-0.5*FN*LOGLIK
C
C      Computes estimates of BETA()
C
      DO 20 I=1,P
        BETA(I)=0.0
        SE(I)=(XD1(I)-ALPHA*XD2(I))/TAU
  20    CONTINUE
      DO 40 I=1,P
        DO 30 J=1,P
  30      BETA(I)=BETA(I)+XXINV(I,J)*SE(J)
  40    CONTINUE
C
C      Computes standard errors for BETA()
C
      U=PSI/((1.0-ALPHA)*TAU)
      DO 50 I=1,P
  50    SE(I)=SQRT(U*XXINV(I,I))
      RETURN
      END


C---- SUBROUTINE ARSUFF ------------------------------------------------
C
C  Computes sufficient statistics for stationary autoregressive
C  normal process Y(1),...,Y(T) of order ORDER.
C  The output values are stored in vector DLIN (linear terms) and
C  DQUAD (quadratic terms, packed upper triangular).
C  A. Azzalini (Dipart. di Statistica, Universita` di Padova) Feb.1982
C-----------------------------------------------------------------------
      SUBROUTINE ARSUFF(Y,T,ORDER,DLIN,DQUAD,IFAULT)
C
C     ALGORITHM AS 246.3  APPL. STATIST. (1989), VOL.38, NO. 2
C
      INTEGER T,ORDER,R,S,SMR
      REAL Y(T),DLIN(ORDER+1),DQUAD((ORDER+1)*(ORDER+2)/2)
      IFAULT=0
      IF(ORDER.LT.0) IFAULT=1
      IF(T.LE.2*ORDER) IFAULT=IFAULT+2
      IF(IFAULT.GT.0) RETURN
      M=ORDER+1
      N=M*(ORDER+2)/2
      SUM1=0.
      SUM2=0.
      DO 10 J=M,T-ORDER
        YJ=Y(J)
        SUM1=SUM1+YJ
        SUM2=SUM2+YJ*YJ
  10    CONTINUE
      DLIN(M)=SUM1
      DQUAD(N)=SUM2
      DO 20 K=1,ORDER
        R=ORDER-K+1
        YJ=Y(R)
        SUM1=SUM1+YJ
        SUM2=SUM2+YJ*YJ
        YJ=Y(T-R+1)
        SUM1=SUM1+YJ
        SUM2=SUM2+YJ*YJ
        DLIN(R)=SUM1
        DQUAD(R*(R+1)/2)=SUM2
  20    CONTINUE
      DO 50 SMR=1,ORDER
        SUM2=0.
        DO 30 J=M,T-ORDER+SMR
          SUM2=SUM2+Y(J)*Y(J-SMR)
  30      CONTINUE
        DQUAD(N-SMR)=SUM2
        DO 40 R=ORDER-SMR,1,-1
          S=R+SMR
          SUM2=SUM2+Y(R)*Y(S)+Y(T-R+1)*Y(T-S+1)
          DQUAD(S*(S-1)/2+R)=SUM2
  40      CONTINUE
  50    CONTINUE
      RETURN
      END

C--------------------------------------------------------------
C
      FUNCTION PROD(XD1,XX,XD2,M)
      REAL A,ZERO,XD1(M),XD2(M),XX(M,M)
      DATA ZERO/0.0/
      PROD=ZERO
      DO 20 I=1,M
        A=ZERO
        DO 10 J=1,M
  10      A=A+XX(I,J)*XD2(J)
  20    PROD=PROD+XD1(I)*A
      RETURN
      END
