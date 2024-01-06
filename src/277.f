      SUBROUTINE MEDIAN(X, Y, N, IWS, XMED, YMED, IFAULT)
C
C        ALGORITHM AS 277.1 APPL.STATIST. (1992) VOL.41, NO.3
C
C        On exit (XMED, YMED) will contain the bivariate median
C        of the points (X(i),Y(i)), i = 1, ... ,N.
C
C     *** WARNING ***
C     The routine as published in the journal contains at least three
C     errors.   This version has been constructed by Alan Miller.   It
C     appears to work correctly but no guarantee is given.
C     N.B. The name of this routine is MEDIAN, not MED as specified in
C     the Structure description given in the journal.
C
C     Auxiliary routines required: M01DAF and M01ZAF from the NAG library.
C
      INTEGER N, IWS(N), IFAULT
      REAL X(N), Y(N), XMED, YMED
C
      EXTERNAL M01DAF, M01ZAF
C
      INTEGER I, J, K, KQ, L, LL, NDIF, NT, NZ, NZERO
      REAL D0, DP, EPS, SMALL, SMALLD, T, TL, TU, W, WW, XKL,
     *      XMAX, XMIN, YKL, YMAX, YMIN
C
      DATA EPS / 1.0E-6 /
C
      IFAULT=0
      CALL M01DAF(X,1,N,'A',IWS,IFAULT)
      CALL M01ZAF(IWS,1,N,IFAULT)
      XMIN=X(IWS(1))
      XMAX=X(IWS(N))
C
      CALL M01DAF(Y,1,N,'A',IWS,IFAULT)
      CALL M01ZAF(IWS,1,N,IFAULT)
      YMIN=Y(IWS(1))
      YMAX=Y(IWS(N))
C
      IFAULT=0
      SMALL=EPS*(XMAX-XMIN)*(YMAX-YMIN)
      SMALLD=SMALL*N
C
      DO 40 NDIF=1,N
	 WW=0.0
	 DO 30 K=1,N-1
	    DO 20 L=K+1,N
	       NT=0
	       XKL=X(K)-X(L)
	       YKL=Y(K)-Y(L)
		 IF(ABS(XKL)+ABS(YKL).LT.SMALL) GO TO 15
	       DO 10 I=1,N
		    W=(Y(I)-Y(L))*XKL-(X(I)-X(L))*YKL
		    WW=WW+ABS(W)
		    IF(W.GT.0.0)NT=NT+1
  10         CONTINUE
  15         CONTINUE
	       IF((ABS(2*NT-N+2).LE.NDIF).AND.WW.GT.SMALLD) GO TO 50
 20        CONTINUE
 30     CONTINUE
	IF(WW.LE.SMALLD) THEN
C
C        The data set is completely collinear.
C
	    IF (MOD(N,2).EQ.0) THEN
		 XMED=(X(IWS(N/2))+X(IWS(N/2+1)))/2
		 YMED=(Y(IWS(N/2))+Y(IWS(N/2+1)))/2
	    ELSE
	       XMED=X(IWS((N+1)/2))
	       YMED=Y(IWS((N+l)/2))
	    END IF
	    RETURN
	 END IF
C
   40 CONTINUE
C
   50 CALL TUTL(X,Y,N,XMIN,YMNI,XMAX,YMAX,K,L,SMALL,TU,TL)
C
C        Start the search along the line through the points (X(k),Y(k))
C        and (X(l),Y(l)), dividing the N-2 other points evenly.
C
   60 KQ=1
C
C        Search for the minimum of the Oja objective function on the
C        line (X(k),Y(k)) to (X(l),Y(l)).
C
      DO 80 I=1,N-1
	   DO 70 J=I+1,N
	   WW=(X(K)-X(L))*(Y(J)-Y(I))-(Y(K)-Y(L))*(X(J)-X(I))
	   IF(ABS(WW).LT.SMALL) GO TO 70
		T=((X(L)-X(J))*(Y(I)-Y(J))-(Y(L)-Y(J))*(X(I)-X(J)))/WW
		IF(T.LT.TU.AND.T.GT.TL)THEN
		   CALL DER(X,Y,N,K,L,T,IWS,SMALL,NZERO,KQ,DP,D0)
c                  A=DP+D0                     ! Redundant instruction
		   IF(D0.GE.ABS(DP)-SMALLD) GO TO 90
		   IF(DP+D0.LT.0.0)TL=T
		   IF(DP+D0.GT.0.0)TU=T
		END IF
   70    CONTINUE
   80 CONTINUE
C
C        The local minimum on the line (X(k),Y(k)) - (X(l),Y(l))
C        was not found.
C
      IFAULT=IFAULT+2
      IF(IFAULT.GT.6) RETURN
      L=N-1
      K=N
   90 KQ=0
      I=K
      J=L
C
      LL=1
      NZ=NZERO
  100 K=INT(IWS(LL)/REAL(N))
      L=MOD(IWS(LL),N)+1
      WW=(X(K)-X(L))*(Y(J)-Y(I))-(Y(K)-Y(L))*(X(J)-X(I))
      IF(ABS(WW).GE.SMALL) THEN
	 T=((X(L)-X(J))*(Y(I)-Y(J))-(Y(L)-Y(J))*(X(I)-X(J)))/WW
	 CALL DER(X,Y,N,K,L,T,IWS,SMALL,NZERO,KQ,DP,D0)
	 IF (D0.LT.ABS(DP)-SMALL) THEN
	     CALL TUTL(X,Y,N,XMIN,YWIN,XMAX,YMAX,K,L,SMALL,TU,TL)
	     IF(DP+D0.LT.0.0)THEN
		TL=T
	    ELSE
		TU=T
	    END IF
	    GO TO 60
	 END IF
	 LL=LL+1
      END IF
C
      IF(LL.EQ.NZ+1)THEN
	 XMED=X(L)+T*(X(K)-X(L))
	 YMED=Y(L)+T*(Y(K)-Y(L))
	 RETURN
      END IF
C
      GO TO 100
      END
C
      SUBROUTINE DER(X,Y,N,K,L,T,IWS,SMALL,NZERO,KQ,DP,D0)
C
C        ALGORITHM AS 277.2 APPL.STATIST. (1992) VOL.41, NO.3
C
C        DP-D0 and DP+D0 are the directional derivatives of the Oja
C         Objective function just before and after the point
C         (X(l)+t*(X(k)-X(l),Y(l)+t*(Y(k)-Y(l))
C
	INTEGER N, K, L, IWS(N), NZERO, KQ
	REAL X(N), Y(N), T, SMALL, DP, D0
C
      REAL DIF, SGN, SMALLD, TT, WW, XKL, YKL
      INTEGER II, JJ
C
      DATA SMALLD / 1E-6 /
C
      DP=0.0
      D0=0.0
      NZERO=0
      XKL=X(K)-X(L)
      YKL=Y(K)-Y(L)
C
      DO 20 II=1,N-1
	  DO 10 JJ=II+1,N
	  WW=(Y(JJ)-Y(II))*XKL-(X(JJ)-X(II))*YKL
	  IF(ABS(WW).LE.SMALL) GO TO 10
	  TT=((X(L)-X(JJ))*(Y(II)-Y(JJ))-(Y(L)-Y(JJ))*(X(II)-X(JJ)))/WW
	  DIF=T-TT
	  IF(ABS(DIF).LE.SMALLD)THEN
	     NZERO=NZERO+1
	     NZERO=MIN(NZERO,N)
	     IF(KQ.NE.0)IWS(NZERO)=N*II+JJ-1
		D0=D0+ABS(WW)
	     ELSE
		DP=DP+ABS(WW)*SGN(DIF)
	     END IF
   10     CONTINUE
   20 CONTINUE
      RETURN
      END
C
      FUNCTION SGN(X)
C
C          ALGORITHM AS 277.3 APPL.STATIST. (1992) VOL.41, NO.3
C
      REAL SGN, X
      SGN=1.0
      IF (X.LT.0.0)SGN=-1.0
      RETURN
      END
C
      SUBROUTINE TUTL(X,Y,N,XMIN,YMIN,XMAX,YMAX,K,L,SMALL,TU,TL)
C
C        ALGORITHM AS 277.4 APPL.STATIST. (1992) VOL.41, NO.3
C
C        This subroutine calculates the upper (TU) and lower limit (TL)
C        for parameter T on the line (X(l)+t*(X(k)-X(l)), Y(l)+t*(Y(k)-Y(l))
C        inside the rectangle with vertices (XMIN,YMIN), (XMIN,YMAX),
C        (XMAX,YMIN) and (XMAX,YMAX).
C
	INTEGER N, K, L
	REAL X(N), Y(N), XMIN, YMIN, XMAX, YMAX, SMALL, TU, TL
	REAL T1, T2, T3, T4, VBIG
      DATA VBIG / 1E38 /
      T1=-VBIG
      T2=-T1
      IF(ABS(X(K)-X(L)).GT.SMALL)THEN
	 T1=-(X(L)-XMIN)/(X(K)-X(L))
	 T2= (XMAX-X(L))/(X(K)-X(L))
      END IF
      T3=-VBIG
      T4=-T3
      IF(ABS(Y(K)-Y(L)).GT.SMALL)THEN
	 T3=-(Y(L)-YMIN)/(Y(K)-Y(L))
	 T4= (YMAX-Y(L))/(Y(K)-Y(L))
      END IF
      TU=MIN(MAX(T1,T2),MAX(T3,T4))+SMALL
      TL=MAX(MIN(T1,T2),MIN(T3,T4))-SMALL
      RETURN
      END

