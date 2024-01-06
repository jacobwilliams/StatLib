      SUBROUTINE INCLUD (NP,NRBAR,WEIGHT,XROW,YELEM,D,RBAR,THETAB,SSERR,
     1 IFAULT)
C 
C     Algorithm AS 75.1 Appl. Statist. (1974) Vol.23, No.3, p448
C 
C     Calling this subroutine updates d, rbar, thetab and sserr
C     by the inclusion of xrow, yelem with specified weight.
C 
      DOUBLE PRECISION XROW(NP),D(NP),RBAR(NRBAR),THETAB(NP),WEIGHT,
     1 SSERR,CBAR,DI,DPI,SBAR,W,XI,XK,Y,ZERO,YELEM
C 
      DATA ZERO/0.0D0/
C 
C        check input parameters
C 
      IFAULT=1
      IF (NP.LT.1.OR.NRBAR.LE.NP*(NP-1)/2) RETURN
      IFAULT=0
C 
      W=WEIGHT
      Y=YELEM
      DO 30 I=1,NP
C 
C     Skip unnecessary transformations.  Test on exact zeros must
C     be used or stability can be destroyed.
C 
        IF (W.EQ.ZERO) RETURN
        IF (XROW(I).EQ.ZERO) GO TO 30
        XI=XROW(I)
        DI=D(I)
        DPI=DI+W*XI*XI
        CBAR=DI/DPI
        SBAR=W*XI/DPI
        W=CBAR*W
        D(I)=DPI
        IF (I.EQ.NP) GO TO 20
        NEXTR=(I-1)*(NP+NP-I)/2+1
        IP=I+1
        DO 10 K=IP,NP
          XK=XROW(K)
          XROW(K)=XK-XI*RBAR(NEXTR)
          RBAR(NEXTR)=CBAR*RBAR(NEXTR)+SBAR*XK
          NEXTR=NEXTR+1
 10     CONTINUE
 20     XK=Y
        Y=XK-XI*THETAB(I)
        THETAB(I)=CBAR*THETAB(I)+SBAR*XK
 30   CONTINUE
      SSERR=SSERR+W*Y*Y
      RETURN
      END
C
C
      SUBROUTINE CONFND (NP,NRBAR,J,RBAR,CONTRA,IFAULT)
C 
C     Algorithm AS 75.2 Appl. Statist. (1974) Vol. 23, No. 3, P448
C 
C     Calling this subroutine obtains the contrast which could not
C     be estimated if D(j) were assumed to be zero, that is, obtains the
C     linear combination of the first j columns which would be zero.  Th
C     obtained by setting the first j-1 elements of contra to the soluti
C     of the triangular system formed by the first j-1 rows and columns
C     rbar with the first j-1 elements of the jth column as right hand
C     side, setting the jth element of contra to -1, and setting the
C     remaining elements of contra to zero.
C 
      DOUBLE PRECISION RBAR(NRBAR),CONTRA(NP),ZERO,ONE
C 
      DATA ZERO/0.0D0/,ONE/1.0D0/
C 
C         check input parameters
C 
      IFAULT=1
      IF (NP.LT.1.OR.NRBAR.LE.NP*(NP-1)/2) RETURN
      IFAULT=0
C 
      JM=J-1
      IF (J.EQ.NP) GO TO 20
      JP=J+1
      DO 10 I=JP,NP
 10     CONTRA(I)=ZERO
 20   CONTRA(J)=-ONE
      IF (J.EQ.1) RETURN
      DO 40 IJ=1,JM
        I=J-IJ
        NEXTR=(I-1)*(NP+NP-I)/2+1
        K=NEXTR+J-I-1
        CONTRA(I)=RBAR(K)
        IF (I.EQ.JM) GO TO 40
        IP=I+1
        DO 30 K=IP,JM
          CONTRA(I)=CONTRA(I)-RBAR(NEXTR)*CONTRA(K)
          NEXTR=NEXTR+1
 30     CONTINUE
 40   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE SSDCMP (NP,D,THETAB,SS,IFAULT)
C 
C     Algorithm AS75.3  Appl. Statist. (1974) Vol.23, No. 3, P448
C 
C     Calling this subroutine computes the np components of the sum
C     of squares decomposition from D and thetab.
C 
      DOUBLE PRECISION D(NP),THETAB(NP),SS(NP)
C 
C        check input parameters
C 
      IFAULT=1
      IF (NP.LT.1) RETURN
      IFAULT=0
C 
      DO 10 I=1,NP
 10     SS(I)=D(I)*THETAB(I)**2
      RETURN
      END
C
C
      SUBROUTINE REGRSS (NP,NRBAR,RBAR,THETAB,BETA,IFAULT)
C 
C     Algorithm AS 75.4  Appl. Statist. (1974), Vol. 23, No. 3, p448
C 
C     Calling this subroutine obtains beta by back-substitution in
C     the triangular system rbar and thetab.
C 
      DOUBLE PRECISION RBAR(NRBAR),THETAB(NP),BETA(NP)
C 
C        check input parameters
C 
      IFAULT=1
      IF (NP.LT.1.OR.NRBAR.LE.NP*(NP-1)/2) RETURN
      IFAULT=0
C 
      DO 20 J=1,NP
        I=NP-J+1
        BETA(I)=THETAB(I)
        NEXTR=(I-1)*(NP+NP-I)/2+1
        IP=I+1
        DO 10 K=IP,NP
          BETA(I)=BETA(I)-RBAR(NEXTR)*BETA(K)
          NEXTR=NEXTR+1
 10     CONTINUE
 20   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE ANALYZ (MAXP,MAXR,NP,NRBAR,INPUT,OUTPUT,X,DD,THETA,R,
     1 IFAULT)
C 
C        Algorithm AS 75.5  Appl. Statist. (1974) Vol. 23, P. 448
C 
C        This subroutine reads data from channel input and
C        produces least squares variance component analyses on
C        channel output.  The input medium is assumed to
C        contain first a value for the integer np, the number of
C        independent variates, and then a sequence of records
C        containing the rows of x and y together with weights.
C        The number zero (read as a weight) will indicate the end of
C        each data set, and upon reading this an analysis will be
C        produced.  Successive analyses will be done until np,
C        the number of independent variates, is read as zero.
C 
C        No need to check the ifault parameter after calling
C        auxiliary routines as parameters have already been
C        checked here.
C 
      INTEGER OUTPUT
      LOGICAL FIRST
      DOUBLE PRECISION X(MAXP),DD(MAXP),THETA(MAXP),R(MAXR),ERROR,W,Y,
     1 ZERO,ONE,EPS,TOL
C 
      DATA ZERO,ONE,EPS,TOL/0.0D0,1.0D0,1.0D-8,1.0D-4/
C 
C        input formats
C 
 10   FORMAT (I3)
 20   FORMAT (F10.6)
 30   FORMAT (F10.6)
 40   FORMAT (I3)
 50   FORMAT (F10.6)
C 
C        output formats
C 
 60   FORMAT (1H1,I5,18H observations read/19H diagonal matrix is/(1X,
     1 5G15.8))
 70   FORMAT (21H confounded contrasts)
 80   FORMAT (1X,5G15.8)
 90   FORMAT (29H sum of squares decomposition/(1X,5G15.8))
 100  FORMAT (22H sum of squares error ,G15.8)
 110  FORMAT (24H regression coefficients,(1X,5G15.8))
C 
C        read np, check its value and initialise arrays
C 
      IFAULT=0
 120  READ (INPUT,10) NP
      IF (NP.EQ.0) RETURN
      NRBAR=NP*(NP-1)/2
      IF (NP.LT.0.OR.NP.GT.MAXP.OR.NRBAR.GT.MAXR) GO TO 240
      DO 130 K=1,NP
        DD(K)=ZERO
        THETA(K)=ZERO
 130  CONTINUE
      DO 140 K=1,NRBAR
 140    R(K)=ZERO
      ERROR=ZERO
      N=0
C 
C        read weight and check its value
C 
 150  READ (INPUT,20) W
      IF (W.EQ.ZERO) GO TO 190
      IF (W.GT.ZERO) N=N+1
      IF (W.LT.ZERO) N=N-1
C 
C        read the y-value and all corresponding x-values
C 
      READ (INPUT,30) Y
      DO 160 K=1,NP
 160    X(K)=ZERO
 170  READ (INPUT,40) K
      IF (K.LT.0.OR.K.GT.NP) GO TO 250
      IF (K.EQ.0) GO TO 180
      READ (INPUT,50) X(K)
      GO TO 170
 180  CALL INCLUD (NP,NRBAR,W,X,Y,DD,R,THETA,ERROR,IFAIL)
      GO TO 150
C 
C        begin output
C 
 190  WRITE (OUTPUT,60) N,(DD(I),I=1,NP)
      FIRST=.TRUE.
      DO 230 J=1,NP
	IF (ABS(DD(J)).GE.EPS) GO TO 230
C 
C        confounding discovered
C 
        IF (.NOT.FIRST) GO TO 200
        FIRST=.FALSE.
        WRITE (OUTPUT,70)
 200    CALL CONFND (NP,NRBAR,J,R,X,IFAIL)
        WRITE (OUTPUT,80) (X(I),I=1,NP)
C 
C        choose resolving constraint
C 
        IF (J.EQ.1) GO TO 220
        JM=J-1
        DO 210 K=1,JM
	  IF (ABS(X(K)).LE.TOL) GO TO 210
          X(K)=ZERO
          GO TO 220
 210    CONTINUE
 220    CALL INCLUD (NP,NRBAR,ONE,X,ZERO,DD,R,THETA,ERROR,IFAIL)
 230  CONTINUE
      CALL SSDCMP (NP,DD,THETA,X,IFAIL)
      WRITE (OUTPUT,90) (X(I),I=1,NP)
      WRITE (OUTPUT,100) ERROR
      CALL REGRSS (NP,NRBAR,R,THETA,X,IFAIL)
      WRITE (OUTPUT,110) (X(I),I=1,NP)
      GO TO 120
C 
C        error returns
C 
 240  IFAULT=2
      RETURN
 250  IFAULT=3
      RETURN
      END
