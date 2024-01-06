! Below are 2 versions of this algorithm, the first as submitted
! to the RSS in Fortran 77; the second is my Fortran 90 version.
! The author has given permission for his version to be submitted
! to the apstat collection.

! As a `title' in the index I suggest:
! Unconstrained variable metric function minimization without derivatives.

! The two versions are separated by a line of !!!!!!!!'s

! Alan
! P.S. I have tested the F90 version.   As it was derived from the F77
! version, that should work too!



C Algorithm AS 319 and test program

C----------------------------------------------------------------------
        PROGRAM  VAR
C----------------------------------------------------------------------
C       A PROGRAM TO IMPLEMENT A QUASI-NEWTON METHOD.
C       USING NEW ALGORITHM VARMET   JULY 1994
C----------------------------------------------------------------------
C
      COMMON /FUNERR/LER
      COMMON /TEST/IG,IFN
      EXTERNAL FUN
      LOGICAL LER
C
      INTEGER N, NMAX
      PARAMETER (N=2, NMAX=50)
      INTEGER IG,IFN,GRADTL,MESS,MAXFN,IER
      PARAMETER (GRADTL = 12, MAXFN = 1000, MESS = 6)
      DOUBLE PRECISION X,XTMP,W,FP
      DIMENSION X(N),XTMP(N),W(225)
C
      WRITE(*,*)'  '
      WRITE(*,*) 'INPUT YOUR STARTING GUESS'
      DO 12 I=1,N
        READ(*,*) XTMP(I)
   12 CONTINUE
      WRITE(*,*)' '
      WRITE(*,*)'INITIALIZATION COMPLETE.'
      WRITE(*,*)'***************************************'
C
        DO 24 J=1,N
  24      X(J)=XTMP(J)
        IFN = 0
        CALL VARME(FUN,N,X,FP,W,GRADTL,MAXFN,MESS,IER)
C
      WRITE(*,*)' '
      WRITE(*,*)' THE NUMBER OF FUNCTION EVALUATIONS IS ',IFN
      WRITE(*,*)' '
      WRITE(*,*)'THE MINIMUM FOUND IS',(X(I),I=1,N)
      WRITE(*,*)' '
      CALL FUN(N,X,FP)
      WRITE(*,*)'THE FUNCTION VALUE IS: ',FP
      STOP
      END
C----------------------------------------------------------------------
      SUBROUTINE FUN(NORD,BP,Q)
C----------------------------------------------------------------------
      COMMON /FUNERR/LER
      COMMON /TEST/IG,IFN
      DIMENSION BP(*)
      LOGICAL LER
      DOUBLE PRECISION BP,Q
      Q=100.*(BP(2)-BP(1)**2)**2+(BP(1)-1.)**2
      IFN = IFN + 1
      LER = .FALSE.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE VARME(FUN,NPAR,B,F0,W,NSIG,MAXFN,IOUT,IER)
C-----------------------------------------------------------------------
C
C     CALLING SUBROUTINE FOR SUBROUTINE VARMET
C
C     ALLOWS FOR SETUP OF DEFAULT PARAMETERS
C     AND EFFICIENT USE OF STORAGE
C     AS WELL AS WRITING OF ERROR MESSAGES
C
C    VERSION 0.1
C    CODED BY JOHN J. KOVAL
C    MARCH 1986
C
C    VERSION 0.2
C    CODED BY JOHN J. KOVAL
C    JULY 1988
C
C    VERSION 0.26
C    CODED BY MURRAY ALEXANDER FOR JOHN J. KOVAL
C    JULY 1989
C
C    VERSION 0.27
C    MODIFIED BY NAZIH HASSAN, JULY 1993
C
C    VERSION 0.28
C    MODIFIED BY JOHN KOVAL, JUNE 1996
C    BECAUSE OF COMMENTS FROM REVIEWER FOR APPLIED STATISTICS
C    CHANGES TO ORDER OF PARAMETERS IN GRAD
C
C    PARAMETERS              MEANING                       DEFAULT
C    ----------              -------                       -------
C
C     FUN         NAME OF FUNCTION TO BE MINIMIZED
C
C     NPAR        ORDER OF PARAMETER VECTOR
C                 (NUMBER OF UNKNOWNS)
C
C     B           ARRAY CONTAINING INITIAL ESTIMATES
C                 ON OUTPUT CONTAINING FINAL ESTIMATES
C
C     F0          VALUE OF FUNCTION AT MINIMUM
C
C     W           WORK ARRAY OF LENGTH FO (NPAR+5)*NPAR
C
C     NSIG        MACHINE ACCURACY AS NEGATIVE POWER       10 OR 5
C                 OF TEN
C
C     MAXFN       MAXIMUM NUMBER OF FUNCTION EVALUATIONS     1000
C                 (DOES INCLUDE EVALUATIONS BY SUBROUTINE
C                 GRAD WHICH CALCULATES APPROXIMATE GRADIENT)
C
C     IOUT        OUTPUT CHANNEL FOR ERROR MESSAGES            0
C                 (IF 0, THEN MESSAGES NOT WRITTEN)
C
C      IER        ERROR INDICATOR                              0
C                 INTEGER
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION W,B,F0,EPD,GRADTL
      EXTERNAL FUN
      DIMENSION B(*),W(*)
      PARAMETER (MAXF = 1000, EPD = 1.0D-05, MSIG = 10)
C
C     INITIALIZE
C
      IER=0
C
      IF(NSIG.EQ.0) NSIG = MSIG
      GRADTL = 1.0/(10.0**(NSIG))
C
      IF(GRADTL.LT.0.0) THEN
        IF(IOUT.GT.0) WRITE(IOUT,300) NSIG, GRADTL
 300    FORMAT(' NSIG VALUE OF ',I3,' CREATES NEGATIVE VALUE OF',
     1' GRADTL, NAMELY, ',G12.5)
        GRADTL = 1.0/(10**(MSIG))
        IF(IOUT.GT.0) WRITE(IOUT,310) MSIG, GRADTL
 310    FORMAT(' PROGRAM SUBSTITUTES NSIG VALUE OF ',I3,' WHICH',
     1' GIVES GRADTL VALUE OF ',G12.5)
      ENDIF
C
      IF(MAXFN.EQ.0) MAXFN = MAXF
C
C      NOW WE ARE READY TO CALL THE MINIMIZATION SUBROUTINE
C
      I1 = NPAR*NPAR + 1
      I2 = I1 + NPAR
      I3 = I2 + NPAR
      I4 = I3 + NPAR
C
      CALL VARMET(FUN,NPAR,B,F0,W(I3),W,W(I1),W(I2),W(I4),
     1 GRADTL,MAXFN,IER)
C
      IF(IER.GT.0.AND.IOUT.GT.0)THEN
        WRITE(IOUT,30) IER
   30   FORMAT(/' SUBROUTINE VARMET ERROR NUMBER ',I3)
        IF(IER.EQ.1) THEN
          WRITE(IOUT,40)
        ELSE IF(IER.EQ.2) THEN
          WRITE(IOUT,60)
        ELSE IF(IER.EQ.3) THEN
          WRITE(IOUT,70)
        ELSE IF(IER.EQ.4) THEN
          WRITE(IOUT,80)
        ENDIF
   40   FORMAT(' FUNCTION UNDEFINED AT INITIAL VALUE         ')
   60   FORMAT(' GRADIENT UNDEFINED IN TOO MANY DIMENSIONS   ')
   70   FORMAT(' FUNCTON NOT MINIMIZED BUT'
     1  /' UNABLE TO FIND MINIMUM IN DIRECTION OF SEARCH')
   80   FORMAT(' TOO MANY FUNCTION EVALUATIONS REQUIRED      ')
C
      ENDIF
C
  200 RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE VARMET(FUN,NPAR,B,F0,G,H,C,D,T,GRADTL,MAXFN,IFAULT)
C
C       ALGORITHM AS 319 APPL.STATIST. (1997), VOL.46, NO.4
C
C               VARIABLE METRIC FUNCTION MINIMISATION
C
      EXTERNAL FUN
      COMMON/FUNERR/LER
      DIMENSION B(NPAR),G(NPAR),H(NPAR,NPAR),C(NPAR),D(NPAR),T(2*NPAR)
      DOUBLE PRECISION B,F0,G,H,C,D,T,GRADTL,W,TOLER,D1,S,CK,F1,D2
      INTEGER IFN,IG
      LOGICAL LER
      PARAMETER (ICMAX=20, TOLER=0.00001, W=0.2)
C
      IG = 0
      IFN = 0
      LER = .FALSE.
      IFAULT = 0
      NP = NPAR + 1
C
      IF (MAXFN.EQ.0) MAXFN = 1000
      IF (GRADTL.EQ.0.0) GRADTL = 1.0D-10
C
      CALL FUN(NPAR,B,F0)
      IF(LER) THEN
         IFAULT = 1
         RETURN
      ENDIF
      IFN = IFN + 1
C
      CALL GRAD(FUN,NPAR,B,F0,G,T(NP),GRADTL,IFAULT)
      IF(IFAULT.GT.0) RETURN
C
      IG = IG + 1
      IFN = IFN + NPAR
      IF(IFN.GT.MAXFN) THEN
         IFAULT = 4
         RETURN
      ENDIF
C
   10 DO 30 K = 1,NPAR
         DO 20 L = 1,NPAR
            H(K,L) = 0.0
   20    CONTINUE
         H(K,K) = 1.00
   30 CONTINUE
      ILAST = IG
C
   40 DO 50 I = 1,NPAR
         D(I) = B(I)
         C(I) = G(I)
   50 CONTINUE
C
      D1 = 0.0
      DO 70 I = 1,NPAR
         S = 0.0
         DO 60 J = 1,NPAR
            S = S - H(I,J)*G(J)
   60    CONTINUE
         T(I) = S
         D1 = D1 - S*G(I)
   70 CONTINUE
C
      IF(D1.LE.0.0) THEN
         IF(ILAST.EQ.IG) THEN
            RETURN
         ENDIF
         GO TO 10
      ELSE
         CK = 1.0
         IC = 0
   90    ICOUNT = 0
         DO 100 I = 1,NPAR
            B(I) = D(I) + CK*T(I)
            IF(B(I).EQ.D(I)) THEN
               ICOUNT = ICOUNT + 1
            ENDIF
  100    CONTINUE
C
         IF(ICOUNT.GE.NPAR) THEN
            IF(ILAST.EQ.IG) THEN
               RETURN
            ENDIF
            GO TO 10
         ELSE
            CALL FUN(NPAR,B,F1)
C
            IFN = IFN + 1
            IF(IFN.GT.MAXFN) THEN
               IFAULT = 4
               RETURN
            ELSE IF(LER) THEN
               CK = W * CK
               IC = IC+1
               IF(IC.GT.ICMAX) THEN
                  IFAULT = 3
                  RETURN
               ENDIF
               GO TO 90
C
            ELSE IF(F1.GE.F0 - D1*CK*TOLER) THEN
               CK = W * CK
               GO TO 90
            ELSE
               F0 = F1
               CALL GRAD(FUN,NPAR,B,F0,G,T(NP),GRADTL,IFAULT)
               IF(IFAULT.GT.0) THEN
                  RETURN
               ENDIF
               IG = IG + 1
               IFN = IFN + NPAR
               IF(IFN.GT.MAXFN) THEN
                  IFAULT = 4
                  RETURN
               ENDIF
C
               D1 = 0.0
               DO 130 I = 1,NPAR
                  T(I) = CK*T(I)
                  C(I) = G(I) - C(I)
                  D1 = D1 + T(I)*C(I)
  130          CONTINUE
C
               IF(D1.LE.0.0) THEN
                  GOTO 10
               ENDIF
C
               D2 = 0.0
               DO 150 I = 1,NPAR
                  S = 0.0
                  DO 140 J = 1,NPAR
                     S = S + H(I,J)*C(J)
  140             CONTINUE
                  D(I) = S
                  D2 = D2 + S*C(I)
  150          CONTINUE
               D2 = 1.0 + D2/D1
C
               DO 170 I = 1,NPAR
                  DO 170 J = 1,NPAR
                     H(I,J) = H(I,J) - (T(I)*D(J) + D(I)*T(J) -
     1               D2*T(I)*T(J))/D1
  170          CONTINUE
            ENDIF
         ENDIF
      ENDIF
      GO TO 40
      END

      SUBROUTINE GRAD(F,NPAR,B,F0,G,SA,ER,IFAULT)
C
C     CALCULATE APPROXIMATE GRADIENT
C
      DIMENSION B(NPAR),G(NPAR),SA(NPAR)
      DOUBLE PRECISION B,F0,G,SA,ER,H,F1
      COMMON/FUNERR/LER
      LOGICAL LER
C
      JCMAX=NPAR-2
      JC = 0
C
      DO 20 I = 1,NPAR
         H =(DABS(B(I)) +DSQRT(ER)) *DSQRT(ER)
         SA(I) = B(I)
         B(I) = B(I) + H
         CALL F(NPAR,B,F1)
         B(I) = SA(I)
C
         IF(LER) THEN
            F1 = F0 + H
            JC = JC + 1
         ENDIF
C
         G(I) = (F1 -F0)/H
   20 CONTINUE
C
      IF(JC.GT.JCMAX) IFAULT = 2
      RETURN
      END