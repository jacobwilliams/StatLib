C UKC NETLIB DISTRIBUTION COPYRIGHT 1990 RSS
C
      FUNCTION CHISQN(X, DF, FL, IFAULT)
C
C<<<<<  Acquired in machine-readable form from 'Applied Statistics'
C<<<<<  algorithms editor, January 1983.
C
C
C        ALGORITHM AS 170  APPL. STATIST. (1981) VOL.30, NO.3
C
C        The non-central chi-squared distribution.
C
C     Auxiliary routines required: GAMMDS = AS147, ALOGAM = CACM 291.
C     See AS245 for an alternative to ALOGAM.
C
      CHISQN = 0.0
      IFAULT = 0
C
C        TEST FOR ADMISSIBILITY OF ARGUMENTS
C
      IF (DF.LE.0.0) IFAULT = 1
      IF (X.LT.0.0) IFAULT = 2
      IF (FL.LT.0.0) IFAULT = 3
      IF (IFAULT.GT.0.OR.X.EQ.0.0) RETURN
C
      DF2 = 0.5*DF
      X2 = 0.5*X
      FXP = GAMMDS(X2,DF2,IFAULT)
      CHISQN = CHI(X2,DF2,FL,FXP)
      RETURN
      END
C
      SUBROUTINE CHISQL(X, DF, FX, FL, IFAULT)
C
C        ALGORITHM AS 170.1  APPL.STATIST. (1981) VOL.30, NO.3
C
C        DEFINE ACCURACY AND INITIALIZE
C
C        N SHOULD BE SPECIFIED SUCH THAT ACC IS GREATER THAN
C        OR EQUAL TO (AU-AL)/2**N
C
      PARAMETER (ACC = 1.0E-6, N = 30)
C
      AL = 0.0
      AINC = 80.0
      AU = 80.0
C
      IFAULT = 0
C
C        TEST FOR ADMISSIBILITY OF ARGUMENTS
C
      IF (DF.LE.0.0) IFAULT = 1
      IF (X.LT.0.0) IFAULT = 2
      IF (FX.LE.0.0) IFAULT = 3
      IF (IFAULT.GT.0) GO TO 4
C
      DF2 = 0.5*DF
      X2 = 0.5*X
      FX1 = GAMMDS(X2,DF2,IFAULT)
    1 APROX = CHI(X2,DF2,AU,FX1)
      IF (FX.GT.APROX) GOTO 2
      IF (FX.LT.APROX) AL = AU
      AU = AU+AINC
      GO TO 1
    2 DO 3 J = 1,N
        FL = 0.5*(AL+AU)
        APROX = CHI(X2, DF2, FL, FX1)
        IF (ABS(FX-APROX).LT.ACC) GO TO 4
        IF (FX.LT.APROX) AL = FL
        IF (FX.GE.APROX) AU = FL
    3 CONTINUE
    4 RETURN
      END
C
      FUNCTION CHI(X, DF, FL, FXC)
C
C        ALGORITHM AS 170.2  APPL. STATIST. (1981) VOL.30, NO.3
C
      PARAMETER (ACC2 = 1.0E-8)
C
      CHI = FXC
      DF1 = DF
      FL2 = 0.5*FL
      C = 1.0
      T = 0.0
    1 T = T+1.0
      C = C*FL2/T
      DF1 = DF1+1.0
      TERM = C*GAMMDS(X, DF1, IFAULT)
      CHI = CHI+TERM
      IF (TERM.GE.ACC2) GO TO 1
      CHI = CHI*EXP(-FL2)
      RETURN
      END
