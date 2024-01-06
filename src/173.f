      SUBROUTINE DESMAT(MAXROW,MAXCOL,NDIM,NSUB,MXDINT,MSUB,
     1LSUB,ISUB,IBEG,MATDES,IFAULT)
C
C<<<<<  Acquired in machine-readable form from 'Applied Statistics'
C<<<<<  algorithms editor, January 1983.
C
C        ALGORITHM AS 173  APPL. STATIST. (1982) VOL.31, NO.1
C
C        GENERATES A DESIGN MATRIX FOR BALANCED FACTORIAL
C        EXPERIMENTS
C
      DIMENSION MATDES(MAXROW,MAXCOL),NSUB(NDIM),MSUB(NDIM),
     1          LSUB(NDIM),ISUB(NDIM),IBEG(NDIM)
C
C     PRIMARY INPUT PARAMETER CHECK
C
      IFAULT = 0
      DO 5 I = 1,NDIM
      IF(NSUB(I) .GE.2) GO TO 5
      IFAULT = 1
      RETURN
    5 CONTINUE
      IF(MXDINT .LE.NDIM) GO TO 8
      IFAULT = 2
      RETURN
C
C     CALCULATE PRODUCT VECTOR MSUB FOR M.E. GENERATION
C     CALCULATE PRODUCT VECTOR LSUB FOR INTRT. GENERATION
C
    8 IBEG(1) = 1
      MSUB(1) = NSUB(1)
      LSUB(1) = NSUB(1)-1
      DO 40 I = 2,NDIM
      IK  = I-1
      IBEG(I) = IBEG(IK)+NSUB(IK)-1
      MSUB(I) = MSUB(IK) * NSUB(I)
      LSUB(I) = LSUB(IK) *(NSUB(I)-1)
   40 CONTINUE
C
C     RESERVE THE TOTAL NUMBER OF OBSERVATIONS
C
      NOBS = MSUB(NDIM)
C
C     SECONDARY INPUT PARAMETER CHECK
C
      IF(NOBS .LE. MAXROW) GO TO 50
      IFAULT = 3
      RETURN
C
C     GENERATE THE MAIN EFFECTS
C
   50 CALL MANEFF(MAXROW,MAXCOL,NDIM,NSUB,MSUB,ISUB,MATDES,JCOL)
C
C     GENERATE ALL POSSIBLE INTERACTIONS
C
      IF(MXDINT.GE.2)CALL INTER(MAXROW,MAXCOL,NDIM,NOBS,NSUB,MXDINT,IBEG
     1,MSUB,LSUB,ISUB,MATDES,JCOL,IFAULT)
C
C        PASS BACK NOBS AND JCOL IRRESPECTIVE OF VALUE OF IFAULT
C
      IBEG(1) = NOBS
      IBEG(2) = JCOL
      RETURN
      END
      SUBROUTINE MANEFF(MAXROW,MAXCOL,NDIM,NSUB,MSUB,ISUB,MATDES,
     1JCOL)
C
C        ALGORITHM AS 173.1  APPL. STATIST. (1982) VOL.31, NO.1
C
C       COMPUTES THOSE COLUMNS OF THE DESIGN MATRIX (MATDES) ASSOCIATED
C       WITH THE MEAN AND THE MAIN EFFECTS
C
      DIMENSION MATDES(MAXROW,MAXCOL),NSUB(NDIM),MSUB(NDIM),ISUB(NDIM)
C
C        SET THE TOTAL NUMBER OF OBSERVATIONS
C
      NOBS = MSUB(NDIM)
C
C     SET UP MUU
C
      DO 10 I = 1,NOBS
      MATDES(I,1) = 1
   10 CONTINUE
C
C     SET UP OTHER MAIN EFFECTS
C
      DO 40 I = 1,NOBS
      CALL SIMDO(.TRUE.,.TRUE.,MSUB,NDIM,I,ISUB,JFALT)
C
C        ROUTINE PROTECTED FROM FAILURE SO FAULT PARAMETER UNTESTED
C
      JCOL = 1
      DO 30 J = 1,NDIM
      IK = NSUB(J)-1
      DO 20 K = 1,IK
      JCOL = JCOL+1
      IF (ISUB(J) .EQ. K) GO TO 15
      IF (ISUB(J) .EQ. NSUB(J))GO TO 18
      MATDES(I,JCOL) = 0
      GO TO 20
   15 MATDES(I,JCOL) = 1
      GO TO 20
   18 MATDES(I,JCOL) = -1
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
      SUBROUTINE INTER(MAXROW,MAXCOL,NDIM,NOBS,NSUB,MXDINT,IBEG,
     1ICOMB,LSUB,IREF,MATDES,JCOL,IFAULT)
C
C        ALGORITHM AS 173.2  APPL. STATIST. (1982) VOL.31, NO.1
C
C        GENERATES THE SPECIFIED INTERACTION TERMS
C
      DIMENSION MATDES(MAXROW,MAXCOL),NSUB(NDIM),IBEG(NDIM),ICOMB(NDIM)
     1,LSUB(NDIM),IREF(NDIM)
      LOGICAL QNXT
C
C     GENERATE ALL POSSIBLE COMBINATIONS OF DIMENSION .LE. MXDINT
C
      QNXT = .FALSE.
      DO 50 ICMB = 2,MXDINT
    5 CALL COMB(NDIM,ICMB,ICOMB,QNXT)
C
C     HERE TO PROCESS CURRENT COMBINATION OF DIMENSION ICMB
C
C     CALCULATE THE ADDITIONAL COLUMNS REQUIRED TO STORE THIS
C     INTERACTION
C
      NPROD = 1
      DO 10 I = 1,ICMB
      IK = ICOMB(I)
      NPROD = NPROD*(NSUB(IK)-1)
   10 CONTINUE
C
C     CHECK IF SUFFICIENT COLUMNS ARE AVAILABLE IN MATDES
C
      IF (JCOL+NPROD .LE. MAXCOL) GO TO 15
      IFAULT = 4
      RETURN
C
C     PROCESS THE CURRENT INTERACTION FOR WHICH THERE IS SPACE
C
   15 DO 40 J = 1,NPROD
      CALL SIMDO(.TRUE.,.TRUE.,LSUB,NDIM,J,IREF,JFALT)
C
C        ROUTINE PROTECTED FROM FAILURE SO FAULT PARAMETER UNTESTED
C
      JCOL = JCOL+1
      DO 30 I = 1,NOBS
      MATDES(I,JCOL) = 1
      DO 20 L = 1,ICMB
      IK = ICOMB(L)
      JREF = IBEG(IK) + IREF(L)
      MATDES(I,JCOL) = MATDES(I,JCOL)*MATDES(I,JREF)
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
C
C     FINISHED ADDING IN THIS INTERACTION - PROCESS NEXT ONE
C
      IF(QNXT)GO TO 5
   50 CONTINUE
      RETURN
      END
      SUBROUTINE COMB(N,K,A,MTC)
C
C        ALGORITHM AS 173.3  APPL. STATIST. (1982) VOL.31, NO.1
C
C        GENERATES ALL POSSIBLE COMBINATIONS OF DIMENSION K
C        FROM N INTEGERS
C
      INTEGER A(K),H
      LOGICAL MTC
C
      IF(.NOT.MTC) GO TO 20
      DO 10 H = 1,K
      M = K + 1 - H
      M2 = A(M)
      IF(M2.NE.N+1-H)GO TO 30
   10 CONTINUE
   20 M2 = 0
      H = K
   30 DO 40 J = 1,H
      M = K + J - H
      A(M) = M2 + J
   40 CONTINUE
      MTC = A(1).NE.N-K+1
      RETURN
      END