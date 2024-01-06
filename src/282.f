C
C THE FOLLOWING ROUTINES ARE WRITTEN IN DOUBLE PRECISION. IT IS STRONGLY
C RECOMMENDED THAT THEY BE RUN THIS WAY. IF IT IS DESIRED TO RUN THEM IN
C SINGLE PRECISION, ALL "IMPLICIT DOUBLE PRECISION (A-H,O-Z)" LINES SHOULD
C BE DELETED, ALL DOUBLE PRECISION CONSTANTS IN DATA STATEMENTS SHOULD BE
C CONVERTED TO REAL CONSTANTS, AND THE VALUE OF "TOLER" SET IN THE DATA
C STATEMENT IN ROUTINE MVELMS SHOULD BE CHANGED TO "1.E-5"
C
      SUBROUTINE MVELMS (X,Y,NCAS,NPRE,IOPTN,MAXTRY,NCLOW,NCHIGH,
     1 COEFFS,EPRMIN,RESID,ROBDSQ,CVEMIN,DATA,IWORK,WORK,NVDIM,
     2 NODIM,IFAULT)
C
C     ALGORITHM AS 282 APPL.STATIST. (1993), VOL.42, NO.2
C
C     HIGH BREAKDOWN REGRESSION AND MULTIVARIATE ESTIMATION
C
C
C ROUTINE TO CALCULATE LEAST MEDIAN OF SQUARES REGRESSION, MINIMUM VOLUME
C ELLIPSOID, AND ASSOCIATED STATISTICS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NODIM,NVDIM), Y(NODIM), COEFFS(NVDIM,*),EPRMIN(*),
     1 RESID(NODIM,*), CVEMIN(*), ROBDSQ(NODIM,*), DATA(0:NVDIM,*),
     2 IWORK(*), WORK(*)
C
C     OPTIONAL ADDITIONAL ARGUMENTS FOR MORE SPECIALIZED USE
C
C     DIMENSION XINV(NVDIM,NVDIM,*),ELSCAL(*),LMSBAS(NVDIM,*),
C    1 MVEBAS(NVDIM,*)
C
      LOGICAL LMS,MVE,ISINT,EXH,EXACT
      DATA TOLER/1.D-9/,BIG/1.D30/,NFRESH/50/,ZERO/0.D0/,HALF/0.5D0/,
     1 ONE/1.D0/,TEN/10.D0/
C
C     EXTRACT THE OPTIONS REQUESTED
C
      ITEMP = IOPTN
      LMS = MOD(ITEMP,2) .EQ. 0
      ITEMP = ITEMP / 2
      MVE = MOD(ITEMP,2) .EQ. 0
      ITEMP = ITEMP / 2
      ISINT = MOD(ITEMP,2) .EQ. 0
      ITEMP = ITEMP / 2
      EXH = MOD(ITEMP,2) .EQ. 0
      IF (ISINT) THEN
        NVAR=NPRE+1
        ADDCON=ONE/FLOAT(NVAR)
        RPP1=SQRT(ADDCON)
        INT=1
      ELSE
        NVAR=NPRE
        ADDCON=ZERO
        RPP1=ONE
        INT=0
      END IF
C
C CHECK FAULT PARAMETERS
C
      IFAULT=0
      IF ((NCLOW.LT.NVAR).OR.(NCLOW.GT.NCHIGH).OR.(NCHIGH.
     1  GT.NCAS)) IFAULT=2
      IF (NCAS.LT.NVAR) IFAULT=1
      IF ((ISINT).AND.(NVAR.EQ.1).AND.(MVE)) IFAULT=IFAULT+4
      IF (IFAULT.GT.0) GO TO 490
      EXACT=.FALSE.
      J2=NCAS
      J3=2*NCAS
      J4=3*NCAS
      I2=NVAR
      I3=2*NVAR
      I4=3*NVAR
      I5=I4+NCAS
      NVR2=2*NVAR
      DETER=RPP1
      POWMED=FLOAT(NPRE)*HALF
      DO 10 I=1,NVAR
 10   IWORK(I)=NCAS+I-NVAR
      IPTR=1
      NCRANG=NCHIGH-NCLOW+1
      DO 20 I=1,NCRANG
      CVEMIN(I)=BIG
      EPRMIN(I)=BIG
 20   CONTINUE
      LFRESH=0
      NCOUNT=0
C
C TARGET VARIABLE, IF PRESENT, IS TREATED IN STANDARDIZED FORM SO THAT EXACT
C FIT CAN BE DETECTED MORE EASILY
C
      IF (LMS) THEN
        IXLO=0
        DO 30 I=1,NCAS
 30     WORK(I)=Y(I)
C
C THE USER-SUPPLIED SUBROUTINE SORTSUB(RA,N,ILOW) SORTS ENTRIES ILOW+1 THROUGH
C ILOW+N IN ASCENDING ORDER FROM A VECTOR RA
C
        CALL SORTSUB(WORK,NCAS,0)
        YMED=(WORK(NCAS/2+1)+WORK((NCAS+1)/2))*HALF
        DO 40 I=1,NCAS
 40     WORK(I)=ABS(Y(I)-YMED)
        CALL SORTSUB(WORK,NCAS,0)
        YMAD=(WORK(NCAS/2+1)+WORK((NCAS+1)/2))*HALF
        IF (YMAD.EQ.ZERO) YMAD=ONE
      ELSE
        IXLO=1
      END IF
      DETADJ=ZERO
      DO 60 J=1,NPRE
      DO 50 I=1,NCAS
 50   WORK(I)=ABS(X(I,J))
      CALL SORTSUB(WORK,NCAS,0)
      WORK(J4+J+INT)=(WORK(NCAS/2+1)+WORK((NCAS+1)/2))*HALF
      IF (WORK(J4+J+INT).EQ.ZERO) WORK(J4+J+INT)=WORK(NCAS)
      DETADJ=DETADJ+LOG10(WORK(J4+J+INT))
 60   CONTINUE
      IF (ISINT) WORK(J4+1)=ONE
C
C DATA IS TRANSFERRED TO WORKAREA; INITIAL SIMPLEX TABLEAU IS SET UP
C
      CALL REFRESH(DATA,NODIM,NVDIM,IXLO,X,NCAS,NPRE,Y,YMAD,NVAR,WORK,
     1  J4,IWORK,I2,I3,I5,0,ISINT,INT,LMS,DETER,LFRESH,TOLER)
C
C INITIAL BASIS IS SET UP. CHECKS ARE MADE THAT INITIAL BASIS MEMBERS ARE IN
C GENERAL POSITION
C
      IF (EXH) THEN
        IWORK(I2+1)=1
        DO 100 I=1,NVAR
 70     J=IWORK(I2+I)
        DO 80 II=I,NVAR
        IF (ABS(DATA(II,J)).GE.TOLER) THEN
          CALL PIVOT(DATA,IXLO,NVAR,NCAS,IWORK,I5,II,J,DETER,
     1      NVDIM)
          IWORK(I3+I)=J
          IF (II.NE.I) CALL SWAP(DATA,NVAR,NCAS,II,I,NVDIM)
          GO TO 90
        END IF
 80     CONTINUE
        IF (J.EQ.NCAS) THEN
          IFAULT=IFAULT+8
          IF (IFAULT.GE.24) IFAULT=IFAULT-16
          GO TO 490
        ELSE
          IFAULT=16
          IWORK(I2+I)=IWORK(I2+I)+1
          GO TO 70
        END IF
 90     IF (I.LT.NVAR) IWORK(I2+I+1)=IWORK(I2+I)+1
100     CONTINUE
        IPTR=NVAR
        GO TO 230
      ELSE
        DO 110 I=1,NCAS
110     IWORK(I4+I)=I
        CALL PERM(IWORK,I4,NCAS,0)
        INXPTR=0
        DO 140 I=1,NVAR
120     INXPTR=INXPTR+1
        J=IWORK(I4+INXPTR)
        DO 130 II=I,NVAR
        IF (ABS(DATA(II,J)).GE.TOLER) THEN
          CALL PIVOT(DATA,IXLO,NVAR,NCAS,IWORK,I5,II,J,DETER,
     1      NVDIM)
          IWORK(I3+I)=J
          IF (II.NE.I) CALL SWAP(DATA,NVAR,NCAS,II,I,NVDIM)
          IWORK(I2+I)=J
          GO TO 140
        END IF
130     CONTINUE
        IF (INXPTR.EQ.NCAS) THEN
          IFAULT=IFAULT+8
          IF (IFAULT.GE.24) IFAULT=IFAULT-16
          GO TO 490
        ELSE
          IFAULT=16
          GO TO 120
        END IF
140     CONTINUE
        GO TO 230
      END IF

C
C MAIN ANALYSIS LOOP. GENERATE ALL SUBSETS (IF EXH) OR A SUBSET (IF
C NOT EXH)
C
150   IF (EXH) THEN
C
C IF EXHAUSTIVE, SUCCESSIVE BASES ARE CONSIDERED
C
160     IWORK(I2+IPTR)=IWORK(I2+IPTR)+1
        IF (IWORK(I2+IPTR).GT.IWORK(IPTR)) THEN
          IPTR=IPTR-1
          IF (IPTR.EQ.0) THEN
            GO TO 440
          ELSE
            GO TO 160
          END IF
        ELSE
C
C IF DATA IS NOT IN GENERAL POSITION, NEXT POSITION IN LIST MUST BE FOUND
C
          IF (ABS(DATA(IPTR,IWORK(I2+IPTR))).GE.TOLER) GO TO 210
          IF (LFRESH.GT.NVR2) THEN
            DETER=RPP1
            CALL REFRESH(DATA,NODIM,NVDIM,IXLO,X,NCAS,NPRE,Y,YMAD,
     1        NVAR,WORK,J4,IWORK,I2,I3,I5,IPTR-1,ISINT,INT,LMS,DETER,
     2        LFRESH,TOLER)
            IF (ABS(DATA(IPTR,IWORK(I2+IPTR))).GE.TOLER) GO TO 210
            IFAULT=16
            IF (IPTR.EQ.NVAR) THEN
              GO TO 160
            ELSE
              GO TO 170
            END IF
          END IF
          IFAULT=16
          IF (IPTR.EQ.NVAR) GO TO 160
          DETER=RPP1
          CALL REFRESH(DATA,NODIM,NVDIM,IXLO,X,NCAS,NPRE,Y,YMAD,
     1      NVAR,WORK,J4,IWORK,I2,I3,I5,IPTR-1,ISINT,INT,LMS,DETER,
     2      LFRESH,TOLER)
170       J=IWORK(I2+IPTR)
          DO 180 II=IPTR,NVAR
          IF (ABS(DATA(II,J)).GE.TOLER) THEN
            CALL PIVOT(DATA,IXLO,NVAR,NCAS,IWORK,I5,II,J,DETER,
     1        NVDIM)
            IWORK(I3+IPTR)=J
            IF (II.NE.IPTR) CALL SWAP(DATA,NVAR,NCAS,II,IPTR,NVDIM)
            GO TO 220
          END IF
180       CONTINUE
          IWORK(I2+IPTR)=IWORK(I2+IPTR)+1
          IF (IWORK(I2+IPTR).GT.IWORK(IPTR)) THEN
            IPTR=IPTR-1
            IF (IPTR.EQ.0) THEN
              GO TO 440
            ELSE
              IWORK(I2+IPTR)=IWORK(I2+IPTR)+1
            END IF
          END IF
          GO TO 170
        END IF
      ELSE
C
C IF NOT EXHAUSTIVE, THE NEXT ENTRY IN RANDOM PERMUTATION VECTOR IS ENTERED
C INTO THE BASIS
C
        IF (NCOUNT.GT.MAXTRY) GO TO 440
        IPTR=MOD(IPTR,NVAR)+1
190     INXPTR=INXPTR+1
        IF (INXPTR.GT.NCAS) THEN
C
C IF COME TO THE END OF THE PERMUTATION VECTOR, GENERATE A NEW ONE
C
          DO 200 I=1,NVAR
          IWORK(I4+NCAS-NVAR+I)=IWORK(I4+I)
          IWORK(I4+I)=IWORK(I2+I)
200       CONTINUE
          CALL PERM(IWORK,I4,NCAS,NVAR)
          INXPTR=NVAR+1
        END IF
        NCOL=IWORK(I4+INXPTR)
        IF (ABS(DATA(IPTR,NCOL)).LT.TOLER) THEN
          IF (LFRESH.GT.NVR2) THEN
            DETER=RPP1
            CALL REFRESH(DATA,NODIM,NVDIM,IXLO,X,NCAS,NPRE,Y,YMAD,
     1        NVAR,WORK,J4,IWORK,I2,I3,I5,IPTR-1,ISINT,INT,LMS,DETER,
     2        LFRESH,TOLER)
            IF (ABS(DATA(IPTR,NCOL)).LT.TOLER) THEN
              IFAULT=16
              GO TO 190
            END IF
          ELSE
            IFAULT=16
            GO TO 190
          END IF
        END IF
        IWORK(I2+IPTR)=NCOL
      END IF
210   NCOUNT=NCOUNT+1
      LFRESH=LFRESH+1
C
C CARRY OUT NEXT PIVOT, THEREBY CREATING A NEW BASIS
C
      CALL PIVOT (DATA,IXLO,NVAR,NCAS,IWORK,I3,IPTR,
     1  IWORK(I2+IPTR),DETER,NVDIM)
C
C RECOMPUTE THE INVERSE BASIS EVERY NFRESH SIMPLEX PIVOTS
C
      IF (LFRESH.GT.NFRESH) THEN
        DETER=RPP1
        CALL REFRESH(DATA,NODIM,NVDIM,IXLO,X,NCAS,NPRE,Y,YMAD,
     1    NVAR,WORK,J4,IWORK,I2,I3,I5,NVAR,ISINT,INT,LMS,DETER,LFRESH,
     2    TOLER)
      END IF
220   IF (EXH.AND.(IPTR.LT.NVAR)) THEN
        IPTR=IPTR+1
        IWORK(I2+IPTR)=IWORK(I2+IPTR-1)
        GO TO 150
      END IF
230   IF (LMS) THEN
C
C CHECK TO SEE IF THIS BASIS HAS A SMALLER LMS CRITERION VALUE THAN PREVIOUS
C BASES. IF THERE IS NO INTERCEPT, A PRELIMINARY COMPARISON IS MADE TO SEE IF
C THIS BASIS COULD BE OPTIMAL.
C
        IF (EXACT) GO TO 340
        IF (ISINT) THEN
          DO 240 J=1,NCAS
240       WORK(J)=DATA(0,J)
        ELSE
          DO 250 J=1,NCAS
250       WORK(J)=ABS(DATA(0,J))
          DO 270 KX=NCLOW,NCHIGH
          KXSTO=KX-NCLOW+1
          TARGET=EPRMIN(KXSTO)
          ICOUNT=0
          DO 260 J=1,NCAS
260       IF (WORK(J).LE.TARGET) ICOUNT=ICOUNT+1
          IF (ICOUNT.GE.KX) GO TO 280
270       CONTINUE
          GO TO 340
        END IF
280     CALL SORTSUB(WORK,NCAS,0)
        DO 330 KX=NCLOW,NCHIGH
        KXSTO=KX-NCLOW+1
        OFFS=ZERO
        IF (ISINT) THEN
C
C CALCULATE PROPER OFFSET FOR INTERCEPT TERM, IF REQUIRED, AND DETERMINE
C CRITERION FOR THIS BASIS
C
          SHORT=BIG
          DO 290 KKX=KX,NCAS
          IF (WORK(KKX)-WORK(KKX-KX+1).LT.SHORT) THEN
            SHORT=WORK(KKX)-WORK(KKX-KX+1)
            CRITER=SHORT*HALF
            OFFS=(WORK(KKX)+WORK(KKX-KX+1))*HALF
          END IF
290       CONTINUE
        ELSE
          CRITER=WORK(KX)
        END IF
C
C UPDATE SUMMARY STATISTICS IF THIS BASIS IS CURRENTLY OPTIMAL
C
        IF (CRITER.LT.EPRMIN(KXSTO)) THEN
          EPRMIN(KXSTO)=CRITER
C         DO 300 J=1,NVAR
C300      LMSBAS(J,KXSTO)=IWORK(I3+J)
          DO 310 J=1,NVAR
310       COEFFS(J,KXSTO)=-DATA(0,NCAS+J)*YMAD
          IF (ISINT) COEFFS(1,KXSTO)=COEFFS(1,KXSTO)+OFFS*YMAD
          DO 320 J=1,NCAS
320       RESID(J,KXSTO)=(DATA(0,J)-OFFS)*YMAD
C
C A CHECK FOR EXACT FIT IS MADE
C
          IF ((CRITER.LT.TOLER).AND.(KX.EQ.NCHIGH)) THEN
            EXACT=.TRUE.
            IF (.NOT.MVE) GO TO 460
          END IF
        END IF
330     CONTINUE
C
C IF THIS IS A UNIVARIATE LOCATION PROBLEM, NO FURTHER ANALYSIS IS NEEDED
C
        IF ((ISINT).AND.(NVAR.EQ.1)) GO TO 460
      END IF
340   IF (MVE) THEN
C
C CHECK TO SEE IF THIS BASIS HAS A SMALLER MVE CRITERION VALUE THAN PREVIOUS
C BASES. A PRELIMINARY COMPARISON IS MADE TO SEE IF THIS BASIS COULD BE
C OPTIMAL.
C
        BASVOL=LOG10(ABS(DETER))
        DO 360 J=1,NCAS
        WORK(J+J2)=-ADDCON
        DO 350 I=1,NVAR
350     WORK(J+J2)=WORK(J+J2)+DATA(I,J)**2
        IF (WORK(J+J2).GT.ZERO) THEN
          WORK(J+J2)=LOG10(WORK(J+J2))
        ELSE
          WORK(J+J2)=-BIG
        END IF
        WORK(J+J3)=WORK(J+J2)
360     CONTINUE
        DO 380 KX=NCLOW,NCHIGH
        KXSTO=KX-NCLOW+1
        TARGET=(CVEMIN(KXSTO)-BASVOL)/POWMED
        ICOUNT=0
        DO 370 J=1,NCAS
370     IF (WORK(J+J3).LE.TARGET) ICOUNT=ICOUNT+1
        IF (ICOUNT.GE.KX) GO TO 390
380     CONTINUE
        GO TO 150
390     CALL SORTSUB(WORK,NCAS,J3)
        DO 430 KX=NCLOW,NCHIGH
        KXSTO=KX-NCLOW+1
        VOLU=BASVOL+POWMED*WORK(KX+J3)
        IF (VOLU.LT.CVEMIN(KXSTO)) THEN
C
C UPDATE SUMMARY STATISTICS IF THIS BASIS IS CURRENTLY OPTIMAL
C
          CVEMIN(KXSTO)=VOLU
C         ELSCAL(KXSTO)=WORK(KX+J3)
C         DO 400 J=1,NVAR
C400      MVEBAS(J,KXSTO)=IWORK(I3+J)
          DO 410 JJ=1,NCAS
410       ROBDSQ(JJ,KXSTO)=WORK(JJ+J2)-WORK(KX+J3)
C         DO 420 J=1,NVAR
C         DO 420 JJ=1,NVAR
C420      XINV(JJ,J,KXSTO)=DATA(J,JJ+NCAS)/WORK(J4+JJ)
        END IF
430     CONTINUE
      END IF
C
C END OF MAIN PROCESSING LOOP. GO TO BEGINNING AND GENERATE ANOTHER BASIS
C
      GO TO 150
440   IF (MVE) THEN
C
C CONVERT ROBUST DISTANCES AND SCALING BACK TO ORIGINAL SCALE
C
        DO 450 K=1,NCRANG
C       ELSCAL(K)=TEN**ELSCAL(K)
        CVEMIN(K)=CVEMIN(K)+DETADJ
        DO 450 J=1,NCAS
        IF (ROBDSQ(J,K).GT.(-BIG)) THEN
          ROBDSQ(J,K)=TEN**ROBDSQ(J,K)
        ELSE
          ROBDSQ(J,K)=ZERO
        END IF
450     CONTINUE
      END IF
460   IF (LMS) THEN
        DO 480 K=1,NCRANG
        EPRMIN(K)=EPRMIN(K)*YMAD
        DO 470 I=1,NPRE
470     COEFFS(I+INT,K)=COEFFS(I+INT,K)/WORK(I+J4+INT)
480     CONTINUE
      END IF
490   RETURN
      END
      SUBROUTINE REFRESH(DATA,NODIM,NVDIM,IXLO,X,NCAS,NPRE,Y,YMAD,
     1 NVAR,WORK,J4,IWORK,I2,I3,I5,IUP,ISINT,INT,LMS,DETER,LFRESH,TOLER)
C
C SUBROUTINE TO REFRESH FIRST IUP ENTRIES OF SIMPLEX BASIS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DATA(0:NVDIM,*), X(NODIM,NVDIM), Y(NODIM),
     1  WORK(*),IWORK(*)
      LOGICAL ISINT,LMS
      DATA ZERO/0.D0/,ONE/1.D0/
      LFRESH=0
      IF (ISINT) THEN
        DO 10 I=1,NCAS
10      DATA(1,I)=ONE
      END IF
      DO 20 I=1,NCAS
      DO 20 J=1,NPRE
20    DATA(J+INT,I)=X(I,J)/WORK(J4+J+INT)
      DO 40 I=1,NVAR
      DO 30 J=1,NVAR
30    DATA(J,NCAS+I)=ZERO
      DATA(I,NCAS+I)=ONE
40    CONTINUE
      IF (LMS) THEN
        DO 50 I=1,NCAS
50      DATA(0,I)=Y(I)/YMAD
        DO 60 I=1,NVAR
60      DATA(0,NCAS+I)=ZERO
      END IF
      DO 80 I=1,IUP
      J=IWORK(I2+I)
      DO 70 II=I,NVAR
      IF (ABS(DATA(II,J)).GE.TOLER) THEN
        CALL PIVOT(DATA,IXLO,NVAR,NCAS,IWORK,I5,II,J,DETER,NVDIM)
        IF (II.NE.I) CALL SWAP(DATA,NVAR,NCAS,II,I,NVDIM)
        IWORK(I3+I)=J
        GO TO 80
      END IF
70    CONTINUE
80    CONTINUE
      RETURN
      END
      SUBROUTINE PIVOT(X,IXLO,NORD,NCAS,IWORK,I3,NROW,NCOL,DETER,
     1  NVDIM)
C
C SUBROUTINE TO PIVOT ENTRY CORRESPONDING TO (NROW,NCOL) INTO SIMPLEX TABLEAU
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:NVDIM,*),IWORK(*)
      DATA ZERO/0.D0/,ONE/1.D0/
      PIVT=X(NROW,NCOL)
      DETER=DETER*PIVT
      NHIGH=NORD+NCAS
      DO 20 J=1,NHIGH
      IF (J.NE.NCOL) THEN
        FMULT=X(NROW,J)/PIVT
        DO 10 I=IXLO,NORD
 10     IF (I.NE.NROW) X(I,J)=X(I,J)-FMULT*X(I,NCOL)
        X(NROW,J)=FMULT
      END IF
 20   CONTINUE
      DO 30 I=IXLO,NORD
 30   X(I,NCOL)=ZERO
      X(NROW,NCOL)=ONE
      IWORK(I3+NROW)=NCOL
      RETURN
      END
      SUBROUTINE PERM(INDEX,IAA,N,IAB)
C
C SUBROUTINE TO RETURN PSEUDORANDOM PERMUTATION OF N - IAB ENTRIES IN
C VECTOR INDEX STARTING AT IAA+1
C
      DIMENSION INDEX(*)
      L=0
      M=N-IAB
C
C GENERATE A RANDOM DIGIT FROM 1 TO M USING THE PSEUDORANDOM U(0,1) VARIATE
C RANDOM() (SUCH AS FROM ALGORITHM AS 183)
C
 10   J=INT(RANDOM()*FLOAT(M))+1
C
C SWAP ENTRIES IN INDEX CORRESPONDING TO J AND M, OFFSET BY IAB
C
      IF (J.NE.M) THEN
        ITEMP=INDEX(IAA+J+IAB)
        INDEX(IAA+J+IAB)=INDEX(IAA+M+IAB)
        INDEX(IAA+M+IAB)=ITEMP
      END IF
      M=M-1
      IF (M.GT.1) GO TO 10
      RETURN
      END
      SUBROUTINE SWAP(DATA,NVAR,NCAS,IR1,IR2,NVDIM)
C
C SUBROUTINE TO SWAP ROWS IR1 AND IR2 OF MATRIX DATA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DATA(0:NVDIM,*)
      NHIGH=NVAR+NCAS
      DO 10 I=1,NHIGH
      TEMP=DATA(IR1,I)
      DATA(IR1,I)=DATA(IR2,I)
      DATA(IR2,I)=TEMP
10    CONTINUE
      RETURN
      END