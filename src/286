C                   EVM PROGRAMS WITH EXAMPLES
C
C
C CONTENTS
C
C This file gives the full computer code, in FORTRAN 77, for the algorithm
C described in "Parameter Estimation in the Error-in- Variables Model", by
C Park M. Reilly, S.E. Keeler, and H.V. Reilly.
C
C It is given in two forms, first for the case where there is any number of
C functional equations in the model, and then for the special case where there
C is only one.  The special case covers a large fraction of practical
C applications and the applicable form of the algorithm is not given in detail
C in the paper.  The special form executes substantially more quickly than the
C more general one.  Also given are the complete input and output files for an
C example using each of the two forms of the algorithm.  The material is
C organised as follows:
C
C           1. EVM.FOR     The algorithm in general form
C                          for any number of functional equations.
C           2. EVMS.FOR    The special form of the algorithm
C                          for only one functional equation.
C           3. EX1.FOR     Example 1. The complete Fortran
C                          code for the example in the paper.
C                          It has three functional equations.
C           4. EX1.DAT     The data required by Example 1 in the
C                          required format.
C           5. EX1.OUT     Output from Example 1.  Shown as
C                          Table 1 in the paper.
C           6. EX2.FOR     The complete program for Example 2
C                          which has only one functional equation.
C           7. EX2.DAT     The data for Example 2 in the required
C                          format.
C           8. EX2.OUT     Output from Example 2.
C
C Park Reilly
C
C
C 1. EVM.FOR, THE GENERAL ALGORITHM
C

      SUBROUTINE EVM(B,BV,DATA,EST,F,G1,G2,G3,H,IFAULT,
     *  NDAT,NEQ,NPAR,NVAR,PHI,Q1,Q2,RESID,S,T,THETA,V,XI,Z)
C
C  PARAMETER ESTIMATION FOR THE ERROR-IN-VARIABLES MODEL
C
      REAL B(NEQ,NVAR), BV(NEQ,NVAR), CRIT1, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), F(NEQ), G1(NPAR,NPAR), G2(NPAR,NPAR),
     * G3(NPAR,NPAR), H(NEQ), ONE, PHI, Q1(NPAR), Q2(NPAR),
     * RESID(NDAT,NVAR), S(NEQ,NEQ), T(NEQ), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, W2, XI(NVAR), Z(NEQ,NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT1=1.E-6, NITER1=50, ZERO=0.E0, ONE=1.E0)
      IFAULT(1) = 0
      IFAULT(2) = 0
      DO 4 ITER = 1,NITER1
C
C  EACH PASS THROUGH THIS LOOP PERFORMS A SINGLE OUTER ITERATION
C
        CALL INNER(B,BV,DATA,EST,F,G1,H,IFAULT,NDAT,
     *  NEQ,NPAR,NVAR,PHI,Q1,RESID,S,T,THETA,V,XI,Z)
C
C  PUT Q1 INTO Q2 AND G1 INTO G2
C
        DO 2 IPAR1 = 1,NPAR
          Q2(IPAR1) = Q1(IPAR1)
          DO 1 IPAR2 = IPAR1,NPAR
1         G2(IPAR1,IPAR2) = G1(IPAR1,IPAR2)
2       CONTINUE
C
C  CALCULATE CORRECTION FOR THETA
C
        CALL CHOL(G2,Q2,NPAR)
        CALL BKSUB(G2,Q2,1,NPAR,1)
C
C CHECK FOR CONVERGENCE
C
        WDIFF = ZERO
        DO 3 IPAR = 1,NPAR
          W1 = Q2(IPAR)
          W2 = THETA(IPAR)
          WDIFF = MAX(WDIFF,ABS(W1/W2))
          THETA(IPAR) = W2-W1
3       CONTINUE
        IF(WDIFF.LE.CRIT1)GO TO 5
4     CONTINUE
      IFAULT(1) = 1
C
C  FILL LOWER TRIANGLE OF G1 AND PREPARE G2
C  FOR THE INVERSION OF G1
C
5     DO 7 IPAR1 = 1,NPAR
        DO 6 IPAR2 = IPAR1,NPAR
          W1 = G1(IPAR1,IPAR2)
          G1(IPAR2,IPAR1) = W1
          G2(IPAR1,IPAR2) = W1
          G3(IPAR1,IPAR2) = ZERO
6       CONTINUE
7     CONTINUE
C
C  INVERT G1 INTO G3
C
      CALL CHOL(G2,Q2,NPAR)
      DO 8 IPAR1 = 1,NPAR
8     G3(IPAR1,IPAR1) = ONE/G2(IPAR1,IPAR1)
      CALL BKSUB(G2,G3,NPAR,NPAR,3)
      RETURN
      END


      SUBROUTINE INNER(B,BV,DATA,EST,F,G1,H,IFAULT,NDAT,NEQ,
     *    NPAR,NVAR,PHI,Q1,RESID,S,T,THETA,V,XI,Z)
      REAL B(NEQ,NVAR), BV(NEQ,NVAR), CRIT2, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), F(NEQ), G1(NPAR,NPAR), H(NEQ), HALF,
     * PHI, Q1(NPAR), RESID(NDAT,NVAR), S(NEQ,NEQ),
     * T(NEQ), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, XI(NVAR), Z(NEQ,NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT2=1.E-6, HALF=.5E0, NITER2=20, ZERO=0.E0)
C
C  INITIALIZE PHI, Q1, AND G1
C
      PHI = ZERO
      DO 2 IPAR1 = 1,NPAR
        Q1(IPAR1) = ZERO
        DO 1 IPAR2 = IPAR1,NPAR
1       G1(IPAR1,IPAR2) = ZERO
2     CONTINUE
C
C  PERFORM INNER ITERATION TO CONVERGENCE FOR EACH IDAT
C
      DO 24 IDAT = 1,NDAT
        DO 3 IVAR1 = 1,NVAR
3       XI(IVAR1) = DATA(IDAT,IVAR1)
        DO 16 ITER = 1,NITER2
C
C  THIS LOOP PERFORMS A SINGLE INNER ITERATION
C
          WDIFF = ZERO
          CALL BF(B,F,THETA,XI)
C
C  CALCULATE BV
C
          DO 6 IEQ1 = 1,NEQ
            DO 5 IVAR1 = 1,NVAR
              W1 = ZERO
              DO 4 IVAR2 = 1,NVAR
4             W1 = W1 + B(IEQ1,IVAR2)*V(IVAR1,IVAR2)
              BV(IEQ1,IVAR1) = W1
5           CONTINUE
6         CONTINUE
C
C  CALCULATE BVB' AND STORE IT IN S
C
          DO 9 IEQ1 = 1,NEQ
            DO 8 IEQ2 = IEQ1,NEQ
              W1 = ZERO
              DO 7 IVAR1 = 1,NVAR
7             W1 = W1 + BV(IEQ1,IVAR1)*B(IEQ2,IVAR1)
              S(IEQ1,IEQ2) = W1
8           CONTINUE
9         CONTINUE
C
C  PUT F+B(X-XI) INTO H
C
          DO 11 IEQ1 = 1,NEQ
            W1 = F(IEQ1)
            DO 10 IVAR1 = 1,NVAR
10          W1 = W1 + B(IEQ1,IVAR1)*(DATA(IDAT,IVAR1)-XI(IVAR1))
            H(IEQ1) = W1
11        CONTINUE
C
C  REPLACE S BY ITS TRIANGULAR FACTORIZATION
C  AND SOLVE ST = H FOR T
C
          CALL CHOL(S,H,NEQ)
          DO 13 IEQ1 = 1,NEQ
            T(IEQ1) = H(IEQ1)
            DO 12 IEQ2 = IEQ1,NEQ
12          S(IEQ2,IEQ1) = S(IEQ1,IEQ2)
13        CONTINUE
          CALL BKSUB(S,T,1,NEQ,1)
C
C  CALCULATE NEW XI AND PUT RELATIVE CHANGE IN XI INTO WDIFF
C
          DO 15 IVAR1 = 1,NVAR
            W1 = DATA(IDAT,IVAR1)
            DO 14 IEQ1 = 1,NEQ
14          W1 = W1 - BV(IEQ1,IVAR1)*T(IEQ1)
            WDIFF = MAX(WDIFF, ABS((W1-XI(IVAR1))/W1))
            XI(IVAR1) = W1
15        CONTINUE
C
C  CHECK FOR CONVERGENCE
C
          IF(WDIFF.LE.CRIT2)GO TO 17
16      CONTINUE
        IFAULT(2) = 1
C
C  ACCUMULATE PHI
C
17      DO 18 IEQ1 = 1,NEQ
18      PHI = PHI + H(IEQ1)**2
C
C  CALCULATE FITTED VALUES AND RESIDUALS
C
        DO 19 IVAR1 = 1,NVAR
          W1 = XI(IVAR1)
          EST(IDAT,IVAR1) = W1
          RESID(IDAT,IVAR1) = DATA(IDAT,IVAR1)-W1
19      CONTINUE
C
C  CALCULATE CONTRIBUTIONS TO Q1 AND G1 AND COMPLETE THE LOOPS
C
        CALL ZED(B,F,THETA,XI,Z)
        CALL BKSUB(S,Z,NPAR,NEQ,2)
        DO 23 IPAR1 = 1,NPAR
          W1 = Q1(IPAR1)
          DO 20 IEQ1 = 1,NEQ
20        W1 = W1 + Z(IEQ1,IPAR1)*H(IEQ1)
          Q1(IPAR1) = W1
          DO 22 IPAR2 = IPAR1,NPAR
            W1 = G1(IPAR1,IPAR2)
            DO 21 IEQ1 = 1,NEQ
21          W1 = W1 + Z(IEQ1,IPAR1)*Z(IEQ1,IPAR2)
            G1(IPAR1,IPAR2) = W1
22        CONTINUE
23      CONTINUE
24    CONTINUE
      PHI = PHI*HALF
      RETURN
      END


      SUBROUTINE CHOL(A,B,N)
C
C  PERFORM CHOLESKY DECOMPOSITION OF A AND PRELIMINARY TREATMENT OF B
C
      REAL A(N,N),B(N),W1,W2
      DO 5 IR1 = 1,N
        W1 = A(IR1,IR1)
        DO 1 IR2 = 1,IR1-1
1       W1 = W1 - A(IR2,IR1)**2
        W1 = SQRT(W1)
        A(IR1,IR1) = W1
        DO 3 IC1 = IR1+1,N
          W2 = A(IR1,IC1)
          DO 2 IR2 = 1,IR1-1
2         W2 = W2 - A(IR2,IR1)*A(IR2,IC1)
          A(IR1,IC1) = W2/W1
3       CONTINUE
        W2 = B(IR1)
        DO 4 IR2 = 1,IR1-1
4       W2 = W2 - A(IR2,IR1)*B(IR2)
        B(IR1) = W2/W1
5     CONTINUE
      RETURN
      END


      SUBROUTINE BKSUB(A,B,M,N,JOB)
C
C  PERFORM BACK-SOLUTION ACCORDING AS JOB = 1, 2, OR 3
C
      REAL A(N,N),B(N,M),W1
      DO 3 IC1 = M,1,-1
        IF(JOB.EQ.2) THEN
          K1 = 1
          K2 = N
          K3 = 1
        ELSE
          K2 = 1
          K3 = -1
          IF(JOB.EQ.1) THEN
            K1 = N
          ELSE
            K1 = IC1
          END IF
        END IF
        DO 2 IR1 = K1,K2,K3
          W1 = B(IR1,IC1)
          IF (JOB.EQ.2) THEN
            K2 = IR1-1
          ELSE
            K1 = N
            K2 = IR1+1
          END IF
          DO 1 IR2 = K1,K2,K3
1         W1 = W1 - A(IR1,IR2)*B(IR2,IC1)
          B(IR1,IC1) = W1/A(IR1,IR1)
          IF(JOB.EQ.3) B(IC1,IR1) = B(IR1,IC1)
2       CONTINUE
3     CONTINUE
      RETURN
      END

C----------------------------------------------------------------------
C
C 2. EVMS.FOR, THE SPECIAL FORM FOR ONLY ONE FUNCTIONAL EQUATION

      SUBROUTINE EVMS(BS,BVS,DATA,EST,G1,G2,G3,IFAULT,
     *    NDAT,NPAR,NVAR,PHI,Q1,Q2,RESID,THETA,V,XI,ZS)
C
C  PARAMETER ESTIMATION FOR THE ERROR-IN-VARIABLES MODEL WITH M = 1
C  91/3/21 PMR
C
      REAL BS(NVAR), BVS(NVAR), BVB, CRIT1, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), FS, G1(NPAR,NPAR), G2(NPAR,NPAR),
     * G3(NPAR,NPAR), HS, ONE, PHI, Q1(NPAR), Q2(NPAR),
     * RESID(NDAT,NVAR), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, W2, XI(NVAR), ZS(NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT1=1.E-6, NITER1=50, ONE=1.E0, ZERO=0.E0)
C
      IFAULT(1) = 0
      IFAULT(2) = 0
      DO 4 ITER = 1,NITER1
C
C  EACH PASS THROUGH THIS LOOP PERFORMS A SINGLE OUTER ITERATION
C
        CALL INNERS(BS,BVS,DATA,EST,G1,IFAULT,NDAT,
     *      NPAR,NVAR,PHI,Q1,RESID,THETA,V,XI,ZS)
C
C  PUT Q1 INTO Q2 AND G1 INTO G2
C
        DO 2 IPAR1 = 1,NPAR
          Q2(IPAR1) = Q1(IPAR1)
          DO 1 IPAR2 = IPAR1,NPAR
1         G2(IPAR1,IPAR2) = G1(IPAR1,IPAR2)
2       CONTINUE
C
C  CALCULATE CORRECTION FOR THETA
C
        CALL CHOL(G2,Q2,NPAR)
        CALL BKSUB(G2,Q2,1,NPAR,1)
C
C CHECK FOR CONVERGENCE
C
        WDIFF = ZERO
        DO 3 IPAR = 1,NPAR
          W1 = Q2(IPAR)
          W2 = THETA(IPAR)
          WDIFF = MAX(WDIFF,ABS(W1/W2))
          THETA(IPAR) = W2-W1
3       CONTINUE
        IF(WDIFF.LE.CRIT1) GO TO 5
4     CONTINUE
      IFAULT(1) = 1
C
C  FILL LOWER TRIANGLE OF G1 AND PREPARE G2 FOR THE INVERSION OF G1
C
5     DO 7 IPAR1 = 1,NPAR
        DO 6 IPAR2 = IPAR1,NPAR
          W1 = G1(IPAR1,IPAR2)
          G1(IPAR2,IPAR1) = W1
          G2(IPAR1,IPAR2) = W1
          G3(IPAR1,IPAR2) = ZERO
6       CONTINUE
7     CONTINUE
C
C  INVERT G1 INTO G3
C
      CALL CHOL(G2,Q2,NPAR)
      DO 8 IPAR1 = 1,NPAR
8     G3(IPAR1,IPAR1) = ONE/G2(IPAR1,IPAR1)
      CALL BKSUB(G2,G3,NPAR,NPAR,3)
      RETURN
      END


      SUBROUTINE INNERS(BS,BVS,DATA,EST,G1,IFAULT,NDAT,
     *    NPAR,NVAR,PHI,Q1,RESID,THETA,V,XI,ZS)
C
C  INNER ITERATION FOR M = 1
C
      REAL BS(NVAR),BVS(NVAR),BVB,CRIT2,DATA(NDAT,NVAR),
     * EST(NDAT,NVAR),FS,HALF,G1(NPAR,NPAR),HS,
     * PHI,Q1(NPAR),RESID(NDAT,NVAR),THETA(NPAR),
     * V(NVAR,NVAR),WDIFF,W1,XI(NVAR),ZS(NPAR),ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT2=1.E-5, HALF=.5E0, NITER2=15, ZERO=0.E0)
C
C  INITIALIZE PHI, Q1, AND G1
C
      PHI = ZERO
      DO 2 IPAR1 = 1,NPAR
        Q1(IPAR1) = ZERO
        DO 1 IPAR2 = IPAR1,NPAR
1       G1(IPAR1,IPAR2) = ZERO
2     CONTINUE
C
C  PERFORM INNER ITERATION TO CONVERGENCE FOR EACH IDAT
C
      DO 14 IDAT = 1,NDAT
C
C  INITIALIZE XI
C
        DO 3 IVAR1 = 1,NVAR
3       XI(IVAR1) = DATA(IDAT,IVAR1)
        DO 9 ITER = 1,NITER2
C
C  THIS LOOP PERFORMS A SINGLE INNER ITERATION
C
          WDIFF = ZERO
          CALL BF(BS,FS,THETA,XI)
C
C  CALCULATE BVS
C
          DO 5 IVAR1 = 1,NVAR
            W1 = ZERO
            DO 4 IVAR2 = 1,NVAR
4           W1 = W1 + BS(IVAR2)*V(IVAR1,IVAR2)

            BVS(IVAR1) = W1
5         CONTINUE
C
C  CALCULATE BVB
C
          BVB = ZERO
          DO 6 IVAR1 = 1,NVAR
6         BVB = BVB + BVS(IVAR1)*BS(IVAR1)
C
C  PUT F+B(X-XI) INTO HS
C
          HS = FS
          DO 7 IVAR1 = 1,NVAR
7         HS = HS + BS(IVAR1)*(DATA(IDAT,IVAR1)-XI(IVAR1))
C
C  CALCULATE NEW XI AND PUT RELATIVE CHANGE IN XI INTO WDIFF
C
          DO 8 IVAR1 = 1,NVAR
            W1 = DATA(IDAT,IVAR1) - BVS(IVAR1)*HS/BVB
            WDIFF = MAX(WDIFF, ABS((W1-XI(IVAR1))/W1))
            XI(IVAR1) = W1
8         CONTINUE
C
C  CHECK FOR CONVERGENCE
C
          IF(WDIFF.LE.CRIT2) GO TO 10
9       CONTINUE
        IFAULT(2) = 1
C
C  ACCUMULATE PHI
C
10      PHI = PHI+HS*HS/BVB
C
C  CALCULATE FITTED VALUES AND RESIDUALS
C
        DO 11 IVAR1 = 1,NVAR
          W1 = XI(IVAR1)
          EST(IDAT,IVAR1) = W1
          RESID(IDAT,IVAR1) = DATA(IDAT,IVAR1) - W1
11      CONTINUE
C
C  CALCULATE CONTRIBUTIONS TO Q1 AND G1 AND COMPLETE THE LOOPS
C
        CALL ZED(BS,FS,THETA,XI,ZS)
        DO 13 IPAR1 = 1,NPAR
          W1 = ZS(IPAR1)
          Q1(IPAR1) = Q1(IPAR1) + W1*HS/BVB
          DO 12 IPAR2 = IPAR1,NPAR
12        G1(IPAR1,IPAR2) = G1(IPAR1,IPAR2) + W1*ZS(IPAR2)/BVB
13      CONTINUE
14    CONTINUE
      PHI = PHI*HALF
      RETURN
      END


      SUBROUTINE CHOL(A,B,N)
C
C  PERFORM CHOLESKY DECOMPOSITION OF A AND PRELIMINARY TREATMENT
C     OF B
C
      REAL A(N,N),B(N),W1,W2
      DO 5 IR1 = 1,N
        W1 = A(IR1,IR1)
        DO 1 IR2 = 1,IR1-1
1       W1 = W1 - A(IR2,IR1)**2
        W1 = SQRT(W1)
        A(IR1,IR1) = W1
        DO 3 IC1 = IR1+1,N
          W2 = A(IR1,IC1)
          DO 2 IR2 = 1,IR1-1
2         W2 = W2 - A(IR2,IR1)*A(IR2,IC1)
          A(IR1,IC1) = W2/W1
3       CONTINUE
        W2 = B(IR1)
        DO 4 IR2 = 1,IR1-1
4       W2 = W2 - A(IR2,IR1)*B(IR2)
        B(IR1) = W2/W1
5     CONTINUE
      RETURN
      END


      SUBROUTINE BKSUB(A,B,M,N,JOB)
C
C  PERFORM BACK-SOLUTION ACCORDING AS JOB = 1, 2, OR 3
C
      REAL A(N,N),B(N,M),W1
      DO 3 IC1 = M,1,-1
        IF(JOB.EQ.2) THEN
          K1 = 1
          K2 = N
          K3 = 1
        ELSE
          K2 = 1
          K3 = -1
          IF(JOB.EQ.1) THEN
            K1 = N
          ELSE
            K1 = IC1
          END IF
        END IF
        DO 2 IR1 = K1,K2,K3
          W1 = B(IR1,IC1)
          IF (JOB.EQ.2) THEN
            K2 = IR1-1
          ELSE
            K1 = N
            K2 = IR1+1
          END IF
          DO 1 IR2 = K1,K2,K3
1         W1 = W1 - A(IR1,IR2)*B(IR2,IC1)
          B(IR1,IC1) = W1/A(IR1,IR1)
          IF(JOB.EQ.3) B(IC1,IR1) = B(IR1,IC1)
2       CONTINUE
3     CONTINUE
      RETURN
      END

C----------------------------------------------------------------------
C
C 3. EX1.FOR, EXAMPLE 1
C
C This example is fully described in the paper. The description is
C repeated here.  The problem is a small artificial one (called the
C New Brunswick problem) in which it is wanted to obtain estimates
C for the elements of the parameter vector ~ in the following three
C functional equations which apply for all i=1,2,...,6
C
C             (see published paper)
C
C The parameters which describe the magnitude of the problem are:
C     the number of data vectors, NDAT=6
C     the number of equations, NEQ=3
C     the number of model parameters to be estimated, NPAR=4
C     the number of measured variables in a data vector, NVAR=5.
C The program parameters CRIT1 and CRIT2 were both set at 1.E-6 while NITER1
C and NITER2 were 50 and 20.  The full Fortran code as used is listed in
C EX1.FOR.  The data in the form required for input to the program are given
C in the file EX1.DAT.  The same data are listed with descriptive headings at
C the beginning of the file EX1.OUT.  The latter part of EX1.OUT shows the
C output from the program.
C
C EX1.FOR

C  TEST PROGRAM FOR EVM USING NEW BRUNSWICK PROBLEM   91/3/27 PMR

      DIMENSION B(3,5), BV(3,5), DATA(6,5), EST(6,5), F(3), G1(4,4),
     * G2(4,4), G3(4,4), H(3), IFAULT(2), Q1(4), Q2(4), RESID(6,5),
     * S(3,3), T(3), THETA(4), V(5,5), XI(5), Z(3,4)
C
      READ(20,*) THETA
      READ(20,*) ((DATA(I,J),J=1,5),I=1,6)
      READ(20,*) V
      WRITE(30,1)
1     FORMAT(36X,'TABLE 1')
      WRITE(30,2)
2     FORMAT(//'DATA')
      WRITE(30,3) ((DATA(I,J),J=1,5),I=1,6)
3     FORMAT(5F9.2)
      WRITE(30,4)
4     FORMAT(//'V')
      WRITE(30,3) V
      WRITE(30,5) THETA
5     FORMAT(//'INITIAL PARAMETER ESTIMATES = ',4F9.2)
      CALL EVM(B,BV,DATA,EST,F,G1,G2,G3,H,IFAULT,
     *    6,3,4,5,PHI,Q1,Q2,RESID,S,T,THETA,V,XI,Z)
      WRITE(30,6) THETA
6     FORMAT(//'FINAL PARAMETER ESTIMATES = ',4F9.6)
      WRITE(30,7) PHI
7     FORMAT(//'FINAL PHI = ',E13.7)
      WRITE(30,8) Q1
8     FORMAT(//'FINAL Q = ',4E16.7)
      WRITE(30,9)
9     FORMAT(//'ESTIMATED TRUE VALUES OF MEASUREMENTS')
      WRITE(30,10) ((EST(I,J),J=1,5),I=1,6)
10    FORMAT(5F9.4)
      WRITE(30,11)
11    FORMAT(//'RESIDUALS')
      WRITE(30,10) ((RESID(I,J),J=1,5),I=1,6)
      WRITE(30,12)
12    FORMAT(//'EXPECTED INFORMATION MATRIX')
      WRITE(30,13) G1
13    FORMAT(4F9.4)
      WRITE(30,14)
14    FORMAT(//'INVERSE OF EXPECTED INFORMATION MATRIX')
      WRITE(30,13) G3
      WRITE(30,15)
15    FORMAT(//'FAULT STATUS')
      WRITE(30,16) IFAULT
16    FORMAT(2I5)
      STOP
      END


      SUBROUTINE BF(B,F,THETA,XI)
      DIMENSION B(3,5), F(3), THETA(4), XI(5)
C
      F(1) = THETA(1)**2 + THETA(2)*XI(1) + THETA(3)*XI(2)**2 - XI(3)
      F(2) = THETA(2)**2 + (THETA(3)*XI(1))**2
     *       + THETA(4)**3*XI(3)**2 - XI(4)
      F(3) = THETA(3)**2 + (THETA(1)*XI(2))**2 - XI(5)
      B(1,1) = THETA(2)
      B(1,2) = 2.E0*THETA(3)*XI(2)
      B(1,3) = -1.E0
      B(1,4) = 0.E0
      B(1,5) = 0.E0
      B(2,1) = 2.E0*THETA(3)**2*XI(1)
      B(2,2) = 0.E0
      B(2,3) = 2.E0*THETA(4)**3*XI(3)
      B(2,4) = -1.E0
      B(2,5) = 0.E0
      B(3,1) = 0.E0
      B(3,2) = 2.E0*THETA(1)**2*XI(2)
      B(3,3) = 0.E0
      B(3,4) = 0.E0
      B(3,5) = -1.E0
      RETURN
      END


      SUBROUTINE ZED(B,F,THETA,XI,Z)
      DIMENSION B(3,5),F(3),THETA(4),XI(5),Z(3,4)
      Z(1,1) = 2.E0*THETA(1)
      Z(1,2) = XI(1)
      Z(1,3) = XI(2)**2
      Z(1,4) = 0.E0
      Z(2,1) = 0.E0
      Z(2,2) = 2.E0*THETA(2)
      Z(2,3) = 2.E0*THETA(3)*XI(1)**2
      Z(2,4) = 3.E0*(THETA(4)*XI(3))**2
      Z(3,1) = 2.E0*THETA(1)*XI(2)**2
      Z(3,2) = 0.E0
      Z(3,3) = 2.E0*THETA(3)
      Z(3,4) = 0.E0
      RETURN
      END


      SUBROUTINE EVM(B,BV,DATA,EST,F,G1,G2,G3,H,IFAULT,
     *    NDAT,NEQ,NPAR,NVAR,PHI,Q1,Q2,RESID,S,T,THETA,V,XI,Z)
C
C  PARAMETER ESTIMATION FOR THE ERROR-IN-VARIABLES MODEL
C
      REAL B(NEQ,NVAR), BV(NEQ,NVAR), CRIT1, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), F(NEQ), G1(NPAR,NPAR), G2(NPAR,NPAR),
     * G3(NPAR,NPAR), H(NEQ), ONE, PHI, Q1(NPAR), Q2(NPAR),
     * RESID(NDAT,NVAR), S(NEQ,NEQ), T(NEQ), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, W2, XI(NVAR), Z(NEQ,NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT1=1.E-6, NITER1=50, ZERO=0.E0, ONE=1.E0)
      IFAULT(1) = 0
      IFAULT(2) = 0
      DO 4 ITER = 1,NITER1
C
C  EACH PASS THROUGH THIS LOOP PERFORMS A SINGLE OUTER ITERATION
C
        CALL INNER(B,BV,DATA,EST,F,G1,H,IFAULT,NDAT,
     *      NEQ,NPAR,NVAR,PHI,Q1,RESID,S,T,THETA,V,XI,Z)
C
C  PUT Q1 INTO Q2 AND G1 INTO G2
C
        DO 2 IPAR1 = 1,NPAR
          Q2(IPAR1) = Q1(IPAR1)
          DO 1 IPAR2 = IPAR1,NPAR
1         G2(IPAR1,IPAR2) = G1(IPAR1,IPAR2)
2       CONTINUE
C
C  CALCULATE CORRECTION FOR THETA
C
        CALL CHOL(G2,Q2,NPAR)
        CALL BKSUB(G2,Q2,1,NPAR,1)
C
C CHECK FOR CONVERGENCE
C
        WDIFF = ZERO
        DO 3 IPAR = 1,NPAR
          W1 = Q2(IPAR)
          W2 = THETA(IPAR)
          WDIFF = MAX(WDIFF,ABS(W1/W2))
          THETA(IPAR) = W2-W1
3       CONTINUE
        IF(WDIFF.LE.CRIT1)GO TO 5
4     CONTINUE
      IFAULT(1) = 1
C
C  FILL LOWER TRIANGLE OF G1 AND PREPARE G2
C  FOR THE INVERSION OF G1
C
5     DO 7 IPAR1 = 1,NPAR
        DO 6 IPAR2 = IPAR1,NPAR
          W1 = G1(IPAR1,IPAR2)
          G1(IPAR2,IPAR1) = W1
          G2(IPAR1,IPAR2) = W1
          G3(IPAR1,IPAR2) = ZERO
6       CONTINUE
7     CONTINUE
C
C  INVERT G1 INTO G3
C
      CALL CHOL(G2,Q2,NPAR)
      DO 8 IPAR1 = 1,NPAR
8     G3(IPAR1,IPAR1) = ONE/G2(IPAR1,IPAR1)
      CALL BKSUB(G2,G3,NPAR,NPAR,3)
      RETURN
      END


      SUBROUTINE INNER(B,BV,DATA,EST,F,G1,H,IFAULT,NDAT,NEQ,
     * NPAR,NVAR,PHI,Q1,RESID,S,T,THETA,V,XI,Z)
      REAL B(NEQ,NVAR), BV(NEQ,NVAR), CRIT2, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), F(NEQ), G1(NPAR,NPAR), H(NEQ), HALF,
     * PHI, Q1(NPAR), RESID(NDAT,NVAR), S(NEQ,NEQ),
     * T(NEQ), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, XI(NVAR), Z(NEQ,NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT2=1.E-6, HALF=.5E0, NITER2=20, ZERO=0.E0)
C
C  INITIALIZE PHI, Q1, AND G1
C
      PHI = ZERO
      DO 2 IPAR1 = 1,NPAR
        Q1(IPAR1) = ZERO
        DO 1 IPAR2 = IPAR1,NPAR
1       G1(IPAR1,IPAR2) = ZERO
2     CONTINUE
C
C  PERFORM INNER ITERATION TO CONVERGENCE FOR EACH IDAT
C
      DO 24 IDAT = 1,NDAT
        DO 3 IVAR1 = 1,NVAR
3       XI(IVAR1) = DATA(IDAT,IVAR1)
        DO 16 ITER = 1,NITER2
C
C  THIS LOOP PERFORMS A SINGLE INNER ITERATION
C
          WDIFF = ZERO
          CALL BF(B,F,THETA,XI)
C
C  CALCULATE BV
C
          DO 6 IEQ1 = 1,NEQ
            DO 5 IVAR1 = 1,NVAR
              W1 = ZERO
              DO 4 IVAR2 = 1,NVAR
4             W1 = W1+B(IEQ1,IVAR2)*V(IVAR1,IVAR2)
              BV(IEQ1,IVAR1) = W1
5           CONTINUE
6         CONTINUE
C
C  CALCULATE BVB' AND STORE IT IN S
C
          DO 9 IEQ1 = 1,NEQ
            DO 8 IEQ2 = IEQ1,NEQ
              W1 = ZERO
              DO 7 IVAR1 = 1,NVAR
7             W1 = W1+BV(IEQ1,IVAR1)*B(IEQ2,IVAR1)
              S(IEQ1,IEQ2) = W1
8           CONTINUE
9         CONTINUE
C
C  PUT F+B(X-XI) INTO H
C
          DO 11 IEQ1 = 1,NEQ
            W1 = F(IEQ1)
            DO 10 IVAR1 = 1,NVAR
10          W1 = W1+B(IEQ1,IVAR1)*(DATA(IDAT,IVAR1)-XI(IVAR1))
            H(IEQ1) = W1
11        CONTINUE
C
C  REPLACE S BY ITS TRIANGULAR FACTORIZATION
C  AND SOLVE ST = H FOR T
C
          CALL CHOL(S,H,NEQ)
          DO 13 IEQ1 = 1,NEQ
            T(IEQ1) = H(IEQ1)
            DO 12 IEQ2 = IEQ1,NEQ
12          S(IEQ2,IEQ1) = S(IEQ1,IEQ2)
13        CONTINUE
          CALL BKSUB(S,T,1,NEQ,1)
C
C  CALCULATE NEW XI AND PUT RELATIVE CHANGE IN XI INTO WDIFF
C
          DO 15 IVAR1 = 1,NVAR
            W1 = DATA(IDAT,IVAR1)
            DO 14 IEQ1 = 1,NEQ
14          W1 = W1 - BV(IEQ1,IVAR1)*T(IEQ1)
            WDIFF = MAX(WDIFF, ABS((W1-XI(IVAR1))/W1))
            XI(IVAR1) = W1
15        CONTINUE
C
C  CHECK FOR CONVERGENCE
C
          IF(WDIFF.LE.CRIT2)GO TO 17
16      CONTINUE
        IFAULT(2) = 1
C
C  ACCUMULATE PHI
C
17      DO 18 IEQ1 = 1,NEQ
18      PHI = PHI + H(IEQ1)**2
C
C  CALCULATE FITTED VALUES AND RESIDUALS
C
        DO 19 IVAR1 = 1,NVAR
          W1 = XI(IVAR1)
          EST(IDAT,IVAR1) = W1
          RESID(IDAT,IVAR1) = DATA(IDAT,IVAR1)-W1
19      CONTINUE
C
C  CALCULATE CONTRIBUTIONS TO Q1 AND G1 AND COMPLETE THE LOOPS
C
        CALL ZED(B,F,THETA,XI,Z)
        CALL BKSUB(S,Z,NPAR,NEQ,2)
        DO 23 IPAR1 = 1,NPAR
          W1 = Q1(IPAR1)
          DO 20 IEQ1 = 1,NEQ
20        W1 = W1 + Z(IEQ1,IPAR1)*H(IEQ1)
          Q1(IPAR1) = W1
          DO 22 IPAR2 = IPAR1,NPAR
            W1 = G1(IPAR1,IPAR2)
            DO 21 IEQ1 = 1,NEQ
21          W1 = W1 + Z(IEQ1,IPAR1)*Z(IEQ1,IPAR2)
            G1(IPAR1,IPAR2) = W1
22        CONTINUE
23      CONTINUE
24    CONTINUE
      PHI = PHI*HALF
      RETURN
      END


      SUBROUTINE CHOL(A,B,N)
C
C  PERFORM CHOLESKY DECOMPOSITION OF A AND PRELIMINARY TREATMENT
C          OF B
C
      REAL A(N,N),B(N),W1,W2
      DO 5 IR1 = 1,N
        W1 = A(IR1,IR1)
        DO 1 IR2 = 1,IR1-1
1       W1 = W1 - A(IR2,IR1)**2
        W1 = SQRT(W1)
        A(IR1,IR1) = W1
        DO 3 IC1 = IR1+1,N
          W2 = A(IR1,IC1)
          DO 2 IR2 = 1,IR1-1
2         W2 = W2 - A(IR2,IR1)*A(IR2,IC1)
          A(IR1,IC1) = W2/W1
3       CONTINUE
        W2 = B(IR1)
        DO 4 IR2 = 1,IR1-1
4       W2 = W2 - A(IR2,IR1)*B(IR2)
        B(IR1) = W2/W1
5     CONTINUE
      RETURN
      END


      SUBROUTINE BKSUB(A,B,M,N,JOB)
C
C  PERFORM BACK-SOLUTION ACCORDING AS JOB = 1, 2, OR 3
C
      REAL A(N,N),B(N,M),W1
      DO 3 IC1 = M,1,-1
        IF(JOB.EQ.2) THEN
          K1 = 1
          K2 = N
          K3 = 1
        ELSE
          K2 = 1
          K3 = -1
          IF(JOB.EQ.1) THEN
            K1 = N
          ELSE
            K1 = IC1
          END IF
        END IF
        DO 2 IR1 = K1,K2,K3
          W1 = B(IR1,IC1)
          IF (JOB.EQ.2) THEN
            K2 = IR1-1
          ELSE
            K1 = N
            K2 = IR1+1
          END IF
          DO 1 IR2 = K1,K2,K3
1         W1 = W1 - A(IR1,IR2)*B(IR2,IC1)
          B(IR1,IC1) = W1/A(IR1,IR1)
          IF(JOB.EQ.3) B(IC1,IR1) = B(IR1,IC1)
2       CONTINUE
3     CONTINUE
      RETURN
      END


C----------------------------------------------------------------------
C
C 4. EX1.DAT

1 1 1 1
-.78 1.39 3.99 11.7 3.44
-1.27 2 5.15 36.58 4.59
2.62 .38 5.77 16.49 .33
2.05 2.92 6.78 75.68 2.74
1.05 2.46 12.85 149.42 7.59
2.96 1.47 3.9 52.11 6.82
1 -.2 .1 .3 0
-.2 1 .2 -.1 0
.1 .2 2 .4 -.2
.3 -.1 .4 5 .4
0 0 -.2 .4 3


C----------------------------------------------------------------------
C
C 5. EX1.OUT
                                    TABLE 1


DATA
    -0.78     1.39     3.99    11.70     3.44
    -1.27     2.00     5.15    36.58     4.59
     2.62     0.38     5.77    16.49     0.33
     2.05     2.92     6.78    75.68     2.74
     1.05     2.46    12.85   149.42     7.59
     2.96     1.47     3.90    52.11     6.82


V
     1.00    -0.20     0.10     0.30     0.00
    -0.20     1.00     0.20    -0.10     0.00
     0.10     0.20     2.00     0.40    -0.20
     0.30    -0.10     0.40     5.00     0.40
     0.00     0.00    -0.20     0.40     3.00


INITIAL PARAMETER ESTIMATES =      1.00     1.00     1.00
1.00


FINAL PARAMETER ESTIMATES =  0.822474 0.445297 1.442967 1.001551


FINAL PHI = 0.7515061E+01


FINAL Q =   -0.4954202E-05   0.1068151E-05   0.1989809E-05
0.2217203E-04


ESTIMATED TRUE VALUES OF MEASUREMENTS
  -0.6633   1.4147   3.2692  11.8515   3.4361
  -1.3739   1.9729   5.6814  36.5576   4.7153
   1.7831   1.0674   3.1146  16.5644   2.8529
   2.9277   1.9753   7.6105  76.2357   4.7217
   0.8652   2.7678  12.1155 149.2269   7.2642
   3.1432   1.5559   5.5692  51.9294   3.7197


RESIDUALS
  -0.1167  -0.0247   0.7208  -0.1515   0.0039
   0.1039   0.0271  -0.5314   0.0224  -0.1253
   0.8369  -0.6874   2.6554  -0.0744  -2.5229
  -0.8777   0.9447  -0.8305  -0.5557  -1.9817
   0.1848  -0.3078   0.7345   0.1931   0.3258
  -0.1832  -0.0859  -1.6692   0.1806   3.1003


EXPECTED INFORMATION MATRIX
  72.8418  -3.7697   4.1724 -58.4874
  -3.7697   2.9379   0.7385   6.4070
   4.1724   0.7385  10.2625  18.3139
 -58.4874   6.4070  18.3139 334.4229


INVERSE OF EXPECTED INFORMATION MATRIX
   0.0187   0.0197  -0.0157   0.0037
   0.0197   0.3781  -0.0315  -0.0021
  -0.0157  -0.0315   0.1218  -0.0088
   0.0037  -0.0021  -0.0088   0.0042


FAULT STATUS
    0    0


C----------------------------------------------------------------------
C
C 6. EX2.FOR, EXAMPLE 2
C
C The data for this example were obtained as digitized observations on the
C outline of the nerve-head in a subject's eye.  There are 32 data (NDAT),
C one functional relation (NEQ), two variables (NVAR), and five parameters
C (NPAR).  The functional relation describes an ellipse with the following
C equation
C
C  F = T3*(XI1-T1)**2 + 2*T4*(XI1-T1)*(XI2-T2) + T5*(XI2-T2)**2 - 1 = 0
C
C where F is the functional value which numerically must be brought to zero,
C T1,T2,T3,T4 and T5 are the functional parameters the values of which are
C to be estimated.  XI1 and XI2 are the estimated true values of the measured
C variables.  The T's are represented by the Greek letter theta in the paper
C and the XI's by xi.  The parameters which describe the magnitude of the
C problem are NDAT = 32, NPAR = 4, NVAR = 2.  The program parameters were
C set at CRIT1 = 1.E-6 and CRIT2 = 1.E-5.  NITER1 was 50 and NITER2 15.  The
C full Fortran code as used is listed under EX2.FOR.  The data in the format
C used by the program is under EX2.DAT.  The output is listed under EX2.OUT.
C
C EX2.FOR


C FITTING OF ELLIPSE BY EVMS  91/3/25 PMR

      DIMENSION BS(2),BVS(2),DATA(32,2),EST(32,2),
     * G1(5,5),G2(5,5),
     * G3(5,5),IFAULT(2),Q1(5),Q2(5),RESID(32,2),THETA(5),V(2,2),
     * XI(2),ZS(5)
      DATA V/1,0,0,1/
C
      READ(21,*)((DATA(I,J),J=1,2),I=1,32)
      DO 2 I = 1,32
        DO 1 J = 1,2
1       DATA(I,J) = .001*DATA(I,J)
2     CONTINUE
      READ(21,*) THETA
      CALL EVMS(BS,BVS,DATA,EST,G1,G2,G3,IFAULT,32,5,2,PHI,
     *    Q1,Q2,RESID,THETA,V,XI,ZS)
      DO 4 I = 1,32
        DO 3 J = 1,2
          DATA(I,J) = 1000.*DATA(I,J)
          EST(I,J) = 1000.*EST(I,J)
          RESID(I,J) = 1000.*RESID(I,J)
3       CONTINUE
4     CONTINUE
      WRITE(31,5)
5     FORMAT(//'FAULT STATUS')
      WRITE(31,6) IFAULT
6     FORMAT(2I5)
      WRITE(31,7)
7     FORMAT(//'PARAMETER ESTIMATES')
      WRITE(31,8) THETA
8     FORMAT(5F9.4)
      WRITE(31,9)
9     FORMAT(//'INFORMATION MATRIX')
      WRITE(31,8) G1
      WRITE(31,10)
10    FORMAT(//'INVERSE OF INFORMATION MATRIX')
      WRITE(31,8) G3
      WRITE(31,11)
11    FORMAT(//6X,'DATA',11X,'ESTIMATES',14X,'RESIDUALS')
      WRITE(31,12)((DATA(I,J),J=1,2),(EST(I,J),J=1,2),
     *  (RESID(I,J),J=1,2),I=1,32)
12    FORMAT(2F7.0,2F11.4,1X,2F10.4)
      STOP
      END


      SUBROUTINE BF(B,F,THETA,XI)
      DIMENSION B(2),THETA(5),XI(2)
C
C CALCULATE F
      W1 = XI(1)-THETA(1)
      W2 = XI(2)-THETA(2)
      F = THETA(3)*W1*W1 + 2.D0*THETA(4)*W1*W2 + THETA(5)*W2*W2 - 1.D0
C
C CALCULATE B
      B(1) = 2.D0*(W1*THETA(3) + W2*THETA(4))
      B(2) = 2.D0*(W1*THETA(4) + W2*THETA(5))
      RETURN
      END


      SUBROUTINE ZED(B,F,THETA,XI,Z)
      DIMENSION B(2),THETA(5),XI(2),Z(5)
C
      W1 = XI(1)-THETA(1)
      W2 = XI(2)-THETA(2)
      Z(1) = -2.D0*(W1*THETA(3) + W2*THETA(4))
      Z(2) = -2.D0*(W1*THETA(4) + W2*THETA(5))
      Z(3) = W1*W1
      Z(4) = 2.D0*W1*W2
      Z(5) = W2*W2
      RETURN
      END


      SUBROUTINE EVMS(BS,BVS,DATA,EST,G1,G2,G3,IFAULT,
     *    NDAT,NPAR,NVAR,PHI,Q1,Q2,RESID,THETA,V,XI,ZS)
C
C  PARAMETER ESTIMATION FOR THE ERROR-IN-VARIABLES MODEL WITH M = 1
C  91/3/21 PMR
C
      REAL BS(NVAR), BVS(NVAR), BVB, CRIT1, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), FS, G1(NPAR,NPAR), G2(NPAR,NPAR),
     * G3(NPAR,NPAR), HS, ONE, PHI, Q1(NPAR), Q2(NPAR),
     * RESID(NDAT,NVAR), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, W2, XI(NVAR), ZS(NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER(CRIT1=1.E-6,NITER1=50,ONE=1.E0,ZERO=0.E0)
C
      IFAULT(1) = 0
      IFAULT(2) = 0
      DO 4 ITER = 1,NITER1
C
C  EACH PASS THROUGH THIS LOOP PERFORMS A SINGLE OUTER ITERATION
C
        CALL INNERS(BS,BVS,DATA,EST,G1,IFAULT,NDAT,
     *      NPAR,NVAR,PHI,Q1,RESID,THETA,V,XI,ZS)
C
C  PUT Q1 INTO Q2 AND G1 INTO G2
C
        DO 2 IPAR1 = 1,NPAR
          Q2(IPAR1) = Q1(IPAR1)
          DO 1 IPAR2 = IPAR1,NPAR
1         G2(IPAR1,IPAR2) = G1(IPAR1,IPAR2)
2       CONTINUE
C
C  CALCULATE CORRECTION FOR THETA
C
        CALL CHOL(G2,Q2,NPAR)
        CALL BKSUB(G2,Q2,1,NPAR,1)
C
C CHECK FOR CONVERGENCE
C
        WDIFF = ZERO
        DO 3 IPAR = 1,NPAR
          W1 = Q2(IPAR)
          W2 = THETA(IPAR)
          WDIFF = MAX(WDIFF,ABS(W1/W2))
          THETA(IPAR) = W2-W1
3       CONTINUE
        IF(WDIFF.LE.CRIT1)GO TO 5
4     CONTINUE
      IFAULT(1) = 1
C
C  FILL LOWER TRIANGLE OF G1 AND PREPARE G2
C  FOR THE INVERSION OF G1
C
5     DO 7 IPAR1 = 1,NPAR
        DO 6 IPAR2 = IPAR1,NPAR
          W1 = G1(IPAR1,IPAR2)
          G1(IPAR2,IPAR1) = W1
          G2(IPAR1,IPAR2) = W1
          G3(IPAR1,IPAR2) = ZERO
6       CONTINUE
7     CONTINUE
C
C  INVERT G1 INTO G3
C
      CALL CHOL(G2,Q2,NPAR)
      DO 8 IPAR1 = 1,NPAR
8     G3(IPAR1,IPAR1) = ONE/G2(IPAR1,IPAR1)
      CALL BKSUB(G2,G3,NPAR,NPAR,3)
      RETURN
      END


      SUBROUTINE INNERS(BS,BVS,DATA,EST,G1,IFAULT,NDAT,
     *    NPAR,NVAR,PHI,Q1,RESID,THETA,V,XI,ZS)
C
C  INNER ITERATION FOR M = 1
C
      REAL BS(NVAR), BVS(NVAR), BVB, CRIT2, DATA(NDAT,NVAR),
     * EST(NDAT,NVAR), FS, HALF, G1(NPAR,NPAR), HS,
     * PHI, Q1(NPAR), RESID(NDAT,NVAR), THETA(NPAR),
     * V(NVAR,NVAR), WDIFF, W1, XI(NVAR), ZS(NPAR), ZERO
      DIMENSION IFAULT(2)
      PARAMETER (CRIT2=1.E-5, HALF=.5E0, NITER2=15, ZERO=0.E0)
C
C  INITIALIZE PHI, Q1, AND G1
C
      PHI = ZERO
      DO 2 IPAR1 = 1,NPAR
        Q1(IPAR1) = ZERO
        DO 1 IPAR2 = IPAR1,NPAR
1       G1(IPAR1,IPAR2) = ZERO
2     CONTINUE
C
C  PERFORM INNER ITERATION TO CONVERGENCE FOR EACH IDAT
C
      DO 14 IDAT = 1,NDAT
C
C  INITIALIZE XI
C
        DO 3 IVAR1 = 1,NVAR
3       XI(IVAR1) = DATA(IDAT,IVAR1)
        DO 9 ITER = 1,NITER2
C
C  THIS LOOP PERFORMS A SINGLE INNER ITERATION
C
          WDIFF = ZERO
          CALL BF(BS,FS,THETA,XI)
C
C  CALCULATE BVS
C
          DO 5 IVAR1 = 1,NVAR
            W1 = ZERO
            DO 4 IVAR2 = 1,NVAR
4           W1 = W1 + BS(IVAR2)*V(IVAR1,IVAR2)

            BVS(IVAR1) = W1
5         CONTINUE
C
C  CALCULATE BVB
C
          BVB = ZERO
          DO 6 IVAR1 = 1,NVAR
6         BVB = BVB + BVS(IVAR1)*BS(IVAR1)
C
C  PUT F+B(X-XI) INTO HS
C
          HS = FS
          DO 7 IVAR1 = 1,NVAR
7         HS = HS + BS(IVAR1)*(DATA(IDAT,IVAR1) - XI(IVAR1))
C
C  CALCULATE NEW XI AND PUT RELATIVE CHANGE IN XI INTO WDIFF
C
          DO 8 IVAR1 = 1,NVAR
            W1 = DATA(IDAT,IVAR1) - BVS(IVAR1)*HS/BVB
            WDIFF = MAX(WDIFF, ABS((W1-XI(IVAR1))/W1))
            XI(IVAR1) = W1
8         CONTINUE
C
C  CHECK FOR CONVERGENCE
C
          IF(WDIFF.LE.CRIT2)GO TO 10
9       CONTINUE
        IFAULT(2) = 1
C
C  ACCUMULATE PHI
C
10      PHI = PHI+HS*HS/BVB
C
C  CALCULATE FITTED VALUES AND RESIDUALS
C
        DO 11 IVAR1 = 1,NVAR
          W1 = XI(IVAR1)
          EST(IDAT,IVAR1) = W1
          RESID(IDAT,IVAR1) = DATA(IDAT,IVAR1) - W1
11      CONTINUE
C
C  CALCULATE CONTRIBUTIONS TO Q1 AND G1 AND COMPLETE THE LOOPS
C
        CALL ZED(BS,FS,THETA,XI,ZS)
        DO 13 IPAR1 = 1,NPAR
          W1 = ZS(IPAR1)
          Q1(IPAR1) = Q1(IPAR1) + W1*HS/BVB
          DO 12 IPAR2 = IPAR1,NPAR
12        G1(IPAR1,IPAR2) = G1(IPAR1,IPAR2) + W1*ZS(IPAR2)/BVB
13      CONTINUE
14    CONTINUE
      PHI = PHI*HALF
      RETURN
      END


      SUBROUTINE CHOL(A,B,N)
C
C  PERFORM CHOLESKY DECOMPOSITION OF A AND PRELIMINARY TREATMENT
C     OF B
C
      REAL A(N,N),B(N),W1,W2
      DO 5 IR1 = 1,N
        W1 = A(IR1,IR1)
        DO 1 IR2 = 1,IR1-1
1       W1 = W1-A(IR2,IR1)**2
        W1 = SQRT(W1)
        A(IR1,IR1) = W1
        DO 3 IC1 = IR1+1,N
          W2 = A(IR1,IC1)
          DO 2 IR2 = 1,IR1-1
2         W2 = W2 - A(IR2,IR1)*A(IR2,IC1)
          A(IR1,IC1) = W2/W1
3       CONTINUE
        W2 = B(IR1)
        DO 4 IR2 = 1,IR1-1
4       W2 = W2 - A(IR2,IR1)*B(IR2)
        B(IR1) = W2/W1
5     CONTINUE
      RETURN
      END


      SUBROUTINE BKSUB(A,B,M,N,JOB)
C
C  PERFORM BACK-SOLUTION ACCORDING AS JOB = 1, 2, OR 3
C
      REAL A(N,N),B(N,M),W1
      DO 3 IC1 = M,1,-1
        IF(JOB.EQ.2) THEN
          K1 = 1
          K2 = N
          K3 = 1
        ELSE
          K2 = 1
          K3 = -1
          IF(JOB.EQ.1) THEN
            K1 = N
          ELSE
            K1 = IC1
          END IF
        END IF
        DO 2 IR1 = K1,K2,K3
          W1 = B(IR1,IC1)
          IF (JOB.EQ.2) THEN
            K2 = IR1-1
          ELSE
            K1 = N
            K2 = IR1+1
          END IF
          DO 1 IR2 = K1,K2,K3
1         W1 = W1-A(IR1,IR2)*B(IR2,IC1)
          B(IR1,IC1) = W1/A(IR1,IR1)
          IF(JOB.EQ.3) B(IC1,IR1) = B(IR1,IC1)
2       CONTINUE
3     CONTINUE
      RETURN
      END

C----------------------------------------------------------------------

7. EX2.DAT

1420 470
1188 476
1015 539
830 655
660 825
543 1041
468 1268
459 1528
510 1804
628 2063
795 2221
1040 2345
1303 2417
1452 2431
1738 2411
1901 2318
2091 2222
2230 2117
2311 1984
2385 1837
2451 1630
2467 1462
2472 1263
2446 1112
2387 963
2255 796
2109 696
1929 597
1846 556
1620 486
1418 437
1259 453
1 1 1  -.001 1


8. EX2.OUT

FAULT STATUS
    0    0


PARAMETER ESTIMATES
   1.4594   1.4454   0.9582   0.0039   1.0186


INFORMATION MATRIX
  14.6768   0.4050  -1.1978   0.2560  -0.2880
   0.4050  17.3231   0.1325  -0.6215   1.2703
  -1.1978   0.1325   3.1070   0.0154   0.9994
   0.2560  -0.6215   0.0154   3.9975   0.1110
  -0.2880   1.2703   0.9994   0.1110   3.1565


INVERSE OF INFORMATION MATRIX
   0.0705  -0.0019   0.0277  -0.0049  -0.0014
  -0.0019   0.0600   0.0051   0.0102  -0.0263
   0.0277   0.0051   0.3698   0.0008  -0.1166
  -0.0049   0.0102   0.0008   0.2524  -0.0137
  -0.0014  -0.0263  -0.1166  -0.0137   0.3647


      DATA           ESTIMATES              RESIDUALS
  1420.   470.  1419.3920   455.4515     0.6078   14.5485
  1188.   476.  1191.7790   490.1515    -3.7791  -14.1515
  1015.   539.  1021.0570   552.0391    -6.0567  -13.0391
   830.   655.   835.8889   662.8408    -5.8889   -7.8408
   660.   825.   663.1567   827.6081    -3.1567   -2.6081
   543.  1041.   531.0313  1035.3420    11.9687    5.6577
   468.  1268.   455.4884  1265.5690    12.5116    2.4308
   459.  1528.   441.1519  1529.4960    17.8481   -1.4961
   510.  1804.   506.2181  1805.5060     3.7819   -1.5061
   628.  2063.   646.4303  2048.4960   -18.4303   14.5044
   795.  2221.   804.9326  2208.6670    -9.9326   12.3329
  1040.  2345.  1038.0070  2349.5790     1.9934   -4.5788
  1303.  2417.  1301.8300  2424.9300     1.1703   -7.9303
  1452.  2431.  1451.9830  2436.1840     0.0166   -5.1842
  1738.  2411.  1734.5620  2398.5040     3.4379   12.4955
  1901.  2318.  1908.4610  2333.5910    -7.4607  -15.5907
  2091.  2222.  2090.8450  2221.7980     0.1550    0.2022
  2230.  2117.  2218.0730  2105.9520    11.9267   11.0476
  2311.  1984.  2313.4760  1985.6700    -2.4755   -1.6704
  2385.  1837.  2394.2520  1841.1940    -9.2518   -4.1943
  2451.  1630.  2461.8330  1632.1880   -10.8328   -2.1882
  2467.  1462.  2480.6950  1462.2970   -13.6950   -0.2967
  2472.  1263.  2464.4650  1264.4130     7.5355   -1.4126
  2446.  1112.  2425.3450  1119.3370    20.6547   -7.3367
  2387.   963.  2361.4770   977.0116    25.5227  -14.0116
  2255.   796.  2242.8230   806.5399    12.1768  -10.5398
  2109.   696.  2117.7400   685.2596    -8.7397   10.7405
  1929.   597.  1942.6510   570.5976   -13.6507   26.4024
  1846.   556.  1856.1900   530.8438   -10.1901   25.1563
  1620.   486.  1622.9630   466.7079    -2.9625   19.2921
  1418.   437.  1418.7830   455.4772    -0.7826  -18.4772
  1259.   453.  1263.0220   473.7592    -4.0221  -20.7591
