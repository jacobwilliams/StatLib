       SUBROUTINE PLYFIT(NN,Y,NINTS,NOBS,DIVBEG,DIVEND,NK,NC,NOTE,
     * AGUESS,BOUNDD,BOUNDU,DBOUND,UBOUND,GROUP,KIND,ITERMX,QUICK,
     * A,SEA,CORR,ALOGL,P,CHISQ,NDF,AINV,AMU,DMUDA,YFN,IFAULT)
C
C         ALGORITHM AS 270.1 APPL.STATIST. (1992), VOL.41, NO.1

C         Subroutine for fitting a parametric 'key' function to data
C         (grouped or ungrouped) and for making polynomial adjustments 
C         to that fit.  An iterative maximum likelihood routine, based
C         on Newton-Raphson search, is employed.


C         Stephen T Buckland, Scottish Agricultural Statistics Service,
C         Aberdeen Unit, MLURI, Craigiebuckler, Aberdeen AB9 2QJ, Scotland

C         January, 1991

C         Auxiliary routines required: LUDCMP and LUBKSB from Numerical
C         Recipes by Press, Flannery, Teukolsky and Vetterling, Cambridge
C         University Press, 1986.

C         To change the key function, edit subroutines STD & KEY

C       IFAULT = 0 if convergence achieved
C              = 1 if convergence not achieved after ITERMX steps
C              = 2 if at most one group contains observations
C              = 3 if DIVBEG and/or DIVEND values are invalid
C              = 4 if both NN=0 (grouped data) and GROUP=.FALSE.
C              = 5 if AGUESS = 0.0 or lies outside specified bounds
C                  for any key parameter
C              = 6 if NK<1 or NK>3
C              = 7 if NC<0 or NC>10
C              = 8 if ITERMX < 10
C              = 9 if QUICK=.TRUE. and NC=0
C              = 10 if NOTE contains unexpected zeros or invalid values
C              = 11 if values of NOTE are not ordered from smallest to
C                   largest or if not all values of NOTE are unique
C              = 12 if at least one correlation between parameter
C                   estimates is very close to unity
C              = 13 if fitted density is negative somewhere in the
C                   range (DIVBEG(1),DIVEND(NINTS))

       PARAMETER (KEYMAX=3,NMAX=13,INTMAX=85,NOBMAX=1000,EPS=1.0E-6,
     * JMAX=10)

C      KEYMAX    maximum number of key parameters
C      NMAX      maximum total number of parameters (key parameters and
C                adjustment terms)
C      INTMAX    maximum number of intervals (groups)
C      NOBMAX    maximum number of observations
C      EPS       determines the level of precision of numerical integration
C      JMAX      maximum number of iterations for numerical integration

       REAL DIVBEG(*),DIVEND(*),Y(*),AGUESS(*),DBOUND(*),UBOUND(*),
     * Y2(NOBMAX),YFN(0:50)

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION DFNDA(NMAX),A(*),A2(NMAX),A3(NMAX),
     * A4(KEYMAX),ALOGL,ALOGL2,AMU,AINFO(NMAX,NMAX),AINV(NMAX,*),
     * DLLDA(NMAX),DPGDA(INTMAX,NMAX),DMUDA(*),DMUDA2(NMAX),
     * SEA(*),CORR(NMAX,*),B(NMAX),P(*),ALL,FNEV,DPDA(INTMAX,NMAX),
     * DLFDA(NMAX),DAREA(NMAX),DAREA2(NMAX),ODAREA(NMAX),SUMD(NMAX),
     * BMAT(NMAX,NMAX),BINV(NMAX,NMAX),DKEYDA(KEYMAX),DXSDA(KEYMAX),
     * BB(NMAX,NMAX)
       INTEGER NOBS(*),NOTE(*),INDX(NMAX)
       LOGICAL GROUP,KIND,QUICK,BOUNDD(*),BOUNDU(*)
 
       DATA SMALL1,SMALL2,P999,BIG1/1.0E-3,1.0E-8,9.99E-1,1.0E10/

C      CALL GETTIM(IHR,IMIN,ISEC,I100TH)
C      TSTRT=(IHR*3600)+(IMIN*60)+ISEC+I100TH/100.0
       IFAULT=0

C      Check for errors in the input parameters

       DO 10 J=1,NINTS
          IF(DIVEND(J).LE.DIVBEG(J)) IFAULT=3

          IF(J.NE.NINTS) THEN
             IF(DIVEND(J+1).LE.DIVEND(J)) IFAULT=3
             IF(DIVBEG(J+1).LE.DIVBEG(J)) IFAULT=3
             IF(DIVEND(J).GT.DIVBEG(J+1)+SMALL2) IFAULT=3
             IF(.NOT.GROUP.AND.DIVEND(J).LT.DIVBEG(J+1)-SMALL2) IFAULT=3
          ENDIF

 10    CONTINUE

       IF(.NOT.GROUP.AND.NN.EQ.0) IFAULT=4
       IF(NK.LT.1.OR.NK.GT.3) IFAULT=6
       IF(NC.LT.0.OR.NC.GT.10) IFAULT=7
       IF(ITERMX.LT.10) IFAULT=8
       IF(NC.EQ.0.AND.QUICK) IFAULT=9

       IF(NC.GT.0) THEN

          DO 20 J=1,NC
             IF(NOTE(J).LE.0.OR.NOTE(J).GT.10) IFAULT=10
             IF(J.GT.1.AND.NOTE(J).LE.NOTE(J-1)) IFAULT=11
 20       CONTINUE

       ENDIF

       DO 30 J=1,NK
          IF(ABS(AGUESS(J)).LT.SMALL2) IFAULT=5
          IF(BOUNDD(J).AND.AGUESS(J).LT.DBOUND(J)) IFAULT=5
          IF(BOUNDU(J).AND.AGUESS(J).GT.UBOUND(J)) IFAULT=5
 30    CONTINUE

       IF(IFAULT.NE.0) RETURN

       IF(NN.GT.0) THEN

C         Calculate group frequencies

          DO 40 J=1,NINTS
 40       NOBS(J)=0

          N=0

          DO 60 I=1,NN

             IF(Y(I).GT.DIVBEG(1)-SMALL2) THEN

                DO 50 J=1,NINTS

                   IF(Y(I).LE.DIVEND(J)+SMALL2) THEN
                      N=N+1
                      Y2(N)=Y(I)
                      NOBS(J)=NOBS(J)+1
                      GOTO 60
                   ENDIF

 50             CONTINUE

             ENDIF

 60       CONTINUE

       ELSE
          N=0

          DO 65 J=1,NINTS
 65       N=N+NOBS(J)

       ENDIF

       NZ=0

       DO 70 J=1,NINTS
          IF(NOBS(J).EQ.0) NZ=NZ+1
 70    CONTINUE

       IF(NZ.GE.NINTS-1) THEN
          IFAULT=2
          RETURN
       ENDIF

       NP=NK+NC
       KOUNT=0
       ITER2=0
       ITER3=0
       ICON=0

       IF(GROUP) THEN

C         Calculate the constant part of the log likelihood function

          CONST=0.0
          K=0

          DO 90 J=2,NINTS
             K=K+NOBS(J-1)

             DO 80 I=1,NOBS(J)
 80          CONST=CONST+LOG(FLOAT(I+K))-LOG(FLOAT(I))

 90       CONTINUE

       ENDIF

       DO 110 J=1,NP

          IF(J.LE.NK) A4(J)=-1.0

          DO 100 K=1,NP
 100      CORR(J,K)=0.0

          B(J)=AGUESS(J)
 110   CONTINUE

       IF(.NOT.QUICK) THEN

C         In MODE 2, the search is conditional on the current key parameter
C         estimates, and only the polynomial coefficients are updated

          MODE=2
       ELSE

C         In MODE 4, only the key parameters are estimated

          MODE=4
       ENDIF

       M=NP
       INIT=0

C      This initialising ensures that erroneous 'convergence'
C      cannot occur at the 1st iteration

       DO 120 J=1,M
 120   A3(J)=-1.0

       ALOGL2=-BIG1

C      If there are no polynomial coefficients to be fitted, only a
C      single set of iterations to fit the key parameters is required

       IF(M.EQ.NK) MODE=3

       IF(QUICK) THEN
          MOLD=M
          M=NK
       ENDIF

       M2=M

 130   ITER=0

 140   ALOGL=ALL(CONST,P,NINTS,DIVBEG,DIVEND,B,DPGDA,AMU,DMUDA,DAREA,
     * DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,DXSDA,M2,NOTE,NK,N,Y2,NOBS,
     * KIND,GROUP,INTMAX,EPS,JMAX)

C      If the log likelihood has decreased, enter the Marquardt
C      procedure to find a point at least as good as the previous one

       IF(ALOGL.LT.ALOGL2) THEN
          CALL MARQUA(B,A2,AINFO,DLLDA,NK,MM,MODE,KOUNT,NMAX,BMAT,BINV,
     *    BB,INDX)
          GOTO 230
       ENDIF

       INIT=1

C      The following subroutine call allows progress in the search to
C      be monitored.  The subroutine is only entered if the likelihood
C      has increased or stayed the same relative to the previous call.
C      If the user wishes to monitor progress, he/she should write a
C      short subroutine MONITR with the arguments listed below that
C      simply prints each argument to a file or to the screen.  Note
C      that, if double precision is used, the first line following the
C      SUBROUTINE statement should be  DOUBLE PRECISION ALOGL,B(*)
C      For single precision, replace DOUBLE PRECISION by REAL

C      ALOGL is the current value of the likelihood
C      MODE is the current mode of search
C      ITER is the number of iterations in the current mode
C      ITER3 is the total number of iterations
C      B(J), J=1,M, is the vector of estimated parameters and 
C      polynomial coefficients

       CALL MONITR(ALOGL,MODE,ITER,ITER3,M,B)

C      CALL GETTIM(IHR,IMIN,ISEC,I100TH)
C      TIME=(IHR*3600)+(IMIN*60)+ISEC+I100TH/100.0
C      PRINT *,' Time =',TIME-TSTRT
       ITER=ITER+1
       ITER3=ITER3+1

       IF(ICON.EQ.1) ITER2=ITER2+1
       IF(ITER.EQ.1) GOTO 210
 
       IF(MODE.EQ.3.OR.QUICK) THEN

C         Has convergence occurred?  If so, GOTO statement 250.
C         If not, return to statement 210

          IF(M.GT.NK) THEN

             DO 150 J=NK+1,M
                IF(ABS(B(J)-A3(J))/(SMALL1+ABS(A3(J))).GT.
     *          SMALL1/SQRT(FLOAT(N))) GOTO 210
 150         CONTINUE

          ENDIF

          IF(MODE.NE.2) THEN

             DO 160 J=1,NK
                IF(ABS(B(J)-A3(J))/(SMALL1+ABS(A3(J))).GT.
     *          SMALL1/SQRT(FLOAT(N))) GOTO 210
 160         CONTINUE

          ENDIF

          GOTO 250
       ENDIF

C      Update parameter estimates

       DO 170 J=1,M
 170   A(J)=B(J)

       IF(MODE.EQ.1) THEN

C         Search alternates between MODE 1 (in which the key parameters
C         are fitted conditional on the current polynomial coefficient
C         estimates) and MODE 2 (the converse).  If MODE 3 has already
C         been entered once, it is automatically returned to after 
C         several cycles, provided convergence seems close.

          IF(ITER2.LE.3.AND.ICON.EQ.1) THEN
             MODE=2
             GOTO 190
          ENDIF

C         If convergence is close, enter MODE 3

          DO 180 J=1,NK

             IF(ABS(B(J)-A4(J))/(SMALL1+ABS(A4(J))).GT.SMALL1) THEN
                MODE=2
                GOTO 190
             ENDIF

 180      CONTINUE

          MODE=3
          ITER=0
          GOTO 210

C         Array A4 is used to test for convergence when the search is
C         conditional on the polynomial coefficient estimates (MODE 1)

 190      DO 200 J=1,NK
 200      A4(J)=B(J)

       ELSEIF(.NOT.QUICK.AND.MODE.EQ.2) THEN

C         After one step in MODE 2, switch to MODE 1

          MODE=1
       ENDIF

 210   IF(ITER3.GT.ITERMX) THEN
          IFAULT=1
          RETURN
       ENDIF

C      If 12 iterations fail to bring convergence in MODE 3, and there is
C      at least one polynomial coefficient to be fitted, revert to
C      alternation between MODEs 1 and 2 for several cycles.

       IF(ITER.GE.12.AND.MODE.EQ.3.AND.NP.GT.NK) THEN
          MODE=2
          ITER2=0
          ICON=1
       ENDIF

C      Initialise variables for Newton-Raphson search

       ALOGL2=ALOGL

       DO 220 J=1,M
 220   A3(J)=B(J)

       KOUNT=0

C      MM is the number of parameters to be fitted at this stage of the search

       MM=NK
       IF(MODE.EQ.2) MM=M-NK
       IF(MODE.EQ.3) MM=M

       CALL INFO(P,NINTS,B,AINFO,DLLDA,DPGDA,AMU,DMUDA,DPDA,DLFDA,
     * DFNDA,DKEYDA,DXSDA,M,MM,NOTE,NK,N,Y2,NOBS,KIND,GROUP,MODE,
     * NMAX,INTMAX)

       CALL NRAPH(B,A2,AINFO,AINV,DLLDA,NK,MM,MODE,BB,INDX,NMAX)

C      If search exceeds bounds on the parameters, enter Marquardt routine

 230   DO 240 J=1,NK

          IF((BOUNDD(J).AND.B(J).LT.DBOUND(J)).OR.(BOUNDU(J).AND.
     *    B(J).GT.UBOUND(J))) THEN
             CALL MARQUA(B,A2,AINFO,DLLDA,NK,MM,MODE,KOUNT,NMAX,
     *       BMAT,BINV,BB,INDX)
             GOTO 230
          ENDIF

 240   CONTINUE

C      Otherwise, continue to the next iteration of Newton-Raphson

       GOTO 140

C      After convergence has occurred:

 250   IF(MODE.EQ.4) THEN

C         Store the results from fitting the key function, then GOTO
C         statement 130 to fit the polynomial coefficients

          DO 270 K=1,M

             DO 260 L=1,M
 260         CORR(K,L)=AINV(K,L)

             IF(AINV(K,K).LT.0.0) AINV(K,K)=1.0E-20
             SEA(K)=SQRT(AINV(K,K))
 270      CONTINUE

          M=MOLD
          M2=M
          MODE=2

          DO 280 J=1,NK
 280      A(J)=B(J)

          GOTO 130
       ENDIF

       IF(QUICK) THEN
          MINUS=NK
       ELSE
          MINUS=0
       ENDIF

       CHISQ=0.0

C      Evaluate chi-square goodness-of-fit statistic;  if analysis is of
C      ungrouped data, first evaluate areas under the estimated probabilty
C      density function corresponding to the divisions defined by the
C      cutpoints DIVEND and DIVBEG

       DO 290 J=1,NINTS

          IF(.NOT.GROUP) THEN

             CALL AREAS(P(J),DMUDA2,DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,
     *       DXSDA,DIVBEG(J),DIVEND(J),B,M,NOTE,NK,KIND,EPS,JMAX)

             P(J)=P(J)/AMU
          ENDIF

          PJN=P(J)*N
          CHISQ=CHISQ+(NOBS(J)-PJN)**2/ABS(PJN)
 290   CONTINUE

C      Degrees of freedom for chi-square test

       NDF=NINTS-1-M

C      Evaluations of the fitted probability density function at 51
C      equally spaced points;  may be used for plotting

       DO 300 I=0,50
          YVAL=DIVBEG(1)+I*(DIVEND(NINTS)-DIVBEG(1))/50.0
          CALL FNEVAL(YVAL,B,M2,NOTE,NK,KIND,FNEV,DFNDA,DKEYDA,DXSDA)
          YFN(I)=FNEV/AMU
          IF(FNEV.LT.0.0) IFAULT=13
 300   CONTINUE

C      Standard errors of the parameter estimates

       DO 310 K=1,M-MINUS
          IF(AINV(K,K).LT.0.0) AINV(K,K)=1.0E-20
          SEA(K+MINUS)=SQRT(AINV(K,K))
 310   CONTINUE

C      Transfer of parameter estimates from array B to A, and corresponding
C      reassignment of the elements of the inverse of the information matrix
C      and of the estimated derivatives of the normalising function (DMUDA)
C      wrt each parameter

       DO 330 K=1,M-MINUS

          DO 320 L=1,M-MINUS
 320      CORR(K+MINUS,L+MINUS)=AINV(K,L)

 330   CONTINUE

       DO 350 K=MINUS+1,M
          A(K)=B(K)

          DO 340 L=MINUS+1,M
 340      AINV(K,L)=CORR(K,L)

 350   CONTINUE

       IF(MINUS.GT.0) THEN

          DO 370 K=1,MINUS

             DO 360 L=1,MINUS
 360         AINV(K,L)=CORR(K,L)

 370      CONTINUE

          DO 390 K=1,MINUS

             DO 380 L=MINUS+1,M
                AINV(K,L)=0.0
                AINV(L,K)=0.0
 380         CONTINUE

 390      CONTINUE

       ENDIF

C      Calculation of correlations between parameter estimates

       IF(MINUS.GT.0) THEN
          CORR(MINUS,MINUS)=1.0

          IF(MINUS.GT.1) THEN

             DO 410 K=1,MINUS-1
                CORR(K,K)=1.0

                DO 400 L=K+1,MINUS
                   CORR(K,L)=AINV(K,L)/SEA(K)/SEA(L)
                   CORR(L,K)=CORR(K,L)
                   IF(ABS(CORR(K,L)).GT.P999) IFAULT=12
 400            CONTINUE

 410         CONTINUE

          ENDIF

       ENDIF

       CORR(M,M)=1.0
       IF(M.LE.MINUS+1) RETURN

       DO 430 K=MINUS+1,M-1
          CORR(K,K)=1.0

          DO 420 L=K+1,M
             CORR(K,L)=AINV(K,L)/SEA(K)/SEA(L)
             CORR(L,K)=CORR(K,L)
             IF(ABS(CORR(K,L)).GT.P999) IFAULT=12
 420      CONTINUE

 430   CONTINUE

       RETURN
       END

C**************************************************************

       SUBROUTINE AREAS(AREA,DAREA,DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,
     * DXSDA,DIV1,DIV2,A,M2,NOTE,NK,KIND,EPS,JMAX)

C      Calculation of area under the fitted curve and under the derivative
C      of the fitted curve wrt each parameter between DIV1 and DIV2
C      using Simpson's rule for numerical integration.  The number of 
C      intervals is increased until a level of precision, determined
C      by EPS, is achieved.  The method is described by Press et al. (1986)

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION AREA,DFNDA(*),A(*),DAREA(*),OAREA,OAREAT,
     * AREA2,DAREA2(*),ODAREA(*),SUM,SUMD(*),FNEV,DKEYDA(*),DXSDA(*)
       INTEGER NOTE(*)
       LOGICAL KIND

       OAREAT=-1.0E30
       OAREA=-1.0E30

       DO 10 K=1,M2
 10    ODAREA(K)=-1.0E30

       DO 100 J=1,JMAX

          IF(J.EQ.1) THEN
             CALL FNEVAL(DIV1,A,M2,NOTE,NK,KIND,FNEV,DFNDA,DKEYDA,
     *       DXSDA)
             AREA2=FNEV

             DO 20 K=1,M2
 20          DAREA2(K)=DFNDA(K)

             CALL FNEVAL(DIV2,A,M2,NOTE,NK,KIND,FNEV,DFNDA,DKEYDA,
     *       DXSDA)
             CON=0.5*(DIV2-DIV1)
             AREA2=CON*(AREA2+FNEV)

             DO 30 K=1,M2
 30          DAREA2(K)=CON*(DAREA2(K)+DFNDA(K))

             IT=1
          ELSE
             TNM=IT
             DIFF=(DIV2-DIV1)/TNM
             ARG=DIV1+0.5*DIFF
             SUM=0.0

             DO 40 K=1,M2
 40          SUMD(K)=0.0

             DO 60 JJ=1,IT
                CALL FNEVAL(ARG,A,M2,NOTE,NK,KIND,FNEV,DFNDA,DKEYDA,
     *          DXSDA)
                SUM=SUM+FNEV

                DO 50 K=1,M2
 50             SUMD(K)=SUMD(K)+DFNDA(K)

                ARG=ARG+DIFF
 60          CONTINUE

             AREA2=0.5*(AREA2+(DIV2-DIV1)*SUM/TNM)

             DO 70 K=1,M2
 70          DAREA2(K)=0.5*(DAREA2(K)+(DIV2-DIV1)*SUMD(K)/TNM)

             IT=2*IT
          ENDIF

          AREA=(4.0*AREA2-OAREAT)/3.0

          DO 80 K=1,M2
 80       DAREA(K)=(4.0*DAREA2(K)-ODAREA(K))/3.0

          IF(ABS(AREA-OAREA).LT.EPS*ABS(AREA+OAREA)) RETURN

          OAREA=AREA
          OAREAT=AREA2

          DO 90 K=1,M2
 90       ODAREA(K)=DAREA2(K)

 100   CONTINUE

       RETURN
       END

C***************************************************************

       SUBROUTINE PROBS(P,NINTS,DIVBEG,DIVEND,A,PTOTAL,DPGDA,DMUDA,
     * DAREA,DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,DXSDA,M2,NOTE,NK,KIND,
     * INTMAX,EPS,JMAX)

C      Evaluation of the area under each section of the fitted
C      probability density function

       REAL DIVBEG(*),DIVEND(*)

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION A(*),DAREA(*),DMUDA(*),AREA,PTOTAL,P(*),
     * DPGDA(INTMAX,*),DFNDA(*),DAREA2(*),ODAREA(*),SUMD(*),
     * DKEYDA(*),DXSDA(*)
       LOGICAL KIND
       INTEGER NOTE(*)

       PTOTAL=0.0

       DO 10 K=1,M2
 10    DMUDA(K)=0.0

       DO 30 J=1,NINTS

C         Calculation of the area under section J of the fitted curve
C         before normalisation

          DIV2=DIVEND(J)
          DIV1=DIVBEG(J)

          CALL AREAS(AREA,DAREA,DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,
     *    DXSDA,DIV1,DIV2,A,M2,NOTE,NK,KIND,EPS,JMAX)

          P(J)=AREA
          PTOTAL=PTOTAL+P(J)

          DO 20 K=1,M2
             DPGDA(J,K)=DAREA(K)
             DMUDA(K)=DMUDA(K)+DPGDA(J,K)
 20       CONTINUE

 30    CONTINUE

C      PTOTAL is now the value that normalises the fitted curve, to make
C      it a valid density function.  DPGDA is the matrix of estimated
C      derivatives of each of the NINTS non-normalised areas wrt the
C      parameters, and DMUDA is the estimated derivative of the
C      normalising value wrt the parameters.  The following DO loop
C      normalises the P's so that they sum to unity.

       DO 40 J=1,NINTS
          IF(PTOTAL.NE.0.0) P(J)=P(J)/PTOTAL
 40    CONTINUE

       RETURN
       END

C***************************************************************

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION FUNCTION ALL(CONST,P,NINTS,DIVBEG,DIVEND,B,
     * DPGDA,AMU,DMUDA,DAREA,DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,DXSDA,M2,
     * NOTE,NK,N,Y,NOBS,KIND,GROUP,INTMAX,EPS,JMAX)

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION P(*),B(*),AMU,DPGDA(INTMAX,*),DMUDA(*),
     * FNEV,DFNDA(*),DAREA2(*),ODAREA(*),SUMD(*),DKEYDA(*),DXSDA(*),
     * DAREA(*)
       REAL DIVBEG(*),DIVEND(*),Y(*)
       INTEGER NOTE(*),NOBS(*)
       LOGICAL KIND,GROUP
       DATA SMALL/1.0E-10/

       IF(GROUP) THEN

C         Grouped data:

          CALL PROBS(P,NINTS,DIVBEG,DIVEND,B,AMU,DPGDA,DMUDA,DAREA,
     *    DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,DXSDA,M2,NOTE,NK,KIND,
     *    INTMAX,EPS,JMAX)

          ALL=CONST

C         Add the part of the log likelihood that varies with
C         the parameter estimates to the constant part

          DO 10 J=1,NINTS

C            If an estimated probability < or = 0, set it equal
C            to a very small value

             IF(P(J).LT.SMALL) P(J)=SMALL
             ALL=ALL+NOBS(J)*LOG(P(J))
 10       CONTINUE

       ELSE

C         Ungrouped data:

          CALL AREAS(AMU,DMUDA,DFNDA,DAREA2,ODAREA,SUMD,DKEYDA,
     *    DXSDA,DIVBEG(1),DIVEND(NINTS),B,M2,NOTE,NK,KIND,EPS,JMAX)

C         Evaluate the log likelihood;  if the estimated fit is
C         impossible, give the log likelihood a large -ve value

          IF(AMU.LE.0.0) AMU=1.0E25
          ALL=-N*LOG(AMU)

          DO 20 J=1,N

             CALL FNEVAL(Y(J),B,M2,NOTE,NK,KIND,FNEV,DFNDA,DKEYDA,
     *       DXSDA)

             IF(FNEV.LT.SMALL) FNEV=SMALL
             ALL=ALL+LOG(FNEV)
 20      CONTINUE

       ENDIF

       RETURN
       END

C***************************************************************

       SUBROUTINE INFO(P,NINTS,B,AINFO,DLLDA,DPGDA,AMU,DMUDA,DPDA,
     * DLFDA,DFNDA,DKEYDA,DXSDA,M,MM,NOTE,NK,N,Y,NOBS,KIND,GROUP,MODE,
     * NMAX,INTMAX)

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION AMU,AINFO(NMAX,*),DPGDA(INTMAX,*),DMUDA(*),
     * DPDA(INTMAX,*),DLFDA(*),DFNDA(*),P(*),B(*),DLLDA(*),FNEV,
     * DKEYDA(*),DXSDA(*)
       INTEGER NOTE(*),NOBS(*)
       REAL Y(*)
       LOGICAL KIND,GROUP

       DO 20 K=1,M
          DLLDA(K)=0.0

          DO 10 L=1,M
 10       AINFO(K,L)=0.0

 20    CONTINUE

       IF(GROUP) THEN

C         Grouped data:

C         Evaluate derivative wrt each parameter to be fitted of the area
C         under each section of the probability density function

          DO 40 K=1,M

             IF((MODE.EQ.1.AND.K.LE.NK).OR.(MODE.EQ.2.AND.K.GT.NK).
     *       OR.MODE.GT.2) THEN

                DO 30 J=1,NINTS
                   DPDA(J,K)=(DPGDA(J,K)-P(J)*DMUDA(K))/AMU
 30             CONTINUE

             ENDIF

 40      CONTINUE

C         Evaluate derivative of the log likelihood function
C         wrt each parameter to be fitted

          DO 60 K=1,M

             IF((MODE.EQ.1.AND.K.LE.NK).OR.(MODE.EQ.2.AND.K.GT.NK).
     *       OR.MODE.GT.2) THEN

                DO 50 J=1,NINTS
 50             DLLDA(K)=DLLDA(K)+NOBS(J)*DPDA(J,K)/P(J)

             ENDIF

 60       CONTINUE

C         Calculate the information matrix

          DO 90 K=1,MM
             KK=K
             IF(MODE.EQ.2) KK=K+NK
 
             DO 80 L=1,MM
                LL=L
                IF(MODE.EQ.2) LL=L+NK

                DO 70 J=1,NINTS
 70             AINFO(K,L)=AINFO(K,L)+N*DPDA(J,KK)*DPDA(J,LL)/P(J)

 80          CONTINUE

 90       CONTINUE

       ELSE

C         Ungrouped data:

          DO 130 J=1,N

C            For each observation, calculate its contribution to the
C            derivatives wrt each parameter to be fitted of the log of
C            the probability density function and of the log likelihood.

             CALL FNEVAL(Y(J),B,M,NOTE,NK,KIND,FNEV,DFNDA,DKEYDA,
     *       DXSDA)

             DO 100 K=1,M

                IF((MODE.EQ.1.AND.K.LE.NK).OR.(MODE.EQ.2.AND.K.GT.NK).
     *          OR.MODE.GT.2) THEN
                   DLFDA(K)=DFNDA(K)/FNEV -DMUDA(K)/AMU
                   DLLDA(K)=DLLDA(K)+DLFDA(K)
                ENDIF

 100         CONTINUE

C            Calculate the contribution of each observation to the
C            information matrix

             DO 120 K=1,MM
                KK=K
                IF(MODE.EQ.2) KK=K+NK

                DO 110 L=1,MM
                   LL=L
                   IF(MODE.EQ.2) LL=L+NK
                   AINFO(K,L)=AINFO(K,L)+DLFDA(KK)*DLFDA(LL)
 110            CONTINUE

 120         CONTINUE

 130      CONTINUE

       ENDIF

       RETURN
       END

C***************************************************************

       SUBROUTINE NRAPH(B,A2,AINFO,AINV,DLLDA,NK,MM,MODE,BB,INDX,
     * NMAX)

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION B(*),A2(*),AINFO(NMAX,*),AINV(NMAX,*),
     * DLLDA(*),BB(NMAX,*)
       INTEGER INDX(*)

C      Invert the information matrix

       IF(MM.EQ.1) THEN
          AINV(1,1)=1.0/AINFO(1,1)
       ELSE
          CALL INVERS(AINFO,AINV,BB,INDX,MM,NMAX)
       ENDIF

C      For each parameter being fitted, update its estimate, using
C      Newton-Raphson

       DO 20 K=1,MM
          KK=K
          IF(MODE.EQ.2) KK=K+NK
          A2(KK)=B(KK)

          DO 10 L=1,MM
             LL=L
             IF(MODE.EQ.2) LL=L+NK
             B(KK)=B(KK)+AINV(K,L)*DLLDA(LL)
 10       CONTINUE

 20    CONTINUE

       RETURN
       END

C***************************************************************

       SUBROUTINE MARQUA(B,A2,AINFO,DLLDA,NK,MM,MODE,KOUNT,NMAX,
     * BMAT,BINV,BB,INDX)

C      The Marquardt procedure, for when Newton-Raphson fails to increase
C      the value of the log likelihood.  The procedure is very similar to
C      Newton-Raphson, but progressively adds larger values to the
C      diagonal elements of the information matrix until a point with a
C      larger log likelihood is found.

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION B(*),A2(*),AINFO(NMAX,*),DLLDA(*),
     * BMAT(NMAX,*),BINV(NMAX,*),BB(NMAX,*)
       INTEGER INDX(*)

       KOUNT=KOUNT+1

       DO 20 K=1,MM

          DO 10 L=1,MM
             BMAT(K,L)=AINFO(K,L)
             IF(K.EQ.L) BMAT(K,L)=BMAT(K,L)*(1.0+2.0**KOUNT)
 10       CONTINUE

 20    CONTINUE

       IF(MM.EQ.1) THEN
          BINV(1,1)=1.0/BMAT(1,1)
       ELSE
          CALL INVERS(BMAT,BINV,BB,INDX,MM,NMAX)
       ENDIF

       DO 40 K=1,MM
          KK=K
          IF(MODE.EQ.2) KK=K+NK
          B(KK)=A2(KK)

          DO 30 L=1,MM
             LL=L
             IF(MODE.EQ.2) LL=L+NK
             B(KK)=B(KK)+BINV(K,L)*DLLDA(LL)
 30       CONTINUE

 40    CONTINUE

       RETURN
       END

C***************************************************************

       SUBROUTINE POLY(XS,H)

C      Simple polynomials up to order 10

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION H(0:10),XS

       H(0)=1.0
       H(1)=XS
       H(2)=XS*XS
       H(3)=XS**3
       H(4)=XS**4
       H(5)=XS**5
       H(6)=XS**6
       H(7)=XS**7
       H(8)=XS**8
       H(9)=XS**9
       H(10)=XS**10

       RETURN
       END

C***************************************************************

       SUBROUTINE HERM(XS,H)

C      Hermite polynomials up to order 10

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION H(0:10),XS

       H(0)=1.0
       H(1)=XS
       H(2)=XS*XS-1.0
       H(3)=XS**3-3.0*XS
       H(4)=XS**4-6.0*XS*XS+3.0
       H(5)=XS**5-10.0*XS**3+15.0*XS
       H(6)=XS**6-15.0*XS**4+45.0*XS*XS-15.0
       H(7)=XS**7-21.0*XS**5+105.0*XS**3-105.0*XS
       H(8)=XS**8-28.0*XS**6+210.0*XS**4-420.0*XS*XS+105.0
       H(9)=XS**9-36.0*XS**7+378.0*XS**5-1260.0*XS**3+945.0*XS
       H(10)=XS**10-45.0*XS**8+630.0*XS**6-3150.0*XS**4+4725.0*XS*XS
     * -945.0

       RETURN
       END

C***************************************************************

       SUBROUTINE FNEVAL(X,A,M2,NOTE,NK,KIND,FNEV,DFNDA,
     * DKEYDA,DXSDA)

C      Evaluates the function and its derivative wrt A (before
C      normalisation) at point X, using parameter estimates A

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION A(*),DFNDA(*),DKEYDA(*),XS,DXSDA(*),
     * VALKEY,H(0:10),FNEV,TOP,DTOP
       INTEGER NOTE(*)
       LOGICAL KIND

       CALL STD(X,A,XS,DXSDA)

       TOP=1.0

C      VALKEY is the value of the key function

       CALL KEY(X,A,VALKEY,DKEYDA)

       IF(M2.GT.NK) THEN

          IF(KIND) THEN
             CALL HERM(XS,H)
          ELSE
             CALL POLY(XS,H)
          ENDIF

C         TOP is the contribution of the polynomial terms to the fitted
C         density, apart from the normalising value, which is a function
C         of the parameters alone

          DO 10 J=NK+1,M2
 10       TOP=TOP+A(J)*H(NOTE(J-NK))

       ENDIF

       FNEV=VALKEY*TOP
 
C      Evaluate the derivative of the curve (before normalisation) at point
C      X, wrt each parameter, using parameter estimates A

       TOP=1.0
       DTOP=0.0

       IF(M2.GT.NK) THEN

          DO 20 J=NK+1,M2
 20       DFNDA(J)=VALKEY*H(NOTE(J-NK))

          DO 30 J=NK+1,M2
             TOP=TOP+A(J)*H(NOTE(J-NK))
             DTOP=DTOP+A(J)*(NOTE(J-NK))*H(NOTE(J-NK)-1)
 30       CONTINUE

       ENDIF

       DO 40 J=1,NK
 40    DFNDA(J)=VALKEY*DTOP*DXSDA(J)+DKEYDA(J)*TOP

       RETURN
       END

C***********************************************************

       SUBROUTINE INVERS(A,Y,B,INDX,N,NMAX)

C      Calculates the inverse of a matrix.
C      See 'Numerical Recipes', Press et al., 1986, p38.

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION A,B,Y
       DIMENSION A(NMAX,*),B(NMAX,*),Y(NMAX,*),INDX(*)

       DO 20 I=1,N

          DO 10 J=1,N
             B(I,J)=A(I,J)
             Y(I,J)=0.0
 10       CONTINUE

          Y(I,I)=1.0
 20    CONTINUE

       CALL LUDCMP(B,N,NMAX,INDX,D)

       DO 30 J=1,N
          CALL LUBKSB(B,N,NMAX,INDX,Y(1,J))
 30    CONTINUE

       RETURN
       END

C****************************************************

       SUBROUTINE STD(X,A,XS,DXSDA)

C      Reduces observations to standard measure, to avoid numerical
C      problems, and differentiates the standardised observation wrt A

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION A(*),XS,DXSDA(*)

       XS=(X-A(1))/A(2)
       DXSDA(1)=-1.0/A(2)
       DXSDA(2)=-(X-A(1))/(A(2)**2)

       RETURN
       END

C****************************************************

       SUBROUTINE KEY(X,A,VALKEY,DKEYDA)

C      Evaluates the fitted key function and its derivative wrt A at X

C      On a 32 bit machine, change REAL to DOUBLE PRECISION below:

       DOUBLE PRECISION A(*),VALKEY,DKEYDA(*)

       DATA SMALL2,SMALL3/1.0E-8,1.0E-10/

       VALKEY=1.0
       D=0.0
       IF(ABS(X-A(1)).GT.SMALL3) D=((X-A(1))/A(2))**2/2.0

       IF(D.GT.20.0) THEN
          VALKEY=0.0
       ELSEIF(D.GT.SMALL2) THEN
          VALKEY=EXP(-D)
       ENDIF

C      Estimate the differential at X of the key function wrt A 

       DKEYDA(1)=(X-A(1))*VALKEY/(A(2)**2)
       DKEYDA(2)=((X-A(1))**2)*VALKEY/(A(2)**3)

       RETURN
       END
