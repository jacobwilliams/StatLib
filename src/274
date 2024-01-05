Notes on AS 274 "Least-squares routines to supplement those of Gentleman"
by Alan Miller.

This algorithm provides a set of high accuracy least squares routines which
expand upon those provided by Gentleman in AS 75.   In particular, it includes
facilities for (weighted) least-squares for a subset of the variables, for
changing the order of the variables, for testing for singularities, and for
calculating an estimated covariance matrix for the regression coefficients.

The algorithm is NOT consistent with those which have been published in the
same journal by Clarke (AS 163), Stirling (AS 164) or Smith (AS 268), in that
the orthogonal reduction is stored in a different way.   It is unfortunate
that these authors have not used the same format as Gentleman.

The basic algorithm is the same as Gentleman's.   As each new case is added,
the orthogonal reduction is updated.   In traditional least-squares notation,
the algorithm goes straight from the X-matrix (the 'design' matrix) to the
Cholesky factorization WITHOUT the intermediate step of forming a sum of
squares and products matrix.   Thus it avoids squaring the condition number.
To be pedantic, it is actually the Banachiewicz factorization which is used.
We use:

		  Q X = sqrt(D) R

where Q is an orthonormal matrix such that Q'Q = I, D is a diagonal matrix
(sqrt(D) is my crude way of indicating a diagonal matrix with elements on the
diagonal which are the square roots of the elements stored in array D in the
routines), and R is an upper triangular matrix with 1's on the diagonal.   The
array R is stored in a 1-dimensional array with elements stored by rows as
shown below (the 1's on the diagonal are NOT stored).

Storage of R:    (1)  1   2   3   4   5
		     (1)  6   7   8   9
		         (1) 10  11  12
			     (1) 13  14
			         (1) 15
				     (1)
N.B. The matrix Q is a product of planar-rotation matrices (Jacobi/Givens
rotations); it is not calculated or stored.
N.B. If you want to fit a constant (intercept) in your model, you can do it
by putting the first element of each row of X equal to 1.0.
N.B. You may like to dimension XROW as XROW(0:K) in your calling program,
where K is the number of variables (other than the constant).   The parameter
NP passed to routine INCLUD (and others) then takes the value K+1.

To perform the calculation of regression coefficients (and nothing else), you
need to call the following routines in the order given:
1. clear           - initializes arrays
2. includ          - call once for each case to update the orthogonal reduction
3. tolset          - calculates tolerances for each variable
4. sing (optional) - tests for singularities
5. regcf           - calculates regression coefficients for the first NREQ
		     variables

You can then add extra data and recalculate the regression coefficients.

Other routines (all require that 1 and 2 above have been used first).
ss    calculates array RSS such that RSS(I) = the residual sum of squares with
      the first I variables in the model
cov   calculates the estimated covariance matrix of the regression coefficients
      (the user must input a value for the residual variance VAR, usually the
      residual sum of squares SSERR divided by NOBS - NP for all variables,
      or RSS(I) divided by NOBS - I for the first I variables).
pcorr calculates partial correlations with the first IN variables in the model.
      If IN = 1 and the first variable is identically equal to 1 then these are
      the usual full correlations.   If IN = 0 it gives the cross products (not
      centered) divided by the square root of the product of the second moments
      (about zero) of the two variables.
vmove moves the variable in position FROM to position TO.
reordr calls vmove repeatedly to re-order the variables from the current order
      in integer array VORDER to that in the integer array LIST.
hdiag calculates the diagonal elements of the `hat' matrix used in various
      diagnostic statistics.

For vmove and reordr, the user must first set up an integer array VORDER which
assigns a unique integer value with each variable.   You may like to associate
the value 0 with the constant in the model.

For cov, vmove, reordr and hdiag you must call tolset first.

*** Warning *** Routine INCLUD (from Gentleman's AS 75) overwrites the elements
in array XROW.

Alan Miller
25 November 1991


Notes on handling singularities
-------------------------------
If there is a singularity amongst your X-variables and routine SING is not
used, you will get wrong answers from most routines.   The orthogonal reduction
can be written as:

	X  =  Q' sqrt(D) R

If we denote the columns of X as X1, X2, ..., Xk, and the columns of Q' as
q1, q2, ..., qk, then we have:

       x1  =  r(11).q1
       x2  =  r(12).q1  +  r(22).q2
       x3  =  r(13).q1  +  r(23).q2  + r(33).q3

etc., where r(ij) = the (i,j)-element of R multiplied by the square root of the
i-th diagonal element of D.   That is x1 is equal to the first orthogonal
direction, q1, multiplied by a scaling factor.   x2 = a projection of length
r(12) in the first direction, q1, plus a projection of length r(22) in a new
direction, q2, which is orthogonal to q1, etc.   If the X-directions are all
orthogonal then all of the r(ij)'s for off-diagonal elements of R will be zero.

Let us suppose that there is a singularity amongst the X-variables, that is
that one of the X-variables is exactly linearly related to some of the others.
For simplicity, let us suppose that x3 = a.x1 + b.x2.   The direction q3 is
formed from x3 by subtracting its projections in directions q1 and q2 and then
scaling so that it has length equal to 1.   When the projections in directions
q1 and q2 are subtracted in this case, there should be nothing left.   In
practice, there will be almost always be small rounding errors.   Direction
q3 will be formed from these rounding errors.   Thus we will have a rogue
direction in the columns of Q'.   The projections of the dependant variable
on this direction can be large.   It is like correlating the Y-variable on a
variable formed by using random numbers or a column of numbers from the phone
directory.

Routine SING sets the projections on this direction equal to zero.

In the case just mentioned, a full rank X could be obtained by removing any
one of the variables x1, x2 or x3.   In AS274, the equivalent to removing a
variable is to move it to the last position and then to only use the first
(k-1) variables in subsequent calculations.   You still tell routines such
as REGCF that there are NP variables, but you set NREQ = NP - 1.

If you just want properties of a subset of variables, then you use either VMOVE
or REORDR to re-order the variables so that the first ones are those which you
are interested in, and you set NREQ equal to the number of those variables.
In most cases, the constant (intercept) will be included in the model and so
will be included in the count of variables, NREQ (N-required).

Features NOT included in AS274
------------------------------
There is no facility for removing or adding variables (columns).   This could
be done but adding variables requires the storage of the Q-matrix, or at least
of the rotations used to form it.

There is no facility for two or more dependant variables, though this is fairly
easy to program.   The extra dependant variables can be treated as X-variables.
Thus if we have dependant variables Y1, Y2 and Y3, then Y2 and Y3 can be added
at the end of the list of X-variables.   When looking at the properties of Y1
related to the X-variables, use only the first (NP-2) columns by setting
NREQ = NP-2.   When you want to look at Y2, just copy column (NP-1) of RBAR
into THETAB - this requires a little thought to work out just where those
elements are stored.   You may want to first copy the original THETAB, or to
copy the column of RBAR into a different array, perhaps THETAB2 say.   The
residual sum of squares, SSERR, for variable Y2 is D(NP-1).

Y3 is a little more difficult.   You must first change the order of the
variables so that Y3 is before Y2, otherwise you will be regressing Y3 on Y2,
and you probably don't want to do that.

The easy way of thinking about the above is by considering each column of R
(after scaling by the square roots of the diagonal elements of D) as the
projections of the corresponding variable on all of those which have preceded
it.   The Y-variable is just another variable being projected upon a set of
orthogonal directions formed by the variables which came before it.   The
vector THETAB is stored separately just for fast access as most users
will have only one Y-variable; it could have been incorporated as the last
column of RBAR.

Alan Miller
2 May 1993


C-----------------------------------------------------------------------

      SUBROUTINE INCLUD(NP, NRBAR, WEIGHT, XROW, YELEM, D, RBAR, THETAB,
     +      SSERR, IER)
C
C     ALGORITHM AS274.1  APPL. STATIST. (1992) VOL 41, NO. 2
C
C     DOUBLE PRECISION VERSION
C
C     Calling this routine updates d, rbar, thetab and sserr by the
C     inclusion of xrow, yelem with the specified weight.
C     This version has been modified to make it slightly faster when the
C     early elements of XROW are not zeroes.
C
C     *** WARNING ***   The elements of XROW are over-written.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION WEIGHT, XROW(NP), YELEM, D(NP), RBAR(*),
     +       THETAB(NP), SSERR
C
C     Local variables
C
      INTEGER I, K, NEXTR
      DOUBLE PRECISION ZERO, W, Y, XI, DI, WXI, DPI, CBAR, SBAR, XK
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      W = WEIGHT
      Y = YELEM
      NEXTR = 1
      DO 30 I = 1, NP
C
C     Skip unnecessary transformations.   Test on exact zeroes must be
C     used or stability can be destroyed.
C
	IF (W .EQ. ZERO) RETURN
	XI = XROW(I)
	IF (XI .EQ. ZERO) THEN
	  NEXTR = NEXTR + NP - I
	  GO TO 30
	END IF
	DI = D(I)
	WXI = W * XI
	DPI = DI + WXI*XI
	CBAR = DI / DPI
	SBAR = WXI / DPI
	W = CBAR * W
	D(I) = DPI
	IF (I .EQ. NP) GO TO 20
	DO 10 K = I+1, NP
	  XK = XROW(K)
	  XROW(K) = XK - XI * RBAR(NEXTR)
	  RBAR(NEXTR) = CBAR * RBAR(NEXTR) + SBAR * XK
	  NEXTR = NEXTR + 1
   10   CONTINUE
   20   XK = Y
	Y = XK - XI * THETAB(I)
	THETAB(I) = CBAR * THETAB(I) + SBAR * XK
   30 CONTINUE
C
C     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
C     residual.
C
      SSERR = SSERR + W * Y * Y
C
      RETURN
      END
C
      SUBROUTINE CLEAR(NP, NRBAR, D, RBAR, THETAB, SSERR, IER)
C
C     ALGORITHM AS274.2  APPL. STATIST. (1992) VOL.41, NO.2
C
C     Sets arrays to zero prior to calling AS75.1
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), SSERR
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      DO 10 I = 1, NP
	D(I) = ZERO
	THETAB(I) = ZERO
   10 CONTINUE
      DO 20 I = 1, NRBAR
   20 RBAR(I) = ZERO
      SSERR = ZERO
      RETURN
      END
C
      SUBROUTINE REGCF(NP, NRBAR, D, RBAR, THETAB, TOL, BETA, NREQ,
     +     IER)
C
C     ALGORITHM AS274.3  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Modified version of AS75.4 to calculate regression coefficients
C     for the first NREQ variables, given an orthogonal reduction from
C     AS75.1.
C
      INTEGER NP, NRBAR, NREQ, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), TOL(NP), BETA(NP)
C
C     Local variables
C
      INTEGER I, J, NEXTR
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (NREQ .LT. 1 .OR. NREQ .GT. NP) IER = IER + 4
      IF (IER .NE. 0) RETURN
C
      DO 20 I = NREQ, 1, -1
	IF (SQRT(D(I)) .LT. TOL(I)) THEN
	  BETA(I) = ZERO
	  D(I) = ZERO
	  GO TO 20
	END IF
	BETA(I) = THETAB(I)
	NEXTR = (I-1) * (NP+NP-I)/2 + 1
	DO 10 J = I+1, NREQ
	  BETA(I) = BETA(I) - RBAR(NEXTR) * BETA(J)
	  NEXTR = NEXTR + 1
   10   CONTINUE
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE TOLSET(NP, NRBAR, D, RBAR, TOL, WORK, IER)
C
C     ALGORITHM AS274.4  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Sets up array TOL for testing for zeroes in an orthogonal
C     reduction formed using AS75.1.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(*), TOL(NP), WORK(NP)
C
C     Local variables.
C
      INTEGER COL, ROW, POS
      DOUBLE PRECISION EPS, SUM
C
C     EPS is a machine-dependent constant.   For compilers which use
C     the IEEE format for floating-point numbers, recommended values
C     are 1.E-06 for single precision and 1.E-15 for double precision.
C
      DATA EPS/1.D-15/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
C     Set TOL(I) = sum of absolute values in column I of RBAR after
C     scaling each element by the square root of its row multiplier.
C
      DO 10 COL = 1, NP
   10 WORK(COL) = SQRT(D(COL))
      DO 30 COL = 1, NP
	POS = COL - 1
	SUM = WORK(COL)
	DO 20 ROW = 1, COL-1
	  SUM = SUM + ABS(RBAR(POS)) * WORK(ROW)
	  POS = POS + NP - ROW - 1
  20    CONTINUE
	TOL(COL) = EPS * SUM
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SING(NP, NRBAR, D, RBAR, THETAB, SSERR, TOL, LINDEP,
     +   WORK, IER)
C
C     ALGORITHM AS274.5  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Checks for singularities, reports, and adjusts orthogonal
C     reductions produced by AS75.1.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), SSERR, TOL(NP),
     +        WORK(NP)
      LOGICAL LINDEP(NP)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, TEMP
      INTEGER COL, POS, ROW, NP2, POS2
C
      DATA ZERO/0.D0/
C
C     Check input parameters
C
      IER = 0
      IF (NP .LE. 0) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      DO 10 COL = 1, NP
   10 WORK(COL) = SQRT(D(COL))
C
      DO 40 COL = 1, NP
C
C     Set elements within RBAR to zero if they are less than TOL(COL) in
C     absolute value after being scaled by the square root of their row
C     multiplier.
C
	TEMP = TOL(COL)
	POS = COL - 1
	DO 30 ROW = 1, COL-1
	  IF (ABS(RBAR(POS)) * WORK(ROW) .LT. TEMP) RBAR(POS) = ZERO
	  POS = POS + NP - ROW - 1
   30   CONTINUE
C
C     If diagonal element is near zero, set it to zero, set appropriate
C     element of LINDEP, and use INCLUD to augment the projections in
C     the lower rows of the orthogonalization.
C
	LINDEP(COL) = .FALSE.
	IF (WORK(COL) .LT. TEMP) THEN
	  LINDEP(COL) = .TRUE.
	  IER = IER - 1
	  IF (COL .LT. NP) THEN
	    NP2 = NP - COL
	    POS2 = POS + NP - COL + 1
	    CALL INCLUD(NP2, NP2*(NP2-1)/2, D(COL), RBAR(POS+1),
     +            THETAB(COL), D(COL+1), RBAR(POS2), THETAB(COL+1),
     +            SSERR, IER)
	  ELSE
	    SSERR = SSERR + D(COL) * THETAB(COL)**2
	  END IF
	  D(COL) = ZERO
	  WORK(COL) = ZERO
	  THETAB(COL) = ZERO
	END IF
   40 CONTINUE
      RETURN
      END
C
      SUBROUTINE SS(NP, D, THETAB, SSERR, RSS, IER)
C
C     ALGORITHM AS274.6  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Calculates partial residual sums of squares from an orthogonal
C     reduction from AS75.1.
C
      INTEGER NP, IER
      DOUBLE PRECISION D(NP), THETAB(NP), SSERR, RSS(NP)
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION SUM
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (IER .NE. 0) RETURN
C
      SUM = SSERR
      RSS(NP) = SSERR
      DO 10 I = NP, 2, -1
	SUM = SUM + D(I) * THETAB(I)**2
	RSS(I-1) = SUM
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE COV(NP, NRBAR, D, RBAR, NREQ, RINV, VAR, COVMAT,
     +      DIMCOV, STERR, IER)
C
C     ALGORITHM AS274.7  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Calculate covariance matrix for regression coefficients for the
C     first NREQ variables, from an orthogonal reduction produced from
C     AS75.1.
C
C     Auxiliary routine called: INV
C
      INTEGER NP, NRBAR, NREQ, DIMCOV, IER
      DOUBLE PRECISION D(NP), RBAR(*), RINV(*), VAR, COVMAT(DIMCOV),
     +       STERR(NP)
C
C     Local variables.
C
      INTEGER POS, ROW, START, POS2, COL, POS1, K
      DOUBLE PRECISION ZERO, ONE, SUM
C
      DATA ZERO/0.D0/, ONE/1.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (DIMCOV .LT. NREQ*(NREQ+1)/2) IER = IER + 4
      DO 10 ROW = 1, NREQ
	IF (D(ROW) .EQ. ZERO) IER = -ROW
   10 CONTINUE
      IF (IER .NE. 0) RETURN
C
      CALL INV(NP, NRBAR, RBAR, NREQ, RINV)
      POS = 1
      START = 1
      DO 40 ROW = 1, NREQ
	POS2 = START
	DO 30 COL = ROW, NREQ
	  POS1 = START + COL - ROW
	  IF (ROW .EQ. COL) THEN
	    SUM = ONE / D(COL)
	  ELSE
	    SUM = RINV(POS1-1) / D(COL)
	  END IF
	  DO 20 K = COL+1, NREQ
	    SUM = SUM + RINV(POS1) * RINV(POS2) / D(K)
	    POS1 = POS1 + 1
	    POS2 = POS2 + 1
   20     CONTINUE
	  COVMAT(POS) = SUM * VAR
	  IF (ROW .EQ. COL) STERR(ROW) = SQRT(COVMAT(POS))
	  POS = POS + 1
   30   CONTINUE
	START = START + NREQ - ROW
   40 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE INV(NP, NRBAR, RBAR, NREQ, RINV)
C
C     ALGORITHM AS274.8  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Invert first NREQ rows and columns of Cholesky factorization
C     produced by AS75.1.
C
      INTEGER NP, NRBAR, NREQ
      DOUBLE PRECISION RBAR(*), RINV(*)
C
C     Local variables.
C
      INTEGER POS, ROW, COL, START, K, POS1, POS2
      DOUBLE PRECISION SUM, ZERO
C
      DATA ZERO/0.D0/
C
C     Invert RBAR ignoring row multipliers, from the bottom up.
C
      POS = NREQ * (NREQ-1)/2
      DO 30 ROW = NREQ-1, 1, -1
	START = (ROW-1) * (NP+NP-ROW)/2 + 1
	DO 20 COL = NREQ, ROW+1, -1
	  POS1 = START
	  POS2 = POS
	  SUM = ZERO
	  DO 10 K = ROW+1, COL-1
	    POS2 = POS2 + NREQ - K
	    SUM = SUM - RBAR(POS1) * RINV(POS2)
	    POS1 = POS1 + 1
   10     CONTINUE
	  RINV(POS) = SUM - RBAR(POS1)
	  POS = POS - 1
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PCORR(NP, NRBAR, D, RBAR, THETAB, SSERR, IN, WORK,
     +      CORMAT, DIMC, YCORR, IER)
C
C     ALGORITHM AS274.9  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Calculate partial correlations after the first IN variables
C     have been forced into the regression.
C
C     Auxiliary routine called: COR
C
      INTEGER NP, NRBAR, IN, DIMC, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), SSERR, WORK(NP),
     +        CORMAT(*), YCORR
C
C     Local variables.
C
      INTEGER START, IN1, I
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IN .LT. 0 .OR. IN .GT. NP-1) IER = IER + 4
      IF (DIMC .LT. (NP-IN)*(NP-IN-1)/2) IER = IER + 8
      IF (IER .NE. 0) RETURN
C
      START = IN * (NP+NP-IN-1)/2 + 1
      IN1 = IN + 1
      CALL COR(NP-IN, D(IN1), RBAR(START), THETAB(IN1), SSERR, WORK,
     +      CORMAT, YCORR)
C
C     Check for zeroes.
C
      DO 10 I = 1, NP-IN
	IF (WORK(I) .LE. ZERO) IER = -I
   10 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE COR(NP, D, RBAR, THETAB, SSERR, WORK, CORMAT, YCORR)
C
C     ALGORITHM AS274.10  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Calculate correlations from an orthogonal reduction.   This
C     routine will usually be called from PCORR, which will have
C     removed the appropriate number of rows at the start.
C
      INTEGER NP
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), SSERR, WORK(NP),
     +      CORMAT(*), YCORR(NP)
C
C     Local variables.
C
      INTEGER ROW, POS, COL1, POS1, COL2, POS2, DIFF
      DOUBLE PRECISION SUMY, SUM, ZERO
C
      DATA ZERO/0.D0/
C
C     Process by columns, including the projections of the dependent
C     variable (THETAB).
C
      SUMY = SSERR
      DO 10 ROW = 1, NP
   10 SUMY = SUMY + D(ROW) * THETAB(ROW)**2
      SUMY = SQRT(SUMY)
      POS = NP*(NP-1)/2
      DO 70 COL1 = NP, 1, -1
C
C     Calculate the length of column COL1.
C
	SUM = D(COL1)
	POS1 = COL1 - 1
	DO 20 ROW = 1, COL1-1
	  SUM = SUM + D(ROW) * RBAR(POS1)**2
	  POS1 = POS1 + NP - ROW - 1
   20   CONTINUE
	WORK(COL1) = SQRT(SUM)
C
C     If SUM = 0, set all correlations with this variable to zero.
C
	IF (SUM .EQ. ZERO) THEN
	  YCORR(COL1) = ZERO
	  DO 30 COL2 = NP, COL1+1, -1
	    CORMAT(POS) = ZERO
	    POS = POS - 1
   30     CONTINUE
	  GO TO 70
	END IF
C
C     Form cross-products, then divide by product of column lengths.
C
	SUM = D(COL1) * THETAB(COL1)
	POS1 = COL1 - 1
	DO 40 ROW = 1, COL1-1
	  SUM = SUM + D(ROW) * RBAR(POS1) * THETAB(ROW)
	  POS1 = POS1 + NP - ROW - 1
   40   CONTINUE
	YCORR(COL1) = SUM / (SUMY * WORK(COL1))
C
	DO 60 COL2 = NP, COL1+1, -1
	  IF (WORK(COL2) .GT. ZERO) THEN
	    POS1 = COL1 - 1
	    POS2 = COL2 - 1
	    DIFF = COL2 - COL1
	    SUM = ZERO
	    DO 50 ROW = 1, COL1-1
	      SUM = SUM + D(ROW) * RBAR(POS1) * RBAR(POS2)
	      POS1 = POS1 + NP - ROW - 1
	      POS2 = POS1 + DIFF
   50       CONTINUE
	    SUM = SUM + D(COL1) * RBAR(POS2)
	    CORMAT(POS) = SUM / (WORK(COL1) * WORK(COL2))
	  ELSE
	    CORMAT(POS) = ZERO
	  END IF
	  POS = POS - 1
   60   CONTINUE
   70 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB, RSS, FROM,
     +    TO, TOL, IER)
C
C     ALGORITHM AS274.11 APPL. STATIST. (1992) VOL 41, NO.2
C
C     Move variable from position FROM to position TO in an
C     orthogonal reduction produced by AS75.1.
C
      INTEGER NP, NRBAR, VORDER(NP), FROM, TO, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), RSS(NP), TOL(NP)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, D1, D2, X, ONE, D1NEW, D2NEW, CBAR, SBAR, Y
      INTEGER M, FIRST, LAST, INC, M1, M2, MP1, COL, POS, ROW
C
      DATA ZERO/0.D0/, ONE/1.D0/
C
C     Check input parameters
C
      IER = 0
      IF (NP .LE. 0) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (FROM .LT. 1 .OR. FROM .GT. NP) IER = IER + 4
      IF (TO .LT. 1 .OR. TO .GT. NP) IER = IER + 8
      IF (IER .NE. 0) RETURN
C
      IF (FROM .EQ. TO) RETURN
C
      IF (FROM .LT. TO) THEN
	FIRST = FROM
	LAST = TO - 1
	INC = 1
      ELSE
	FIRST = FROM - 1
	LAST = TO
	INC = -1
      END IF
      DO 70 M = FIRST, LAST, INC
C
C     Find addresses of first elements of RBAR in rows M and (M+1).
C
	M1 = (M-1)*(NP+NP-M)/2 + 1
	M2 = M1 + NP - M
	MP1 = M + 1
	D1 = D(M)
	D2 = D(MP1)
C
C     Special cases.
C
	IF (D1 .EQ. ZERO .AND. D2 .EQ. ZERO) GO TO 40
	X = RBAR(M1)
	IF (ABS(X) * SQRT(D1) .LT. TOL(MP1)) THEN
	  X = ZERO
	END IF
	IF (D1 .EQ. ZERO .OR. X .EQ. ZERO) THEN
	  D(M) = D2
	  D(MP1) = D1
	  RBAR(M1) = ZERO
	  DO 10 COL = M+2, NP
	    M1 = M1 + 1
	    X = RBAR(M1)
	    RBAR(M1) = RBAR(M2)
	    RBAR(M2) = X
	    M2 = M2 + 1
   10     CONTINUE
	  X = THETAB(M)
	  THETAB(M) = THETAB(MP1)
	  THETAB(MP1) = X
	  GO TO 40
	ELSE IF (D2 .EQ. ZERO) THEN
	  D(M) = D1 * X**2
	  RBAR(M1) = ONE / X
	  DO 20 COL = M+2, NP
	    M1 = M1 + 1
	    RBAR(M1) = RBAR(M1) / X
   20     CONTINUE
	  THETAB(M) = THETAB(M) / X
	  GO TO 40
	END IF
C
C     Planar rotation in regular case.
C
	D1NEW = D2 + D1*X**2
	CBAR = D2 / D1NEW
	SBAR = X * D1 / D1NEW
	D2NEW = D1 * CBAR
	D(M) = D1NEW
	D(MP1) = D2NEW
	RBAR(M1) = SBAR
	DO 30 COL = M+2, NP
	  M1 = M1 + 1
	  Y = RBAR(M1)
	  RBAR(M1) = CBAR*RBAR(M2) + SBAR*Y
	  RBAR(M2) = Y - X*RBAR(M2)
	  M2 = M2 + 1
   30   CONTINUE
	Y = THETAB(M)
	THETAB(M) = CBAR*THETAB(MP1) + SBAR*Y
	THETAB(MP1) = Y - X*THETAB(MP1)
C
C     Swap columns M and (M+1) down to row (M-1).
C
   40   IF (M .EQ. 1) GO TO 60
	POS = M
	DO 50 ROW = 1, M-1
	  X = RBAR(POS)
	  RBAR(POS) = RBAR(POS-1)
	  RBAR(POS-1) = X
	  POS = POS + NP - ROW - 1
   50   CONTINUE
C
C     Adjust variable order (VORDER), the tolerances (TOL) and
C     the vector of residual sums of squares (RSS).
C
   60   M1 = VORDER(M)
	VORDER(M) = VORDER(MP1)
	VORDER(MP1) = M1
	X = TOL(M)
	TOL(M) = TOL(MP1)
	TOL(MP1) = X
	RSS(M) = RSS(MP1) + D(MP1) * THETAB(MP1)**2
   70 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE REORDR(NP, NRBAR, VORDER, D, RBAR, THETAB, RSS, TOL,
     +      LIST, N, POS1, IER)
C
C     ALGORITHM AS274.12  APPL. STATIST. (1992) VOL 41, NO.2
C
C     Re-order the variables in an orthogonal reduction produced by
C     AS75.1 so that the N variables in LIST start at position POS1,
C     though will not necessarily be in the same order as in LIST.
C     Any variables in VORDER before position POS1 are not moved.
C
C     Auxiliary routine called: VMOVE
C
      INTEGER NP, NRBAR, VORDER(NP), N, LIST(N), POS1, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), RSS(NP), TOL(NP)
C
C     Local variables.
C
      INTEGER NEXT, I, L, J
C
C     Check N.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (N .LT. 1 .OR. N .GE. NP+1-POS1) IER = IER + 4
      IF (IER .NE. 0) RETURN
C
C     Work through VORDER finding variables which are in LIST.
C
      NEXT = POS1
      I = POS1
   10 L = VORDER(I)
      DO 20 J = 1, N
	IF (L .EQ. LIST(J)) GO TO 40
   20 CONTINUE
   30 I = I + 1
      IF (I .LE. NP) GO TO 10
C
C     If this point is reached, one or more variables in LIST has not
C     been found.
C
      IER = 8
      RETURN
C
C     Variable L is in LIST; move it up to position NEXT if it is not
C     already there.
C
   40 IF (I .GT. NEXT) CALL VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB,
     +      RSS, I, NEXT, TOL, IER)
      NEXT = NEXT + 1
      IF (NEXT .LT. N+POS1) GO TO 30
C
      RETURN
      END
C
      SUBROUTINE HDIAG(XROW, NP, NRBAR, D, RBAR, TOL, NREQ, HII, WK,
     +     IFAULT)
C
C     ALGORITHM AS274.13  APPL. STATIST. (1992) VOL.41, NO.2
C
      INTEGER NP, NRBAR, NREQ, IFAULT
      DOUBLE PRECISION XROW(NP), D(NP), RBAR(*), TOL(NP), HII, WK(NP)
C
C     Local variables
C
      INTEGER COL, ROW, POS
      DOUBLE PRECISION ZERO, SUM
C
      DATA ZERO /0.0D0/
C
C     Some checks
C
      IFAULT = 0
      IF (NP .LT. 1) IFAULT = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IFAULT = IFAULT + 2
      IF (NREQ .GT. NP) IFAULT = IFAULT + 4
      IF (IFAULT .NE. 0) RETURN
C
C     The elements of XROW.inv(RBAR).sqrt(D) are calculated and stored
C     in WK.
C
      HII = ZERO
      DO 20 COL = 1, NREQ
	IF (SQRT(D(COL)) .LE. TOL(COL)) THEN
	  WK(COL) = ZERO
	  GO TO 20
	END IF
	POS = COL - 1
	SUM = XROW(COL)
	DO 10 ROW = 1, COL-1
	  SUM = SUM - WK(ROW)*RBAR(POS)
	  POS = POS + NP - ROW - 1
   10   CONTINUE
	WK(COL) = SUM
	HII = HII + SUM**2 / D(COL)
   20 CONTINUE
C
      RETURN
      END
