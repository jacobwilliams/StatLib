(File: LSQ.DOC)
The module LSQ is for unconstrained least-squares fitting.   It is based
upon Applied Statistics algorithm AS 274 (see comments at the start of
the module).    A planar-rotation algorithm is used to update the QR-
factorization.   This makes it suitable for updating regressions as more
data become available.

By taking advantage of the MODULE facility, it has been possible to remove
most of the arguments to routines.   Apart from the new function VARPRD,
and a back-substitution routine BKSUB2 which it calls, the routines behave
as in AS 274.

This release, version 1.00 dated 24 May 1995, has only been tested using
one Fortran 90 compiler, Lahey's LF90 versions up to 1.10c.   It would be
appreciated if users could report their findings with other compilers, or
with the various translators to Fortran 77 and C.

The package as posted comprises:

      LSQ.DOC        this brief note
      LSQ.F90        the least-squares module
      DEMO.F90       a simple demonstration program
      FUELCONS.DAT   a file of data to be read by the demo program
      TEST1.F90      another program which fits a cubic polynomial and has
                     a deliberate singularity

Reference:
Miller, A.J. (1992). Algorithm AS 274: Least squares routines to supplement
      those of Gentleman.  Appl. Statist., vol.41(2), 458-478.

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

(File: LSQ.F90)
MODULE lsq

!     Module for unconstrained linear least-squares calculations.
!     The algorithm is suitable for updating LS calculations as more
!     data are added.   This is sometimes called recursive estimation.
!     Only one dependent variable is allowed.
!     Based upon Applied Statistics algorithm AS 274.
!     Translation from Fortran 77 to Fortran 90 by Alan Miller.
!     A subroutine, VARPRD, has been added for calculating the variances
!     of predicted values, and this uses a subroutine BKSUB2.

!     Version 1.00, 25 May 1995
!     Author: Alan Miller
!             CSIRO Division of Mathematics & Statistics
!             Private Bag 10, Rosebank MDC
!             Clayton 3169, Victoria, Australia
!     Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
!     e-mail: Alan.Miller @ mel.dms.csiro.au

!     The PUBLIC variables are:
!     lsq_kind = a KIND parameter for the floating-point quantities calculated
!                in this module.   See the more detailed explanation below.
!                This KIND parameter should be used for all floating-point
!                arguments passed to routines in this module.

!     nobs    = the number of observations processed to date.
!     ncol    = the total number of variables, including one for the constant,
!               if a constant is being fitted.
!     r_dim   = the dimension of array r = ncol*(ncol-1)/2
!     vorder  = an integer vector storing the current order of the variables
!               in the QR-factorization.   The initial order is 0, 1, 2, ...
!               if a constant is being fitted, or 1, 2, ... otherwise.
!     initialized = a logical variable which indicates whether space has
!                   been allocated for various arrays.
!     tol_set = a logical variable which is set when subroutine TOLSET has
!               been called to calculate tolerances for use in testing for
!               singularities.
!     rss_set = a logical variable indicating whether residual sums of squares
!               are available and usable.
!     d()     = array of row multipliers for the Cholesky factorization.
!               The factorization is X = Q.sqrt(D).R where Q is an ortho-
!               normal matrix which is NOT stored, D is a diagonal matrix
!               whose diagonal elements are stored in array d, and R is an
!               upper-triangular matrix with 1's as its diagonal elements.
!     rhs()   = vector of RHS projections (after scaling by sqrt(D)).
!               Thus Q'y = sqrt(D).rhs
!     r()     = the upper-triangular matrix R.   The upper triangle only,
!               excluding the implicit 1's on the diagonal, are stored by
!               rows.
!     tol()   = array of tolerances used in testing for singularities.
!     rss()   = array of residual sums of squares.   rss(i) is the residual
!               sum of squares with the first r variables in the model.
!               By changing the order of variables, the residual sums of
!               squares can be found for all possible subsets of the variables.
!               The residual sum of squares with NO variables in the model,
!               that is the total sum of squares of the y-values, can be
!               calculated as rss(1) + d(1)*rhs(1)^2.   If the first variable
!               is a constant, then rss(1) is the sum of squares of
!               (y - ybar) where ybar is the average value of y.
!     sserr   = residual sum of squares with all of the variables included.
!
!--------------------------------------------------------------------------

!     General declarations

INTEGER                      :: nobs, ncol, r_dim
INTEGER, ALLOCATABLE         :: vorder(:)
LOGICAL                      :: initialized = .false.,                  &
                                tol_set = .false., rss_set = .false.

!     Note. lsq_kind is being set to give at least 10 decimal digit
!           representation of floating point numbers.   This should be
!           adequate for most problems except the fitting of polynomials.
!           lsq_kind is being set so that the same code can be run on PCs
!           and Unix systems, which will usually represent floating-point
!           numbers in `double precision', and other systems with larger
!           word lengths which will give similar accuracy in `single
!           precision'.

INTEGER, PARAMETER           :: lsq_kind = SELECTED_REAL_KIND(10,70)
REAL (lsq_kind), ALLOCATABLE :: d(:), rhs(:), r(:), tol(:), rss(:)
REAL (lsq_kind)              :: sserr, zero = 0.D0, one = 1.D0,         &
                                vsmall = 1.d-69

PUBLIC                       :: lsq_kind, nobs, ncol, r_dim, vorder,    &
                                initialized, tol_set, rss_set,          &
                                d, rhs, r, tol, rss, sserr
PRIVATE                      :: zero, one, vsmall

SAVE                         ! Save everything

CONTAINS

SUBROUTINE startup(nvar, fit_const)

!     Allocates dimensions for arrays and initializes to zero
!     The calling program must set nvar = the number of variables, and
!     fit_const = .true. if a constant is to be included in the model,
!     otherwise fit_const = .false.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)            :: nvar
LOGICAL, INTENT(IN)            :: fit_const

!     Local variable
INTEGER                        :: i

nobs = 0
IF (fit_const) THEN
  ncol = nvar + 1
ELSE
  ncol = nvar
END IF

IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, vorder)
r_dim = ncol * (ncol - 1)/2
ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol) )

d = zero
rhs = zero
r = zero
sserr = zero

IF (fit_const) THEN
  DO i = 1, ncol
    vorder(i) = i-1
  END DO
ELSE
  DO i = 1, ncol
    vorder(i) = i
  END DO
END IF ! (fit_const)

initialized = .true.
tol_set = .false.
rss_set = .false.

END SUBROUTINE startup




SUBROUTINE includ(weight, xrow, yelem)

!     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

!     Calling this routine updates D, R, RHS and SSERR by the
!     inclusion of xrow, yelem with the specified weight.

!     *** WARNING  Array XROW is overwritten.

!     N.B. As this routine will be called many times in most applications,
!          checks have been eliminated.
!
!--------------------------------------------------------------------------


IMPLICIT NONE
REAL (lsq_kind),INTENT(IN)    :: weight, yelem
REAL (lsq_kind),INTENT(INOUT) :: xrow(ncol)

!     Local variables

INTEGER           :: i, k, nextr
REAL (lsq_kind)   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

nobs = nobs + 1
w = weight
y = yelem
rss_set = .false.
nextr = 1
DO i = 1, ncol

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (ABS(w) < vsmall) RETURN
  xi = xrow(i)
  IF (ABS(xi) < vsmall) THEN
    nextr = nextr + ncol - i
  ELSE
    di = d(i)
    wxi = w * xi
    dpi = di + wxi*xi
    cbar = di / dpi
    sbar = wxi / dpi
    w = cbar * w
    d(i) = dpi
    DO k = i+1, ncol
      xk = xrow(k)
      xrow(k) = xk - xi * r(nextr)
      r(nextr) = cbar * r(nextr) + sbar * xk
      nextr = nextr + 1
    END DO
    xk = y
    y = xk - xi * rhs(i)
    rhs(i) = cbar * rhs(i) + sbar * xk
  END IF
END DO ! i = 1, ncol

!     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
!     residual.

sserr = sserr + w * y * y

RETURN
END SUBROUTINE includ



SUBROUTINE regcf(beta, nreq, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients
!     for the first NREQ variables, given an orthogonal reduction from
!     AS75.1.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)            :: nreq
INTEGER, INTENT(OUT)           :: ifault
REAL (lsq_kind) ,INTENT(OUT)   :: beta(nreq)

!     Local variables

INTEGER             :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq .LT. 1 .OR. nreq .GT. ncol) ifault = ifault + 4
IF (ifault .NE. 0) RETURN

IF (.NOT. tol_set) CALL tolset

DO i = nreq, 1, -1
  IF (SQRT(d(i)) .LT. tol(i)) THEN
    beta(i) = zero
    d(i) = zero
  ELSE
    beta(i) = rhs(i)
    nextr = (i-1) * (ncol+ncol-i)/2 + 1
    DO j = i+1, nreq
      beta(i) = beta(i) - r(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
  END IF
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE regcf



SUBROUTINE tolset

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Sets up array TOL for testing for zeroes in an orthogonal
!     reduction formed using AS75.1.

!     Local variables.
!
!--------------------------------------------------------------------------

IMPLICIT NONE

!     Local variables

INTEGER           :: col, row, pos
REAL (lsq_kind)   :: eps, ten = 10.0, total, work(ncol)

!     EPS is a machine-dependent constant.

eps = ten * EPSILON(ten)

!     Set tol(i) = sum of absolute values in column I of R after
!     scaling each element by the square root of its row multiplier,
!     multiplied by EPS.

work = SQRT(d)
DO col = 1, ncol
  pos = col - 1
  total = work(col)
  DO row = 1, col-1
    total = total + ABS(r(pos)) * work(row)
    pos = pos + ncol - row - 1
  END DO
  tol(col) = eps * total
END DO

tol_set = .TRUE.
RETURN
END SUBROUTINE tolset




SUBROUTINE sing(lindep, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Checks for singularities, reports, and adjusts orthogonal
!     reductions produced by AS75.1.

!     Auxiliary routine called: INCLUD
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(OUT)        :: ifault
LOGICAL, INTENT(OUT)        :: lindep(ncol)

!     Local variables

REAL (lsq_kind)   :: temp, x(ncol), work(ncol), y, weight
INTEGER           :: col, pos, row, pos2

ifault = 0

work = SQRT(d)

DO col = 1, ncol

!     Set elements within R to zero if they are less than tol(col) in
!     absolute value after being scaled by the square root of their row
!     multiplier.

  temp = tol(col)
  pos = col - 1
  DO row = 1, col-1
    IF (ABS(r(pos)) * work(row) .LT. temp) r(pos) = zero
    pos = pos + ncol - row - 1
  END DO

!     If diagonal element is near zero, set it to zero, set appropriate
!     element of LINDEP, and use INCLUD to augment the projections in
!     the lower rows of the orthogonalization.

  lindep(col) = .FALSE.
  IF (work(col) .LE. temp) THEN
    lindep(col) = .TRUE.
    ifault = ifault - 1
    IF (col .LT. ncol) THEN
      pos2 = pos + ncol - col + 1
      x = zero
      x(col+1:ncol) = r(pos+1:pos2-1)
      y = rhs(col)
      weight = d(col)
      r(pos+1:pos2-1) = zero
      d(col) = zero
      rhs(col) = zero
      CALL includ(weight, x, y)
    ELSE
      sserr = sserr + d(col) * rhs(col)**2
    END IF ! (col .LT. ncol)
  END IF ! (work(col) .LE. temp)
END DO ! col = 1, ncol
RETURN
END SUBROUTINE sing



SUBROUTINE ss

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculates partial residual sums of squares from an orthogonal
!     reduction from AS75.1.
!
!--------------------------------------------------------------------------

!     Local variables

IMPLICIT NONE
INTEGER           :: i
REAL (lsq_kind)   :: total

total = sserr
rss(ncol) = sserr
DO i = ncol, 2, -1
  total = total + d(i) * rhs(i)**2
  rss(i-1) = total
END DO

rss_set = .TRUE.
RETURN
END SUBROUTINE ss



SUBROUTINE cov(nreq, var, covmat, dimcov, sterr, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate covariance matrix for regression coefficients for the
!     first nreq variables, from an orthogonal reduction produced from
!     AS75.1.

!     Auxiliary routine called: INV
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)          :: nreq, dimcov
INTEGER, INTENT(OUT)         :: ifault
REAL (lsq_kind), INTENT(OUT) :: var, covmat(dimcov), sterr(nreq)

!     Local variables.

INTEGER                      :: dim_rinv, pos, row, start, pos2, col, pos1, k
REAL (lsq_kind)              :: total
REAL (lsq_kind), ALLOCATABLE :: rinv(:)

!     Check that dimension of array covmat is adequate.

IF (dimcov < nreq*(nreq+1)/2) THEN
  ifault = 1
  RETURN
END IF

!     Check for small or zero multipliers on the diagonal.

ifault = 0
DO row = 1, nreq
  IF (ABS(d(row)) .LT. vsmall) ifault = -row
END DO
IF (ifault .NE. 0) RETURN

!     Calculate estimate of the residual variance.

IF (nobs > nreq) THEN
  var = rss(nreq) / (nobs - nreq)
ELSE
  ifault = 2
  RETURN
END IF

dim_rinv = nreq*(nreq-1)/2
ALLOCATE ( rinv(dim_rinv) )

CALL INV(nreq, rinv)
pos = 1
start = 1
DO row = 1, nreq
  pos2 = start
  DO col = row, nreq
    pos1 = start + col - row
    IF (row .EQ. col) THEN
      total = one / d(col)
    ELSE
      total = rinv(pos1-1) / d(col)
    END IF
    DO K = col+1, nreq
      total = total + rinv(pos1) * rinv(pos2) / d(k)
      pos1 = pos1 + 1
      pos2 = pos2 + 1
    END DO ! K = col+1, nreq
    covmat(pos) = total * var
    IF (row .EQ. col) sterr(row) = SQRT(covmat(pos))
    pos = pos + 1
  END DO ! col = row, nreq
  start = start + nreq - row
END DO ! row = 1, nreq

RETURN
END SUBROUTINE cov



SUBROUTINE inv(nreq, rinv)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Invert first nreq rows and columns of Cholesky factorization
!     produced by AS 75.1.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)            :: nreq
REAL (lsq_kind), INTENT(OUT)   :: rinv(*)

!     Local variables.

INTEGER           :: pos, row, col, start, K, pos1, pos2
REAL (lsq_kind)   :: total

!     Invert R ignoring row multipliers, from the bottom up.

pos = nreq * (nreq-1)/2
DO row = nreq-1, 1, -1
  start = (row-1) * (ncol+ncol-row)/2 + 1
  DO col = nreq, row+1, -1
    pos1 = start
    pos2 = pos
    total = zero
    DO K = row+1, col-1
      pos2 = pos2 + nreq - K
      total = total - r(pos1) * rinv(pos2)
      pos1 = pos1 + 1
    END DO ! K = row+1, col-1
    rinv(pos) = total - r(pos1)
    pos = pos - 1
  END DO ! col = nreq, row+1, -1
END DO ! row = nreq-1, 1, -1

RETURN
END SUBROUTINE inv



SUBROUTINE partial_corr(in, cormat, dimc, ycorr, ifault)

!     Replaces subroutines PCORR and COR of:
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate partial correlations after the variables in rows
!     1, 2, ..., IN have been forced into the regression.
!     If IN = 1, and the first row of R represents a constant in the
!     model, then the usual simple correlations are returned.

!     If IN = 0, the value returned in array CORMAT for the correlation
!     of variables Xi & Xj is:
!       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )

!     On return, array CORMAT contains the upper triangle of the matrix of
!     partial correlations stored by rows, excluding the 1's on the diagonal.
!     e.g. if IN = 2, the consecutive elements returned are:
!     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
!     Array YCORR stores the partial correlations with the Y-variable
!     starting with YCORR(IN+1) = partial correlation with the variable in
!     position (IN+1).
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)          :: in, dimc
INTEGER, INTENT(OUT)         :: ifault
REAL (lsq_kind), INTENT(OUT) :: cormat(dimc), ycorr(ncol)

!     Local variables.

INTEGER           :: base_pos, pos, row, col, col1, col2, pos1, pos2
REAL (lsq_kind)   :: rms(in+1:ncol), sumxx, sumxy, sumyy, work(in+1:ncol)

!     Some checks.

ifault = 0
IF (in .LT. 0 .OR. in .GT. ncol-1) ifault = ifault + 4
IF (dimc .LT. (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
IF (ifault .NE. 0) RETURN

!     Base position for calculating positions of elements in row (IN+1) of R.

base_pos = in*ncol - (in+1)*(in+2)/2

!     Calculate 1/RMS of elements in columns from IN to (ncol-1).

IF (d(in+1) > zero) rms(in+1) = one / SQRT(d(in+1))
DO col = in+2, ncol
  pos = base_pos + col
  sumxx = d(col)
  DO row = in+1, col-1
    sumxx = sumxx + d(row) * r(pos)**2
    pos = pos + ncol - row - 1
  END DO ! row = in+1, col-1
  IF (sumxx > zero) THEN
    rms(col) = one / SQRT(sumxx)
  ELSE
    rms(col) = zero
    ifault = -col
  END IF ! (sumxx > zero)
END DO ! col = in+1, ncol-1

!     Calculate 1/RMS for the Y-variable

sumyy = sserr
DO row = in+1, ncol
  sumyy = sumyy + d(row) * rhs(row)**2
END DO ! row = in+1, ncol
IF (sumyy > zero) sumyy = one / SQRT(sumyy)

!     Calculate sums of cross-products.
!     These are obtained by taking dot products of pairs of columns of R,
!     but with the product for each row multiplied by the row multiplier
!     in array D.

pos = 1
DO col1 = in+1, ncol
  sumxy = zero
  work(col1+1:ncol) = zero
  pos1 = base_pos + col1
  DO row = in+1, col1-1
    pos2 = pos1 + 1
    DO col2 = col1+1, ncol
      work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
      pos2 = pos2 + 1
    END DO ! col2 = col1+1, ncol
    sumxy = sumxy + d(row) * r(pos1) * rhs(row)
    pos1 = pos1 + ncol - row - 1
  END DO ! row = in+1, col1-1

!     Row COL1 has an implicit 1 as its first element (in column COL1)

  pos2 = pos1 + 1
  DO col2 = col1+1, ncol
    work(col2) = work(col2) + d(col1) * r(pos2)
    pos2 = pos2 + 1
    cormat(pos) = work(col2) * rms(col1) * rms(col2)
    pos = pos + 1
  END DO ! col2 = col1+1, ncol
  sumxy = sumxy + d(col1) * rhs(col1)
  ycorr(col1) = sumxy * rms(col1) * sumyy
END DO ! col1 = in+1, ncol-1

ycorr(1:in) = zero

RETURN
END SUBROUTINE partial_corr




SUBROUTINE vmove(from, to, ifault)

!     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2

!     Move variable from position FROM to position TO in an
!     orthogonal reduction produced by AS75.1.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)            :: from, to
INTEGER, INTENT(OUT)           :: ifault

!     Local variables

REAL (lsq_kind)   :: d1, d2, X, d1new, d2new, cbar, sbar, Y
INTEGER           :: m, first, last, inc, m1, m2, mp1, col, pos, row

!     Check input parameters

ifault = 0
IF (from .LT. 1 .OR. from .GT. ncol) ifault = ifault + 4
IF (to .LT. 1 .OR. to .GT. ncol) ifault = ifault + 8
IF (ifault .NE. 0) RETURN

IF (from .EQ. to) RETURN

IF (.NOT. rss_set) CALL ss

IF (from .LT. to) THEN
  first = from
  last = to - 1
  inc = 1
ELSE
  first = from - 1
  last = to
  inc = -1
END IF

   DO 70 m = first, last, inc

!     Find addresses of first elements of R in rows M and (M+1).

     m1 = (m-1)*(ncol+ncol-m)/2 + 1
     m2 = m1 + ncol - m
     mp1 = m + 1
     d1 = d(m)
     d2 = d(mp1)

!     Special cases.

     IF (d1 .LT. vsmall .AND. d2 .LT. vsmall) GO TO 40
     X = r(m1)
     IF (ABS(x) * SQRT(d1) .LT. tol(mp1)) THEN
       X = zero
     END IF
     IF (d1 .LT. vsmall .OR. ABS(x) .LT. vsmall) THEN
       d(m) = d2
       d(mp1) = d1
       r(m1) = zero
       DO col = m+2, ncol
         m1 = m1 + 1
         X = r(m1)
         r(m1) = r(m2)
         r(m2) = x
         m2 = m2 + 1
       END DO ! col = m+2, ncol
       X = rhs(m)
       rhs(m) = rhs(mp1)
       rhs(mp1) = x
       GO TO 40
     ELSE IF (d2 .LT. vsmall) THEN
       d(m) = d1 * x**2
       r(m1) = one / x
       r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
       rhs(m) = rhs(m) / x
       GO TO 40
     END IF
!
!     Planar rotation in regular case.
!
     d1new = d2 + d1*x**2
     cbar = d2 / d1new
     sbar = x * d1 / d1new
     d2new = d1 * cbar
     d(m) = d1new
     d(mp1) = d2new
     r(m1) = sbar
     DO col = m+2, ncol
       m1 = m1 + 1
       y = r(m1)
       r(m1) = cbar*r(m2) + sbar*y
       r(m2) = y - x*r(m2)
       m2 = m2 + 1
     END DO ! col = m+2, ncol
     y = rhs(m)
     rhs(m) = cbar*rhs(mp1) + sbar*y
     rhs(mp1) = y - x*rhs(mp1)

!     Swap columns M and (M+1) down to row (M-1).

40   IF (m .EQ. 1) GO TO 60
     pos = m
     DO row = 1, m-1
       x = r(pos)
       r(pos) = r(pos-1)
       r(pos-1) = x
       pos = pos + ncol - row - 1
     END DO ! row = 1, m-1

!     Adjust variable order (VORDER), the tolerances (TOL) and
!     the vector of residual sums of squares (RSS).

60   m1 = vorder(m)
     vorder(m) = vorder(mp1)
     vorder(mp1) = m1
     x = tol(m)
     tol(m) = tol(mp1)
     tol(mp1) = x
     rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
70 CONTINUE

RETURN
END SUBROUTINE vmove



SUBROUTINE reordr(list, n, pos1, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Re-order the variables in an orthogonal reduction produced by
!     AS75.1 so that the N variables in LIST start at position POS1,
!     though will not necessarily be in the same order as in LIST.
!     Any variables in VORDER before position POS1 are not moved.

!     Auxiliary routine called: VMOVE
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)         :: n, list(n), pos1
INTEGER, INTENT(OUT)        :: ifault

!     Local variables.

INTEGER        :: next, i, l, j

!     Check N.

   ifault = 0
   IF (n .LT. 1 .OR. n .GT. ncol+1-pos1) ifault = ifault + 4
   IF (ifault .NE. 0) RETURN

!     Work through VORDER finding variables which are in LIST.

   next = pos1
   i = pos1
10 l = vorder(i)
   DO j = 1, n
     IF (l .EQ. list(j)) GO TO 40
   END DO
30 i = i + 1
   IF (i .LE. ncol) GO TO 10

!     If this point is reached, one or more variables in LIST has not
!     been found.

   ifault = 8
   RETURN

!     Variable L is in LIST; move it up to position NEXT if it is not
!     already there.

40 IF (i .GT. next) CALL vmove(i, next, ifault)
   next = next + 1
   IF (next .LT. n+pos1) GO TO 30

   RETURN
END SUBROUTINE reordr



SUBROUTINE hdiag(xrow, nreq, hii, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)          :: nreq
INTEGER, INTENT(OUT)         :: ifault
REAL (lsq_kind), INTENT(IN)  :: xrow(ncol)
REAL (lsq_kind), INTENT(OUT) :: hii

!     Local variables

INTEGER             :: col, row, pos
REAL (lsq_kind)     :: total, wk(ncol)

!     Some checks

ifault = 0
IF (nreq .GT. ncol) ifault = ifault + 4
IF (ifault .NE. 0) RETURN

!     The elements of xrow.inv(R).sqrt(D) are calculated and stored
!     in WK.

hii = zero
DO col = 1, nreq
  IF (SQRT(d(col)) .LE. tol(col)) THEN
    wk(col) = zero
  ELSE
    pos = col - 1
    total = xrow(col)
    DO row = 1, col-1
      total = total - wk(row)*r(pos)
      pos = pos + ncol - row - 1
    END DO ! row = 1, col-1
    wk(col) = total
    hii = hii + total**2 / d(col)
  END IF
END DO ! col = 1, nreq

RETURN
END SUBROUTINE hdiag



REAL (lsq_kind) FUNCTION varprd(x, nreq, var, ifault)

!     Calculate the variance of x'b where b consists of the first nreq
!     least-squares regression coefficients.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)          :: nreq
INTEGER, INTENT(OUT)         :: ifault
REAL (lsq_kind), INTENT(IN)  :: x(*)
REAL (lsq_kind), INTENT(OUT) :: var

!     Local variables

INTEGER           :: row
REAL (lsq_kind)   :: wk(nreq)

!     Check input parameter values

varprd = zero
ifault = 0
IF (nreq .LT. 1 .OR. nreq .GT. ncol) ifault = ifault + 4
IF (nobs .LE. nreq) ifault = ifault + 8
IF (ifault .NE. 0) RETURN

!     Calculate the residual variance estimate.

var = sserr / (nobs - nreq)

!     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
!     First call BKSUB2 to calculate (inv R')x by back-substitution.

CALL BKSUB2(x, wk, nreq)
DO row = 1, nreq
  IF(d(row) .GT. tol(row)) varprd = varprd + wk(row)**2 / d(row)
END DO

varprd = varprd * var

RETURN
END FUNCTION varprd



SUBROUTINE bksub2(x, b, nreq)

!     Solve x = R'b for b given x, using only the first nreq rows and
!     columns of R, and only the first nreq elements of R.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)          :: nreq
REAL (lsq_kind), INTENT(IN)  :: x(nreq)
REAL (lsq_kind), INTENT(OUT) :: b(nreq)

!     Local variables

INTEGER           :: pos, row, col
REAL (lsq_kind)   :: temp

!     Solve by back-substitution, starting from the top.

DO row = 1, nreq
  pos = row - 1
  temp = x(row)
  DO col = 1, row-1
    temp = temp - r(pos)*b(col)
    pos = pos + ncol - col - 1
  END DO
  b(row) = temp
END DO

RETURN
END SUBROUTINE bksub2


END MODULE lsq

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

(File: DEMO.F90)
PROGRAM demo

!     Program to demonstrate the use of module LSQ for unconstrained
!     linear least-squares regression.
!     The data on fuel consumption are from:
!     Sanford Weisberg `Applied Linear Regression', 2nd edition, 1985,
!     pages 35-36.   Publisher: Wiley

!     The program is designed so that users can easily edit it for their
!     problems.   If you are likely to have more than 500 cases then increase
!     the value of maxcases below.   Similarly, if you may have more than
!     30 predictor variables, change maxvar below.

!     The subscript 0 is used to relate to the constant in the model.
!     Module LSQ treats the constant in the model, if there is one, just as
!     any other variable.   The numbering of all arrays starts from 1 in LSQ.

USE lsq                      ! Must be the first declaration in the program
IMPLICIT NONE

INTEGER, PARAMETER :: maxcases = 500, maxvar = 30,                            &
                      max_cdim = maxvar*(maxvar+1)/2, lout = 6
INTEGER            :: case, nvar, iostatus, nreq, ifault, i, list(maxvar), j, &
                      in
REAL ( lsq_kind )  :: xx(0:maxvar), yy, wt, one = 1.D0, beta(0:maxvar),       &
                      var, covmat(max_cdim), sterr(0:maxvar), hii,            &
                      cormat(max_cdim), ycorr(maxvar)
REAL (KIND(0.E0))  :: x(maxcases, maxvar), y(maxcases), t(0:maxvar), tmin,    &
                      r2, resid(maxcases), std_resid(maxcases), std_err_pred, &
                      fitted, stdev_res
CHARACTER          :: heading*80, vname(0:maxvar)*8, state(50)*2, y_name*8, key
LOGICAL            :: fit_const, lindep(0:maxvar)

WRITE(*, *)'The data on fuel consumption are from:'
WRITE(*, *)'Sanford Weisberg "Applied Linear Regression", 2nd edition, 1985,'
WRITE(*, *)'pages 35-36.   Publisher: Wiley   ISBN: 0-471-87957-6'
WRITE(*, *)

OPEN(8, file='fuelcons.dat', status='old')

!     The first line of this file contains the names of the variables.

READ(8, '(a)') heading

vname(0) = 'Constant'

!     Read the heading to get the variable names & the number of variables

CALL separate_text(heading, vname(1), nvar)

nvar = nvar - 2              ! nvar is the number of variables.
                             ! 1 is subtracted as the first field in this
                             ! file contains the abbreviated state name,
                             ! & 1 is subtracted for the Y-variable.
vname(1:nvar) = vname(2:nvar+1)
y_name = vname(nvar+2)

WRITE(*, *)'No. of variables =', nvar
WRITE(*, *)'Predictor variables are:'
WRITE(*, '(1x, 9a9)') vname(1:nvar)
WRITE(*, *)'Dependent variable is: ', y_name

fit_const = .true.           ! Change to .false. if fitting a model without
                             ! a constant.

CALL startup(nvar, fit_const)          ! Initializes the QR-factorization

!     Read in the data, one line at a time, and progressively update the
!     QR-factorization.

wt = one
case = 1

DO
  READ(8, *, IOSTAT=iostatus) state(case), x(case, 1:nvar), y(case)
  IF (iostatus > 0) CYCLE              ! Error in data
  IF (iostatus < 0) EXIT               ! End of file

  xx(0) = one                          ! A one is inserted as the first
                                       ! variable if a constant is being fitted.
  xx(1:nvar) = x(case, 1:nvar)         ! New variables and transformed variables
                                       ! will often be generated here.
  yy = y(case)
  CALL includ(wt, xx, yy)
  case = case + 1
END DO

WRITE(*, *)'No. of observations =', nobs

CALL sing(lindep, ifault)              ! Checks for singularities

IF (ifault == 0) THEN
  WRITE(*, *)'QR-factorization is not singular'
ELSE
  DO i = 1, nvar
    IF (lindep(i)) THEN
      WRITE(*, *) vname(i), ' is exactly linearly related to earlier variables'
    END IF
  END DO ! i = 1, nvar
END IF ! (ifault == 0)
WRITE(*, *)

!     Show correlations (IN = 1 for `usual' correlations)

in = 1
CALL partial_corr(in, cormat, max_cdim, ycorr, ifault)
CALL printc(in, cormat, max_cdim, ycorr, vname, y_name, 1, lout, ifault)
WRITE(*, *)
WRITE(*, *)'Press ENTER to continue'
READ(*, '(a)') key

!     Weisberg only uses variables TAX, INC, ROAD and DLIC.   These are
!     currently in positions 2, 4, 5 & 7.
!     We could have set NVAR = 4 and copied only these variables from X to
!     XX, but we will show how regressions can be performed for a subset
!     of the variables in the QR-factorization.   Routine REORDR will be used
!     to re-order the variables.

list(1:4) = (/ 2, 4, 5, 7/)
CALL reordr(list, 4, 2, ifault)        ! Re-order so that the first 4 variables
                                       ! appear in positions 2,3,4 & 5.
                                       ! N.B. Though variables # 2, 4, 5 & 7
                                       !      will occupy positions 2-5, they
                                       !      may be in any order.
                                       ! The constant remains in position 1.
WRITE(*, 910) (vname(vorder(1:ncol)))
910 FORMAT(1x, 'Current order of variables:'/1x, 8a9/)

!     Calculate regression coefficients of Y against the first variable, which
!     was just a constant = 1.0, and the next 4 predictors.

CALL tolset                            ! Calculate tolerances before calling
                                       ! subroutine regcf.
nreq = 5                               ! i.e. Const, TAX, INC, ROAD, DLIC
CALL regcf(beta, nreq, ifault)

CALL ss                                ! Calculate residual sums of squares

!     Calculate covariance matrix of the regression coefficients & their
!     standard errors.

CALL cov(nreq, var, covmat, max_cdim, sterr, ifault)

!     Calculate t-values

t(0:nreq-1) = beta(0:nreq-1) / sterr(0:nreq-1)

!     Output regression table, residual sums of squares, and R-squared.

WRITE(*, *)
WRITE(*, *)'Variable   Regn.coeff.   Std.error  t-value   Res.sum of sq.'
DO i = 0, nreq-1
  WRITE(*, 900) vname(vorder(i+1)), beta(i), sterr(i), t(i), rss(i+1)
  900 FORMAT(1x, a8, 2x, g12.4, 2x, g11.4, 1x, f7.2, 2x, g14.6)
END DO
WRITE(*, *)

!     Now delete the variable with the smallest t-value by moving it
!     to position 5 and then repeating the calculations for the constant
!     and the next 3 variables.

j = 1
tmin = ABS(t(1))
DO i = 2, nreq-1
  IF (ABS(t(i)) < tmin) THEN
    j = i
    tmin = ABS(t(i))
  END IF
END DO
j = j + 1                    ! Add 1 as the t-array started at subscript 0

WRITE(*, *)'Removing variable in position', j
CALL vmove(j, nreq, ifault)
nreq = nreq - 1
CALL regcf(beta, nreq, ifault)
CALL ss
CALL cov(nreq, var, covmat, max_cdim, sterr, ifault)
t(0:nreq-1) = beta(0:nreq-1) / sterr(0:nreq-1)

!     Output regression table, residual sums of squares, and R-squared.

WRITE(*, *)
WRITE(*, *)'Variable   Regn.coeff.   Std.error  t-value   Res.sum of sq.'
DO i = 0, nreq-1
  WRITE(*, 900) vname(vorder(i+1)), beta(i), sterr(i), t(i), rss(i+1)
END DO
WRITE(*, *)
stdev_res = SQRT( rss(nreq)/(nobs - nreq) )
r2 = one - rss(nreq) / rss(1)          ! RSS(1) is the rss if only the constant
                                       ! is fitted; RSS(nreq) is the rss after
                                       ! fitting the requested NREQ variables.
WRITE(*, '(1x, "R^2 =", f8.4, 5x, "Std. devn. of residuals =", g12.4/)')      &
      r2, stdev_res
WRITE(*, *)'N.B. Some statistical packages wrongly call the standard deviation'
WRITE(*, *)'     of the residuals the standard error of prediction'
WRITE(*, *)

!     Calculate residuals, hii, standardized residuals & standard errors of
!     prediction.

WRITE(*, *)'Press ENTER to continue'
READ(*, '(a)') key

WRITE(*, *)'State     Actual    Fitted  Residual  Std.resid.  SE(prediction)'
DO i = 1, nobs
  xx(0) = one
  DO j = 1, nreq-1
    xx(j) = x(i, vorder(j+1))          ! N.B. Regression coefficient j is for
                                       !      the variable vorder(j+1)
  END DO
  fitted = DOT_PRODUCT( beta(0:nreq-1), xx(0:nreq-1) )
  resid(i) = y(i) - fitted
  CALL hdiag(xx, nreq, hii, ifault)
  std_resid(i) = resid(i) / SQRT(var*(one - hii))
  std_err_pred = varprd(xx, nreq, var, ifault)
  WRITE(*, 920) state(i), y(i), fitted, resid(i), std_resid(i), std_err_pred
  920 FORMAT(3x, a2, 2x, 3f10.1, f9.2, 5x, f10.0)

  IF (i .EQ. 24) THEN
    WRITE(*, *)'Press ENTER to continue'
    READ(*, '(a)') key
    WRITE(*, *)'State     Actual    Fitted  Residual  Std.resid.  SE(prediction)'
  END IF
END DO ! i = 1, nobs

END PROGRAM demo



SUBROUTINE separate_text(text, name, number)

!     Takes the character string in `text' and separates it into `names'.
!     This version only allows spaces as delimiters.

IMPLICIT NONE
CHARACTER (LEN = *), INTENT(INOUT)  :: text
CHARACTER (LEN = *), INTENT(OUT)    :: name(*)
INTEGER, INTENT(OUT)                :: number

!     Local variables

INTEGER  :: length, len_name, pos1, pos2, nchars

length = LEN_TRIM(text)
number = 0
IF (length == 0) RETURN

len_name = LEN(name(1))
IF (len_name < 1) RETURN

pos1 = 1
DO
  DO                                   ! Remove any leading blanks
    IF (text(pos1:pos1) == ' ') THEN
      text(pos1:) = text(pos1+1:)
      length = length - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO
  pos2 = INDEX(text(pos1:length), ' ') ! Find position of next blank
  IF (pos2 > 0) THEN
    pos2 = pos2 + pos1 - 1
  ELSE                                 ! Last name, no blank found
    pos2 = length + 1
  END IF
  number = number + 1
  nchars = MIN(len_name, pos2-pos1)
  name(number) = text(pos1:pos1+nchars-1)
  pos1 = pos2 + 1                      ! Move past the blank
  IF (pos1 > length) RETURN
END DO

END SUBROUTINE separate_text



SUBROUTINE printc(in, cormat, dimc, ycorr, vname, yname, iopt, lout, ier)

!     Print (partial) correlations calculated using partial_corr to unit lout.
!     If IOPT = 0, print correlations with the Y-variable only.

USE lsq

IMPLICIT NONE
INTEGER, INTENT(IN)           :: in, dimc, iopt, lout
INTEGER, INTENT(OUT)          :: ier
REAL ( lsq_kind ), INTENT(IN) :: cormat(dimc), ycorr(ncol)
CHARACTER, INTENT(IN)         :: vname(0:ncol-1)*8, yname*8

!     Local variables.

INTEGER            :: nrows, j1, j2, i1, i2, row, upos, tpos, last
CHARACTER          :: text*74, empty*65 = ' ', char1*9 = ' 1.0'

!     Check validity of arguments

ier = 0
IF (in .GE. ncol) ier = 1
IF (ncol .LE. 1) ier = ier + 2
nrows = ncol - in
IF (dimc .LE. nrows*(nrows-1)/2) ier = ier + 4
IF (ier .NE. 0) RETURN

!     If iopt.NE.0 output heading

IF (iopt .EQ. 0) GO TO 30
WRITE(lout, 900)
900 FORMAT(/5x, 'Correlation matrix')
j1 = in + 1

DO
  j2 = MIN(j1+6, ncol)
  i1 = j1 - in
  i2 = j2 - in
  WRITE(lout, 910) vname(vorder(j1:j2))
  910 FORMAT(11X, 7(a8, 1x))

!     Print correlations for rows 1 to i2, columns i1 to i2.

  DO row = 1, i2
    text = ' ' // vname(vorder(row+in)) // empty
    IF (i1 .GT. row) THEN
      upos = (row-1) * (nrows+nrows-row) /2 + (i1-row)
      last = upos + i2 - i1
      WRITE(text(12:74), '(7(F8.5, 1x))') cormat(upos:last)
    ELSE
      upos = (row-1) * (nrows+nrows-row) /2 + 1
      tpos = 12 + 9*(row-i1)
      text(tpos:tpos+8) = char1
      last = upos + i2 - row - 1
      IF (row .LT. i2) WRITE(text(tpos+9:74), '(6(F8.5, 1X))') cormat(upos:last)
    END IF
    WRITE(lout, '(a)') text
  END DO ! row = 1, i2

!     Move onto the next block of columns.

  j1 = j2 + 1
  IF (j1 .GT. ncol) EXIT
END DO

!     Correlations with the Y-variable.

30 WRITE(lout, 920) yname
   920 FORMAT(/5x, 'Correlations with the dependent variable: ', a)
   j1 = in + 1
DO
  j2 = MIN(j1+7, ncol)
  WRITE(lout, 930) vname(vorder(j1:j2))
  930 FORMAT(/1X, 8(A8, 1X))
  WRITE(lout, 940) ycorr(j1:j2)
  940 FORMAT(1X, 8(F8.5, 1X))
  j1 = j2 + 1
  IF (j1 .GT. ncol) EXIT
END DO

END SUBROUTINE printc

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

(File: FUELCONS.DAT)
State Popn     Tax     NLic  Income  RoadMls FuelCon   DLic  Fuel_Pop
ME    1029.    9.00    540.   3.571   1.976    557.    52.5    541.
NH     771.    9.00    441.   4.092   1.250    404.    57.2    524.
VT     462.    9.00    268.   3.865   1.586    259.    58.0    561.
MA    5787.    7.50   3060.   4.870   2.351   2396.    52.9    414.
RI     968.    8.00    527.   4.399    .431    397.    54.4    410.
CN    3082.   10.00   1760.   5.342   1.333   1408.    57.1    457.
NY   18366.    8.00   8278.   5.319  11.868   6312.    45.1    344.
NJ    7367.    8.00   4074.   5.126   2.138   3439.    55.3    467.
PA   11926.    8.00   6312.   4.447   8.577   5528.    52.9    464.
OH   10783.    7.00   5948.   4.512   8.507   5375.    55.2    498.
IN    5291.    8.00   2804.   4.391   5.939   3068.    53.0    580.
IL   11251.    7.50   5903.   5.126  14.186   5301.    52.5    471.
MI    9082.    7.00   5213.   4.817   6.930   4768.    57.4    525.
WI    4520.    7.00   2465.   4.207   6.580   2294.    54.5    508.
MN    3896.    7.00   2368.   4.332   8.159   2204.    60.8    566.
IA    2883.    7.00   1689.   4.318  10.340   1830.    58.6    635.
MO    4753.    7.00   2719.   4.206   8.508   2865.    57.2    603.
ND     632.    7.00    341.   3.718   4.725    451.    54.0    714.
SD     579.    7.00    419.   4.716   5.915    501.    72.4    865.
NE    1525.    8.50   1033.   4.341   6.010    976.    67.7    640.
KS    2258.    7.00   1496.   4.593   7.834   1466.    66.3    649.
DE     565.    8.00    340.   4.983    .602    305.    60.2    540.
MD    4056.    9.00   2073.   4.897   2.449   1883.    51.1    464.
VA    4764.    9.00   2463.   4.258   4.686   2604.    51.7    547.
WV    1781.    8.50    982.   4.574   2.619    819.    55.1    460.
NC    5214.    9.00   2835.   3.721   4.746   2953.    54.4    566.
SC    2665.    8.00   1460.   3.448   5.399   1537.    54.8    577.
GA    4720.    7.50   2731.   3.846   9.061   2979.    57.9    631.
FL    7259.    8.00   4084.   4.188   5.975   4169.    56.3    574.
KY    3299.    9.00   1626.   3.601   4.650   1761.    49.3    534.
TN    4031.    7.00   2088.   3.640   6.905   2301.    51.8    571.
AL    3510.    7.00   1801.   3.333   6.594   1946.    51.3    554.
MS    2263.    8.00   1309.   3.063   6.524   1306.    57.8    577.
AR    1978.    7.50   1081.   3.357   4.121   1242.    54.7    628.
LA    3720.    8.00   1813.   3.528   3.495   1812.    48.7    487.
OK    2634.    6.58   1657.   3.802   7.834   1695.    62.9    644.
TX   11649.    5.00   6595.   4.045  17.782   7451.    56.6    640.
MT     719.    7.00    421.   3.897   6.385    506.    58.6    704.
ID     756.    8.50    501.   3.635   3.274    490.    66.3    648.
WY     345.    7.00    232.   4.345   3.905    334.    67.2    968.
CO    2357.    7.00   1475.   4.449   4.639   1384.    62.6    587.
NM    1065.    7.00    600.   3.656   3.985    744.    56.3    699.
AZ    1945.    7.00   1173.   4.300   3.635   1230.    60.3    632.
UT    1126.    7.00    572.   3.745   2.611    666.    50.8    591.
NV     527.    6.00    354.   5.215   2.302    412.    67.2    782.
WN    3443.    9.00   1966.   4.476   3.942   1757.    57.1    510.
OR    2182.    7.00   1360.   4.296   4.083   1331.    62.3    610.
CA   20468.    7.00  12130.   5.002   9.794  10730.    59.3    524.


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

(File: TEST1.F90)
PROGRAM test1

!     Test vmove.
!     The LS fit to the data is:  Y = 1 + X + X**2 + X**3
!     i.e. all of the regression coefficients should be exactly 1.0.
!     An extra variable equal to X + X**2 is inserted to give a singularity.

USE lsq
IMPLICIT NONE

REAL ( lsq_kind ) :: x(11,4), y(11), xrow(5), yelem, beta(5), one = 1.D0
INTEGER           :: nvar, i, j, ifault
LOGICAL           :: fit_const, lindep(5)
CHARACTER (LEN=1) :: key

DATA y/65647., 70638., 75889., 81399., 87169., 93202., 99503., &
       106079., 112939., 120094., 127557./

WRITE(*, *)'Fitting nasty cubic'
WRITE(*, *)'1st 4 regression coefficients should equal 1.0'
WRITE(*, *)'Last variable = X + X^2, to introduce a deliberate singularity'
WRITE(*, *)

DO i = 1, 11
  x(i,1) = i + 39
  x(i,2) = x(i,1)**2
  x(i,3) = x(i,1) * x(i,2)
  x(i,4) = x(i,1) + x(i,2)
END DO

!     Use includ to form the orthogonal reduction.

nvar = 4
fit_const = .true.
CALL startup(nvar, fit_const)

DO i = 1, 11
  xrow(1) = one
  DO j = 1, 4
    xrow(j+1) = x(i,j)
  END DO
  yelem = y(i)
  CALL includ(one, xrow, yelem)
END DO

!     CALL tolset to set tolerances, then sing to detect singularities.

CALL tolset
CALL sing(lindep, ifault)
WRITE(*, *)'From routine SING, IFAULT =', ifault
WRITE(*, *)'Array LINDEP:', lindep
WRITE(*, *)
WRITE(*, *)'sserr = ', sserr, '   Should be 286.000'

!     Calculate residual sums of squares (RSS).

      CALL ss

!     Calculate regression coefficients, using vmove to cycle through
!     the order of the variables.

    DO i = 1, 6
      CALL regcf(beta, 5, ifault)
      WRITE(*, 920) vorder
920   FORMAT(1x, 'Variable order:'/ 1x, i10, 4i15)
      WRITE(*, 900) beta
900   FORMAT(' Regn. coeffs.:'/ 4x, 5g15.7)
      WRITE(*, 910) d, r, rhs, rss
910   FORMAT(' d: ', 5g15.7/ ' r: ', 4g15.7/19x,3g15.7/34x,2g15.7/ &
          49x,g15.7/ ' rhs: '/4x, 5g15.7/ ' rss: '/4x, 5g15.7//)

WRITE(*, *)'Press ENTER to continue'
READ(*, '(a)') key

!     Move variable in 1st position to the last.

      CALL vmove(1, 5, ifault)
    END DO
END PROGRAM test1
