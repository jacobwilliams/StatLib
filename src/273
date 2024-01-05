Notes on AS 273 "Comparing subsets of regressor variables" by Alan Miller

This algorithm allows two regression subsets to be compared.   It uses the
method of Spj0tvoll (Ann. Math. Statist. vol.43, 1076-1088, 1972).   The
algorithm will typically be used after the use of some subset selection
procedure, to make comparisons between subsets.   It is a Scheffe-type
procedure which gives a confidence level which applies simultaneously to
all pairs of subsets from the same data set.   The subsets may be of different
numbers of variables, and there is no requirement that some or all of the
variables in one subset must be in the other.   For more details of the
algorithm see chapter 4 of `Subset selection in regression' by Alan Miller,
Chapman & Hall, 1990 (ISBN 0 412 35380 6).

The algorithm calculates upper and lower confidence limits for the difference
in regression sums of squares for the two subsets.   If both limits are
positive then the first subset has a significantly larger regression sum of
squares than the first.   If the range includes zero, the subsets do not
differ significantly.  

The user must enter F-values for the confidence levels required.   The number
of degrees of freedom for the numerator is equal to the total number of
variables less any (such as the optional constant) forced into all subsets.
The number of degrees of freedom for the denominator is equal to the number of
degrees of freedom of the residual variance, or NCASES - NP.

Before calling this algorithm, the user must have first performed an
orthogonal reduction of the data using Gentleman's algorithm (AS 75 or 274),
or an equivalent algorithm.

The calling arguments are (real/dp means real or double precision):
NP     Input, integer   No. of variables in the orthogonal reduction
NRBAR  Input, integer   No. of elements in array RBAR (at least NP(NP-1)/2).
D      Input, real/dp   Row multipliers from AS 75 or AS 274
              array
RBAR   Input, real/dp   Upper triangle of Cholesky factor without the diagonal
              array        elements (which are implicit 1's)
THETAB Input, real/dp   The least-squares projections of the dependent variable
              array        on the NP orthogonal directions
RSS    Input, real/dp   RSS(I) is the residual sum of squares with the first I
              array        variables in the model (output from SS of AS 274)
TOL    Input, real/dp   Array of tolerances from TOLSET in AS 274
              array
NCASES Input, integer   The number of cases on which the orthogonal reduction
                        was based
VORDER Input, integer   An array of integer identifiers, one per variable
              array
NFORCE Input, integer   The number of variables which were forced into all
                        subsets in the subset selection procedure (often = 1
                        for the constant term in all models)
LIST1  Input, integer   Array of identifying integers (see VORDER above) for
              array        the variables in the first subset to be compared.
                        Include variables which are forced in.
N1     Input, integer   The number of variables in LIST1.
LIST2  Input, integer   The list of variables in the second subset.
              array
N2     Input, integer   The number of variables in the second subset.
FVAL   Input, integer   Array of F-values for the required confidence levels.
              array
NF     Input, integer   The number of confidence levels required, usually only
                        1 or 2.
A1     Output, real/dp  Array of lower confidence limits corresponding to the
              array                input F values.
A2     Output, real/dp  Array of the corresponding upper confidence limits.
              array
WK     --, real/dp      Workspace, must be dimensioned at least 2N(N+1) where
              array        N = the number of variables which are in one subset
                        but not both.
DIMWK  Input, integer   The dimension of WK in the calling program.
IWK    --, integer      Workspace, must be dimensioned at least N.
              array
DIMIWK Input, integer   The dimension of IWK in the calling program.
IFAULT Output, integer  Error indicator
                        = 0, no error detected
                        = -NEED, insufficient workspace provided for WK,
                                 increase to NEED.
                        = 1  NP < 2
                        = 2  NRBAR < NP(NP-1)/2
                        = 4  NCASES <= NP
                        = 8  either N1 or N2 exceeds NP
                        = 16 either N1 or N2 is negative
                        = 32 NF < 1
                        = 64 either NFORCE < 0 or NFORCE >= NP
                        = 128 LIST1 = LIST2, i.e. identical subsets
                        = 129 error in both lists, e.g. a repeated number or
                              a number which is not in array VORDER
                        = 130 as for error 129 but only in LIST1
                        = 131 as for error 129 but only in LIST2
                        = 132 inadequate integer work space

Auxiliary routines required: TRED2 and TQL2 from EISPACK, or AS 60 to calculate
the eigenvalues and vectors for a symmetric matrix.   REORDR and VMOVE from
algorithm AS 274.

Alan Miller
25 November 1991

Correction to the value of L for address in WK passed to routines TRED2 and
TQL2 made 1 April 2003.   Also IER changed to IFAULT throughout.
Thanks to Ron Prigeon for pointing these errors out.

C----------------------------------------------------------------------

      SUBROUTINE CONREG(NP, NRBAR, D, RBAR, THETAB, RSS, TOL, NCASES,
     +    VORDER, NFORCE, LIST1, N1, LIST2, N2, FVAL, NF, A1, A2, WK,
     +    DIMWK, IWK, DIMIWK, IFAULT)
C
C     ALGORITHM AS273  APPL. STATIST. (1992) VOL. 41, NO. 2
C
C     Calculate the Spj0tvoll confidence limits for the difference in
C     regression sums of squares of two subsets of variables.
C
C     Auxiliary routines required: TRED2, TQL2 (from EISPACK or AS 60),
C     and REORDR, VMOVE (from AS 274).
C
      INTEGER NP, NRBAR, NCASES, VORDER(NP), NFORCE, N1, LIST1(N1), N2,
     +    LIST2(N2), NF, DIMWK, DIMIWK, IWK(DIMIWK)
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), RSS(NP), TOL(NP),
     +    FVAL(NF), A1(NF), A2(NF), WK(DIMWK)
C
C     Local variables
C
      INTEGER N0, I, L, J, K, N10, N20, N10N20, NEED, POS2, INC, POS1,
     +    ROW, COL, ID, ITHETA, IRSS, ITOL, ROWPN0, ROWL1, NVAR, IZ,
     +    M, IVAL, IMAX, IMIN
      DOUBLE PRECISION ZERO, ONE, HALF, SCALE, TEMP, SUM, DIFF, SUME,
     +    EMAX, EMIN, VAR, TEMPS, X, X1, X2, XTOL, A, B, GAMMA2
C
      DATA ZERO/0.D0/, ONE/1.D0/, HALF/0.5D0/, XTOL/0.001D0/
C
C     Some checks.
C
      IFAULT = 0
      IF (NP .LT. 2) IFAULT = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IFAULT = IFAULT + 2
      IF (NCASES .LE. NP) IFAULT = IFAULT + 4
      IF (N1 .GT. NP .OR. N2 .GT. NP) IFAULT = IFAULT + 8
      IF (N1 .LT. 0 .OR. N2 .LT. 0) IFAULT = IFAULT + 16
      IF (NF .LT. 1) IFAULT = IFAULT + 32
      IF (NFORCE .LT. 0 .OR. NFORCE .GE. NP) IFAULT = IFAULT + 64
      IF (IFAULT .NE. 0) RETURN
C
C     Find which variable numbers are on both lists and move them to
C     the start.
C
      N0 = 0
      IF (N1 .EQ. 0 .OR. N2 .EQ. 0) GO TO 70
      DO 60 I = 1, N1
        L = LIST1(I)
        DO 10 J = N0+1, N2
          IF (L .EQ. LIST2(J)) GO TO 20
   10   CONTINUE
        GO TO 60
   20   N0 = N0 + 1
        IF (I .EQ. N0) GO TO 40
        K = I
   30   LIST1(K) = LIST1(K-1)
        K = K - 1
        IF (K .GT. N0) GO TO 30
        LIST1(N0) = L
C
   40   IF (J .EQ. N0) GO TO 60
        K = J
   50   LIST2(K) = LIST2(K-1)
        K = K - 1
        IF (K .GT. N0) GO TO 50
        LIST2(N0) = L
   60 CONTINUE
C
C     Exit if LIST1 = LIST2.
C
   70 N10 = N1 - N0
      N20 = N2 - N0
      N10N20 = N10 + N20
      IF (N10N20 .EQ. 0) THEN
        IFAULT = 128
        RETURN
      END IF
C
C     Re-order the orthogonal reduction so that the N0 variables common
C     to both lists are at the start.
C
      IF (N0 .GT. 0) CALL REORDR(NP, NRBAR, VORDER, D, RBAR, THETAB,
     +         RSS, TOL, LIST1, N0, 1, IFAULT)
      IF (IFAULT .NE. 0) THEN
        IFAULT = 129
        RETURN
      END IF
C
C     Follow with the rest of the variables, if any, in LIST1, then any
C     any remaining from LIST2.
C
      IF (N10 .GT. 0) CALL REORDR(NP, NRBAR, VORDER, D, RBAR, THETAB,
     +         RSS, TOL, LIST1(N0+1), N10, N0+1, IFAULT)
      IF (IFAULT .NE. 0) THEN
        IFAULT = 130
        RETURN
      END IF
      IF (N20 .GT. 0) CALL REORDR(NP, NRBAR, VORDER, D, RBAR, THETAB,
     +         RSS, TOL, LIST2(N0+1), N20, N1+1, IFAULT)
      IF (IFAULT .NE. 0) THEN
        IFAULT = 131
        RETURN
      END IF
C
C     Estimate residual variance.
C
      VAR = RSS(NP) / (NCASES - NP)
C
C     If N10 > 0 and N20 > 0, calculate the product of planar rotations
C     needed to place the N10 variables from LIST1 after the N20 from
C     LIST2 instead of before them.
C
      IF (N10 .EQ. 0 .OR. N20 .EQ. 0) GO TO 400
C
C     Check that adequate workspace has been provided.
C
      NEED = 2 * N10N20 * (N10N20 + 1)
      IF (DIMWK .LT. NEED) THEN
        IFAULT = -NEED
        RETURN
      END IF
      IF (DIMIWK .LT. N10N20) THEN
        IFAULT = 132
        RETURN
      END IF
C
C     Copy the triangle of RBAR for the N10N20 rows and columns of
C     interest into WK with an N10N20 x N10N20 identity matrix added
C     as additional columns to store the product.   Thus the contents
C     of WK look like (for N10N20 = 4):
C              X  X  X  1  0  0  0
C                 X  X  0  1  0  0
C                    X  0  0  1  0
C                       0  0  0  1
C
C     As row multipliers are used by the square-root free method, the
C     reciprocals of the square roots of these multipliers are stored
C     on the diagonal instead of 1's.
C     POS2 = position within RBAR; INC = number of columns of RBAR to
C     be skipped in each row.
C
      POS2 = N0 * (NP + NP - N0 - 1)/2 + 1
      INC = NP - N0 - N10N20
      POS1 = 1
      DO 100 ROW = 1, N10N20
        DO 80 COL = ROW+1, N10N20
          WK(POS1) = RBAR(POS2)
          POS1 = POS1 + 1
          POS2 = POS2 + 1
   80   CONTINUE
        POS2 = POS2 + INC
        DO 90 COL = 1, N10N20
          IF (COL .EQ. ROW) THEN
            WK(POS1) = ONE / SQRT(D(ROW+N0))
          ELSE
            WK(POS1) = ZERO
          END IF
          POS1 = POS1 + 1
   90   CONTINUE
  100 CONTINUE
C
C     Set up dummy arrays in WK for D, THETAB, RSS and TOL, and in IWK
C     for VORDER.
C
      ID = POS1
      ITHETA = ID + N10N20
      IRSS = ITHETA + N10N20
      ITOL = IRSS + N10N20
      DO 110 ROW = 1, N10N20
        ROWPN0 = ROW + N0
        ROWL1 = ROW - 1
        WK(ID + ROWL1) = D(ROWPN0)
        WK(ITHETA + ROWL1) = ZERO
        WK(IRSS + ROWL1) = RSS(ROWPN0)
        WK(ITOL + ROWL1) = TOL(ROWPN0)
        IWK(ROW) = VORDER(ROWPN0)
  110 CONTINUE
C
C     Swap positions of the N10 and N20 variables in WK.
C
      NVAR = 2 * N10N20
      DO 120 ROW = N10, 1, -1
  120 CALL VMOVE(NVAR, DIMWK, IWK, WK(ID), WK, WK(ITHETA), WK(IRSS),
     +    ROW, N10N20+ROW-N10, WK(ITOL), IFAULT)
C
C     Scale the first N20 rows of the product and place at the start of
C     WK.   This is (P1, P2).
C
      POS2 = N10N20
      POS1 = 1
      DO 140 ROW = 1, N20
        SCALE = SQRT(WK(ID))
        ID = ID + 1
        DO 130 COL = 1, N10N20
          WK(POS1) = SCALE * WK(POS2)
          POS1 = POS1 + 1
          POS2 = POS2 + 1
  130   CONTINUE
        POS2 = POS2 + N10N20 - ROW - 1
  140 CONTINUE
C
C     Construct the matrix, Z, in WK from position IZ.
C           Z  =  ( I  0 )  -  ( P1' ) ( P1  P2 )
C                 ( 0  0 )     ( P2' )
C
      IZ = N10N20**2 + 1
      POS1 = IZ
      POS2 = IZ
      DO 170 ROW = 1, N10N20
        TEMP = ZERO
        IF (ROW .LE. N10) TEMP = ONE
        DO 160 COL = ROW, N10N20
          L = ROW
          J = COL
          DO 150 I = 1, N20
            TEMP = TEMP - WK(L)*WK(J)
            L = L + N10N20
            J = J + N10N20
  150     CONTINUE
          WK(POS1) = TEMP
          WK(POS2) = TEMP
          POS1 = POS1 + 1
          POS2 = POS2 + N10N20
          TEMP = ZERO
  160   CONTINUE
        POS1 = POS1 + ROW
        POS2 = POS1
  170 CONTINUE
C
C     Put a copy of Z at the start of WK.
C
      POS2 = IZ
      DO 180 I = 1, IZ-1
        WK(I) = WK(POS2)
        POS2 = POS2 + 1
  180 CONTINUE
C
C     Calculate the eigenvalues (stored from location IVAL in WK), and
C     the normalized eigenvectors (stored from location IZ).
C
      IVAL = 2*IZ - 1
      L = IVAL + N10N20
      CALL TRED2(N10N20, N10N20, WK, WK(IVAL), WK(L), WK(IZ))
      CALL TQL2(N10N20, N10N20, WK(IVAL), WK(L), WK(IZ), IFAULT)
C
C     Calculate the transformed projections (Spj0tvoll's gammas) and
C     store at the start of WK.   Put the sum of squares of gamma(i).
C     lambda(i) into SUME.
C
      L = IZ
      SUME = ZERO
      DO 250 I = 1, N10N20
        M = N0 + 1
        SUM = ZERO
        DO 240 J = 1, N10N20
          SUM = SUM + WK(L)*THETAB(M)*SQRT(D(M))
          L = L + 1
          M = M + 1
  240   CONTINUE
        WK(I) = SUM
        SUME = SUME + (SUM*WK(IVAL+I-1))**2
  250 CONTINUE
C
C     Find the max. and min. non-zero eigenvalues.
C
      EMAX = WK(IVAL)
      EMIN = EMAX
      IMAX = 1
      IMIN = 1
      L = IVAL + 1
      DO 260 I = 2, N10N20
        TEMP = WK(L)
        L = L + 1
        IF (TEMP .GT. EMAX) THEN
          EMAX = TEMP
          IMAX = I
        ELSE IF (TEMP .LT. EMIN) THEN
          EMIN = TEMP
          IMIN = I
        END IF
  260 CONTINUE
C
C     Cycle through the NF confidence levels.
C
      DO 310 K = 1, NF
        TEMP = VAR * (NP - NFORCE) * FVAL(K)
        TEMPS = SQRT(TEMP)
C
C     Find largest and smallest roots of Spj0tvoll's equation (3.2).
C
        X1 = EMIN - ABS(EMIN * WK(IMIN)) / TEMPS
        X2 = EMIN - SQRT(SUME / TEMP)
        DO 290 I = 1, 2
  270     X = HALF * (X1 + X2)
          L = IVAL
          SUM = -TEMP
          DO 280 J = 1, N10N20
            SUM = SUM + (WK(L)*WK(J)/(WK(L) - X))**2
            L = L + 1
  280     CONTINUE
          IF (SUM .LE. ZERO) THEN
            X2 = X
          ELSE
            X1 = X
          END IF
          IF (ABS(X1 - X2) .GT. XTOL) GO TO 270
          X = HALF * (X1 + X2)
          IF (I .EQ. 1) THEN
            A = - X
            X1 = EMAX + ABS(EMAX * WK(IMAX)) / TEMPS
            X2 = EMAX + SQRT(SUME / TEMP)
          ELSE
            B =  X
          END IF
  290   CONTINUE
C
C     Calculate A1 and A2.
C
        A1(K) = -TEMP
        A2(K) = TEMP
        L = IVAL
        DO 300 I = 1, N10N20
          GAMMA2 = WK(I)**2
          TEMP = WK(L)
          L = L + 1
          A1(K) = A1(K) + TEMP * GAMMA2 / (A + TEMP)
          A2(K) = A2(K) + TEMP * GAMMA2 / (B - TEMP)
  300   CONTINUE
        A1(K) = A * A1(K)
        A2(K) = B * A2(K)
  310 CONTINUE
C
      RETURN
C----------------------------------------------------------------------
C     Special cases; either N10 or N20 = 0.
C
  400 GAMMA2 = RSS(N0) - RSS(N1+N20)
      DO 410 K = 1, NF
        TEMP = SQRT(VAR * (NP - NFORCE) * FVAL(K))
        TEMPS = SQRT(GAMMA2)
        SUM = (TEMPS + TEMP)**2
        DIFF = MAX(ZERO, TEMPS-TEMP)
        DIFF = DIFF**2
        IF (N10 .EQ. 0) THEN
          A1(K) = -SUM
          A2(K) = -DIFF
        ELSE
          A1(K) = DIFF
          A2(K) = SUM
        END IF
  410 CONTINUE
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
	 T = 4.0D0 + R
	 IF (T .EQ. 4.0D0) GO TO 20
	 S = R/T
	 U = 1.0D0 + 2.0D0*S
	 P = U*P
	 R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end



      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2

   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10

   20 PYTHAG = P
      RETURN
      END



      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end



      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
