	SUBROUTINE BTREE(FIND, LOC, NODE, LLINK, RLINK, BF, MAXT, NTREE,
     +		ROOT, LAST, IFAULT)
C
C	ALGORITHM AS219  APPL. STATIST. (1986) VOL. 35, NO. 2
C
C	BTREE returns the location (LOC) of an item (FIND) in a
C	height-balanced tree composed of NODE, LLINK, RLINK  and BF.
C	The tree is rooted at ROOT and contains NTREE nodes.   LAST
C	points to the next available free location in the tree.
C
C***	N.B. The user is expected to provided a function ICOMP to replace
C	     the example provided here. ***
C
	INTEGER LLINK(MAXT), RLINK(MAXT), BF(MAXT)
	INTEGER NTREE, T, S, PNTR, R, A, MAXT, K, ROOT, ICOMP
	REAL	NODE(MAXT), FIND
C
	IFAULT = 0
	IF (NTREE .GE. 1) GO TO 20
C
C	First item in tree
C
	LAST = LAST + 1
	IF (LAST .GT. MAXT) THEN
	  IFAULT = 1
	  RETURN
	END IF
	NODE(LAST) = FIND
	LLINK(LAST) = 0
	RLINK(LAST) = 0
	BF(LAST) = 0
	ROOT = LAST
	NTREE = 1
	LOC = LAST
	RETURN
C
C	Search through tree for current item.
C	If found, return.    PNTR moves through the tree.   S points to
C	the place from which re-balancing may be needed, and T points
C	to the parent of S.
C
   20	T = 0
	S = ROOT
	PNTR = S
   30	K = ICOMP(FIND, NODE(PNTR))
	IF (K .NE. 0) GO TO 40
	LOC = PNTR
	RETURN
   40	IF (K .GT. 0) THEN
	  LOC = RLINK(PNTR)
	ELSE IF (K .LT. 0) THEN
	  LOC = LLINK(PNTR)
	END IF
	IF (LOC .NE. 0) GO TO 60
	NTREE = NTREE + 1
	LAST = LAST + 1
	IF (LAST .LE. MAXT) GO TO 50
	IFAULT = 1
	RETURN
C
   50	LOC = LAST
	IF (K .GT. 0) THEN
	  RLINK(PNTR) = LOC
	ELSE IF (K .LT. 0) THEN
	  LLINK(PNTR) = LOC
	END IF
	GO TO 100
C
   60	IF (BF(LOC) .EQ. 0) GO TO 70
	T = PNTR
	S = LOC
   70	PNTR = LOC
	GO TO 30
C
C	Add new node to tree and adjust balance factors
C
  100	NODE(LOC) = FIND
	LLINK(LOC) = 0
	RLINK(LOC) = 0
	BF(LOC) = 0
	A = ICOMP(FIND, NODE(S))
	IF (A .GT. 0) THEN
	  PNTR = RLINK(S)
	ELSE IF (A .LT. 0) THEN
	  PNTR = LLINK(S)
	END IF
	R = PNTR
C
  110	IF (PNTR .EQ. LOC) GO TO 130
	K = ICOMP(FIND, NODE(PNTR))
	IF (K .GT. 0) GO TO 120
	BF(PNTR) = -1
	PNTR = LLINK(PNTR)
	GO TO 110
  120	BF(PNTR) = 1
	PNTR = RLINK(PNTR)
	GO TO 110
C
  130	IF (BF(S) .NE. 0) GO TO 140
	BF(S) = A
	RETURN
  140	IF (BF(S) .NE. -A) GO TO 150
	BF(S) = 0
	RETURN
  150	IF (BF(R) .NE. A) GO TO 190
C
C	Single rotation
C
	PNTR = R
	IF (A .NE. -1) GO TO 170
	LLINK(S) = RLINK(R)
	RLINK(R) = S
	GO TO 180
  170	RLINK(S) = LLINK(R)
	LLINK(R) = S
  180	BF(S) = 0
	BF(R) = 0
	GO TO 240
C
C	Double rotation
C
  190	IF (A .NE. -1) GO TO 200
	PNTR = RLINK(R)
	RLINK(R) = LLINK(PNTR)
	LLINK(PNTR) = R
	LLINK(S) = RLINK(PNTR)
	RLINK(PNTR) = S
	GO TO 210
  200	PNTR = LLINK(R)
	LLINK(R) = RLINK(PNTR)
	RLINK(PNTR) = R
	RLINK(S) = LLINK(PNTR)
	LLINK(PNTR) = S
  210	BF(S) = 0
	BF(R) = 0
	IF (BF(PNTR) .NE. A) GO TO 220
	BF(S) = -A
	GO TO 230
  220	IF (BF(PNTR) .EQ. -A) BF(R) = A
  230	BF(PNTR) = 0
  240	IF (T .NE. 0) GO TO 250
	ROOT = PNTR
	RETURN
  250	IF (S .NE. RLINK(T)) THEN
	  LLINK(T) = PNTR
	ELSE
	  RLINK(T) = PNTR
	END IF
C
	RETURN
	END
C
C
C
	SUBROUTINE TREEP(PERM, LLINK, RLINK, MAXT, NTREE, ROOT, STACK,
     +		MAXSTA, IFAULT)
	INTEGER NTREE, PERM(NTREE), LLINK(MAXT), RLINK(MAXT)
	INTEGER PNTR, STACK(MAXSTA), SPTR, ROOT
C
C	Returns a permutation vector (PERM) to the tree contained in
C	NODE, LLINK, RLINK and rooted at ROOT
C
	IFAULT = 0
	IF (NTREE .LT. 1) RETURN
	SPTR = 0
	PNTR = ROOT
	I = 1
   10	IF (PNTR .EQ. 0) GO TO 20
C
C	Traverse left subtree
C
	SPTR = SPTR + 1
	IF (SPTR .GT. MAXSTA) GO TO 900
	STACK(SPTR) = PNTR
	PNTR = LLINK(PNTR)
	GO TO 10
   20	IF (SPTR .LT. 1) RETURN
C
C	Process NODE and traverse right subtree
C
	PNTR = STACK(SPTR)
	SPTR = SPTR - 1
	PERM(I) = PNTR
	I = I + 1
	PNTR = RLINK(PNTR)
	GO TO 10
C
  900	IFAULT = 1
	RETURN
	END
C
C
C
	INTEGER FUNCTION SEARCH(FIND, NODE, LLINK, RLINK, MAXT, NTREE,
     +		ROOT)
	INTEGER LLINK(MAXT), RLINK(MAXT)
	INTEGER NTREE, PNTR, MAXT, ROOT, ICOMP
	REAL	NODE(MAXT), FIND
C
C	Returns the location in the tree of the item FIND.
C	Returns 0 if FIND is not in the tree.
C
	SEARCH = 0
	IF (NTREE .LT. 0) RETURN
	PNTR = ROOT
   10	K = ICOMP(FIND, NODE(PNTR))
	IF (K .EQ. 0) GO TO 20
	IF (K .LT. 0) THEN
	  PNTR = LLINK(PNTR)
	ELSE
	  PNTR = RLINK(PNTR)
	END IF
	IF (PNTR .NE. 0) GO TO 10
C
   20	SEARCH = PNTR
	RETURN
	END
C
C
C
	INTEGER FUNCTION ICOMP(X, Y)
	REAL	X, Y, EPS, Z
C
C	The body of ICOMP will be supplied by the user.
C	ICOMP will return 0 if X and Y are equal to the accuracy
C	needed, and return -1 if X < Y, and +1 if X > Y.
C	The example given tests whether X and Y differ by more than
C	0.001 of the larger in magnitude.
C
	DATA ZERO, ONE, EPS /0.0, 1.0, 0.001/
C
C	Avoid comparison of very different X and Y
C
	IF (X .LT. ONE .OR. Y .GT. -ONE) GO TO 10
	ICOMP = 1
	RETURN
   10	IF (X .GT. -ONE .OR. Y .LT. ONE) GO TO 20
	ICOMP = -1
	RETURN
C
   20	Z = MAX(ABS(X), ABS(Y))
	ICOMP = 0
	IF (Z .EQ. ZERO) RETURN
	Z = (X - Y) / Z
	IF (ABS(Z) .LT. EPS) RETURN
	ICOMP = 1
	IF (Z .GT. ZERO) RETURN
	ICOMP = -1
	RETURN
	END
