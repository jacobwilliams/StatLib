      SUBROUTINE MSTUPD(KN, KI, DK, DIST, NUM, IREF, IPUSH)
C
C     ALGORITHM AS 40  APPL. STATIST. (1971) VOL.20, NO.2
C
C     This subroutine updates the minimal spanning tree of KN points,
C     defined by KI(J) being the point nearest to J on the chain leading
C     from J to the point KN and DK(J) being the length link J to KI(J),
C     with the new point (KN+1) whose distances from the points of the
C     original M.S.T. are given in the array DIST.
C     The update version of the M.S.T. is left in the arrays KI and DK
C     suitable for re-entry to the subroutine with the next point.
C     A working space of three arrays NUM, IREF and IPUSH, each at least
C     of size KN, is required.
C
      INTEGER KN, KI(KN), NUM(KN), IREF(KN), IPUSH(KN)
      REAL DK(KN), DIST(KN)
C
      KT = KN - 1
      IP = 0
C
C     If J = IREF(I) > 0 it refers to link J to KI(J), while if < 0 it
C     refers to the link from point (KN+1) to (-J).
C
      DO 1 I = 1, KN
	IREF(I) = -I
    1 CONTINUE
C
C     NUM(I) holds the number of links pointing to I.   Hence I is an end
C     point if NUM(I) = 0.
C
      DO 2 I = 1, KT
	K = KI(I)
	NUM(K) = NUM(K) + 1
    2 CONTINUE
C
C     Start of algorithm
C     NUM(I) = 1 if the link I to KI(I) is included in the new M.S.T.
C     otherwise it is 0 or < 0.
C
      DO 3 NI = 1, KT
	IF (NUM(NI) .NE. 0) GO TO 3
	I = NI
   31   J = KI(I)
	IF (DIST(I) .GT. DK(I)) GO TO 34
	IF (DK(I) .GE. DIST(J)) GO TO 32
	IREF(J) = I
	DIST(J) = DK(I)
   32   IF (IREF(I) .GT. 0) GO TO 33
	IP = IP + 1
        IPUSH(IP) = -IREF(I)
	GO TO 36
   33   K = IREF(I)
	NUM(K) = 1
	GO TO 36
   34   IF (DIST(I) .GE. DIST(J)) GO TO 35
	IREF(J) = IREF(I)
	DIST(J) = DIST(I)
   35   NUM(I) = 1
   36   NUM(J) = NUM(J) - 1
	IF (NUM(J) .NE. 0) GO TO 3
	I = J
	NUM(I) = -1
	GO TO 31
    3 CONTINUE
C
      K = -IREF(KN)
      IF (K .GE. 0) GO TO 42
      K = -K
      NUM(K) = 1
C
C     KI is updated so that all the chains point towards (KN+1).
C
   41 K = IPUSH(IP)
      IP = IP - 1
   42 KT = KN + 1
      X = DIST(K)
   43 KP = KI(K)
      KI(K) = KT
      Y = DK(K)
      DK(K) = X
      X = Y
      KT = K
      K = KP
      IF (NUM(KT) .GT. 0) GO TO 43
      IF (IP .NE. 0) GO TO 41
C
      RETURN
      END
