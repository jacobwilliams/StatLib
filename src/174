      SUBROUTINE NONPAR(ISTAT, IT, NPOP, IN, IP, Y, NAR, ST, IDF,
     *   VM, AR, W, VP, VR, VY, VV, IFAULT)
C
C     ALGORITHM AS 174  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Computes the distribution-free test statistics for either the
C     multivariate multi-sample rank sum test (MMRST) or the multi-
C     variate multi-sample median test (MMMT)
C
      INTEGER IN(NPOP)
      REAL Y(IT,IP), VM(IP,IP), AR(NAR), W(IP), VP(NPOP,IP), VR(IP),
     *     VY(IT), VV(NAR)
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/
C
      IFAULT = 3
      IF (ISTAT .NE. 0 .AND. ISTAT .NE. 1) GO TO 22
      IFAULT = 4
      IF (NPOP .LT. 2) GO TO 22
      IFAULT = 5
      IF (IP .LT. 1) GO TO 22
      IFAULT = 6
      ISUM = 0
      DO 1 I = 1, NPOP
    1 ISUM = ISUM + IN(I)
      IF (ISUM .NE. IT) GO TO 22
      IFAULT = 7
      JP = IP * (IP+1) / 2
      IF (NAR .NE. JP) GO TO 22
      IFAULT = 0
C
      DO 6 J = 1, IP
	DO 5 I = 1, IT
	  R = HALF
	  DO 4 II = 1, IT
	    IF (Y(I,J) - Y(II,J)) 4, 2, 3
    2       R = R + HALF
	    GO TO 4
    3       R = R + ONE
    4     CONTINUE
	  VY(I) = R
    5   CONTINUE
	DO 6 I = 1, IT
	  Y(I,J) = VY(I)
    6 CONTINUE
C
C     Matrix Y now has the ranks replacing the observed values
C
      B = IT
      IF (ISTAT .EQ. 1) GO TO 9
      A = HALF * B
      DO 8 J = 1, IP
	DO 8 IJ = 1, IT
	  IF (Y(IJ,J) .LE. A) GO TO 7
	  Y(IJ,J) = ZERO
	  GO TO 8
    7     Y(IJ,J) = ONE
    8 CONTINUE
C
C     If MMMT test is selected (ISTAT = 0), the Y-matrix will now
C     contain 1's and 0's.   1 indicates an observation <= the median;
C     0 if greater than the median.
C
    9 DO 13 I = 1, IP
	IA = 0
	DO 11 IB = 1, NPOP
	  IST = IN(IB)
	  CST = ZERO
	  DO 10 IC = 1, IST
	    IA = IA + 1
	    CST = CST + Y(IA,I)
   10     CONTINUE
	  VP(IB,I) = CST
   11   CONTINUE
	DO 13 IR = 1, IP
	  V = ZERO
	  DO 12 JR = 1, IT
   12     V = V + Y(JR,I) * Y(JR,IR)
	  VM(I,IR) = V
   13 CONTINUE
C
C     Matrix VP contains the sum of the Y-matrix elements for each
C     population and each multivariate response.
C
C     Matrix VM contains the sum of squares and cross-products of
C     the columns of Y, that is VM = Y'Y.
C
      DO 14 I = 1, IP
   14 VR(I) = ZERO
      DO 15 I = 1, IP
	DO 15 IPR = 1, NPOP
	  VR(I) = VR(I) + VP(IPR,I)
   15 CONTINUE
C
C     Vector VR contains the sum of elements in each column of Y.
C
      KPR = 0
      DO 16 J = 1, IP
	DO 16 I = 1, J
	  KPR = KPR + 1
	  AR(KPR) = VM(I,J)/B - VR(I)*VR(J)/B**2
   16 CONTINUE
C
C     Vector AR contains the upper triangular portion of VM adjusted for
C     the overall means read in columnwise.   This is the upper triangle
C     of the dispersion matrix required for subroutine SYMINV.
C
      DO 17 I = 1, NPOP
	AA = IN(I)
	DO 17 J = 1, IP
	  VP(I,J) = VP(I,J)/AA - VR(J)/B
   17 CONTINUE
C
C     Matrix VP now contains the sum of the Y-matrix elements for each
C     population and multivariate response adjusted for the overall
C     population response mean.
C
      CALL SYMINV(AR, IP, VV, W, NULLTY, IFAULT)
      IF (IFAULT .NE. 0) GO TO 22
      K = 0
      DO 18 J = 1, IP
	DO 18 I = 1, J
	  K = K + 1
	  VM(I,J) = VV(K)
	  VM(J,I) = VV(K)
   18 CONTINUE
C
C     Matrix VM now is the inverse of the dispersion matrix.
C
      ST = ZERO
      DO 21 IPR = 1, NPOP
	TSTAT = ZERO
	DO 20 JPR = 1, IP
	  ACC = ZERO
	  DO 19 KPR = 1, IP
   19     ACC = ACC + VM(KPR,JPR) * VP(IPR,KPR)
	  TSTAT = TSTAT + ACC * VP(IPR,JPR) * IN(IPR)
   20   CONTINUE
	ST = ST + TSTAT
   21 CONTINUE
      IDF = IP * (NPOP - 1)
C
   22 RETURN
      END
