      SUBROUTINE DENEST(DT, NDT, DLO, DHI, WINDOW, FT, SMOOTH,
     *	NFT, ICAL, IFAULT)
      REAL DT(NDT), FT(NFT), SMOOTH(NFT)
C
C     ALGORITHM AS 176  APPL. STATIST. (1982) VOL.31, NO.1
C     Modified using AS R50 (Appl. Statist. (1984))
C
C     Find density estimate by kernel method using Gaussian kernel.
C     The interval on which the estimate is evaluated has end points
C     DLO and DHI.   If ICAL is not zero then it is assumed that the
C     routine has been called before with the same data and end points
C     and that the array FT has not been altered.
C
C     Auxiliary routines called: FORRT & REVRT from AS 97
C
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/, SIX/6.0/, THIR2/32.0/
      DATA BIG/30.0/, KFTLO/5/, KFTHI/11/
C
C     The constant BIG is set so that exp(-BIG) can be calculated
C     without causing underflow problems and can be considered = 0.
C
C     Initialize and check for valid parameter values.
C
      IF (WINDOW .LE. ZERO) GO TO 92
      IF (DLO .GE. DHI) GO TO 93
      II = 2**KFTLO
      DO 1 K = KFTLO, KFTHI
	IF (II .EQ. NFT) GO TO 2
	II = II + II
    1 CONTINUE
      IFAULT = 1
      RETURN
    2 STEP = (DHI - DLO) / FLOAT(NFT)
      AINC = ONE / (NDT * STEP)
      NFT2 = NFT / 2
      HW = WINDOW / STEP
      FAC1 = THIR2 * (ATAN(ONE) * HW / NFT) ** 2
      IF (ICAL .NE. 0) GO TO 10
C
C     Discretize the data
C
      DLO1 = DLO - STEP * HALF
      DO 3 J = 1, NFT
    3 FT(J) = ZERO
      DO 4 I = 1, NDT
	WT = (DT(I) - DLO1) / STEP
	JJ = INT(WT)
	IF (JJ .LT. 1 .OR. JJ .GT. NFT) GO TO 4
	WT = WT - JJ
	WINC = WT * AINC
	KK = JJ + 1
	IF (JJ .EQ. NFT) KK = 1
	FT(JJ) = FT(JJ) + AINC - WINC
	FT(KK) = FT(KK) + WINC
    4 CONTINUE
C
C     Transform to find FT.
C
      CALL FORRT(FT, NFT)
C
C     Find transform of density estimate.
C
   10 JHI = SQRT(BIG / FAC1)
      JMAX = MIN(NFT2 - 1, JHI)
      SMOOTH(1) = FT(1)
      RJ = ZERO
      DO 11 J = 1, JMAX
	RJ = RJ + ONE
	RJFAC = RJ * RJ * FAC1
	BC = ONE - RJFAC / (HW * HW * SIX)
	FAC = EXP(-RJFAC) / BC
	J1 = J + 1
	J2 = J1 + NFT2
	SMOOTH(J1) = FAC * FT(J1)
	SMOOTH(J2) = FAC * FT(J2)
   11 CONTINUE
C
C     Cope with underflow by setting tail of transform to zero.
C
      IF (JHI + 1 - NFT2) 21, 23, 20
   20 SMOOTH(NFT2 + 1) = EXP(-FAC1 * FLOAT(NFT2)**2) * FT(NFT2 + 1)
      GO TO 24
   21 J2LO = JHI + 2
      DO 22 J1 = J2LO, NFT2
	J2 = J1 + NFT2
	SMOOTH(J1) = ZERO
	SMOOTH(J2) = ZERO
   22 CONTINUE
   23 SMOOTH(NFT2 + 1) = ZERO
C
C     Invert Fourier transform of SMOOTH to get estimate and eliminate
C     negative density values.
C
   24 CALL REVRT(SMOOTH, NFT)
      DO 25 J = 1, NFT
   25 IF (SMOOTH(J) .LT. ZERO) SMOOTH(J) = ZERO
      IFAULT = 0
      RETURN
C
   92 IFAULT = 2
      RETURN
   93 IFAULT = 3
      RETURN
      END

