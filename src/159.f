      SUBROUTINE RCONT2(NROW, NCOL, NROWT, NCOLT, JWORK, MATRIX,
     *    KEY, IFAULT)
C
C     ALGORITHM AS 159  APPL. STATIST. (1981) VOL.30, NO.1
C
C     Generate random two-way table with given marginal totals
C
C     N.B. The call to the NAG random number function G05AAF has been
C     replaced by the use of the random number generator from AS 183.
C
C     N.B. See the journal article for the significance of the lines
C     which start with C*
C
      INTEGER NROWT(NROW), NCOLT(NCOL), MATRIX(NROW, NCOL),
     *     JWORK(NCOL)
      REAL FACT(5001), DUMMY, ZERO, ONE, HALF
      LOGICAL KEY
      LOGICAL LSP, LSM
      COMMON /B/ NTOTAL, NROWM, NCOLM, FACT
C
C*      COMMON /TEMPRY/ HOP
C
      DATA MAXTOT /5000/, ZERO/0.0/, ONE/1.0/, HALF/0.5/
C
      IFAULT = 0
      IF (KEY) GO TO 103
C
C     Set KEY for subsequent calls
C
      KEY = .TRUE.
C
C     Check for faults and prepare for future calls
C
      IF (NROW .LE. 1) GO TO 212
      IF (NCOL .LE. 1) GO TO 213
      NROWM = NROW - 1
      NCOLM = NCOL - 1
      DO 100 I = 1, NROW
	IF (NROWT(I) .LE. 0) GO TO 214
  100 CONTINUE
      NTOTAL = 0
      DO 101 J = 1, NCOL
	IF (NCOLT(J) .LE. 0) GO TO 215
	NTOTAL = NTOTAL + NCOLT(J)
  101 CONTINUE
      IF (NTOTAL .GT. MAXTOT) GO TO 216
C
C     Calculate log-factorials
C
      X = ZERO
      FACT(1) = ZERO
      DO 102 I = 1, NTOTAL
	X = X + LOG(FLOAT(I))
	FACT(I+1) = X
  102 CONTINUE
C
C     -----------------------
C     Construct random matrix
C     -----------------------
C
  103 DO 105 J = 1, NCOLM
  105 JWORK(J) = NCOLT(J)
      JC = NTOTAL
C
C*      HOP = ONE
C
      DO 190 L = 1, NROWM
	NROWTL = NROWT(L)
	IA = NROWTL
	IC = JC
	JC = JC - NROWTL
	DO 180 M = 1, NCOLM
	  ID = JWORK(M)
	  IE = IC
	  IC = IC - ID
	  IB = IE - IA
	  II = IB - ID
C
C     Test for zero entries in MATRIX
C
	  IF (IE .NE. 0) GO TO 130
	  DO 121 J = M, NCOL
  121     MATRIX(L,J) = 0
	  GO TO 190
C
C     Generate pseudo-random number
C
  130	  DUMMY = RAND()
C
C     Compute conditional expected value of MATRIX(L, M)
C
  131     NLM = IA * ID / FLOAT(IE) + HALF
	  IAP = IA + 1
	  IDP = ID + 1
	  IGP = IDP - NLM
	  IHP = IAP - NLM
	  NLMP = NLM + 1
	  IIP = II + NLMP
          X = EXP(FACT(IAP) + FACT(IB+1) + FACT(IC+1) + FACT(IDP) -
     *      FACT(IE+1) - FACT(NLMP) - FACT(IGP) - FACT(IHP) - FACT(IIP))
          IF (X .GE. DUMMY) GO TO 160
	  SUMPRB = X
	  Y = X
	  NLL = NLM
	  LSP = .FALSE.
	  LSM = .FALSE.
C
C     Increment entry in row L, column M
C
  140     J = (ID - NLM) * (IA - NLM)
	  IF (J .EQ. 0) GO TO 156
	  NLM = NLM + 1
	  X = X * J / FLOAT(NLM * (II + NLM))
	  SUMPRB = SUMPRB + X
	  IF (SUMPRB .GE. DUMMY) GO TO 160
  150     IF (LSM) GO TO 155
C
C     Decrement entry in row L, column M
C
	  J = NLL * (II + NLL)
	  IF (J .EQ. 0) GO TO 154
	  NLL = NLL - 1
	  Y = Y * J / FLOAT((ID - NLL) * (IA - NLL))
	  SUMPRB = SUMPRB + Y
	  IF (SUMPRB .GE. DUMMY) GO TO 159
	  IF (.NOT. LSP) GO TO 140
	  GO TO 150
  154     LSM = .TRUE.
  155     IF (.NOT. LSP) GO TO 140
	  DUMMY = SUMPRB * RAND()
	  GO TO 131
  156     LSP = .TRUE.
	  GO TO 150
  159     NLM = NLL
C
C*          HOP = HOP * Y
C*          GO TO 161
C*160       HOP = HOP * X
C*161       MATRIX(L,M) = NLM
C
  160     MATRIX(L,M) = NLM
	  IA = IA - NLM
	  JWORK(M) = JWORK(M) - NLM
  180   CONTINUE
	MATRIX(L,NCOL) = IA
  190 CONTINUE
C
C     Compute entries in last row of MATRIX
C
      DO 192 M = 1, NCOLM
  192 MATRIX(NROW,M) = JWORK(M)
      MATRIX(NROW,NCOL) = IB - MATRIX(NROW,NCOLM)
      RETURN
C
C     Set faults
C
  212 IFAULT = 1
      RETURN
  213 IFAULT = 2
      RETURN
  214 IFAULT = 3
      RETURN
  215 IFAULT = 4
      RETURN
  216 IFAULT = 5
      RETURN
      END
