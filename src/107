      SUBROUTINE SEQUAL(BDRY, LBDRY, VEC, LVEC, NSTRD, LNSTRD, NDECNS,
     *    NSTOP, PROBS, IPROB1, IPROB2, TSTEPS, INCS, UNITT, RESULT,
     *    IFAULT)
C
C     ALGORITHM AS 107  APPL.STATIST. (1977), VOL.26, NO.1
C
C     Calculates the first (M-1) operating characteristics and the
C     average sampling number, given the parameter, for a general closed
C     sequential sampling plan with M terminal decisions.
C
      INTEGER LBDRY, LVEC, LNSTRD, NSTRD(LNSTRD), NDECNS, NSTOP, IPROB1,
     *    IPROB2, INCS(IPROB2), IFAULT
      REAL BDRY(LBDRY), VEC(LVEC), PROBS(IPROB2), TSTEPS(IPROB2), UNITT,
     *    RESULT(NDECNS)
C
C     Local variables
C
      INTEGER NL, NPAIRS, NBDRYS, JJ, N, IBASEC, IBASEL, NPC, NPL, KPC,
     *    L1, L2, J, IT, ITLIM, IBC, KKPC, KNPC, IR, IB, KKPL, KPL,
     *    KNPL, L, K, NS, LL, JPOS
      LOGICAL FLAG
      REAL U2, Z, T, TLIM, PN, PXGTH
      REAL SMALL, ZERO, HALF, ONE
C
C     SMALL is a machine-dependent constant used to test that TSTEPS
C     values divided by UNITT equal INCS values.
C
      DATA SMALL /1.E-06/, ZERO /0.0/, HALF /0.5/, ONE /1.0/
C
      NL = (LVEC / 2) / NDECNS
      NPAIRS = NDECNS - 1
      NBDRYS = 2 * NPAIRS
C
C     Test for admissibility of parameters.
C
      IFAULT = 1
      IF (NDECNS .LE. 1) RETURN
      IFAULT = 2
      IF (NSTOP .LE. 0) RETURN
      IFAULT = 3
      IF (IPROB1 .LE. 0) RETURN
      IFAULT = 4
      IF (IPROB2 .LE. IPROB1) RETURN
      IFAULT = 5
      IF (LVEC .LE. 0) RETURN
      IFAULT = 6
      IF (LNSTRD .LT. NBDRYS) RETURN
      IFAULT = 7
      IF (LBDRY .LT. NBDRYS * NSTOP) RETURN
C
      U2 = HALF * UNITT
      Z = UNITT * U2
      JJ = 2
C
C     Test that sampling plan is closed at observation count NSTOP.
C
      IFAULT = 10
      DO 10 N = 1, NPAIRS
	IF (BDRY(JJ) - BDRY(JJ-1) .GT. Z) GO TO 190
	JJ = JJ + 2
   10 CONTINUE
C
C     Test that TSTEPS values divided by UNITT equal INCS values.
C
      IFAULT = 100
      DO 20 N = IPROB1, IPROB2
	T = INCS(N)
	IF (ABS(TSTEPS(N) / UNITT - T) .GT. SMALL) GO TO 190
   20 CONTINUE
C
      IFAULT = 0
      IBASEC = 0
      IBASEL = 0
      NPC = 0
      NPL = NPAIRS
      N = NSTOP - 1
      FLAG = .FALSE.
      KPC = 0
      L1 = NPL + 1
      L2 = NPL + NPAIRS
      DO 30 J = L1, L2
   30 NSTRD(J) = 0
C
C     Enter backward instruction cycle.
C
   40 KPL = KPC
      IF (N .GT. 0) GO TO 50
      FLAG = .TRUE.
      T = ZERO
      GO TO 80
   50 KPC = KPC + NBDRYS
      IT = IBASEC
      ITLIM = IT + NL
      IBC = IBASEC
      J = 1
      KKPC = KPC
      KNPC = NPC
C
C     Enter loop through boundary pairs for sample size N.
C
   60 KKPC = KKPC + 2
      TLIM = BDRY(KKPC) - U2
      T = BDRY(KKPC-1)
C
C     Enter cycle to compute results for each T in continuation region
C     at sample size N.
C
   70 T = T + UNITT
      IF (T .GT. TLIM) GO TO 170
      IT = IT + 1
C
C     Test that number of continuation states does not exceed the limit.
C
      IF (IT .GT. ITLIM) GO TO 180
C
   80 DO 90 IR = 1, NPAIRS
   90 RESULT(IR) = ZERO
      RESULT(NDECNS) = ONE
      IB = IBASEL
      KKPL = KPL
      KNPL = NPL
      L1 = IPROB1
C
C     Loop through reachable continuation and stopping regions at
C     sample size (N+1).
C
      DO 140 K = 1, NPAIRS
	KKPL = KKPL + 2
	KNPL = KNPL + 1
	NS = NSTRD(KNPL)
	Z = T - BDRY(KKPL-1)
	JJ = Z / UNITT + SIGN(HALF, Z - U2)
	JJ = JJ + IB - 1
	Z = U2 - Z
	PN = ZERO
C
C     Loop through stop region for decision K.
C
	DO 100 L = L1, IPROB2
	  IF (TSTEPS(L) .GT. Z) GO TO 110
	  PN = PN + PROBS(L)
  100   CONTINUE
	RESULT(K) = RESULT(K) + PN
	GO TO 150
  110   L1 = L
	RESULT(K) = RESULT(K) + PN
	IF (NS .EQ. 0) GO TO 140
	Z = BDRY(KKPL) - T - U2
C
C     Loop through continuation region between decisions K and (K+1).
C
	DO 120 LL = L1, IPROB2
	  IF (TSTEPS(LL) .GT. Z) GO TO 130
	  JPOS = (JJ + INCS(LL)) * NDECNS
	  PXGTH = PROBS(LL)
	  DO 120 IR = 1, NDECNS
	    JPOS = JPOS + 1
	    RESULT(IR) = RESULT(IR) + PXGTH * VEC(JPOS)
  120   CONTINUE
	GO TO 150
  130   IB = IB + NS
	L1 = LL
  140 CONTINUE
C
C     Test for completion of calculations.
C
  150 IF (FLAG) RETURN
C
C     Store results for (N, T).
C
      JPOS = (IT - 1) * NDECNS
      DO 160 IR = 1, NDECNS
	JPOS = JPOS + 1
	VEC(JPOS) = RESULT(IR)
  160 CONTINUE
      GO TO 70
C
  170 KNPC = KNPC + 1
      NSTRD(KNPC) = IT - IBC
      IBC = IT
      J = J + 1
      IF (J .LE. NPAIRS) GO TO 60
C
C     End of calculations for sample size N.
C
      IBASEL = IBASEC
      IBASEC = NL - IBASEC
      NPL = NPC
      N = N - 1
      GO TO 40
C
  180 IFAULT = 1000
  190 IFAULT = IFAULT + N
      RETURN
      END

