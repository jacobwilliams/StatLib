      SUBROUTINE IRRED(ITBTYP, NTAXA, NTESTS, NDIST, LSTAB, LGROUP,
     *                 GROUP, NBI, ISETS, LSETS, LTBP, IOUT, IFAULT)
C
C        ALGORITHM AS 263.1  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Constructs irredundant test sets
C
      INTEGER IFAULT, IOUT, ITBTYP, LSETS, LTBP, NBI, NDIST, NTAXA,
     *        NTESTS
      INTEGER GROUP(NTAXA), ISETS(LSETS)
      LOGICAL LGROUP, LSTAB
      INTEGER I, ITAX1, ITAX2, IT1LIM, JTLIST, KPOWR2, LDBP, NIBP, NB0
C
C        Check parameters
C
      IF ((NTAXA .LE. 0) .OR. (NTESTS .LE. 0) .OR. (NDIST .LE. 0) .OR.
     *    (NBI .LE. 0) .OR. (IOUT .LT. 0)) GO TO 40
      IFAULT = 0
C
C        Set number of integers to store each irredundant test set
C
      NIBP = (NTESTS - 1) / NBI + 1
C
C        Set number of unused bits within each bit pattern
C
      NB0 = NIBP * NBI - NTESTS
C
C        Set up Powers of 2 (2**0,1,2...) in
C        ISETS(KPOWR2...KPOWR2+NBI-1)
C
      KPOWR2 = LSETS - NBI + 1
      ISETS(KPOWR2) = 1
      DO 10 I = KPOWR2 + 1, LSETS
         ISETS(I) = ISETS(I - 1) * 2
   10 CONTINUE
C
C        Initialise taxon indices and limits
C
      ITAX1 = ITBTYP
      IT1LIM = NTAXA
      IF (ITBTYP .GT. 0) THEN
         IF (LGROUP) THEN
C
C        Find 1st taxon of group ITBTYP (to be distinguished from
C        other groups)
C
            DO 20 ITAX1 = 1, NTAXA
               IF (GROUP(ITAX1) .EQ. ITBTYP) GO TO 30
   20       CONTINUE
C
C        No taxa in group
C
            IFAULT = 3
            RETURN
         ELSE
C
C        Test sets to distinguish an individual taxon
C
            IF ((ITBTYP .LT. 0) .OR. (ITBTYP .GT. NTAXA)) GO TO 40
            IT1LIM = ITAX1
         END IF
      END IF
C
C        Form irredundant sets
C
   30 LDBP = NIBP
C
C        Initialise b.p. (i.e. bit pattern) of essential tests
C
      CALL SETVAL(ISETS(1), 0, NIBP)
      I = ITAX1
      ITAX2 = 0
      CALL IRRDTB(ITBTYP, NTAXA, NTESTS, NDIST, LSTAB, LGROUP, GROUP,
     *            NBI, ISETS, LSETS, IOUT, IFAULT, ITAX1, ITAX2, IT1LIM,
     *            NIBP, NB0, KPOWR2, 0, NIBP, LDBP, JTLIST, KPOWR2 - 1,
     *            .FALSE.)
      IF (IFAULT .NE. 0) RETURN
      ITAX1 = I
      ITAX2 = 0
      CALL IRRMLT(ITBTYP, NTAXA, NTESTS, NDIST, LSTAB, LGROUP, GROUP,
     *            NBI, ISETS, LSETS, LTBP, IOUT, IFAULT, ITAX1, ITAX2,
     *            IT1LIM, NIBP, NB0, KPOWR2, 0, NIBP, LDBP, JTLIST,
     *            KPOWR2 - 1)
      RETURN
   40 IFAULT = 1
      RETURN
      END
      SUBROUTINE IRRDTB(ITBTYP, NTAXA, NTESTS, NDIST, LSTAB, LGROUP,
     *                  GROUP, NBI, ISETS, LSETS, IOUT, IFAULT, ITAX1,
     *                  ITAX2, IT1LIM, NIBP, NB0, KPOWR2, KBPEST, KDBP,
     *                  LDBP, JTLIST, LTLIST, LNEXTE)
C
C        ALGORITHM AS 263.2  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Forms entries of the taxon x taxon table of sets of tests
C        that distinguish each pair of taxa.  LNEXTE is .TRUE. on entry
C        if required to form tests in next non-null entry (by IRRMLT
C        when taxon x taxon table not stored)
C
      INTEGER IFAULT, IOUT, ITAX1, ITAX2, ITBTYP, IT1LIM, JTLIST,
     *        KBPEST, KDBP, KPOWR2, LDBP, LSETS, LTLIST, NBI, NB0,
     *        NDIST, NIBP, NTAXA, NTESTS
      INTEGER GROUP(NTAXA), ISETS(LSETS)
      LOGICAL LGROUP, LNEXTE, LSTAB
      INTEGER IBT, JBT, J1, J2, K1, K2, K3, NTD, NTE
      INTEGER IAND, ICNT, IOR
      LOGICAL LTDIST
      EXTERNAL IAND, ICNT, IOR, LTDIST
C
      IF (ITAX1 .GT. 0) GO TO 20
C
C        Increment ITAX1/2 (scan rows/cols of table)
C
   10 ITAX1 = ITAX1 + 1
      IF (ITAX1 .GT. IT1LIM) RETURN
      IF (LGROUP .AND. (ITBTYP .GT. 0)) THEN
         IF (GROUP(ITAX1) .NE. ITBTYP) GO TO 10
         ITAX2 = 0
      ELSE
C
C        Initialise ITAX2 (to scan across next line of table)
         ITAX2 = ITAX1
      END IF
C
C        Increment ITAX2
   20 ITAX2 = ITAX2 + 1
      IF (ITAX2 .GT. NTAXA) GO TO 10
      IF (LGROUP) THEN
         IF (GROUP(ITAX1) .EQ. GROUP(ITAX2)) GO TO 20
      ELSE
         IF (ITAX1 .EQ. ITAX2) GO TO 20
      END IF
C
C        Form list of tests that distinguish between ITAX1 and ITAX2
C
      JTLIST = LTLIST + 1
      NTE = 0
      NTD = 0
      DO 30 J1 = 1, NTESTS
C
C        Exit to end of loop if test J1 does not distinguish ITAX1 and
C        ITAX2
C
         IF ( .NOT. LTDIST(J1, ITAX1, ITAX2)) GO TO 30
         NTD = NTD + 1
C
C        If NDIST tests already distinguish ITAX1 and ITAX2, and the
C        aim is merely to form the b.p. of essential tests, take the
C        next pair of taxa
C
         IF ((NTD .GT. NDIST) .AND. ( .NOT. (LNEXTE .OR.
     *       LSTAB))) GO TO 20
C
C        Check whether J1 is an essential test
C
         JBT = (J1 + NB0 - 1) / NBI + 1
         IBT = ISETS(KPOWR2 + JBT * NBI - (J1 + NB0))
         IF (IAND(IBT, ISETS(KBPEST + JBT)) .GT. 0) THEN
C
C        Check whether sufficient essential tests already separate
C        ITAX1/2
C
            NTE = NTE + 1
            IF (NTE .GE. NDIST) GO TO 20
         END IF
C
C        Enter J1 in the list of tests that distinguish ITAX1 and ITAX2
C
         JTLIST = JTLIST - 2
         IF (JTLIST .LE. LDBP) GO TO 80
         ISETS(JTLIST) = JBT
         ISETS(JTLIST + 1) = IBT
   30 CONTINUE
      IF (NTD .LE. NDIST) THEN
C
C        Taxa ITAX1 and ITAX2 distinguished by </= NDIST tests
C
         IF (LNEXTE) GO TO 20
         IF (NTD .EQ. 0) THEN
C
C        Taxa ITAX1 and ITAX2 cannot be distinguished
C
            WRITE (IOUT,'(A,I6,A,I6,A)') ' WARNING: TAXA', ITAX1,
     *                         ' AND', ITAX2,' CANNOT BE DISTINGUISHED'
            GO TO 20
         END IF
C
C        No more than NDIST tests distinguish ITAX1 and ITAX2
C
         IF (NTD .LT. NDIST) WRITE (IOUT,'(A,I6,A,I6,A,I6,A)')
     *      ' WARNING: TAXA', ITAX1, ' AND', ITAX2,
     *      ' DISTINGUISHED BY ONLY', NTD, ' TEST(S)'
C
C        Add to set of essential tests
C
         DO 40 J1 = JTLIST, LTLIST, 2
            ISETS(KBPEST + ISETS(J1)) = IOR(ISETS(KBPEST + JBT),
     *                                      ISETS(J1 + 1))
   40    CONTINUE
C
C        Delete any sets containing >/= NDIST essential tests
C
         IF (KDBP .GE. LDBP) GO TO 20
         K2 = KDBP
         K3 = LDBP
         LDBP = KDBP
         DO 60 K1 = KDBP, K3, NIBP
            IF (K1 .LT. K3) THEN
               IF (NDIST .EQ. 1) THEN
C
C        Check whether set contains new essential test (IBT/JBT from
C        loop 30)
C
                  IF (IAND(ISETS(K1 + JBT), IBT) .EQ. 0) GO TO 60
               ELSE
C
C        Check whether set contains >/= NDIST essential tests
C
                  NTE = 0
                  DO 50 J2 = K1 + 1, K1 + NIBP
                     NTE = NTE + ICNT(IAND(ISETS(J2),
     *                     ISETS(J2 + KBPEST - K1)))
   50             CONTINUE
                  IF (NTE .LT. NDIST) GO TO 60
               END IF
            END IF
C
C        All sets checked or current set deleted, compact the list
C
            CALL CPSETS(ISETS, LSETS, LDBP, K2, K1 - K2)
            LDBP = LDBP + K1 - K2
            K2 = K1 + NIBP
   60    CONTINUE
      ELSE
C
C        >NDIST tests distinguish ITAX1 and ITAX2
C
         IF (LNEXTE) RETURN
C
C        Form bit pattern for new set at LDBP
C
         IF (LDBP + NIBP .GE. JTLIST) GO TO 80
         CALL SETVAL(ISETS(LDBP + 1), 0, NIBP)
         DO 70 J1 = JTLIST, LTLIST, 2
            ISETS(LDBP + ISETS(J1)) = ISETS(LDBP + ISETS(J1)) +
     *                                ISETS(J1 + 1)
   70    CONTINUE
C
C        Merge new set into list KDBP...LDBP
C
         J1 = LDBP
         CALL IRRMRG(KDBP, J1, LDBP, NIBP, ISETS, LSETS)
      END IF
      GO TO 20
C
C        Insufficient workspace
C
   80 IFAULT = 2
      RETURN
      END
      SUBROUTINE IRRMLT(ITBTYP, NTAXA, NTESTS, NDIST, LSTAB, LGROUP,
     *                  GROUP, NBI, ISETS, LSETS, LTBP, IOUT, IFAULT,
     *                  ITAX1, ITAX2, IT1LIM, NIBP, NB0, KPOWR2, KBPEST,
     *                  KDBP, LDBP, JTLIST, LTLIST)
C
C        ALGORITHM AS 263.3  APPL.STATIST. (1991), VOL.40, NO.1
C
C        'multiplies' entries in the taxon x taxon table of (sets of)
C        distinguishing tests, thus forming irredundant test sets
C
      INTEGER IFAULT, IOUT, ITAX1, ITAX2, ITBTYP, IT1LIM, JTLIST,
     *        KBPEST, KDBP, KPOWR2, LTBP, LDBP, LSETS, LTLIST, NBI,
     *        NIBP, NB0, NDIST, NTAXA
      INTEGER NTESTS, GROUP(NTAXA), ISETS(LSETS)
      LOGICAL LGROUP, LSTAB
      INTEGER IBT, INTD, JBT, J1, J2, J3, J4, KBPENL, KTBCAB, KTBMBE,
     *        KTBP, KTBP2, K1, K2, K3, K4, LTBCAB, LTBMBE, LTBMRG,
     *        LTBNAB, LTIAB, NTD
      INTEGER IAND, ICNT, IOR
      LOGICAL IRRWSP
      EXTERNAL IAND, ICNT, IOR, IRRWSP
C
C        Initial set contains just the essential tests
C
      KTBP = LDBP
      CALL CPSETS(ISETS, LSETS, KTBP, KBPEST, NIBP)
      LTBP = KTBP + NIBP
C
C        Get the next entry from the taxon x taxon table
C
   10 KTBP2 = KTBP
      IF (LSTAB) THEN
C
C        Generate test list from stored bit pattern
C
         IF (KTBP .LE. KDBP) RETURN
         JTLIST = LTLIST + 1
         KTBP = KTBP - NIBP
         K2 = KPOWR2 + NBI - NB0 - 1
         DO 40 JBT = 1, NIBP
            K1 = ISETS(KTBP + JBT)
            DO 20 J1 = K2, KPOWR2, -1
               IBT = ISETS(J1)
               IF (IAND(K1, IBT) .EQ. 0) GO TO 20
               JTLIST = JTLIST - 2
               IF (JTLIST .LE. LTBP + NIBP) GO TO 180
               ISETS(JTLIST) = JBT
               ISETS(JTLIST + 1) = IBT
               K1 = K1 - IBT
               IF (K1 .LE. 0) GO TO 30
   20       CONTINUE
   30       K2 = KPOWR2 + NBI - 1
   40    CONTINUE
      ELSE
         CALL IRRDTB(ITBTYP, NTAXA, NTESTS, NDIST, LSTAB, LGROUP, GROUP,
     *               NBI, ISETS, LSETS, IOUT, IFAULT, ITAX1, ITAX2,
     *               IT1LIM, NIBP, NB0, KPOWR2, KBPEST, KDBP, LDBP,
     *               JTLIST, LTLIST, .TRUE.)
         IF ((ITAX1 .GT. IT1LIM) .OR. (IFAULT .NE. 0)) RETURN
      END IF
C
C        Sort sets into (a) pass-cannot-absorb (in KTBP...LTBNAB),
C        (b) pass-can-absorb (in KTBCAB (initially LTBP)...LTBCAB),
C        (c/d) must-be-enlarged (in KTBMBE...LTBMBE); also, if NDIST=1,
C        sort tests into contained-in-can-absorb (JTLIST...LTIAB-1)
C        and not-contained-in-can-absorb (LTIAB...LTLIST)
C
      LTBNAB = KTBP
      KTBCAB = LTBP
      LTBCAB = LTBP
      LTBMBE = JTLIST - 1
C
C        Leave space above LTBMBE for sets with none of the new tests(d)
C
      IF (NDIST .GT. 1) LTBMBE = LTBCAB + (JTLIST - LTBCAB) / 2
      KTBMBE = LTBMBE
      LTIAB = JTLIST
      DO 60 K1 = KTBP2, LTBP - 1, NIBP
C
C        Count number of tests from list already in each set
C
         NTD = 0
         DO 50 J1 = JTLIST, LTLIST, 2
            IF (IAND(ISETS(K1 + ISETS(J1)), ISETS(J1 + 1)) .EQ.
     *          0) GO TO 50
            NTD = NTD + 1
            K3 = J1
   50    CONTINUE
         IF (NDIST .LT. NTD) THEN
C
C        Pass, cannot absorb
C
            CALL CPSETS(ISETS, LSETS, LTBNAB, K1, NIBP)
            LTBNAB = LTBNAB + NIBP
         ELSE
            IF (NDIST .EQ. NTD) THEN
C
C        Pass, can absorb
C
               LTBCAB = LTBCAB + NIBP
               IF (LTBCAB .GT. KTBMBE) THEN
                  IF (IRRWSP(LTBCAB, KTBMBE, LTBMBE, JTLIST, ISETS,
     *                LSETS)) GO TO 180
               END IF
               CALL CPSETS(ISETS, LSETS, LTBCAB - NIBP, K1, NIBP)
               IF ((NDIST .EQ. 1) .AND. (K3 .GE. LTIAB)) THEN
                  JBT = ISETS(K3)
                  IBT = ISETS(K3 + 1)
                  ISETS(K3) = ISETS(LTIAB)
                  ISETS(K3 + 1) = ISETS(LTIAB + 1)
                  ISETS(LTIAB) = JBT
                  ISETS(LTIAB + 1) = IBT
                  LTIAB = LTIAB + 2
               END IF
            ELSE
C
C        Must be enlarged
C
               IF ((NDIST .EQ. 1) .OR. (NTD .GT. 0)) THEN
                  KTBMBE = KTBMBE - NIBP
                  IF (LTBCAB .GE. KTBMBE) THEN
                     IF (IRRWSP(LTBCAB, KTBMBE, LTBMBE, JTLIST, ISETS,
     *                   LSETS)) GO TO 180
                  END IF
                  CALL CPSETS(ISETS, LSETS, KTBMBE, K1, NIBP)
                  IF (NDIST .EQ. 1) GO TO 60
C
C        (NDIST>1) store number of new tests already in set
C
                  ISETS(KTBMBE) = NTD
                  KTBMBE = KTBMBE - 1
               ELSE
                  LTBMBE = LTBMBE + NIBP + 1
                  IF (LTBMBE .GE. JTLIST) THEN
                     IF (IRRWSP(LTBCAB, KTBMBE, LTBMBE, JTLIST, ISETS,
     *                   LSETS)) GO TO 180
                  END IF
                  CALL CPSETS(ISETS, LSETS, LTBMBE - NIBP, K1, NIBP)
                  ISETS(LTBMBE - NIBP) = 0
               END IF
            END IF
         END IF
   60 CONTINUE
C
C        Copy down pass-can-absorb sets
C
      J3 = LTBCAB - KTBCAB
      CALL CPSETS(ISETS, LSETS, LTBNAB, KTBCAB, J3)
      KTBCAB = LTBNAB
      LTBCAB = KTBCAB + J3
      LTBP = LTBCAB
      IF (LTBMBE .LE. KTBMBE) GO TO 10
C
C        Form and check new b.p.'s from must-be-enlarged list
C
      IF (NDIST .EQ. 1) THEN
         NTD = 0
         K2 = JTLIST
         K3 = LTLIST
         K4 = 2
      ELSE
C
C        Calc number of 'augmenting' patterns (no. combinations of
C        NDIST tests)
C
         J1 = (LTLIST - JTLIST + 1) / 2
         J3 = 1
         J4 = 1
         DO 70 J2 = 1, NDIST
            J3 = J3 * J2
            J4 = J4 * J1
            J1 = J1 - 1
   70    CONTINUE
C
C        Reserve space for 'augmenting' patterns
C
         KBPENL = JTLIST - 1 - J4 * NIBP / J3
C
C        Copy must-be-enlarged patterns to below 'augmenting' patterns
C
         J1 = LTBMBE - KTBMBE
         KTBMBE = KBPENL - J1
         CALL CPSETS(ISETS, LSETS, KTBMBE, LTBMBE - J1, J1)
         LTBMBE = KTBMBE + J1
C
C        Form 'augmenting' patterns (all combinations of NDIST tests)
C
         K1 = KBPENL
         K4 = LTBP + NDIST
         IF (K4 .GT. KTBMBE) GO TO 180
         K3 = LTBP + 1
         K2 = LTBP
         J1 = JTLIST - 2
   80    K2 = K2 + 1
   90    J1 = J1 + 2
         ISETS(K2) = J1
         IF (K2 .LT. K4) GO TO 80
         CALL SETVAL(ISETS(K1 + 1), 0, NIBP)
         DO 100 J1 = K3, K4
            J2 = ISETS(J1)
            ISETS(K1 + ISETS(J2)) = ISETS(K1 + ISETS(J2)) +
     *                              ISETS(J2 + 1)
  100    CONTINUE
         K1 = K1 + NIBP
  110    J1 = ISETS(K2)
         IF (LTLIST - J1 .GT. (K4 - K2 + 1) * 2) GO TO 90
         K2 = K2 - 1
         IF (K2 .GE. K3) GO TO 110
         K2 = KBPENL
         K3 = K1 - 1
         K4 = NIBP
         INTD = 1
      END IF
C
C        Enlarge b.p.'s (and check for absorption)
C
  120 IF (NDIST .GT. 1) THEN
         KTBMBE = KTBMBE + 1
         NTD = ISETS(KTBMBE)
         LTBMRG = LTBP
         IF ((NTD .EQ. 0) .AND. (INTD .EQ. 1)) THEN
C
C        Sets initially with 0 new tests may be absorbed by sets that
C        had >0
C
            LTBCAB = LTBP
            INTD = 0
         END IF
      END IF
      DO 170 K1 = K2, K3, K4
C
C        Check space
C
         IF (LTBP + NIBP .GT. KTBMBE) GO TO 180
         CALL CPSETS(ISETS, LSETS, LTBP, KTBMBE, NIBP)
         IF (NDIST .EQ. 1) THEN
C
C        Add 1 test
C
            ISETS(LTBP + ISETS(K1)) = IOR(ISETS(LTBP + ISETS(K1)),
     *                                ISETS(K1 + 1))
            IF (K1 .GE. LTIAB) GO TO 160
         ELSE
C
C        Add >1 tests
C
            J2 = 0
            DO 130 J1 = 1, NIBP
               J2 = J2 + ICNT(IAND(ISETS(LTBP + J1), ISETS(K1 + J1)))
               ISETS(LTBP + J1) = IOR(ISETS(LTBP + J1), ISETS(K1 + J1))
  130       CONTINUE
C
C        Sets augmented with >NDIST-NTD tests will be absorbed
C
            IF (J2 .LT. NTD) GO TO 170
         END IF
C
C        Check whether absorbed by pass-can-absorb b.p. (b)
C
         DO 150 J2 = KTBCAB, LTBCAB - 1, NIBP
            DO 140 J3 = 1, NIBP
               IF (IAND(ISETS(LTBP + J3), ISETS(J2 + J3)) .NE.
     *             ISETS(J2 + J3)) GO TO 150
  140       CONTINUE
            GO TO 170
  150    CONTINUE
  160    IF (NTD .EQ. 0) THEN
            LTBP = LTBP + NIBP
         ELSE
C
C        Check whether contains/contained in b.p. formed from other
C        sets in (c)
C
            CALL IRRMRG(LTBCAB, LTBMRG, LTBP, NIBP, ISETS, LSETS)
         END IF
  170 CONTINUE
C
C        Next set to be enlarged
C
      KTBMBE = KTBMBE + NIBP
      IF (KTBMBE - LTBMBE) 120, 10, 10
C
C        Insufficient workspace
C
  180 IFAULT = 2
      RETURN
      END
      SUBROUTINE IRRMRG(KBP, LBP, KNBP, NIBP, ISETS, LSETS)
C
C        ALGORITHM AS 263.4  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Checks new b.p. ISETS(KNBP+1...) vs list in ISETS(KBP+1...LBP):
C        deletes it if iT contains any b.p. in the list; otherwise,
C        deletes b.p.'s that contain it, adjusts LBP & sets KNBP to
C        length of the full list
C
      INTEGER KBP, KNBP, LBP, LSETS, NIBP
      INTEGER ISETS(LSETS)
      INTEGER IBPINT, IBPL, IBPN, INCL, J1, K1, K2, K3
      INTEGER IAND
      EXTERNAL IAND
C
      K2 = KBP
      K3 = KBP
      DO 20 K1 = KBP, LBP - 1, NIBP
C
C        Check whether new set contains/is contained in set at K1
C        (INCL=-/+1)
C
         INCL = 0
         DO 10 J1 = 1, NIBP
            IBPL = ISETS(K1 + J1)
            IBPN = ISETS(KNBP + J1)
            IF (IBPL .EQ. IBPN) GO TO 10
            IBPINT = IAND(IBPL, IBPN)
            IF (IBPINT .EQ. IBPL) THEN
               IF (INCL .GT. 0) GO TO 20
               INCL = -1
            ELSE
               IF ((IBPINT .NE. IBPN) .OR. (INCL .LT. 0)) GO TO 20
               INCL = 1
            END IF
   10    CONTINUE
C
C        Return if new set contains/is same as an already stored set
C
         IF (INCL .LE. 0) RETURN
         CALL CPSETS(ISETS, LSETS, K2, K3, K1 - K3)
C
C        K2 is length of compacted list
         K2 = K2 + K1 - K3
C
C        K3=origin of section to be copied down when next containing
C        set found
         K3 = K1 + NIBP
   20 CONTINUE
      IF (K2 .GT. KBP) THEN
C
C        Copy final section, intermediate b.p.'s (LBP...KNBP) and
C        new b.p.
C
         CALL CPSETS(ISETS, LSETS, K2, K3, KNBP + NIBP - K3)
         LBP = LBP - K3 + K2
         KNBP = K2 + KNBP + NIBP - K3
      ELSE
C
C        No contained/containing sets found
C
         KNBP = KNBP + NIBP
      END IF
      RETURN
      END
      LOGICAL FUNCTION IRRWSP(KWSP, KUSED, LUSED, LWSP, ISETS, LSETS)
C
C        ALGORITHM AS 263.5  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Reorganises workspace to give room before and after
C        KUSED...LUSED; returns value .TRUE. if insufficient space
C
      INTEGER KUSED, KWSP, LSETS, LUSED, LWSP
      INTEGER ISETS(LSETS)
      INTEGER NUSED
C
      NUSED = LUSED - KUSED
      KUSED = KWSP + (KUSED - KWSP + LWSP - LUSED) / 2
      CALL CPSETS(ISETS, LSETS, KUSED, LUSED - NUSED, NUSED)
      LUSED = KUSED + NUSED
      IRRWSP = (KUSED .LE. KWSP) .OR. (LUSED .GE. LWSP)
      RETURN
      END
      SUBROUTINE CPSETS(IVAL, LVAL, KTO, KFROM, NVAL)
C
C        ALGORITHM AS 263.6  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Copies IVAL(KFROM+1...+NVAL) to IVAL(KTO+1...)
C
      INTEGER KFROM, KTO, LVAL, NVAL
      INTEGER IVAL(LVAL)
      INTEGER FROMTO, INC, ITO, JTO, LTO
C
C        Check whether new location above, same as, or below old
C
      FROMTO = KFROM - KTO
      IF (FROMTO) 10, 50, 20
C
C        Values may overlap, copy in reverse order
C
   10 INC = -1
      JTO = KTO + NVAL
      LTO = KTO + 1
      GO TO 30
C
C        Set/copy values from origin upwards
C
   20 INC = 1
      JTO = KTO + 1
      LTO = KTO + NVAL
C
C        Copy values
C
   30 DO 40 ITO = JTO, LTO, INC
         IVAL(ITO) = IVAL(ITO + FROMTO)
   40 CONTINUE
   50 RETURN
      END
      SUBROUTINE SETVAL(IVAL, ICON, NVAL)
C
C        ALGORITHM AS 263.7  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Sets IVAL(1...NVAL) to ICON
C
      INTEGER ICON, NVAL
      INTEGER IVAL(NVAL)
      INTEGER I
      DO 10 I = 1, NVAL
         IVAL(I) = ICON
   10 CONTINUE
      RETURN
      END
