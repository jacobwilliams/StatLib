      SUBROUTINE PRINBP(NBI, NFACT, NBP, NPRINT, IBPS, LIBPS, IOUT,
     *                  IFAULT)
C
C        ALGORITHM AS 264  APPL.STATIST. (1991), VOL.40, NO.1
C
C        Prints bit patterns on unit IOUT
C
      INTEGER IOUT, IFAULT, LIBPS, NBI, NPRINT, NBP, NFACT
      INTEGER IBPS(LIBPS)
      INTEGER JP2IT1, JTLIST, J1, J2, J3, J4, J5, KPOWR2, K1, K2, K3,
     *        K4, K5, LTLIST, N, NB0, NIBP, NMIN, NMAX, NSETP
      INTEGER IAND, ICNT
      EXTERNAL IAND, ICNT
C
      IF ((NFACT .LE. 0) .OR. (NBP .LE. 0)) GO TO 160
C
C        Calculate number of integers to store each bit pattern
C
      NIBP = (NFACT - 1) / NBI + 1
C
C        Set up powers of 2 (2**0,1,2...) in IBPS(KPOWR2...LIBPS)
C
      KPOWR2 = LIBPS - NBI + 1
      IBPS(KPOWR2) = 1
      DO 10 J1 = KPOWR2 + 1, LIBPS
         IBPS(J1) = IBPS(J1 - 1) * 2
   10 CONTINUE
C
C        Address of the power of 2 used to represent the first item
C
      JP2IT1 = KPOWR2 + NFACT - NBI * (NIBP - 1) - 1
      LTLIST = KPOWR2 - 1
      IF (NBP + NBP / NIBP .GE. LTLIST) GO TO 150
      NSETP = 0
      NMIN = NFACT
      NMAX = 1
C
C        Store number of factors in each set in IBPS(NBP+1...K1)
C
      K1 = NBP
      DO 30 K3 = 1, NBP, NIBP
         J2 = 0
         DO 20 J1 = K3, K3 + NIBP - 1
            J2 = J2 + ICNT(IBPS(J1))
   20    CONTINUE
         NMIN = MIN(NMIN, J2)
         NMAX = MAX(NMAX, J2)
         K1 = K1 + 1
         IBPS(K1) = J2
   30 CONTINUE
      IF (NPRINT .LT. 0) NMAX = MIN(NMAX, NMIN - NPRINT - 1)
      IF (NPRINT .GT. 0) NMAX = MIN(NMAX, NPRINT)
C
      DO 140 N = NMIN, NMAX
         JTLIST = LTLIST - N + 1
C
C        Store addresses of sets with N factors in IBPS(K1+1...K2)
C
         K2 = K1
         J1 = NBP
         DO 90 K3 = 1, NBP, NIBP
            J1 = J1 + 1
            IF (IBPS(J1) .NE. N) GO TO 90
C
C        Place address of set into list, according to
C        lexicographic order
C
            DO 70 J2 = K1 + 1, K2
               J3 = IBPS(J2)
               DO 60 J4 = K3, K3 + NIBP - 1
                  J3 = J3 + 1
                  IF (IBPS(J3) - IBPS(J4)) 40, 60, 70
C
C        Move other sets up and insert current set
C
   40             DO 50 J5 = K2, J2, -1
                     IBPS(J5 + 1) = IBPS(J5)
   50             CONTINUE
                  IBPS(J2) = K3 - 1
                  GO TO 80
   60          CONTINUE
   70       CONTINUE
C
C        Insert set at end of list
C
            IBPS(K2 + 1) = K3 - 1
   80       K2 = K2 + 1
            IF (K2 .GE. JTLIST) GO TO 150
   90    CONTINUE
C
C        Decode each bit pattern, store factors in
C        IBPS(JTLIST...LTLIST), then print
C
         DO 130 J2 = K1 + 1, K2
            K3 = 0
            K4 = JTLIST
            K5 = JP2IT1
            DO 110 J3 = IBPS(J2) + 1, IBPS(J2) + NIBP
               J4 = IBPS(J3)
               DO 100 J5 = K5, KPOWR2, -1
                  K3 = K3 + 1
                  IF (IAND(IBPS(J5), J4) .EQ. 0) GO TO 100
                  IBPS(K4) = K3
                  IF (K4 .EQ. LTLIST) GO TO 120
                  K4 = K4 + 1
  100          CONTINUE
               K5 = LIBPS
  110       CONTINUE
C
C        Print contents of bit pattern (now stored
C        in IBPS(JTLIST...LTLIST))
C
  120       NSETP = NSETP + 1
C
C ***** Lines below can be modified to change the style of output *****
C
            WRITE (IOUT, 9000) NSETP, N, (IBPS(J3), J3=JTLIST, LTLIST)
 9000       FORMAT(1X, I6, ')', I5, '  FACTOR(S):  ',20I4 / (27X, 20I4))
  130    CONTINUE
  140 CONTINUE
      RETURN
C
  150 IFAULT = IFAULT + 1
  160 IFAULT = IFAULT + 1
      RETURN
      END
