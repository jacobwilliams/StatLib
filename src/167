      SUBROUTINE FIT(A, NA, NTERMS, MSIZE, IE, INC, NINC, NADD,
     *  GVAR, KEY, EFF, EPS, IFAULT)
C
C        ALGORITHM AS 167  APPL. STATIST. (1981)  VOL.30, NO.3
C
C        CALCULATES EFFICIENCIES OF ESTIMATIOM OF TERMS IN
C        THE CURRENT MODEL AND ITS GENERALIZED VARIANCE
C
C     *** Warning: This text has been read using a scanner and may
C                  contain errors
C
      REAL A(NA), EFF(NTERMS)
      INTEGER INC(NINC), IE(NTERMS), KEY(NTERMS)
      DATA ONE /1.0/
      IFAULT = 0
      IA = NTERMS * (NTERMS + 1) / 2
      IF (NA .LT. IA) IFAULT = 2
      IF (MSIZE .LT. 0 .OR. NADD .LE. 0 .OR. NADD .GT. NINC)
     +    IFAULT = 3
      IF (MSIZE .NE. 0 .OR. IFAULT .NE. 0) GOTO 5
      ONEPLS = ONE + EPS
      ONEMIN = ONE - EPS
      IA = 0
      DO 3 J1 = 1, NTERMS
      DO 1 J2 = 1, J1
      IA = IA + 1
      R = A(IA)
      IF (R .LT. -ONEPLS .OR. R .GT. ONEPLS) IFAULT = 1
    1 CONTINUE
      IF (R .LT. ONEMIN) IFAULT = 1
      IE(J1) = -1
    3 CONTINUE
      GVAR = ONE
    5 IF (IFAULT .NE. 0) RETURN
      DO 20 J = 1, NADD
      JP = INC(J)
      IF (JP .LE. 0 .OR. JP .GT. NTERMS) GOTO 50
      IF (IABS(IE(JP)) .NE. 1) GOTO 50
      IA = JP * (JP + 1) / 2
      IF (ABS(A(IA)) .LE. EPS) GOTO 19
      IF (IE(JP) .LT. 0) GOTO 16
C
C        DROP TERM FROM MODEL
C
      K = MSIZE - 1
      DO 12 J1 = 1, K
      IF (JP .NE. KEY(J1)) GOTO 12
      DO 8 J2 = J1, K
    8 KEY(J2) = KEY(J2 + 1)
      GOTO 13
   12 CONTINUE
   13 MSIZE = K
      GVAR = -GVAR / A(IA)
      CALL PIVOT(A, NA, NTERMS, JP, IE)
      GOTO 20
C
C        ADD TERM TO MODEL
C
   16 MSIZE = MSIZE + 1
      KEY(MSIZE) = JP
      CALL PIVOT(A, NA, NTERMS, JP, IE)
      GVAR = -GVAR * A(IA)
      GOTO 20
C
C        LEAVE TERM IN OR OUT OF MODEL AS PIVOT TOO SMALL
C
   19 IFAULT = -JP
   20 CONTINUE
C
C        CALCULATE EFFICIENCIES OF EACH TERM IN MODEL
C
      DO 40 J = 1, MSIZE
      IA = IA * (IA + 1) / 2
      EFF(J) = -ONE / A(IA)
   40 CONTINUE
      RETURN
   50 IFAULT = 4
      RETURN
      END
C
      SUBROUTINE PIVOT(A, NA, K, IP, IE)
C
C        ALGORITHM AS 167.1  APPL. STATIST. (1981) VOL.30, NO.3
C
C        ROW BY ROW INVERSION IN SITU OF A SYMMETRIC MATRIX STORED AS
C        LOWCR TRIANGLE. A FORTRAN TRANSLATION OF ALGORITHM AS 37
C        (GARSIDE, 1971)
      REAL A(NA)
      INTEGER IE(K)
C
      DATA ONE /1.0/
      JP = IABS(IP)
      JPP = JP * (JP + 1) / 2
      JOP = JPP - JP
      AA = ONE / A(JPP)
      A(JPP) = -AA
      NIJ = 0
      IF (JP .EQ. 1) GOTO 100
      JP1 = JP - 1
      DO 80 I = I, JP1
      NIP = JOP + I
      AIP = -A(NIP) * FLOAT(IE(JP))
      A(NIP) = AIP * AA
      DO 60 J = 1, I
      NIJ = NIJ + 1
      NJP = JOP + J
      A(NIJ) = A(NIJ) - AIP * A(NJP)
   60 CONTINUE
   80 CONTINUE
  100 IF (JP .EQ. K) GOTO 190
      NIJ = NIJ + JP
      JP2 = JP + 1
      NIP = JPP
      DO 150 I = JP2, K
      NIP = NIP + I - 1
      AIP = -A(NIP) * FLOAT(IE(JP))
      A(NIP) = AIP * AA
      IF (JP .EQ. 1) GOTO 130
      DO 110 J = 1, JP1
      NIJ = NIJ + 1
      NJP = JOP + J
      A(NIJ) = A(NIJ) - AIP * A(NJP)
  110 CONTINUE
  130 NIJ = NIJ + 1
      NJP = JPP + JP
      DO 140 J = JP2, I
      NIJ = NIJ + 1
      A(NIJ) = A(NIJ) - AIP * A(NJP)
      NJP = NJP + J
  140 CONTINUE
  150 CONTINUE
  190 IE(JP) = -IE(JP)
      RETURN
      END

