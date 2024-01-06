      SUBROUTINE BXPLT(NG, RVALS, IVALS, SCALEX, SCALEY, NY, NX, XTYPE,
     *   WIDTH, NOTCH, IWRITE, OUTARR, INTARR, IVALS2, IFAULT)
C
C     ALGORITHM AS 272.1 APPL.STATIST. (1992), VOL.41, NO.1
C
C     PRODUCES A BOX-PLOT USING SUPPLIED QUARTILES
C
      REAL RVALS(6, NG), SCALEX(2), SCALEY(2), OUTARR(*), VALS(20)
      INTEGER IVALS(6, NG), WIDTH, WID, NX, IVALS2(8, NG),
     *   INTARR(*), POINTS, PNTR
      CHARACTER INTCH(0:9), COLON, DASH, DOT
      CHARACTER*19 IFORM1
      CHARACTER*20 IFORM2
      CHARACTER*120 IOUT
      LOGICAL LWIDOK, NOTCH, XTYPE
      EXTERNAL AXIS, SCALE
C
      DATA MPVX / 10 /, MPVY / 5 /, MAXHT / 66 /, MAXWID / 130 /
      DATA IFORM1 /'(1H ,F8.0,1X,120A1)' /
      DATA IFORM2 /'(1H ,5X,16(F8.0,2X))' /
      DATA COLON / ':' / , DASH / '-' / ,  DOT / '.' /
      DATA INTCH / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
C
C CHECK PARAMETERS
C
      IFAULT = 0
      IF (NG .LT. 1) GOTO 900
      XMIN = SCALEX(1)
      XMAX = SCALEX(2)
      YMIN = SCALEY(1)
      YMAX = SCALEY(2)
C
C CHECK IF NLX NLY ARE IN LEGAL RANGE
C
      NLX = NX
      NLY = NY
      IF (NLY + 5 .GT. MAXHT) NLY = MAXHT - 5
      IF (NLY .LT. MPVY) NLY = MPVY
      IF (NLX + 11 .GT. MAXWID) NLX = MAXWID - 11
      IF (NLX .LT. MPVX) NLX = MPVX
      IF (XTYPE) NLX = (NG + 1) * (NLX / (NG + 1))
      NOUT = 0
      DO  10 I = 1, NG
         NOUT = NOUT + IVALS(1, I) + IVALS(2, I)
   10 CONTINUE
      WID = WIDTH
      IF (WIDTH .EQ. 0) WID = NLX / (2 * NG + 1) + 1
C
C CALCULATE WIDTH EQUAL MARK SPACE RATIO
C
      IF (XTYPE) THEN
C
C QUALITATIVE X SIMULATE CALL  OF SCALE
C
         XMIN = 0.0
         XSTEP = 1.0 / (NLX / (NG + 1))
C
C MUST CALCULATE SCALE, FIRST CHOOSE XMIN AND XMAX
C
      ELSE IF (XMIN .GE. XMAX) THEN
         XMIN = RVALS(1, 1)
         XMAX = XMIN
         DO 20 I = 2, NG
         IF (RVALS(1, I) .LT. XMIN) THEN
            XMIN = RVALS(1, I)
         ELSE IF (RVALS(1, I) .GT. XMAX) THEN
            XMAX = RVALS(1, I)
         ENDIF
   20    CONTINUE
         XMN = XMIN
         XMX = XMAX
   30    CALL SCALE(XMN, XMX, NLX, MPVX, TEMP, XVSTEP, NXVALS, IRX,
     *        IFAIL)
         IF (IFAIL .GT. 0) GOTO 901
         XSTEP = XVSTEP / MPVX
         IF ((XMIN - TEMP) / XSTEP .LT. WID) THEN
C
C OFFPLOT ON LEFT
C
            XMN = XMN - XSTEP
            GOTO  30
         ENDIF
         IF ((XMAX - TEMP) / XSTEP .GT. (NLX - WID)) THEN
C
C OFFPLOT ON RIGHT
C
            XMX = XMX + XSTEP
            GOTO  30
         ENDIF
         XMIN = TEMP
      ELSE
         CALL SCALE(XMIN, XMAX, NLX, MPVX, TEMP, XVSTEP, NXVALS, IRX,
     *        IFAIL)
         IF (IFAIL .GT. 0) GOTO 901
         XMIN = TEMP
         XSTEP = XVSTEP / MPVX
      ENDIF
      IF (YMIN .GE. YMAX) THEN
C
C CHOOSE YMIN AND YMAX
C
         YMAX = RVALS(2, 1)
         YMIN = RVALS(3, 1)
         DO 40 I = 2, NG
         IF (RVALS(2, I) .GT. YMAX) THEN
            YMAX = RVALS(2, I)
         ELSE IF (RVALS(3, I) .LT. YMIN) THEN
            YMIN = RVALS(3, I)
         ENDIF
   40    CONTINUE
C
C CHECK IN OUTLIER ARRAY, ADJUST YMIN AND YMAX
C
         DO 50 I = 1, NOUT
         IF (OUTARR(I) .GT. YMAX) THEN
            YMAX = OUTARR(I)
         ELSE IF (OUTARR(I) .LT. YMIN) THEN
            YMIN = OUTARR(I)
         ENDIF
   50    CONTINUE
      ENDIF
      CALL SCALE(YMIN, YMAX, NLY, MPVY, TEMP, YVSTEP, NYVALS, IRY,
     * IFAIL)
      IF (IFAIL .GT. 0) GOTO 902
      YMIN = TEMP
      YSTEP = YVSTEP / MPVY
C
C CONVERT TO PLOT UNITS
C
      DO 70 I = 1, NG
      DO 60 J = 2, 6
      IVALS2(J, I) = (RVALS(J, I) - YMIN) / YSTEP + 0.5
   60 CONTINUE
C
C NOTCHES AT MEDIAN +- 1.57 * IQR/SQRT(N)
C
      IF(NOTCH .AND. IVALS(5, I) .GT. 0) THEN
         QR = (RVALS(4, I) - RVALS(5, I)) * 1.57/SQRT(REAL(IVALS(5, I)))
         YDIFF = RVALS(6, I) - YMIN
         IVALS2(7, I) = (YDIFF + QR)/YSTEP + 0.5
         IVALS2(8, I) = (YDIFF - QR)/YSTEP + 0.5
      ENDIF
   70 CONTINUE
      DO 80 I = 1, NG
         IF(XTYPE) THEN
            IVALS2(1, I) = (I - XMIN) / XSTEP + 1.5
         ELSE
            IVALS2(1, I) = (RVALS(1, I) - XMIN) / XSTEP + 1.5
         ENDIF
   80 CONTINUE
      DO 90 I = 1, NOUT
      INTARR(I) = (OUTARR(I) - YMIN) / YSTEP + 0.5
   90 CONTINUE
C
C CHECK OVERPRINTING ONLY IF WIDTH CALCULATED
C
      IF (WIDTH .EQ. 0) THEN
  100    IF (WID .GT. 1) THEN
            LWIDOK = (IVALS2(1, 1) .GE. WID / 2 + 1 .AND.
     *                IVALS2(1, NG) .LE. NLX - (WID / 2 + 1))
            IF (.NOT. LWIDOK) GOTO 903
            DO 110 I = 2, NG
C
C CHECK GAPS BETWEEN BOXES
C
            IF (IVALS2(1, I) - IVALS2(1, I - 1) .LT. WID)
     *          LWIDOK = .FALSE.
  110       CONTINUE
            IF (LWIDOK) GOTO 120
            WID = WID - 2
            GOTO 100
         ENDIF
  120    CONTINUE
C
C LIMIT MUST PROVIDE ROOM FOR NOTCHES
C
         LIMIT = 3
         IF (NOTCH) LIMIT = 5
         IF (WID .LT. LIMIT) GOTO 903
      ENDIF
      DO 130 I = 1, NG
      IF (WIDTH .GE. 0) IVALS(6, I) = WID
C
C WIDTH MUST BE ODD EVEN IF EVEN VALUE SUPPLIED
C
      IF (MOD(IVALS(6, I), 2) .EQ. 0) IVALS(6, I) = IVALS(6, I) + 1
      IF (IVALS(6, I) .LT. 1) IVALS(6, I) = 1
  130 CONTINUE
      CALL AXIS(YMIN, YVSTEP, NYVALS, 6, IRY, IRPR, OFFSET, IFACT,
     *       VALS, 20, IFAIL)
      IF (IFAIL .GT. 0) GOTO 904
      IFORM1(9:9) = INTCH(IRPR)
      IF (IFACT .NE. 0) WRITE(IWRITE, 2) IFACT
      IF (OFFSET .NE. 0.0) WRITE(IWRITE, 3) OFFSET
      DO 170 I = 1, NLY
C
C MAIN LOOP TO BUILD PLOT STARTS HERE
C
      IY = NLY - I
      IOUT(1:NLX) = ' '
C
C FIND CORRECT CODE
C
      DO 160 J = 1, NG
      IF (IVALS(5, J) .GE. 1) THEN
         IX = IVALS2(1, J)
         ICODE = 0
         IF (IY .LE. IVALS2(2, J) .AND. IY .GE. IVALS2(3, J)) ICODE = 1
         POINTS = 0
         IF (IVALS(1, J) .GT. 0) THEN
            PNTR = IVALS(3, J)
            DO 140 K = 1, IVALS(1, J)
            IF (INTARR(PNTR) .EQ. IY) POINTS = POINTS + 1
            PNTR = PNTR + 1
  140       CONTINUE
         ENDIF
         IF (IVALS(2, J) .GT. 0) THEN
            PNTR = IVALS(4, J)
            DO 150 K = 1, IVALS(2, J)
            IF (INTARR(PNTR) .EQ. IY) POINTS = POINTS + 1
            PNTR = PNTR + 1
  150       CONTINUE
         ENDIF
         IF (POINTS .EQ. 1) ICODE = 2
         IF (POINTS .GT. 1) ICODE = 3
         IF (IY .EQ. IVALS2(4, J) .OR. IY .EQ. IVALS2(5, J)) ICODE = 4
         IF (IY .LT. IVALS2(4, J) .AND. IY .GT. IVALS2(5, J)) THEN
            ICODE = 5
            IF (NOTCH) THEN
               IF ((IY .LE. IVALS2(7, J) .AND. IY .GE. IVALS2(8, J)))
     *              ICODE = 6
            ENDIF
         ENDIF
         IF (IY .EQ. IVALS2(6, J)) ICODE = 7
         CALL SETGR(IOUT, IVALS(6, J), ICODE, IX, NOTCH, NLX, IFAIL)
         IF (IFAIL .EQ. 1) IFAULT =  - 1
      ENDIF
C
C UNABLE TO FIT ON PLOT
C
  160 CONTINUE
      IF (MOD(IY, MPVY) .EQ. 0) THEN
         WRITE(IWRITE, IFORM1) VALS(NYVALS), DASH, COLON,
     *                        (IOUT(IX:IX), IX = 1, NLX)
         NYVALS = NYVALS - 1
      ELSE
         WRITE(IWRITE, 1)(IOUT(IX:IX), IX = 1, NLX)
      ENDIF
  170 CONTINUE
      WRITE(IWRITE, 1)(DOT, I = 1, NLX)
      IF ( .NOT. XTYPE) THEN
         CALL AXIS(XMIN, XVSTEP, NXVALS, 6, IRX, IRPR, OFFSET, IFACT,
     *             VALS, 20, IFAIL)
         IF (IFAIL .GT. 0) GOTO 905
         IFORM2(15:15) = INTCH(IRPR)
         WRITE(IWRITE, 6)(COLON, I = 1, NXVALS)
         WRITE(IWRITE, IFORM2)(VALS(I), I = 1, NXVALS)
         IF (IFACT .NE. 0) WRITE(IWRITE, 4) IFACT
         IF (OFFSET .NE. 0.0) WRITE(IWRITE, 5) OFFSET
      ENDIF
C
C CHECK POINTS AND PARTS OF BOXES LYING BEYOND Y LIMITS
C
      IFAIL = 0
      DO 180 I = 1, NG
      IF (IVALS2(2, I) .GT. NLY - 1 .OR. IVALS2(3, I) .LT. 0) GOTO 190
  180 CONTINUE
      GOTO 200
  190 IFAULT = IFAULT - 2
  200 DO 210 I = 1, NOUT
      IF (INTARR(I) .GT. NLY - 1 .OR. INTARR(I) .LT. 0) GOTO 220
  210 CONTINUE
      RETURN
  220 IFAULT = IFAULT - 4
      RETURN
1     FORMAT(11X, ':', 120A1)
2     FORMAT(' TIMES 10**', I3)
3     FORMAT(' OFFSET', F10.0)
4     FORMAT(15X, 'TIMES 10**', I3)
5     FORMAT(15X, 'OFFSET', F10.0)
6     FORMAT(3X, 16(9X, A1))
900   IFAULT = 1
      RETURN
901   IFAULT = 2
      RETURN
902   IFAULT = 3
      RETURN
903   IFAULT = 4
      RETURN
904   IFAULT = 5
      RETURN
905   IFAULT = 6
      RETURN
      END
C
C
      SUBROUTINE SETGR(IOUT, WID, ICODE, IX, NOTCH, NLX, IFAIL)
C
C     ALGORITHM AS 272.2 APPL.STATIST. (1992), VOL.41, NO.1
C
C PLACES CHARACTERS INTO IOUT AT POSITION IX
C CURRENT PART OF BOX CODED BY ICODE
C
      INTEGER WID, ICODE, IX, IFAIL
      CHARACTER VERT, STAR, PLUS, HORIZ, DASH
      CHARACTER*120 IOUT
      LOGICAL NOTCH
      DATA VERT / ':' / , STAR / 'X' / , PLUS / '+' / , HORIZ / '-' / ,
     * DASH / '=' /
      IFAIL = 1
      IF (IX .LT. 1 .OR. IX .GT. NLX) RETURN
      IFAIL = 0
      IW = WID / 2
      IF (ICODE .EQ. 0) RETURN
      IF (ICODE .EQ. 1) IOUT(IX:IX) = VERT
      IF (ICODE .EQ. 2) IOUT(IX:IX) = STAR
      IF (ICODE .EQ. 3) IOUT(IX:IX) = PLUS
      IF (ICODE .LT. 4) RETURN
      IXMIN = MAX(IX - IW, 1)
      IXPLUS = MIN(IX + IW, NLX)
      IF (ICODE .EQ. 4) THEN
         DO 10 I = IXMIN, IXPLUS
         IOUT(I:I) = HORIZ
   10    CONTINUE
         RETURN
      ENDIF
      IF (ICODE .LT. 7) THEN
C
C IF WID LT 3 NO SIDES OF BOX OR NOTCHES
C
         IF (IW .LE. 0) RETURN
         IF (ICODE .EQ. 5) THEN
            IF (IX - IW .GE. 1) IOUT(IXMIN:IXMIN) = VERT
            IF (IX + IW .LE. NLX) IOUT(IXPLUS:IXPLUS) = VERT
         ENDIF
         IF (ICODE .EQ. 6) THEN
C
C IF ROOM MOVE SIDES OF BOX IN TO FORM NOTCH
C
            IF (WID .GT. 3) IXPLUS = IXPLUS - 1
            IF (IX + IW .LE. NLX) IOUT(IXPLUS:IXPLUS) = VERT
            IF (WID .GT. 3) IXMIN = IXMIN + 1
            IF (IX - IW .GE. 1) IOUT(IXMIN:IXMIN) = VERT
         ENDIF
         RETURN
      ENDIF
      IF (NOTCH .AND. WID .GT. 3) THEN
C
C MEDIAN SHORTER IF NOTCHED
C
         IXPLUS = IXPLUS - 1
         IXMIN = IXMIN + 1
      ENDIF
      DO 20 I = IXMIN, IXPLUS
      IOUT(I:I) = DASH
   20 CONTINUE
      IF (IW .GT. 0) THEN
C
C IF WID LT 3 DO NOT PUT SIDES OF BOX AT MEDIAN
C
         IF (IX - IW .GE. 1) IOUT(IXMIN:IXMIN) = VERT
         IF (IX + IW .LE. NLX) IOUT(IXPLUS:IXPLUS) = VERT
      ENDIF
      RETURN
      END