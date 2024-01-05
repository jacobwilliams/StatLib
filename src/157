      SUBROUTINE UDRUNS(X, N, UV, DV, IFAULT)
C
C        ALGORITHM AS 157 APPL. STATIST. (1981), VOL. 30, NO. 1
C
C        THE RUNS-UP AND FUNS-DOWN TEST 
C
      INTEGER UCOUNT, DCOUNT, RU, RD
      DIMENSION X(N), A(6, 6), B(6), UCOUNT(6), DCOUNT(6)
C
C        SET-UP THE A AND B MATRICES
C
      DATA
     *  A(1, 1), A(1, 2), A(1, 3), A(1, 4), A(1, 5), A(1, 6), A(2, 2), 
     *  A(2, 3), A(2, 4), A(2, 5), A(2, 6), A(3, 3), A(3, 4), A(3, 5), 
     *  A(3, 6), A(4, 4), A(4, 5), A(4, 6), A(5, 5), A(5, 6), A(6, 6)
     *  /4529.4, 9044.9, 13568.0, 18091.0, 22615.0, 27892.0, 18097.0,
     *  27139.0, 36187.0, 45234.0, 55789.0, 40721.0, 54281.0, 67852.0,
     *  83685.0, 72414.0, 90470.0, 111580.0, 113262.0, 139476.0,
     *  172860.0/
      IFAULT = 0
      IF (N .LT. 4000) GOTO 500
      DO 1 J = 2, 6
        J1 = J - 1
        DO 1 I = 1, J1
          A(J, I) = A(I, J)
    1 CONTINUE
      B(1) = 1.0 / 6.0
      B(2) = 5.0 / 24.0
      B(3) = 11.0 / 120.0
      B(4) = 19.0 / 720.0
      B(5) = 29.0 / 5040.0
      B(6) = 1.0 / 840.0
C
      DO 100 I = 1, 6
        UCOUNT(I) = 0
        DCOUNT(I) = 0
  100 CONTINUE
C
C        THE LOOP THAT ENDS AT LINE 300 DETERMINES THE NUMBER OF 
C        RUNS-UP AND RUNS-DOWN OF LENGTH I FOR I=1(1)5 AND THE NUMBER
C        OF RUNS-UP AND RUNS-DOWN OF LENGTH GREATER OR EQUAL TO 6
C
      RU = 1
      RD = 1
      DO 300 J = 2, N
C
C        THE FOLLOWING STATEMENT TESTS FOR BOTH RUNS-UP AND
C        RUNS-DOWN BREAK-POINTS.  IF A RUN-DOWN BREAK-POINT IS
C        DETECTED - GOTO 200, OTHERWISE - GOTO 150.
C        ( RU AND RD ACT AS LOCAL COUNTERS FOR THE NUMBER OF
C        RUNS-UP AND RUNS-DOWN RESPECTIVELY.)
C        A TEST IS ALSO MADE FOR DATA TIES BEWTEEN ADJACENT
C        ELEMENTS.  IF A TIE IS DETECTED - GOTO 600.
C
        IF (X(J) - X(J - 1)) 150, 600, 200
  150   UCOUNT(RU) = UCOUNT(RU) + 1
        RU = 1
        IF (RD .LT. 6) RD = RD + 1
        GOTO 300
  200   DCOUNT(RD) = DCOUNT(RD) + 1
        RD = 1
        IF (RU .LT. 6) RU = RU + 1
  300 CONTINUE
      UCOUNT(RU) = UCOUNT(RU) + 1
      DCOUNT(RD) = DCOUNT(RD) + 1
C
C         CALCULATE THE TEST STATISTICS UV AND DV.
C
      UV = 0.0
      DV = 0.0
      RN = FLOAT(N)
      DO 400 I = 1, 6
        DO 400 J = 1, 6
          UV = UV + (FLOAT(UCOUNT(I)) - RN * B(I)) *
     *       (FLOAT(UCOUNT(J)) - RN * B(J)) * A(I, J)
          DV = DV + (FLOAT(DCOUNT(I)) - RN * B(I)) *
     *       (FLOAT(DCOUNT(J)) - RN * B(J)) * A(I, J)
  400 CONTINUE
      UV = UV / RN
      DV = DV / RN
      GOTO 700
  500 IFAULT = N
      GOTO 700
  600 IFAULT = 1
  700 RETURN
      END
