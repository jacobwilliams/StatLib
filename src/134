      SUBROUTINE BET1GT(BETAR)
C
C        ALGORITHM AS 134.1 APPL. STATIST. (1979) VOL.28 NO.1
C
C        GENERATES BETA RANDOM VARIABLES
C        WITH 0.0 .LT. ALPHA .LT. 1.0 AND BETA  .GT.  1.0
C        BY THE SWITCHING ALGORITHM OF ATKINSON AND WHITTAKER
C
C     The random number call has been changed from RANF(..) to RAND()
C     so that the Wichmann & Hill random number generator, AS 183
C     can be used.   To conform with Fortran-77, SAVE instructions have
C     been added.
C
      REAL ONE
      DATA ONE/1.0/
      COMMON /BETCOM/ AM1, BM1, ARECIP, BRECIP, T, R
      SAVE /BETCOM/
C
    1 UVAR = RAND()
      CALL FNE(EVAR)
      IF(UVAR .GT. R) GOTO 3
C
C     WARNING.  With some compilers, underflow can occur for small
C               values of alpha in executing the next instruction.
C
      BETAR = T*((UVAR/R)**ARECIP)
      IF(-BM1*LOG(ONE-BETAR) .GT. EVAR) GOTO 1
      RETURN
    3 BETAR = ONE - (ONE-T)*(((ONE-UVAR)/(ONE-R))**BRECIP)
      IF(-AM1*LOG(BETAR/T) .GT. EVAR) GOTO 1
      RETURN
      END
C
      SUBROUTINE TOPT(ALPHA, BETA, IFAULT)
C
C        ALGORITHM AS 134.2 APPL. STATIST. (1979) VOL.28 NO.1
C
C        SETS CONSTANTS FOR BETA VARIABLE GENERATION.
C        TOL IS TOLERANCE FOR ACCEPTING T AS SOLUTION OF GX = 0
C
      COMMON /BETCOM/ AM1, BM1, ARECIP, BRECIP, T, R
      SAVE /BETCOM/
      REAL ZERO, ONE, TWO, TOL
      DATA ZERO/0.0/, ONE/1.0/, TWO/2.0/, TOL/1.0E-5/
C
C        TEST PARAMETERS AND CALCULATE CONSTANTS
C
      IF(ALPHA .GT. ZERO .AND. ALPHA .LT. ONE) GOTO 1
      IFAULT = 1
      RETURN
    1 IF(BETA .GT. ONE) GOTO 2
      IFAULT = 2
      RETURN
    2 AM1 = ALPHA - ONE
      TTILDE = AM1/(AM1-BETA)
      IF(ONE-TTILDE*TTILDE .LT. ONE) GOTO 3
      IFAULT = 3
      RETURN
    3 BM1 = BETA - ONE
      ARECIP = ONE/ALPHA
      BRECIP = ONE/BETA
C
C        INITIAL VALUES FOR FALSE POSITION
C
      XA = ZERO
      GA = AM1
      XB = ONE
      GB = BETA
C
C        FALSE POSITION
C
    4 X = (XA*GB - XB*GA)/(GB - GA)
      AA = ONE - X
      BT = BETA*X
      GX = BT + AA**BM1*(AM1*AA - BT)
      IF(ABS(GX) .LT. TOL) GOTO 6
      IF(GX*GA .LT. ZERO) GOTO 5
C
C        MODIFICATION
C
      XA = X
      GA = GX
      GB = GB/TWO
      GOTO 4
C
C        MODIFICATION
C
    5 XB = XA
      GB = GA
      XA = X
      GA = GX
      GOTO 4
C
C        CALCULATE REMAINING CONSTANTS
C
    6 T = X
      R = BT/(BT + ALPHA*AA**BETA)
      RETURN
      END
C
      SUBROUTINE FNE(REX)
C
C        ALGORITHM AS134.3 APPL. STATIST. (1979) VOL.28 NO.1
C
C        GENERATES EXPONENTIAL RANDOM VARIABLES
C        BY THE METHOD OF VON NEUMANN
C
      REAL ZERO, ONE
      DATA ZERO/0.0/, ONE/1.0/
      A = ZERO
    1 U = RAND()
      UO = U
    2 USTAR = RAND()
      IF(U .LT. USTAR) GOTO 3
      U = RAND()
      IF(U .LT. USTAR) GOTO 2
      A = A + ONE
      GOTO 1
    3 REX = A + UO
      RETURN
      END

