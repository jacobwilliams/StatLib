      SUBROUTINE HALPRN(PARAM, NP, NCOLS, NOBS, NP1, NP2, P1LO, P1HI,
     +         P2LO, P2HI, CONF, NCONF, F, GRID, DIM1, NGRID1, NGRID2,
     +         DERIV, WK, DIMWK, LINDEP, IFAULT)
C
C     ALGORITHM AS290  APPL. STATIST. (1994) VOL. 43, NO. 1
C
C     Generate a rectangular 2-D grid of variance ratios from which
C     to plot confidence regions for two parameters using Halperin's
C     method.   If the model is linear in all parameters other than
C     the two selected, the confidence regions are exact; otherwise
C     they are approximate and the user should test the sensitivity of
C     the confidence regions to variation in the other parameters.
C
C     Auxiliary routines required: CLEAR, INCLUD, SS, TOLSET and SING,
C     from AS 274 and a routine DERIV supplied by the user to calculate 
C     first derivatives.
C
C Arguments:
C PARAM(NP)  Parameter values
C NP         Number of parameters
C NCOLS      Number of first derivatives, usually the same as NP
C NOBS       Number of observations (needed to calculate deg. of freedom
C NP1, NP2   The numbers (positions) of the two parameters for which
C            the confidence regions are needed
C P1LO, P1HI Upper and lower limits for parameter number NP1
C P2LO, P2HI Same for the other parameter
C CONF(NCONF) The confidence levels in percent
C NCONF      The number of confidence levels
C F          (Output) F ratios corresponding to the confidence levels
C GRID(DIM1, NGRID2) (Output) 2-D array of variance ratios on a regular
C            grid of values of the two parameters
C DIM1       1st dimension of array GRID in the calling program
C NGRID1, NGRID2 The numbers of grid points required
C DERIV      The name of the user's routine to calculate 1st derivatives
C WK(DIMWK)  Workspace
C DIMWK      Must be at least NCOLS(NCOLS + 9)/2
C LINDEP(NCOLS) (Output) Logical array.   LINDEP(I) = .TRUE. if the Ith
C            column of derivatives is linearly related to previous columns
C IFAULT     (Output) Error indicator
C            = 0 if no errors detected
C            = 1 if NCOLS < 2
C            = 2 if NOBS < NCOLS
C            = 4 if NP1 outside the range 1 to NP
C            = 8 ditto for NP2
C            = 16 if NP1 = NP2
C            = 32 if NGRID1 < 2 or NGRID2 < 2 or NGRID1 > DIM1
C            = 64 if DIMWK is too small
C            = -I if the Ith confidence level is not between 0 and 100
C N.B. IFAULT may equal a sum of some of the above values if there is more
C      than one error.
C
      INTEGER NP, NCOLS, NOBS, NP1, NP2, NCONF, DIM1, NGRID1, NGRID2,
     +          DIMWK, IFAULT
      REAL PARAM(NP), P1LO, P1HI, P2LO, P2HI, CONF(*), F(NCONF),
     +          GRID(DIM1, NGRID2), WK(DIMWK)
      LOGICAL LINDEP(NCOLS)
      EXTERNAL DERIV
C
C     Local variables
C
      INTEGER I, NDF, NRBAR, XPTR, DPTR, THPTR, RSSPTR, TOLPTR, I1, I2,
     +        IOBS
      REAL ZERO, ONE, HUNDRD, HALF, TWO, Q, RESID, SSERR, STEP1, STEP2,
     +          WT, P1SAVE, P2SAVE, SSQ1
      DATA ZERO/0.0/, ONE/1.0/, HUNDRD/100.0/, HALF/0.5/, TWO/2.0/,
     +          WT/1.0/
C
C     Some checks
C
      IFAULT = 0
      IF (NCOLS .LT. 2) IFAULT = 1
      IF (NOBS .LE. NCOLS) IFAULT = IFAULT + 2
      IF (NP1 .LT. 1 .OR. NP1 .GT. NP) IFAULT = IFAULT + 4
      IF (NP2 .LT. 1 .OR. NP2 .GT. NP) IFAULT = IFAULT + 8
      IF (NP1 .EQ. NP2) IFAULT = IFAULT + 16
      IF (NGRID1 .LT. 2 .OR. NGRID2 .LT. 2 .OR. NGRID1 .GT. DIM1)
     +           IFAULT = IFAULT + 32
      IF (IFAULT .GT. 0) RETURN
C
C     Set up pointers to the workspace.
C
      NRBAR = NCOLS * (NCOLS - 1) / 2
      XPTR = NRBAR + 1
      DPTR = XPTR + NCOLS
      THPTR = DPTR + NCOLS
      RSSPTR = THPTR + NCOLS
      TOLPTR = RSSPTR + NCOLS
      IF (TOLPTR + NCOLS - 1 .GT. DIMWK) THEN
        IFAULT = 64
        RETURN
      END IF
C
C     Calculate step sizes.
C
      STEP1 = (P1HI - P1LO) / (NGRID1 - 1)
      STEP2 = (P2HI - P2LO) / (NGRID2 - 1)
      P1SAVE = PARAM(NP1)
      P2SAVE = PARAM(NP2)
C----------------------------------------------------------------------
C
C     Start of cycle through the grid.
C
      PARAM(NP1) = P1LO
      DO 30 I1 = 1, NGRID1
        PARAM(NP2) = P2LO
        DO 20 I2 = 1, NGRID2
C
C     Initialize orthogonal reduction.
C
          CALL CLEAR(NCOLS, NRBAR, WK(DPTR), WK, WK(THPTR), SSERR,
     *               IFAULT)
          DO 10 IOBS = 1, NOBS
            CALL DERIV(PARAM, NP, NCOLS, IOBS, WK(XPTR), RESID, WT)
C
C     Re-order the derivatives so that those for parameters NP1 & NP2
C     are at the end.
C
            IF (NP1 .LT. NCOLS-1) THEN
              IF (NP2 .NE. NCOLS-1) THEN
                Q = WK(XPTR + NCOLS - 2)
                WK(XPTR + NCOLS - 2) = WK(XPTR + NP1 - 1)
                WK(XPTR + NP1 - 1) = Q
              ELSE
                Q = WK(XPTR + NCOLS - 1)
                WK(XPTR + NCOLS - 1) = WK(XPTR + NP1 - 1)
                WK(XPTR + NP1 - 1) = Q
              END IF
            END IF
            IF (NP2 .LT. NCOLS-1) THEN
              Q = WK(XPTR + NCOLS - 1)
              WK(XPTR + NCOLS - 1) = WK(XPTR + NP2 - 1)
              WK(XPTR + NP2 - 1) = Q
            END IF
C
C     Update QR factorization.
C
            CALL INCLUD(NCOLS, NRBAR, WT, WK(XPTR), RESID, WK(DPTR), WK,
     +              WK(THPTR), SSERR, IFAULT)
   10     CONTINUE
C
C----------------------------------------------------------------------
C
C     Test for singularities
C
      CALL TOLSET(NCOLS, NRBAR, WK(DPTR), WK, WK(TOLPTR), WK(XPTR),
     +            IFAULT)
      CALL SING(NCOLS, NRBAR, WK(DPTR), WK, WK(THPTR), SSERR,
     +          WK(TOLPTR), LINDEP, WK(XPTR), IFAULT)
C
C     If IFAULT < 0, it contains the rank deficiency.
C
      NDF = NOBS - NCOLS - IFAULT
C
C     Calculate variance ratio.
C
          CALL SS(NCOLS, WK(DPTR), WK(THPTR), SSERR, WK(RSSPTR), IFAULT)
          IF (NCOLS .GT. 2) THEN
            SSQ1 = WK(RSSPTR + NCOLS - 3)
          ELSE
            SSQ1 = WK(RSSPTR) + WK(DPTR)*WK(THPTR)**2
          END IF
          GRID(I1, I2) = HALF * (SSQ1 - SSERR) / (SSERR / NDF)
          PARAM(NP2) = PARAM(NP2) + STEP2
   20   CONTINUE
        PARAM(NP1) = PARAM(NP1) + STEP1
   30 CONTINUE
C
      PARAM(NP1) = P1SAVE
      PARAM(NP2) = P2SAVE
C
C     Convert confidence levels in % to F-values, if NCONF > 0.
C
      DO 40 I = 1, NCONF
        Q = ONE - CONF(I) / HUNDRD
        IF (Q .LE. ZERO .OR. Q .GT. ONE) THEN
          IFAULT = -I
          RETURN
        END IF
        F(I) = HALF * NDF * (Q**(-TWO/NDF) - ONE)
   40 CONTINUE
C
      RETURN
      END
