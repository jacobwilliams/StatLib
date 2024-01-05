      REAL FUNCTION PSI2(X, P, Q, A2, DELTA, MAXITR, IFAULT)
C
C        ALGORITHM AS 278.1 APPL.STATIST. (1992) VOL.41, NO.3
C
C        Calculates the probability that a random variable distributed
C        according to the psi-square distribution with P and Q degrees
C        of freedom and A2 eccentricity parameter, is less than or
C        equal to X
C
      REAL X, P, Q, A2, DELTA
      INTEGER MAXITR, IFAULT
C
      REAL ALNGAM, BETAIX
      EXTERNAL ALNGAM, BETAIX
C
      REAL AQAL, BL, CL, DAX, DJ, DJ1, ERP, ERR, EXPL, G,
     *     GCJ, GCL, HALF, ONE, PQ2, PREC, P2, QQAL, Q2, R,
     *     RL, SUM, SUMC, SUMCP, TWO, V, W, XLOW, XP, Y, YYL,
     *     Z, ZERO, ZM, ZN
      INTEGER IOK, J, JJ, JJJ, JOK, KOK
      LOGICAL LALF
C
      PARAMETER (ZERO=0.0E0, HALF=0.5E0, ONE=1.0E0, TWO=2.0E0)
C
C        machine-dependent constants
C
      PARAMETER (PREC=1.0E-6, XLOW=1.2E-38, EXPL=-87.315E0)
      PARAMETER (XP=XLOW/PREC)
C
      IFAULT = 0
C
C        test for valid input arguments
C
      IF (X .LT. ZERO .OR. A2 .LT. ZERO
     *   .OR. P .LE. ZERO .OR. Q .LE. ZERO
     *   .OR. DELTA .GE. ONE .OR. DELTA .LE. PREC) THEN
         IFAULT = 1
         PSI2 = -ONE
         RETURN
      ENDIF
C
C        define usefull parameters
C
      P2  = P * HALF
      Q2  = Q * HALF
      PQ2 = P2 + Q2
C
C        case A2 = 0
C
      IF (A2 .LT. PREC) THEN
         PSI2 = BETAIX( P*X/(P*X+Q), P2, Q2, IOK )
         IF (IOK .NE. 0) THEN
            IFAULT = 4
            PSI2 = -ONE
         ENDIF
         RETURN
      ENDIF
C
C        case P = 1
C
      IF (ABS(P-ONE) .LT. PREC) THEN
         DAX  = TWO * SQRT( A2*X )
         PSI2 = HALF * ( SIGN( ONE, X-A2 ) *
     *          BETAIX( (A2+X-DAX)/(A2+X-DAX+Q),HALF,Q2,IOK ) +
     *          BETAIX( (A2+X+DAX)/(A2+X+DAX+Q),HALF,Q2,JOK ) )
         IF (IOK .NE. 0 .OR. JOK .NE. 0) THEN
            IFAULT = 4
            PSI2 = -ONE
         ENDIF
         RETURN
      ENDIF
C
      Y = P*X / (A2+Q+P*X)
C
C        case P = Q and Y = 0.5
C
      IF (ABS(Y-HALF) .LT. PREC .AND. ABS(P-Q) .LT. PREC) THEN
         PSI2 = HALF
         RETURN
      ENDIF
C
C        calculate 1-F(X) or F(X) (LALF is the indicator)
C
      IF (Y .GT. HALF .OR. (Y .EQ. HALF .AND. P .GT. Q) ) THEN
         LALF = .TRUE.
         Z    = ONE - Y
         ZM   = Q2
         ZN   = P2
        ELSE
         LALF = .FALSE.
         Z    = Y
         ZM   = P2
         ZN   = Q2
      ENDIF
C
C        Y near 0 or 1 ?
C
      IF (Z .LT. PREC) THEN
         IF( LALF )THEN
            PSI2 = ONE
           ELSE
            PSI2 = ZERO
         ENDIF
         RETURN
      ENDIF
C
C        General case (iterations)
C        G's are decreasing for J >= JJ = -V/W + 1
C        CL's  are decreasing for J >= JJJ = A2/2 + 1
C        GCJ, SUMC, SUMCP are only used for stopping rule
C        logs are used to avoid underflows
C
      ERR  = DELTA * HALF
      ERP  = ERR / PREC
      YYL  = LOG( Y*(ONE-Y) )
      QQAL = Q2 * LOG( Q/(Q+A2) )
      AQAL = LOG( A2/(Q+A2) )
      V    = PQ2 * Z - ZN
      W    = Z + Z - ONE
C
C        initialize
C
      BL = ZERO
      CL = QQAL
      G  = BETAIX( Z, ZM, ZN, IOK )
      IF (IOK .NE. 0) THEN
         IFAULT = 4
         PSI2 = -ONE
         RETURN
      ENDIF
      RL = P2 * LOG(Y)  +  Q2 * LOG(ONE-Y)  +  ALNGAM(PQ2,IOK)
     *     - ALNGAM(P2,JOK) - ALNGAM(Q2,KOK) - LOG(P2) - LOG(Q2)
      IF (IOK+JOK+KOK .NE. 0) THEN
         IFAULT = 5
         PSI2 = -ONE
         RETURN
      ENDIF
      IF (RL .GE. EXPL) THEN
         R = EXP( RL )
        ELSE
         R = ZERO
      ENDIF
      SUM   = ZERO
      GCJ   = ZERO
      SUMC  = ZERO
      SUMCP = SUMC
      DJ    = ZERO
C
C        define minimum numbers of iterations (JJ and JJJ)
C
      IF (ZM .GT. ZN .AND. W .NE. ZERO) THEN
         JJ = INT( -V/W ) + 1
        ELSE
         JJ = 0
      ENDIF
      JJJ = INT( A2*HALF ) + 1
C
C        iteration loop
C
      DO 10 J = 0, MAXITR
C
         IF (G .GT. ZERO) THEN
            GCL = LOG(G) + CL
            IF (GCL .GE. EXPL) THEN
               GCL = EXP( GCL )
               SUM = SUM + GCL
               GCJ = GCJ + GCL*DJ
C                 check loss of accuracy
               IF (GCJ .GE. ERP) THEN
                  IFAULT = 3
                  GOTO 20
               ENDIF
            ENDIF
         ENDIF
C
         IF (CL .GE. EXPL) SUMC = SUMC + EXP(CL)
C
C           check accuracy (stopping rule)
C
         IF (J .GE. JJ) THEN
C              XP is used to prevent possible underflow
            IF (GCJ .GE. XP) THEN
               IF (PREC*GCJ+G*(ONE-SUMC) .LT. ERR) GOTO 20
              ELSE
               IF (         G*(ONE-SUMC) .LT. ERR) GOTO 20
            ENDIF
            IF (J .GE. JJJ .AND. SUMC .EQ. SUMCP .AND.
     *         ABS(ONE-SUMC) .GE. PREC) THEN
               IFAULT = 3
               GOTO 20
            ENDIF
         ENDIF
C
C           prepare next iteration
C
         SUMCP = SUMC
         DJ1   = DJ + ONE
         BL    = BL + LOG( (Q2+DJ)/DJ1 )
         CL    = BL + QQAL + DJ1*AQAL
         G     = G  + (V+W*DJ)*R
         RL    = RL + YYL +
     *           LOG( ((DJ+DJ+PQ2)/(DJ1+P2)) * ((DJ+DJ1+PQ2)/(DJ1+Q2)) )
         IF (RL .GE. EXPL) THEN
            R = EXP( RL )
           ELSE
            R = ZERO
         ENDIF
         DJ = DJ1

   10 CONTINUE
C
C        we get here if maximum number of iterations is reached
C
      IFAULT = 2
C
C        the end
C
   20 CONTINUE
      IF (LALF) THEN
         PSI2 = ONE - SUM
        ELSE
         PSI2 = SUM
      ENDIF
C
      END
      REAL FUNCTION BETAIX(X, P, Q, IFAULT)
C
C        ALGORITHM AS 278.2 APPL.STATIST. (1992) VOL.41, NO.3
C
      REAL X, P, Q
      INTEGER IFAULT
C
      REAL ALNGAM, BETAIN
      EXTERNAL ALNGAM, BETAIN
C
      REAL BETA, EXPL, ZERO
      INTEGER IOK, JOK, KOK
C
      PARAMETER (ZERO=0.0E0)
C
C        machine-dependent constant
C
      PARAMETER (EXPL=-87.315E0)
C
      BETAIX = ZERO
C
C        calculate beta parameter for betain
C
      BETA = ALNGAM(P,JOK) + ALNGAM(Q,KOK) - ALNGAM(P+Q,IOK)
      IF (IOK+JOK+KOK .NE. 0) THEN
         IFAULT = 5
         RETURN
      ENDIF
      IF (BETA .GE. EXPL) THEN
         BETA =  EXP(BETA)
        ELSE
         IFAULT = 4
         RETURN
      ENDIF
      BETAIX = BETAIN(X, P, Q, BETA, IFAULT)
      END

