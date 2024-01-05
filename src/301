      SUBROUTINE LOGF1(CELL, DELTA, EPS, INQCOL, INQROW, LQ, R, RS,
     *     S, UQ, COL, IFAULT, LFCELL, LFCOL, LFROW, LOGEF1, LOGPH,
     *     LOGPHP, MU1, N, QCELL, QCOL, QROW, ROW)
C
C        ALGORITHM AS 301.1 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes the logarithms of F1, P(H) and P(H')
C
      INTEGER R, RS, S, CELL(RS), COL(S), IFAULT, N, ROW(R)   
      REAL DELTA, EPS, INQCOL(S), INQROW(R), LQ, UQ, LOGEF1, LOGPH, 
     *     LOGPHP, MU1, QCELL(RS), QCOL(S), QROW(R) 
      DOUBLE PRECISION LFCELL(RS), LFCOL(S), LFROW(R) 
C
      INTEGER I, J, K 
      DOUBLE PRECISION LFM
      REAL C1, C2, INTOFG, LGPHI1, LGPHI2, LGPHI3, LOGFY,
     *     MLOGG, RPI, XMAX
C
      EXTERNAL MRGINS, KAPPA0, PARRHO, LIMITS, MXLOGG
C
      IFAULT = 0
      CALL MRGINS(CELL, INQCOL, INQROW, R, RS, S, COL, IFAULT,
     *            LFCELL, LFCOL, LFM, LFROW, N, ROW)
      IF (IFAULT .EQ. 0) THEN
         CALL KAPPA0(N, DELTA, EPS, LFM, LFROW,
     *               ROW, INQROW, R, IFAULT, QROW)
         CALL KAPPA0(N, DELTA, EPS, LFM, LFCOL,
     *               COL, INQCOL, S, IFAULT, QCOL)
         CALL PARRHO(LQ, UQ, IFAULT, MU1, RPI)
         IF (IFAULT .NE. 1) THEN
C
C        Compute log(Fisher-Yates probability) and the cell
C        probabilities
C
            LOGFY = -SNGL(LFM)
            K = 0
            DO 10 I = 1, R
               LOGFY = LOGFY + SNGL(LFROW(I))
               DO 5 J = 1, S
                  K = K + 1
                  QCELL(K) = QROW(I)*QCOL(J)
                  LOGFY = LOGFY - SNGL(LFCELL(K))
    5          CONTINUE
   10       CONTINUE
            DO 15 J = 1, S
   15       LOGFY = LOGFY + SNGL(LFCOL(J))
C
C        Compute log PHI1 = log INTOFG((N(I, J)), EXP(X)(Q(I.)Q(.J)))
C
            CALL LIMITS(N, LFM, LFCELL, CELL, QCELL, RS,    C1, C2)
            CALL MXLOGG(N, C1, C2, LFM, LFCELL, CELL,
     *                  QCELL, RS, 1, 1, MLOGG, XMAX)
            LGPHI1 = ALOG(RPI*INTOFG(N, C1, C2, DELTA, EPS, LFM, LFCELL,
     *              CELL, MLOGG, QCELL, RS, 1, XMAX, IFAULT)) + MLOGG
C
C        Compute log PHI2 = log INTOFG((N(I.)), EXP(X)(Q(I.)))
C
            CALL LIMITS(N, LFM, LFROW, ROW, QROW, R,    C1, C2)
            CALL MXLOGG(N, C1, C2, LFM, LFROW, ROW,
     *                  QROW, R, 1, 1, MLOGG, XMAX)
            LGPHI2 = ALOG(RPI*INTOFG(N, C1, C2, DELTA, EPS, LFM,
     *          LFROW, ROW, MLOGG, QROW, R, 1, XMAX, IFAULT)) + MLOGG
C
C        Compute log PHI3 = log INTOFG((N(.J)), EXP(X)(Q(.J)))
C
            CALL LIMITS(N, LFM, LFCOL, COL, QCOL, S,    C1, C2)
            CALL MXLOGG(N, C1, C2, LFM, LFCOL, COL,
     *                  QCOL, S, 1, 1, MLOGG, XMAX)
            LGPHI3 = ALOG(RPI*INTOFG(N, C1, C2, DELTA, EPS, LFM,
     *          LFCOL, COL, MLOGG, QCOL, S, 1, XMAX, IFAULT)) + MLOGG
C
C        Compute the logarithms of F1, P(H) and P(H')
C
            LOGEF1 = LGPHI1 - (LOGFY + LGPHI2 + LGPHI3)
            LOGPH = LOGFY + LGPHI3
            LOGPHP = LGPHI1 - LGPHI2
         ENDIF
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE MRGINS(CELL, INQCOL, INQROW, R, RS, S, COL, IFAULT,
     *                  LFCELL, LFCOL, LFM, LFROW, N, ROW)
C
C        ALGORITHM AS 301.2 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes the table total N, the row and column totals, N(I.)
C        and N(.J), log N(I.)!, log N(.J)!, and checks whether (Q(I.))
C        and (Q(.J)) sum to 1, CELL(I) < 0 and N < 1 or N > 10E+09
C
      INTEGER R, RS, S, CELL(RS), COL(S), IFAULT, N, ROW(R) 
      REAL INQCOL(S), INQROW(R) 
      DOUBLE PRECISION LFCELL(RS), LFCOL(S), LFM, LFROW(R)
C
      DOUBLE PRECISION LFACT 
      REAL APPONE, NINES, SUMQ, ZERO
      INTEGER I, J, K, L
C
      DATA APPONE, NINES, ZERO /1.000001, 0.999999, 0.0/
C
      DO 5 I = 1, RS
         IF (CELL(I) .LT. 0) IFAULT = 2
    5 CONTINUE
      IF (IFAULT .NE. 2) THEN
         SUMQ = ZERO
         N = 0
         DO 15 I = 1, R
            SUMQ = SUMQ + INQROW(I)
            ROW(I) = 0
            L = (I - 1)*S
            DO 10 J = 1, S
               K = L + J
               ROW(I) = ROW(I) + CELL(K)
               LFCELL(K) = LFACT(CELL(K))
   10       CONTINUE
            LFROW(I) = LFACT(ROW(I))
            N = N + ROW(I)
   15    CONTINUE
         IF (N .GE. 1 .AND. N .LE. 1000000000) THEN
            IF (SUMQ .LE. NINES .OR. SUMQ .GE. APPONE) THEN
               IFAULT = 3
            ELSE
               SUMQ = ZERO
               LFM = LFACT(N)
               DO 25 J = 1, S
                  SUMQ = SUMQ + INQCOL(J)
                  COL(J) = 0
                  K = J
                  DO 20 I = 1, R
                     COL(J) = COL(J) + CELL(K)
                     K = K + S
   20             CONTINUE
                  LFCOL(J) = LFACT(COL(J))
   25          CONTINUE
               IF (SUMQ .LE. NINES .OR. SUMQ .GE. APPONE) IFAULT = 3
            ENDIF
         ELSE
            IFAULT = 4
         ENDIF
      ENDIF
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION LFACT(N)
C
C        ALGORITHM AS 301.3 APPL. STATIST. (1996), VOL.45, NO.2
C
C       Computes log Gamma(N + 1) = log N!
C
      INTEGER N
C
      DOUBLE PRECISION C1, C2, LF(0:4), L2PID2, ONE, ONED2, U
C
      DATA C1, C2 /0.08333333333333333D0, 0.0027777777777777778D0/,
     *     L2PID2, LF(0), LF(1) /0.9189385332046728D0, 0.0D0, 0.0D0/,
     *     LF(2), LF(3) /0.6931471805599454D0, 1.791759469228055D0/,
     *     LF(4), ONE, ONED2 /3.178053830347946D0, 1.0D0, 0.5D0/
C
      IF (N .LE. 4) THEN
         LFACT = LF(N)
      ELSE
         U = DBLE(N + 1)
      LFACT = (U - ONED2)*DLOG(U) - U + L2PID2 + C1/U - C2*(ONE/U)**3
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE KAPPA0(BIGM, DELTA, EPS, LFM, LOGFAC,
     *                  M, PRQ, T, IFAULT, Q)
C
C        ALGORITHM AS 301.4 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes K0 and then updates (Q(I.)) and (Q(.J))
C
      INTEGER BIGM, T, M(T), IFAULT
      REAL DELTA, EPS, PRQ(T), Q(T)
      DOUBLE PRECISION LFM, LOGFAC(T)
C
      REAL C1, C2, DUMMY1, DUMMY2, FIFTY, INTOFG, K0,
     *     LGNUM, LGDNOM, MLOGG, TEN, XMAX, ZERO
      INTEGER I
C
      DATA FIFTY, TEN, ZERO /50.0, 10.0, 0.0/
C
      EXTERNAL PARRHO, LIMITS, MXLOGG
C
      CALL PARRHO(TEN, FIFTY, IFAULT, DUMMY1, DUMMY2)
      CALL LIMITS(BIGM, LFM, LOGFAC, M, PRQ, T, C1, C2)
      CALL MXLOGG(BIGM, ZERO, C2, LFM, LOGFAC, M,
     *            PRQ, T, 2, 2, MLOGG, XMAX)
      LGNUM = ALOG(INTOFG(BIGM, ZERO, C2, DELTA, EPS, LFM,
     *            LOGFAC, M, MLOGG, PRQ, T, 2, XMAX, IFAULT)) + MLOGG
      CALL MXLOGG(BIGM, C1/BIGM, ZERO, LFM, LOGFAC, M,
     *            PRQ, T, 3, 3, MLOGG, XMAX)
      LGDNOM = ALOG(INTOFG(BIGM, C1/BIGM, ZERO, DELTA, EPS, LFM,
     *            LOGFAC, M, MLOGG, PRQ, T, 3, XMAX, IFAULT)) + MLOGG
      K0 = EXP(LGNUM - LGDNOM)
      DO 5 I = 1, T
         Q(I) = (FLOAT(M(I)) + PRQ(I)*K0)/(BIGM + K0)
    5 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE LIMITS(BIGM, LFM, LOGFAC, M, Q, T, C1, C2)
C
C        ALGORITHM AS 301.5 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes the logarithms of the limits C1 and C2
C
      INTEGER BIGM, T, M(T)
      DOUBLE PRECISION LFM, LOGFAC(T)
      REAL Q(T), C1, C2 
C
      DOUBLE PRECISION DPC2
      REAL ZERO
      INTEGER I
C
      DATA ZERO /0.0/
C
      C1 = ZERO
      DPC2 = LFM
      DO 5 I = 1, T
         DPC2 = DPC2 + DBLE(M(I))*DLOG(DBLE(Q(I))) - LOGFAC(I)
         IF (M(I) .EQ. BIGM) C1 = ALOG(Q(I))
    5 CONTINUE
      C2 = SNGL(DPC2)
C
      RETURN
      END
C
      SUBROUTINE MXLOGG(BIGM, C1, C2, LFM, LOGFAC, M,
     *                  Q, T, TAU, ZETA, MLOGG, XMAX)
C
C        ALGORITHM AS 301.6 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Finds the approximate maximum of log G(X), 
C        log G(X)/(1 + M/EXP(X)) or log G(X)/(M + EXP(X)) where 
C        M denotes BIGM.
C
      INTEGER BIGM, T, M(T), TAU, ZETA
      REAL C1, C2, Q(T), MLOGG, XMAX
      DOUBLE PRECISION LFM, LOGFAC(T)
C
      DOUBLE PRECISION D1, D2, DPZERO, EPS, LOGDM
      REAL CONDTN, KMAX, LMB, LMB2, LMU1, LXLMU1, ONE,
     *     ONED2, RHO, RLMB, SAVEK, THRTY2, U, Y, ZERO
      INTEGER I
C
      COMMON /RHOPAR/ LMB, LMB2, LMU1, RLMB
C
      DATA DPZERO, ONE, ONED2, EPS, THRTY2, ZERO
     *     /0.0D0, 1.0, 0.5, 0.001D0, 32.0, 0.0/
C
      RHO(U) = LMB/(LMB2 + (U - LMU1)**2)
      LXLMU1 = SNGL(LOGDM(BIGM, LFM, LOGFAC, M, Q, T, TAU, LMU1)) +
     *         ALOG(RHO(LMU1))
      IF (C1 .NE. ZERO) THEN
         IF (C1 + ALOG(RHO(-THRTY2)) .LT. LXLMU1) THEN
            XMAX = LMU1
            MLOGG = LXLMU1
         ELSE
            XMAX = -THRTY2
            MLOGG = C1 + ALOG(RHO(-THRTY2))
         ENDIF
      ELSE
         CONDTN = ONE
         DO 5 I = 1, T
            Y = Q(I)*BIGM
            CONDTN = CONDTN + ((FLOAT(M(I)) - Y)**2 - FLOAT(M(I)))/Y
    5    CONTINUE
         IF (CONDTN .LE. ZERO .AND. C2 .NE. ZERO) THEN
            IF (C2 + ALOG(RHO(THRTY2)) .LT. LXLMU1) THEN
               XMAX = LMU1
               MLOGG = LXLMU1
            ELSE
               XMAX = THRTY2
               MLOGG = C2 + ALOG(RHO(THRTY2))
            ENDIF
         ELSE
C
C        Find the zero of d(log F(X))/dx, 
C        d(log F(X)/(1 + M/EXP(X)))/dx or d(log F(X)/(M + EXP(X)))/dx
C        by Newton-Raphson's method where M denotes BIGM.
C
            KMAX = ONE
            SAVEK = KMAX
            CALL DERIV(BIGM, M, Q, T, KMAX, ZETA,    D1, D2)
   10       IF (DABS(D1) .LT. EPS .OR. ABS(ALOG(KMAX)) .GT .THRTY2)
     *          GO TO 15
               IF (D2 .LT. DPZERO) THEN
                  KMAX = KMAX - D1/D2
                  IF (KMAX .LT. ZERO) THEN
                     SAVEK = ONED2*SAVEK
                     KMAX = SAVEK
                  ENDIF
               ELSE
                  SAVEK = ONED2*SAVEK
                  KMAX = SAVEK
               ENDIF
               CALL DERIV(BIGM, M, Q, T, KMAX, ZETA,    D1, D2)
               GO TO 10
   15       XMAX = ALOG(KMAX)
            IF (XMAX .GT. THRTY2) XMAX = THRTY2
            MLOGG = SNGL(LOGDM(BIGM, LFM, LOGFAC, M, Q, T, TAU, XMAX))
     *              + ALOG(RHO(XMAX))
            IF (MLOGG .LT. LXLMU1) THEN
               XMAX = LMU1
               MLOGG = LXLMU1
            ENDIF
            IF (TAU .LT. 3 .AND. MLOGG .LT. C2 + ALOG(RHO(THRTY2))) THEN
               XMAX = THRTY2
               MLOGG = C2 + ALOG(RHO(THRTY2))
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE DERIV(BIGM, M, Q, T, K, ZETA, D1, D2)
C
C        ALGORITHM AS 301.7 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes the derivatives of d(log F(X))/dx,
C        d(log F(X)/(1 + M/EXP(X)))/dx or d(log F(X)/(M + EXP(X)))/dx
C        where M denotes BIGM.
C
      INTEGER BIGM, T, M(T), ZETA
      REAL Q(T), K
      DOUBLE PRECISION D1, D2
C
      INTEGER ALPHA, I
      DOUBLE PRECISION C1, C2, ONE, ONED2, PSIAPP, PSIPAP,
     *                 TEN, TWO, U, Z, ZERO
C
      DATA C1, C2 /0.08333333333333333D0, 0.1666666666666667D0/,
     * ONE, ONED2, TEN, TWO, ZERO /1.0D0, 0.5D0, 10.0D0, 2.0D0, 0.0D0/
C
      PSIAPP(U) = DLOG(U) - ONED2/U - C1*(ONE/U)**2
      PSIPAP(U) = ONE/U + ONED2*(ONE/U)**2 + C2*(ONE/U)**3
      D1 = ZERO
      D2 = ZERO
      DO 15 ALPHA = 1, T
         IF (M(ALPHA) .GT. 0) THEN
            Z = DBLE(Q(ALPHA)*K)
            I = 0
    5       IF (I + Z .GT. TEN .OR. I .EQ. M(ALPHA)) GO TO 10
               D1 = D1 + Q(ALPHA)/(I + Z)
               D2 = D2 - (Q(ALPHA)/(I + Z))**2
               I = I + 1
               GO TO 5
   10       IF (I .LT. M(ALPHA)) THEN
               D1 = D1 + Q(ALPHA)*(PSIAPP(M(ALPHA) + Z) - PSIAPP(I + Z))
               D2 = D2 + Q(ALPHA)**2*(PSIPAP(M(ALPHA) + Z) -
     *              PSIPAP(I + Z))
            ENDIF
         ENDIF
   15 CONTINUE
      I = 0
   20 IF (DBLE(I + K) .GT. TEN .OR. I .EQ. BIGM) GO TO 25
         D1 = D1 - ONE/(I + K)
         D2 = D2 + (ONE/(I + K))**2
         I = I + 1
         GO TO 20
   25 IF (I .LT. BIGM) THEN
         D1 = D1 - PSIAPP(DBLE(BIGM + K)) + PSIAPP(DBLE(I + K))
         D2 = D2 - PSIPAP(DBLE(BIGM + K)) + PSIPAP(DBLE(I + K))
      ENDIF
      IF (ZETA .EQ. 2) THEN
         D1 = D1 + ONE/(K*(ONE + K/BIGM))
         D2 = D2 - (ONE + TWO*K/BIGM)*(ONE/(K*(ONE + K/BIGM)))**2
      ELSE IF (ZETA .EQ. 3) THEN
         D1 = D1 - ONE/(BIGM + K)
         D2 = D2 + (ONE/(BIGM + K))**2
      ENDIF
C
      RETURN
      END
C
      REAL FUNCTION INTOFG(BIGM, C1, C2, DELTA, EPS, LFM, LGF,
     *                     M, MLOGG, Q, T, TAU, XMAX, IFAULT)
C
C        ALGORITHM AS 301.8 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes the integral of G(X), G(X)/(1 + M/EXP(X)) or
C        G(X)/(M + EXP(X)) where M denotes BIGM.
C
      INTEGER BIGM, T, M(T), TAU, IFAULT
      REAL C1, C2, DELTA, EPS, MLOGG, Q(T), XMAX
      DOUBLE PRECISION LFM, LGF(T)
C
      REAL C, CONST, ONE,
     *     QG8, THRTY2, XLIMIT, YL, YU, ZERO
      INTEGER CVERG, I
C
      DATA CONST, ONE, THRTY2, ZERO /-23.0, 1.0, 32.0, 0.0/
C
      EXTERNAL CHKCVG
C
      INTOFG = ZERO
C
C        Integrate to the left from XMAX and then to the right
C        from XMAX
C
      XLIMIT = -THRTY2
      C = C1
      YL = -(ONE - ONE/DELTA)
      YU = ZERO
      DO 15 I = 1, 2
         CVERG = 0
         IF (XMAX .EQ. XLIMIT) CALL CHKCVG(BIGM, C, MLOGG + CONST,
     *     DELTA, EPS, LFM, LGF, M, MLOGG, Q, I, T, TAU, XMAX, CVERG,
     *     IFAULT, INTOFG, ZERO, ZERO)
    5    IF (CVERG .EQ. 1) GO TO 10
            INTOFG = INTOFG + QG8(BIGM, LFM, LGF, M, MLOGG, Q, T,
     *                            TAU, XMAX + YL, XMAX + YU)
            IF (YL .LT. ZERO) THEN
            CALL CHKCVG(BIGM, C, MLOGG + CONST, DELTA, EPS, LFM, LGF,
     *         M, MLOGG, Q, I, T, TAU, XMAX, CVERG, IFAULT, INTOFG,
     *         YU, YL)
            ELSE
            CALL CHKCVG(BIGM, C, MLOGG + CONST, DELTA, EPS, LFM, LGF,
     *         M, MLOGG, Q, I, T, TAU, XMAX, CVERG, IFAULT, INTOFG,
     *         YL, YU)
            ENDIF
            GO TO 5
   10    XLIMIT = THRTY2
         C = C2
         YL = ZERO
         YU = ONE - ONE/DELTA
   15 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CHKCVG(BIGM, C, CRITER, DELTA, EPS, LFM, LOGFAC, M,
     *   MLOGG, Q, SIDE, T, TAU, XMAX, CVERG, IFAULT, INTOFG, YL, YU)
C
C        ALGORITHM AS 301.9 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Checks for the convergence to the limit C and updates the
C        end-points of the integration intervals
C
      INTEGER BIGM, T, M(T), SIDE, TAU, CVERG, IFAULT
      REAL C, CRITER, DELTA, EPS, MLOGG, Q(T), XMAX, INTOFG, YL, YU
      DOUBLE PRECISION LFM, LOGFAC(T) 
C
      DOUBLE PRECISION LOGDM
      REAL LDIRMT, LMB, LMB2,
     *     LMU1, LTERM, ONE, PID2, RHO, RLMB,
     *     THRTY2, U, V, ZERO
C
      COMMON /RHOPAR/ LMB, LMB2, LMU1, RLMB
C
      DATA ONE, PID2, THRTY2, ZERO /1.0, 1.5707963, 32.0, 0.0/
C
      RHO(U) = LMB/(LMB2 + (U - LMU1)**2)
      LDIRMT = SNGL(LOGDM(BIGM, LFM, LOGFAC, M, Q, T, TAU, XMAX + YU))
      IF (C .EQ. ZERO) THEN
         IF (LDIRMT + ALOG(RHO(XMAX + YU)) .LT. CRITER) CVERG = 1
      ELSE
         IF (ABS(LDIRMT - C) .LT. ALOG(ONE + EPS)) CVERG = 1
      ENDIF
      IF (CVERG .EQ. 0) THEN
         IF (ABS(XMAX + DELTA*YU) .GE. THRTY2) THEN
            CVERG = 1
            IF (TAU .EQ. 1) IFAULT = 5
         ELSE
            YL = YU
            YU = DELTA*YU
         ENDIF
      ENDIF
      IF (CVERG .EQ. 1) THEN
         IF (C .NE. ZERO) THEN
C
C        LTERM is assigned the value of the integral of DX
C
            V = ATAN((XMAX + YU - LMU1)*RLMB)
            LTERM = ZERO
            IF (SIDE .EQ. 1) THEN
               IF (V + PID2 .GT. ZERO) LTERM = V + PID2
            ELSE
               IF (PID2 - V .GT. ZERO) LTERM = PID2 - V
            ENDIF
            IF (LTERM .GT. ZERO) INTOFG = INTOFG +
     *                                EXP(C + ALOG(LTERM) - MLOGG)
         ENDIF
      ENDIF
C
      RETURN
      END
C
      REAL FUNCTION QG8(BIGM, LFM, LOGFAC, M, MLOGG, Q, T, TAU, XL, XU)
C
C        ALGORITHM AS 301.10 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Integrates from XL to XU by eight-point Gaussian quadrature
C
      INTEGER BIGM, T, M(T), TAU 
      DOUBLE PRECISION LFM, LOGFAC(T)
      REAL MLOGG, Q(T), XL, XU 
C
      DOUBLE PRECISION LOGDM
      REAL C(8), CONST, G, LMB, LMB2, LMU1, ONED2, RHO, RLMB, U, X,
     *     Y, Z, ZERO
      INTEGER I, J
C
      COMMON /RHOPAR/ LMB, LMB2, LMU1, RLMB
C
      DATA C /0.48014493, 0.050614268, 0.39833324, 0.11119052,
     *        0.26276620, 0.15685332, 0.091717321, 0.18134189/,
     *     CONST, ONED2, ZERO /-23.0, 0.5, 0.0/
C
      RHO(U) = LMB/(LMB2 + (U - LMU1)**2)
      Y = ONED2*(XL + XU)
      QG8 = ZERO
      DO 10 I = 2, 8, 2
         Z = C(I - 1)*(XU - XL)
         X = Y - Z
         DO 5 J = 1, 2
            G = SNGL(LOGDM(BIGM, LFM, LOGFAC, M, Q, T, TAU, X))
     *          + ALOG(RHO(X)) - MLOGG
            IF (G .GT. CONST) QG8 = QG8 + C(I)*EXP(G)
            X = Y + Z
    5    CONTINUE
   10 CONTINUE
      QG8 = (XU - XL)*QG8
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION LOGDM(BIGM, LFM, LOGFAC, M, Q, T,
     *                                TAU, X)
C
C        ALGORITHM AS 301.11 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes log F(X), log F(X)/(1 + M/EXP(X)) or
C        log F(X)/(M + EXP(X)) where M denotes BIGM.
C
      INTEGER BIGM, T, M(T), TAU, X
      DOUBLE PRECISION LFM, LOGFAC(T)
      REAL Q(T)
C
      DOUBLE PRECISION C, K, LGAPP, L2PID2, ONE, ONED2, TEN,
     *                 U, Z
      INTEGER ALPHA, I
C
      DATA C /0.08333333333333333D0/, L2PID2, ONE, ONED2, TEN
     * /0.9189385332046728D0, 1.0D0, 0.5D0, 10.0D0/
C
      LGAPP(U) = (U - ONED2)*DLOG(U) - U + L2PID2 + C/U
      K = DEXP(DBLE(X))
      LOGDM = LFM
      DO 15 ALPHA = 1, T
         IF (M(ALPHA) .GT. 0) THEN
            LOGDM = LOGDM - LOGFAC(ALPHA)
            Z = Q(ALPHA)*K
            I = 0
    5       IF (I + Z .GT. TEN .OR. I .EQ. M(ALPHA)) GO TO 10
               LOGDM = LOGDM + DLOG(I + Z)
               I = I + 1
               GO TO 5
   10       IF (I .LT. M(ALPHA)) LOGDM = LOGDM +
     *                           LGAPP(M(ALPHA) + Z) - LGAPP(I + Z)
         ENDIF
   15 CONTINUE
      I = 0
   20 IF (I + K .GT. TEN .OR. I .EQ. BIGM) GO TO 25
         LOGDM = LOGDM - DLOG(I + K)
         I = I + 1
         GO TO 20
   25 IF (I .LT. BIGM) LOGDM = LOGDM - LGAPP(BIGM + K) + LGAPP(I + K)
      IF (TAU .EQ. 2) THEN
         LOGDM = LOGDM - DLOG(ONE + BIGM/K)
      ELSE IF (TAU .EQ. 3) THEN
         LOGDM = LOGDM - DLOG(BIGM + K)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PARRHO(LQ, UQ, IFAULT, MU1, RPI)
C
C        ALGORITHM AS 301.12 APPL. STATIST. (1996), VOL.45, NO.2
C
C        Computes the parameter values of the log-Cauchy density
C
      REAL LQ, UQ, MU1, RPI
      INTEGER IFAULT
C
      REAL LMB, LMB2, LMU1, ONE, PI, RLMB, TENE9, TWO, ZERO
C
      COMMON /RHOPAR/ LMB, LMB2, LMU1, RLMB
C
      DATA ONE, PI, TENE9, TWO, ZERO
     *     /1.0, 3.1415927, 1.0E+09, 2.0, 0.0/
C
      IF (LQ .LE. ZERO .OR. LQ .GE. UQ .OR. UQ .GE. TENE9) THEN
         IFAULT = 1
      ELSE
         MU1 = SQRT(LQ*UQ)
         LMU1 = ALOG(MU1)
         LMB = ALOG(UQ/LQ)/TWO
         LMB2 = LMB**2
         RLMB = ONE/LMB
         RPI = ONE/PI
      ENDIF
C
      RETURN
      END
