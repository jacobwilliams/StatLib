      SUBROUTINE LSINF(IDIST, ITYPE, ZL, ZR, F11, F12, F22, IFAULT)
C
C        ALGORITHM AS 292.1 APPL.STATIST. (1994), VOL.43, NO.3
C
C        Computes Fisher information matrix elements for time (type I)
C        or failure (type II) censored units from the smallest extreme
C        value (sev), largest extreme value (lev), normal or logistic
C        distribution
C
      INTEGER IDIST, ITYPE, IFAULT
      REAL ZL, ZR, F11, F12, F22
C
      INTEGER NDIST
C
      PARAMETER (NDIST = 4)
C
      REAL BIG, BOUNDS(2, NDIST), CDF, DEN, ETA, PDF, SUR, THET0L,
     *     THET1L, THET2L, THET0R, THET1R, THET2R, ZERO, ZMAX, ZMIN,
     *     ZU
C
      DATA BOUNDS / 3.6E + 00, -34.0E + 00,
     *             34.0E + 00,  -3.6E + 00,
     *              8.0E + 00,  -8.0E + 00,
     *             34.0E + 00, -34.0E + 00 /
C
      PARAMETER (BIG = 1.0E + 10, ZERO = 0.0E + 00)
C
      EXTERNAL LSINT
C
C        Check for illegal idist
C
      IFAULT = 1
      IF (IDIST .LT. 1 .OR. IDIST .GT. NDIST) RETURN
C
C        Check for illegal itype
C
      IFAULT = 2
      IF (ITYPE .LT. 1 .OR. ITYPE .GT. 4) RETURN
C
C        Check for illegal arguments
C
      IFAULT = 3
      IF (ITYPE .EQ. 4 .AND. ZL .GT. ZR) RETURN
      IFAULT = 0
      ZMAX = BOUNDS(1, IDIST)
      ZMIN = BOUNDS(2, IDIST)
C
C        Identify the censoring type
C
      GOTO (10, 20, 30, 40), ITYPE
C
C        Uncensored case
C
   10 CALL LSINT(IDIST, ZMIN, ZMAX, BIG, F11, F12, F22, IFAULT)
      RETURN
C
C        Right censored at zr
C
   20 CALL LSINT(IDIST, ZMIN, ZMAX, ZR, F11, F12, F22, IFAULT)
      RETURN
C
C        Left censored at zl
C
   30 ZU = BIG
C
C        Left and right censored data
C
   40 IF (ITYPE .EQ. 4) ZU = ZR
      CALL LSINT(IDIST, ZMIN, ZMAX, ZU, THET0R, THET1R, THET2R, IFAULT)
      IF (IFAULT .EQ. 4 .OR. IFAULT .EQ. 5) RETURN
      CALL LSINT(IDIST, ZMIN, ZMAX, ZL, THET0L, THET1L, THET2L, IFAULT)
      IF (IFAULT .EQ. 4 .OR. IFAULT .EQ. 5) RETURN
C
      ETA = ZERO
      IF (ZL .GT. ZMIN .AND. ZL .LT. ZMAX) THEN
         CALL PCSFUN(IDIST, ZL, PDF, CDF, SUR)
         DEN = CDF * SUR
         IF (DEN .GT. ZERO) ETA = PDF * PDF / DEN
      ENDIF
C
C        Left or left and right censored data - Fisher matrix elements
C
      F11 = THET0R - THET0L + ETA
      F12 = THET1R - THET1L + ZL * ETA
      F22 = THET2R - THET2L + ZL * ZL * ETA
      RETURN
      END
C
      SUBROUTINE LSINT(IDIST, ZMIN, ZMAX, Z, THETA0, THETA1, THETA2,
     *                 IFAULT)
C
C        ALGORITHM AS 292.2 APPL.STATIST. (1994), VOL.43, NO.3
C
C        Computes theta0, theta1, theta2 for the sev, lev, normal and
C        logistic distributions
C
      INTEGER IDIST, IFAULT
      REAL ZMIN, ZMAX, Z, THETA0, THETA1, THETA2
C
      INTEGER NDIST
      PARAMETER (NDIST = 4)
C
      INTEGER IDISEV, J, JMAXLO, JTERM, N
      REAL ASYM(3, NDIST), ANUZ, CDF, CONST0, DEN, DUM, EMABSZ, ETASEV,
     *     HALF, ONE, PDF, SUR, S3, S3OLD, THREE, TMP0, TMP1, TOL, TWO,
     *     T0SEV, T1SEV, T2SEV, WKSP(40), X1, X2, X3, X4, ZERO, ZSEV
C
      DATA ASYM /
     * 1.0E + 00, 0.4227843350984671E + 00, 0.1823680660852879E + 01,
     * 1.0E + 00,-0.4227843350984671E + 00, 0.1823680660852879E + 01,
     * 1.0E + 00, 0.0E + 00, 2.0E + 00,
     * 0.3333333333333333E + 00, 0.0E + 00, 0.14299560445654842E + 01 /
C
      PARAMETER (CONST0 = 0.1644934066848226E + 01, HALF = 0.5E + 00,
     *           IDISEV = 1, JMAXLO = 40, ONE = 1.0E + 00,
     *           THREE = 3.0E + 00, TOL = 1.0E - 11, TWO = 2.0E + 00,
     *           ZERO = 0.0E + 00)
C
      EXTERNAL INTEGR, PCSFUN
C
      IFAULT = 0
C
C        Routing computations
C
      IF (Z .LE. ZMIN) THEN
C
C        Integrals are negligible
C
         THETA0 = ZERO
         THETA1 = ZERO
         THETA2 = ZERO
      ELSE IF (Z .GE. ZMAX) THEN
C
C        Integrals take limiting values
C
         THETA0 = ASYM(1, IDIST)
         THETA1 = ASYM(2, IDIST)
         THETA2 = ASYM(3, IDIST)
      ELSE
C
C        Computed integrals - z is in interval (zmin, zmax)
C        First select the distribution
C
         GOTO (10, 20, 30, 40), IDIST
C
C        Smallest extreme value (sev) distribution
C
   10    CALL INTEGR(Z, THETA0, THETA1, THETA2, IFAULT)
         IF (IFAULT .NE. 0) IFAULT = 4
         RETURN
C
C        Largest extreme value (lev) distribution
C        First compute thetas and eta for the sev at -z
C
   20    ZSEV = -Z
         CALL INTEGR(ZSEV, T0SEV, T1SEV, T2SEV, IFAULT)
         IF (IFAULT .NE. 0) THEN
            IFAULT = 4
            RETURN
         ENDIF
         ETASEV = ZERO
         CALL PCSFUN(IDISEV, ZSEV, PDF, CDF, SUR)
         DEN = CDF * SUR
         IF (DEN .GT. ZERO) ETASEV = PDF * PDF / DEN
C
C        Compute thetas for the lev
C
         THETA0 = ONE - T0SEV + ETASEV
         THETA1 = ASYM(2, IDIST) + T1SEV + Z * ETASEV
         THETA2 = ASYM(3, IDIST) - T2SEV + Z * Z * ETASEV
         RETURN
C
C        Normal distribution
C
   30    CALL PCSFUN(IDIST, Z, PDF, CDF, SUR)
         THETA0 = CDF - Z * PDF
         IF (SUR .GT. ZERO) THETA0 = THETA0 + PDF * PDF / SUR
         THETA1 = Z * THETA0 - Z * CDF - PDF
         THETA2 = Z * THETA1 + TWO * CDF
         RETURN
C
C        Logistic distribution
C
   40    CALL PCSFUN(IDIST, Z, PDF, CDF, SUR)
         THETA0 = (ONE - SUR ** 3) / THREE
         TMP0 = LOG(ONE + EXP(Z))
         THETA1 = Z * THETA0 + (PDF - TMP0) / THREE
C
C        Compute theta2 for the logistic
C        S3 is computed using a power series expansion with
C        accelerated convergence using Euler transformation
C
         EMABSZ = EXP(-ABS(Z))
         X1 = EMABSZ
         S3 = ZERO
C
C        Euler partial sum based on the first term
C
         N = 1
         WKSP(1) = X1
         S3 = HALF * X1
         S3OLD = S3
         DO 60 JTERM = 2, JMAXLO
            X2 = REAL(JTERM) ** 2
            X3 = REAL(JTERM - 1) ** 2
            X1 = -X1 * EMABSZ * X3 / X2
C
C        Euler partial sum based on two or more terms
C
            TMP1 = WKSP(1)
            WKSP(1) = X1
            DO 50 J = 1, N
               IF (J .LT. N) DUM = WKSP(J + 1)
               WKSP(J + 1) = HALF * (WKSP(J) + TMP1)
               IF (J .LT. N) TMP1 = DUM
   50       CONTINUE
C
C        Euler improved partial sums
C
            IF (ABS(WKSP(N + 1)) .GT. ABS(WKSP(N))) THEN
               S3 = S3 + WKSP(N + 1)
            ELSE
               S3 = S3 + HALF * WKSP(N + 1)
               N = N + 1
            ENDIF
C
C        Tests for convergence of the series - a fault is declared if
C        convergence is not reached in a maximum of jmaxlo terms
C
            X4 = ABS(S3OLD - S3)
            S3OLD = S3
            IF (X4 .LT. TOL) THEN
C
C        Add terms to obtain the integral
C
               ANUZ = S3
               IF (Z .GT. ZERO) ANUZ = CONST0 + Z * Z / TWO - S3
               THETA2 = Z * THETA1
     *               + (CDF - Z * TMP0 + Z * PDF + TWO * ANUZ) / THREE
               RETURN
            ENDIF
   60    CONTINUE
         IFAULT = 5
      ENDIF
      RETURN
      END
C
      SUBROUTINE PCSFUN(IDIST, Z, PDF, CDF, SUR)
C
C        ALGORITHM AS 292.3 APPL.STATIST. (1994), VOL.43, NO.3
C
C        Computes the p.d.f., c.d.f. and survival functions for the
C        sev, lev, normal and logistic distributions
C
      INTEGER IDIST
      REAL Z, PDF, CDF, SUR
C
      REAL ALNORM, CVAL, ENZ, EPZ, HALF, ONE
C
      PARAMETER (CVAL = 0.3989422804014326E + 00, HALF = 0.5E + 00,
     *           ONE = 1.0E + 00)
C
      EXTERNAL ALNORM
C
C        Select the distribution
C
      GOTO (10, 20, 30, 40), IDIST
C
C        sev
C
   10 EPZ = EXP(Z)
      PDF = EXP(Z - EPZ)
      CDF = ONE - EXP(-EPZ)
      SUR = EXP(-EPZ)
      RETURN
C
C        lev
C
   20 ENZ = EXP(-Z)
      PDF = EXP(-Z - ENZ)
      CDF = EXP(-ENZ)
      SUR = ONE - EXP(-ENZ)
      RETURN
C
C        Normal
C
   30 PDF = CVAL * EXP(-HALF * Z * Z)
      CDF = ALNORM(Z, .FALSE.)
      SUR = ALNORM(Z, .TRUE.)
      RETURN
C
C        Logistic
C
   40 ENZ = EXP(-Z)
      CDF = ONE / (ONE + ENZ)
      SUR = ENZ / (ONE + ENZ)
      PDF = CDF * SUR
      RETURN
C
      END
