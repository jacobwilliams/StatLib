      SUBROUTINE INTBD(G, N, E, ET, DELTA, ITAIL)
C
C        ALGORITHM AS R96 APPL. STATIST. (1996) VOL.45, NO.2
C
C        Interval for trapezoidal rule integration to limit error to
C        E using moment generating function bound for weighted sum of
C        chi-squared(1) random variables
C
      INTEGER N, ITAIL
      REAL DELTA, G(N), E, ET
C
      INTEGER I
      REAL ZERO, TWO, POINT1, POINT5, TWOPI, ULIM, PSI, DPSI, D2PSI,
     *     HB, U, U1, U2, ELOG, EPS, HM, HSD, SGN, ZABS, ZEXP, ZFLOAT,
     *     ZLOG, ZMAX1, ZMIN1, ZSQRT
C
      EXTERNAL LMGFEV
C
      DATA ZERO, TWO, POINT1, POINT5 /0.0E0, 2.0E0, 0.1E0, 0.5E0/
      DATA TWOPI /6.28318530717959E0/
      DATA EPS /0.09531E0/
C
      ZABS(U) = ABS(U)
      ZEXP(U) = EXP(U)
      ZFLOAT(I) = FLOAT(I)
      ZLOG(U) = ALOG(U)
      ZMAX1(U, U1) = AMAX1(U, U1)
      ZMIN1(U, U1) = AMIN1(U, U1)
      ZSQRT(U) = SQRT(U)
C
C        Tail check
C
      SGN = ZFLOAT(ITAIL)
      IF (ITAIL .EQ. 1) THEN
         IF (G(N) .LE. ZERO) RETURN
         ULIM = POINT5 / G(N)
      ELSE
         IF (G(1) .GE. ZERO) RETURN
         ULIM = POINT5 / G(1)
      ENDIF
      HM = ZERO
      HSD = ZERO
      DO 10 I = 1, N
         HM = HM + G(I)
         HSD = HSD + G(I) ** 2
   10 CONTINUE
C
C        Newton search for zero log(m.g.f.) derivative
C
      HSD = ZSQRT(HSD)
      U1 = ZERO
      U2 = ZERO
      IF (SGN * HM .LT. ZERO .AND. -SGN * HM / HSD .GT. TWO) THEN
         DPSI = HM
   20    IF (DPSI * SGN .LT. ZERO) THEN
            U1 = U2
            U2 = POINT5 * (U2 + ULIM)
            CALL LMGFEV (U2, G, N, PSI, DPSI, D2PSI)
            GOTO 20
         ENDIF
   25    IF (DPSI ** 2 .GE. POINT1 * D2PSI) THEN
            U2 = U2 - DPSI / D2PSI
            CALL LMGFEV(U2, G, N, PSI, DPSI, D2PSI)
         ENDIF
         IF (ZEXP(PSI) .LE. E) RETURN
      ENDIF
C
C        Bound equality
C
      ELOG = - ZLOG(ET)
      CALL LMGFEV(U2, G, N, PSI, DPSI, D2PSI)
      HB = PSI - U2 * DPSI + ELOG
   30 IF (HB .GT. ZERO) THEN
         U1 = U2
         U2 = POINT5 * (U2 + ULIM)
         CALL LMGFEV(U2, G, N, PSI, DPSI, D2PSI)
         HB = PSI - U2 * DPSI + ELOG
         GOTO 30
      ENDIF
C
C        Newton/binary search solution
C
      U = U2
   40 IF (ZABS(HB) .GE. EPS) THEN
         U = U + HB / (U * D2PSI)
         IF (U .LE. ZMIN1(U1, U2) .OR. U .GT. ZMAX1(U1, U2))
     *      U = POINT5 * (U1 + U2)
         CALL LMGFEV(U, G, N, PSI, DPSI, D2PSI)
         HB = PSI - U * DPSI + ELOG
         IF (HB .GT. ZERO) THEN
            U1 = U
         ELSE
            U2 = U
         ENDIF
         GOTO 40
      ENDIF
      IF (ITAIL .EQ. 1 .AND. DPSI .GT. ZERO) THEN
         DELTA = TWOPI / DPSI
         ITAIL = 0
      ELSEIF (ITAIL .EQ. -1 .AND. DPSI .LT. ZERO) THEN
         DELTA = - TWOPI / DPSI
         ITAIL = 0
      ENDIF
C
      RETURN
      END
