      SUBROUTINE SCALE(FMN, FMX, N, VALMIN, STEP, VALMAX, IFAULT)
C
C     ALGORITHM AS 96  APPL. STATIST. (1976) VOL.25, NO.1
C
C     Given extreme values FMN, FMX, and the need for a scale with N
C     marks, calculates value for the lowest scale mark (VALMIN) and
C     step length (STEP) and highest scale mark (VALMAX).
C
      REAL FMN, FMX, VALMIN, STEP, VALMAX
      INTEGER N, IFAULT
C
C     Units for step lengths
C
      REAL UNIT(11)
C
C     Local variables
C
      INTEGER NUNIT, I, J
      REAL TOL, ZERO, HALF, ONE, TEN, BIAS, FMAX, FMIN, RN, X, S, RANGE
C
C     Array length unit()
C
      DATA NUNIT /11/
C
C     Local constant, defining effective equality of values.
C
      DATA TOL /5.0E-6/
      DATA ZERO /0.0/, HALF /0.5/, ONE /1.0/, TEN /10.0/
      DATA BIAS /1.0E-4/
      DATA UNIT /1.0, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,
     *        10.0/
C
      FMAX = FMX
      FMIN = FMN
      IFAULT = 1
C
C     Test for valid parameter values
C
      IF (FMAX .LT. FMIN .OR. N .LE. 1) RETURN
      IFAULT = 0
      RN = N - 1
      X = ABS(FMAX)
      IF (X .EQ. ZERO) X = ONE
      IF ((FMAX - FMIN) / X .GT. TOL) GO TO 20
C
C     All values effectively equal
C
      IF (FMAX .LT. ZERO) THEN
        FMAX = ZERO
      ELSE IF (FMAX .EQ. ZERO) THEN
        FMAX = ONE
      ELSE
        FMIN = ZERO
      END IF
C
   20 STEP = (FMAX - FMIN) / RN
      S = STEP
C
C     Find power of 10
C
   25 IF (S .GE. ONE) GO TO 30
      S = S * TEN
      GO TO 25
   30 IF (S .LT. TEN) GO TO 35
      S = S / TEN
      GO TO 30
C
C     Calculate STEP
C
   35 X = S - BIAS
      DO 40 I = 1, NUNIT
        IF (X .LE. UNIT(I)) GO TO 45
   40 CONTINUE
   45 STEP = STEP * UNIT(I) / S
      RANGE = STEP * RN
C
C     Make first estimate of VALMIN
C
      X = HALF * (ONE + (FMIN + FMAX - RANGE) / STEP)
      J = X - BIAS
      IF (X .LT. ZERO) J = J - 1
      VALMIN = STEP * FLOAT(J)
C
C     Test if VALMIN could be zero
C
      IF (FMIN .GE. ZERO .AND. RANGE .GE. FMAX) VALMIN = ZERO
      VALMAX = VALMIN + RANGE
C
C     Test if VALMAX could be zero
C
      IF (FMAX .GT. ZERO .OR. RANGE .LT. -FMIN) RETURN
      VALMAX = ZERO
      VALMIN = -RANGE
      RETURN
      END

