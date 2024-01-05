      REAL FUNCTION STUDNT (T, DOFF, IFAULT)
C
C     ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
C
C     Calculate the upper tail area under Student's t-distribution
C
C     Translated from Algol by Alan Miller
C
      INTEGER IFAULT
      REAL T, DOFF
C
C     Local variables
C
      REAL V, X, TT, TWO, FOUR, ONE, ZERO, HALF
      REAL A1, A2, A3, A4, A5, B1, B2,
     *     C1, C2, C3, C4, C5, D1, D2,
     *     E1, E2, E3, E4, E5, F1, F2,
     *     G1, G2, G3, G4, G5, H1, H2,
     *     I1, I2, I3, I4, I5, J1, J2
      LOGICAL POS
      DATA TWO /2.0/, FOUR /4.0/, ONE /1.0/, ZERO /0.0/, HALF /0.5/
      DATA A1, A2, A3, A4, A5 /0.09979441, -0.581821, 1.390993,
     *     -1.222452, 2.151185/, B1, B2 /5.537409, 11.42343/
      DATA C1, C2, C3, C4, C5 /0.04431742, -0.2206018, -0.03317253,
     *     5.679969, -12.96519/, D1, D2 /5.166733, 13.49862/
      DATA E1, E2, E3, E4, E5 /0.009694901, -0.1408854, 1.88993,
     *     -12.75532, 25.77532/, F1, F2 /4.233736, 14.3963/
      DATA G1, G2, G3, G4, G5 /-9.187228E-5, 0.03789901, -1.280346,
     *     9.249528, -19.08115/, H1, H2 /2.777816, 16.46132/
      DATA I1, I2, I3, I4, I5 /5.79602E-4, -0.02763334, 0.4517029,
     *     -2.657697, 5.127212/, J1, J2 /0.5657187, 21.83269/
C
C     Check that number of degrees of freedom > 4.
C
      IF (DOFF .LT. TWO) THEN
	IFAULT = 1
	STUDNT = - ONE
	RETURN
      END IF
C
      IF (DOFF .LE. FOUR) THEN
	IFAULT = DOFF
      ELSE
	IFAULT = 0
      END IF
C
C     Evaluate series.
C
      V = ONE / DOFF
      POS = (T .GE. ZERO)
      TT = ABS(T)
      X = HALF * (ONE +
     *    TT * (((A1 + V * (A2 + V * (A3 + V * (A4 + V * A5)))) /
     *        (ONE - V * (B1 - V * B2))) +
     *    TT * (((C1 + V * (C2 + V * (C3 + V * (C4 + V * C5)))) /
     *        (ONE - V * (D1 - V * D2))) +
     *    TT * (((E1 + V * (E2 + V * (E3 + V * (E4 + V * E5)))) /
     *        (ONE - V * (F1 - V * F2))) +
     *    TT * (((G1 + V * (G2 + V * (G3 + V * (G4 + V * G5)))) /
     *        (ONE - V * (H1 - V * H2))) +
     *    TT * ((I1 + V * (I2 + V * (I3 + V * (I4 + V * I5)))) /
     *        (ONE - V * (J1 - V * J2))) ))))) ** (-8)
      IF (POS) THEN
	STUDNT = X
      ELSE
	STUDNT = ONE - X
      END IF
C
      RETURN
      END
