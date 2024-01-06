CSTART OF AS 5
      REAL FUNCTION PRNCST(ST, IDF, D, IFAULT)
C
C        ALGORITHM AS 5  APPL. STATIST. (1968) VOL.17, P.193
C
C        COMPUTES LOWER TAIL AREA OF NON-CENTRAL T-DISTRIBUTION
C
      REAL ST, D, G1, G2, G3, ZERO, ONE, TWO, HALF, EPS, EMIN, F,
     $  A, B, RB, DA, DRB, FMKM1, FMKM2, SUM, AK, FK, FKM1,
     $  ALNORM, TFN, ALOGAM, ZSQRT, ZEXP
C
C        CONSTANTS - G1 IS 1.0 / SQRT(2.0 * PI)
C                    G2 IS 1.0 / (2.0 * PI)
C                    G3 IS SQRT(2.0 * PI)
C
      DATA G1, G2, G3 /0.3989422804, 0.1591549431, 2.5066282746/
      DATA ZERO, ONE, TWO, HALF,    EPS, EMIN
     $     /0.0, 1.0, 2.0,  0.5, 1.0E-6, 12.5/
C
      ZSQRT(A) = SQRT(A)
      ZEXP(A) = EXP(A)
C
      F = IDF
      IF (IDF .GT. 100) GOTO 50
      IFAULT = 0
      IOE = MOD(IDF, 2)
      A = ST / ZSQRT(F)
      B = F / (F + ST ** 2)
      RB = ZSQRT(B)
      DA = D * A
      DRB = D * RB
      SUM = ZERO
      IF (IDF .EQ. 1) GOTO 30
      FMKM2 = ZERO
      IF (ABS(DRB) .LT. EMIN) FMKM2 = A * RB * ZEXP(-HALF * DRB ** 2)
     $  * ALNORM(A * DRB, .FALSE.) * G1
      FMKM1 = B * DA * FMKM2
      IF (ABS(D) .LT. EMIN)
     $  FMKM1 = FMKM1 + B * A * G2 * ZEXP(-HALF * D ** 2)
      IF (IOE .EQ. 0) SUM = FMKM2
      IF (IOE .EQ. 1) SUM = FMKM1
      IF (IDF .LT. 4) GOTO 20
      IFM2 = IDF - 2
      AK = ONE
      FK = TWO
      DO 10 K = 2, IFM2, 2
      FKM1 = FK - ONE
      FMKM2 = B * (DA * AK * FMKM1 + FMKM2) * FKM1 / FK
      AK = ONE / (AK * FKM1)
      FMKM1 = B * (DA * AK * FMKM2 + FMKM1) * FK / (FK + ONE)
      IF (IOE .EQ. 0) SUM = SUM + FMKM2
      IF (IOE .EQ. 1) SUM = SUM + FMKM1
      AK = ONE / (AK * FK)
      FK = FK + TWO
   10 CONTINUE
   20 IF (IOE .EQ. 0) GOTO 40
   30 PRNCST = ALNORM(DRB, .TRUE.) + TWO * (SUM + TFN(DRB, A))
      RETURN
   40 PRNCST = ALNORM(D, .TRUE.) + SUM * G3
      RETURN
C
C        NORMAL APPROXIMATION - K IS NOT TESTED AFTER THE TWO CALLS
C        OF ALOGAM, BECAUSE A FAULT IS IMPOSSIBLE WHEN F EXCEEDS 100
C
   50 IFAULT = 1
      A = ZSQRT(HALF * F) * ZEXP(ALOGAM(HALF * (F - ONE), K)
     $  - ALOGAM(HALF * F, K)) * D
      PRNCST = ALNORM((ST - A) / ZSQRT(F * (ONE + D ** 2)
     $  / (F - TWO) - A ** 2), .FALSE.)
      RETURN
      END
CEND OF AS 5
