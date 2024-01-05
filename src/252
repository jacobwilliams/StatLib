      SUBROUTINE GESAT(NVAR, NP, ISET, NUM, SATMOD, IFAULT)
C
C        ALGORITHM AS 252.1  APPL.STATIST. (1990), VOL.39, NO.1
C
C        Obtains the description of all interactions
C        in the saturated model
C
C     Auxiliary routine required: COMBO from AS 160 (included here)
C
      INTEGER NVAR, NP, IFAULT, ISET(NVAR), NUM(NVAR), SATMOD(NVAR, NP)
      INTEGER NVMAX, I, J, K, M, N
      LOGICAL LAST
      DATA NVMAX / 15 /
C
      IFAULT = 1
      IF (NVAR .LT. 1 .OR. NVAR .GT. NVMAX) RETURN
      IFAULT = 2
      I = 2 ** NVAR - 1
      IF (NP .LT. I) RETURN
      IFAULT = 0
      DO 20 I = 1, NVAR
         DO 10 J = 1, NP
            SATMOD(I, J) = 0
   10    CONTINUE
   20 CONTINUE
      M = 0
      N = 0
      DO 50 K = 1, NVAR
         I = NVAR - K + 1
         LAST = .TRUE.
   30    CALL COMBO(ISET, NVAR, I, LAST)
         IF ( .NOT. LAST) THEN
            N = N + 1
            DO 40 J = 1, I
               SATMOD(J, N) = ISET(J)
   40       CONTINUE
            GO TO 30
         END IF
         M = M + 1
         NUM(M) = N
   50 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CLAGEN(NVAR, NP, NUM, SATMOD, IPAR, DES, NG, INTG,
     *                  IGEN, ISET, JSET, IDES, JDES, NCON, CONFIG,
     *                  IFAULT)
C
C        ALGORITHM AS 252.2  APPL.STATIST. (1990), VOL.39, NO.1
C
C        Obtains the generating class or dual generating class
C        of the model defined in terms of interactions which
C        must be present/absent or which must be added/omitted
C        to/from given model
C
      INTEGER NVAR, NP, IPAR, NG, IGEN, NCON, IFAULT, NUM(NVAR),
     *        SATMOD(NVAR, NP), ISET(NVAR), JSET(NVAR), INTG(NVAR, NG),
     *        CONFIG(NVAR, NP)
      LOGICAL DES(NP), IDES(NP), JDES(NP)
      INTEGER I, IP, J, J1, J2, K, KZ, LF, M, M1, M2, M3, M4, N, NI,
     *        NR, NR1, NR2, NST, NVMAX
      LOGICAL IND, LG, LW, STER
      EXTERNAL INSET
      LOGICAL INSET
C
C        NVMAX:  Set to 15 for 16 bit integers - see text.
C
      DATA NVMAX / 15 /
C
C        Check input parameters
C
      IFAULT = 1
      IF (NVAR .LT. 1 .OR. NVAR .GT. NVMAX) RETURN
      IFAULT = 2
      I = 2 ** NVAR - 1
      IF (NP .LT. I) RETURN
      IFAULT = 3
      IF (IPAR .LT. 1 .OR. IPAR .GT. 4) RETURN
      STER = IPAR .EQ. 1 .OR. IPAR .EQ. 3
      IFAULT = 4
      IF (NG .LT. 1) RETURN
      IFAULT = 5
      DO 20 I = 1, NG
         DO 10 J = 1, NVAR
            K = INTG(J, I)
            IF (K .LT. 0 .OR. K .GT. NVAR) RETURN
   10    CONTINUE
   20 CONTINUE
      IFAULT = 6
      IF (IGEN .NE. 1 .AND. IGEN .NE. 2) RETURN
      IFAULT = 0
C
C        Initialize workspaces and variables
C
      LG = IPAR .EQ. 3 .OR. IPAR .EQ. 4
      IF ( .NOT. LG) THEN
         DO 30 I = 1, NP
            DES(I) = .NOT. STER
   30    CONTINUE
      END IF
      DO 40 I = 1, NP
         LW = DES(I)
         IDES(I) = LW
         JDES(I) = LW
   40 CONTINUE
      IF (STER) THEN
         K = -1
         KZ = 1
      ELSE
         K = 1
         KZ = NVAR
      END IF
C
C        Find full description of the model
C
      DO 140 NI = 1, NG
         DO 50 I = 1, NVAR
            JSET(I) = 0
   50    CONTINUE
         NR = 0
         DO 60 I = 1, NVAR
            J = INTG(I, NI)
            IF (J .EQ. 0) GO TO 70
            NR = NR + 1
            JSET(NR) = J
   60    CONTINUE
   70    N = NR
         IF (N .EQ. 0) THEN
            IF (STER) GO TO 140
            DO 80 I = 1, NP
               IDES(I) = .FALSE.
   80       CONTINUE
            GO TO 150
         END IF
   90    IF (N .EQ. NVAR) THEN
            M3 = 1
            M4 = 1
         ELSE
            J = NVAR - N
            M3 = NUM(J) + 1
            M4 = NUM(J + 1)
         END IF
         DO 120 I = M3, M4
            IDES(I) = STER
            IF (JDES(I) .NEQV. STER) THEN
               IP = 0
               DO 100 J = 1, NVAR
                  IF (SATMOD(J, I) .EQ. 0) GO TO 110
                  IP = IP + 1
                  ISET(IP) = SATMOD(J, I)
  100          CONTINUE
  110          IF (STER) THEN
                  IND = INSET(NR, JSET, IP, ISET)
               ELSE
                  IND = INSET(IP, ISET, NR, JSET)
               END IF
               IDES(I) = JDES(I)
               IF (IND) IDES(I) = STER
            END IF
  120    CONTINUE
         N = N + K
         IF (STER) THEN
            IF (N .GE. KZ) GO TO 90
         ELSE
            IF (N .LE. KZ) GO TO 90
         END IF
         DO 130 J = 1, NP
            JDES(J) = IDES(J)
  130    CONTINUE
  140 CONTINUE
  150 IF (LG) GO TO 170
C
C        If the model is defined by interactions which must be
C        present/absent, store its full description in DES
C
      DO 160 I = 1, NP
         DES(I) = IDES(I)
  160 CONTINUE
C
C        Find the (dual) generating class if the model is saturated
C
  170 STER = IGEN .EQ. 1
      NCON = 1
      LW = IDES(1)
      LG = LW .EQV. STER
      DO 180 I = 1, NVAR
         K = 0
         IF (LG) K = I
         CONFIG(I, 1) = K
  180 CONTINUE
      IF (LW) RETURN
C
C        Find the (dual) generating class if the model is unsaturated
C
      LW = .TRUE.
      NCON = 0
      NR1 = NVAR - 1
      IF (STER) THEN
         NR2 = NR1
         NR1 = 1
         NST = 1
      ELSE
         NR2 = 1
         NST = -1
      END IF
      DO 240 M = NR1, NR2, NST
         I = NVAR - M
         LF = 0
         M1 = NUM(M) + 1
         M2 = NUM(M + 1)
         DO 230 J = M1, M2
            IF (STER .EQV. IDES(J)) THEN
               IF ( .NOT. LW) THEN
                  DO 190 J1 = 1, I
                     ISET(J1) = SATMOD(J1, J)
  190             CONTINUE
                  IF (I .NE. NR2) THEN
                     M3 = NUM(M - NST) + 1
                     M4 = NUM(M - NST + 1)
                     DO 210 J1 = M3, M4
                        IF (IDES(J1) .EQV. STER) THEN
                           IP = I + NST
                           DO 200 J2 = 1, IP
                              JSET(J2) = SATMOD(J2, J1)
  200                      CONTINUE
                           IF ( .NOT. STER) THEN
                              LG = INSET(I, ISET, IP, JSET)
                           ELSE
                              LG = INSET(IP, JSET, I, ISET)
                           END IF
                           IF (LG) GO TO 230
                        END IF
  210                CONTINUE
                  END IF
               END IF
               NCON = NCON + 1
               DO 220 J1 = 1, NVAR
                  CONFIG(J1, NCON) = SATMOD(J1, J)
  220          CONTINUE
            ELSE
               LF = LF + 1
            END IF
  230    CONTINUE
         IF (LF .EQ. 0) RETURN
         J = M2 - M1 + 1
         LW = LF .EQ. J
  240 CONTINUE
      IF (NCON .EQ. 0) NCON = 1
      RETURN
      END
C
C
      LOGICAL FUNCTION INSET(N, ISET, M, JSET)
C
C        ALGORITHM AS 252.3  APPL.STATIST. (1990), VOL.39, NO.1
C
      INTEGER N, M, ISET(N), JSET(M)
      INTEGER I, J
C
      INSET = .FALSE.
      DO 20 I = 1, M
         DO 10 J = 1, N
            IF (JSET(I) .EQ. ISET(J)) GO TO 20
   10    CONTINUE
         RETURN
   20 CONTINUE
      INSET = .TRUE.
      RETURN
      END
C
C
      SUBROUTINE COMBO(ISET, N, M, LAST)
C
C     ALGORITHM AS 160.2  APPL. STATIST. (1981) VOL.30, NO.1
C
C     Subroutine to generate all possible combinations of M of the
C     integers from 1 to N in a stepwise fashion.   Prior to the first
C     call, LAST should be set to .FALSE.   Thereafter, as long as LAST
C     is returned .FALSE., a new valid combination has been generated.
C     When LAST goes .TRUE., there are no more combinations.
C
      LOGICAL LAST
      INTEGER N, M, ISET(M)
C
      IF (LAST) GO TO 110
C
C     Get next element to increment
C
      K = M
100   L = ISET(K) + 1
      IF (L + M - K .LE. N) GO TO 150
      K = K - 1
C
C     See if we are done
C
      IF (K .LE. 0) GO TO 130
      GO TO 100
C
C     Initialize first combination
C
110   DO 120 I = 1, M
120   ISET(I) = I
130   LAST = .NOT. LAST
      RETURN
C
C     Fill in remainder of combination.
C
150   DO 160 I = K, M
	ISET(I) = L
	L = L + 1
160   CONTINUE
C
      RETURN
      END
