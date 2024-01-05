      SUBROUTINE ALLRQR(R, NR, N, W1, NW1, W2, IW1, NW2, W3, NW3, W4,
     *                  IW2, NW4, IW3, NW5, IORD, IFIN, OUTR, IFAULT)
C
C        ALGORITHM AS 268.1  APPL.STATIST. (1991), VOL.40, NO.3
C
C        Calculates statistics for all possible subset regressions using
C        the R* matrix from a QR decomposition of (X|Y)
C
C	 Auxiliary routine required: AS 164 must be called first to form
C        the orthogonal reduction.
C
      INTEGER NR, N, NW1, NW2, IW1(NW2), NW3, IW2(NW4), NW4,
     *       IW3(NW5, 2), NW5, IORD, IFIN, IFAULT
      REAL R(NR), W1(NW1), W2(NW2, 2), W3(NW3), W4(NW4, 3)
      INTEGER I, J, J1, J2, J3, K, L, N1, N2, N3, N4, NB, NC, NF, NI,
     *        NL, NS1, NS2, NS3
      REAL TSS, GC, GS, E1, E2, E3, E4, E5, ONE
      EXTERNAL OUTR
      DATA ONE / 1.E+0 /
C
      IFAULT = 0
C
C        Checks on input arguments
C
      I = N * (N + 1) / 2
      IF (I .NE. NR) THEN
         IFAULT = 1
         RETURN
      END IF
      N1 = N - 1
      IF (N1 .LT. 2) THEN
         IFAULT = 2
         RETURN
      END IF
      N2 = I - N
      IF (N2 .GT. NW3) THEN
         IFAULT = 3
         RETURN
      END IF
      IF (N1 .GT. NW4) THEN
         IFAULT = 4
         RETURN
      END IF
      IF (N .LT. 5) THEN
         J1 = 1
         J2 = 1
         J3 = 1
      ELSE
         J1 = N - 4
         J2 = J1 * (3 + NR / 3)
         J3 = I - 2 * (N + 1)
      END IF
      IF (J1 .GT. NW5) THEN
         IFAULT = 5
         RETURN
      END IF
      IF (J2 .GT. NW1) THEN
         IFAULT = 6
         RETURN
      END IF
      IF (J3 .GT. NW2) THEN
         IFAULT = 7
         RETURN
      END IF
      IF (IORD .LT. 1 .OR. IORD .GT. 2) THEN
         IFAULT = 8
         RETURN
      END IF
      IF (IFIN .LT. 1 .OR. IFIN .GT. (N1 - 1)) THEN
         IFAULT = 9
         RETURN
      END IF
C
C        Calculate RSS's & TSS from original R* matrix
C
      DO 10 I = 1, N1
         IW2(I) = I
   10 CONTINUE
C*******************************************************************************
      TSS = ONE/R(NR)                                                   *  1A  *
C*******************************************************************************
      CALL OUTR(IW2, TSS, R, N2, R(N2 + 1), N1)
      J = N2
      DO 20 I = N1, 1, -1
C*******************************************************************************
          W4(I, 3) = R(N2 + I) * R(N2 + I) / R(J)                       *  2A  *
C*******************************************************************************
         TSS = TSS + W4(I, 3)
         J = J - I
         IF (I .GE. IFIN .AND. I .GT. 1) THEN
            CALL OUTR(IW2, TSS, R, J, R(N2 + 1), (I - 1))
         END IF
   20 CONTINUE
      DO 30 I = 2, N1
         W4(I, 3) = W4(I, 3) + W4(I - 1, 3)
   30 CONTINUE
C
C        Proceed dropping each column in turn
C
      IF (IORD .EQ. 1) THEN
         NB = 1
         NI = 1
         NF = N1 - IFIN
         NL = 0
      ELSE
         NB = N1 - IFIN
         NI = -1
         NF = 1
      END IF
   40 NC = NB
      IF (IORD .EQ. 2) THEN
         NL = NC
      END IF
C
C        Initialise temporary storage markers
C
      NS1 = 0
      NS2 = 0
      NS3 = 0
C
C        Start with original matrix.
C
      DO 50 I = 1, N2
         W3(I) = R(I)
   50 CONTINUE
      DO 60 I = 1, N1
         W4(I, 1) = R(N2 + I)
         W4(I, 2) = W4(I, 3)
         IW2(I) = I
   60 CONTINUE
C
C        Calculate givens rotation factors, perform rotations & fill
C        in elements of work matrix dropping the appropriate column.
C
      N4 = N - NC
   70 N3 = (N4 - 1) * (N4 - 2) / 2
   80 K = N3
      DO 110 L = 1, NC
         J = N4 + L - 2
         K = K + J
         J1 = N4 - NI * NC - NL - 2
C
C        Shift unchanged elements of R matrix
C
         IF (J1 .GE. 1) THEN
            J2 = K - J
            DO 90 I = 1, J1
               W3(J2 + I) = W3(K + I)
   90       CONTINUE
         END IF
C
C        Calculate first diagonal element
C
         J2 = K + J
         GC = W3(K)
         GS = W3(J2 + 1)
         IF (L .EQ. 1) THEN
            E1 = W3(J2)
C*******************************************************************************
            E5 = E1 * E1 * GS + GC                                      *  3A  *
C*******************************************************************************
         ELSE
C*******************************************************************************
            E5 = GS + GC                                                *  4A  *
C*******************************************************************************
         END IF
C*******************************************************************************
         W3(K) = GC * GS / E5                                           *      *
         GC = GC / E5                                                   *  5A  *
         GS = GS / E5                                                   *      *
C*******************************************************************************
         IF (L .EQ. 1) THEN
C*******************************************************************************
            GS = E1 * GS                                                *  6A  *
C*******************************************************************************
         END IF
C
C        Calculate second diagonal element & first off diagonal element
C
         IF (L .LT. NC) THEN
            J1 = J2 + J + 1
            E3 = W3(J1 + 1)
            IF (L .EQ. 1) THEN
               E2 = W3(J1)
               E4 = -E2 + E1 * E3
            ELSE
               E2 = W3(J2)
               E4 = -E2 + E3
            END IF
C*******************************************************************************
            W3(J2) = E2 * GS + E3 * GC                                  *  7A  *
            W3(J2 + 1) = E5 / (E4 * E4)                                 *      *
C*******************************************************************************
         END IF
C
C        Calculate  elements of rest of R matrix
C
         IF (L .LT. (NC - 1)) THEN
            J2 = K + 2 * J + 1
            DO 100 I = 1, NC - L - 1
               J1 = J2 + J + I + 1
               E3 = W3(J1 + 1)
               IF (L .EQ. 1) THEN
                  E2 = W3(J1)
                  W3(J2 + 1) = (-E2 + E1 * E3) / E4
               ELSE
                  E2 = W3(J2)
                  W3(J2 + 1) = (-E2 + E3) / E4
               END IF
C*******************************************************************************
               W3(J2) = E2 * GS + E3 * GC                               *  8A  *
C*******************************************************************************
               J2 = J1
  100       CONTINUE
         END IF
C
C        Apply rotation to vector of C's
C
         E2 = W4(J, 1)
         E3 = W4(J + 1, 1)
         IF (L .LT. NC) THEN
            IF (L .EQ. 1) THEN
               W4(J + 1, 1) = (-E2 + E1 * E3) / E4
            ELSE
               W4(J + 1, 1) = (-E2 + E3) / E4
            END IF
         END IF
C*******************************************************************************
         W4(J, 1) = E2 * GS + E3 * GC                                   *  9A  *
         W4(J, 2) = W4(J, 1) * W4(J, 1) / W3(K)                         *      *
C*******************************************************************************
         IF (J .GT. 1) W4(J, 2) = W4(J, 2) + W4(J - 1, 2)
  110 CONTINUE
C
C        Remove deleted variable from list of variables in model
C
      DO 120 L = N4 - 1, J
         IW2(L) = IW2(L + 1)
  120 CONTINUE
C
C        Print RSS's etc.
C
C*******************************************************************************
      CALL OUTR(IW2, (TSS - W4(J, 2)), W3, K, W4(1, 1), J)              *  1B  *
C*******************************************************************************
      IF (NC .GT. 1) THEN
C*******************************************************************************
C                                                                       *  2B  *
C*******************************************************************************
         J2 = K
         J1 = J + 1
         DO 130 L = 1, NC - 1
            I = J - L
            J2 = J2 - J1 + L
C*******************************************************************************
            CALL OUTR(IW2, (TSS - W4(I, 2)), W3, J2, W4(1, 1), I)       *  3B  *
C*******************************************************************************
  130    CONTINUE
      END IF
C
C        Storing  required R matrices in temporary storage
C
      IF (NC .LE. 2) GO TO 160
      DO 140 J1 = 1, K
         NS1 = NS1 + 1
         W1(NS1) = W3(J1)
  140 CONTINUE
      DO 150 J1 = 1, J
         NS2 = NS2 + 1
         W2(NS2, 1) = W4(J1, 1)
         W2(NS2, 2) = W4(J1, 2)
         IW1(NS2) = IW2(J1)
  150 CONTINUE
      NS3 = NS3 + 1
      IF (IORD .EQ. 1) THEN
         IW3(NS3, 1) = 2
         N4 = N1 - NS3
         IW3(NS3, 2) = NC
         NC = 2
      ELSE
         IW3(NS3, 1) = NC - 2
         IW3(NS3, 2) = 1
      END IF
  160 NC = NC - 1
      IF (IORD .EQ. 1 .AND. NC .EQ. 1) THEN
         GO TO 70
      END IF
      IF (IORD .EQ. 2 .AND. NC .NE. 0) THEN
         GO TO 80
      END IF
C
C        Working through R matrices in temporary storage
C
  170 IF (NS3 .EQ. 0) GO TO 200
      J = N1 - NS3
      J1 = NS2
      DO 180 L = J, 1, -1
         W4(L, 1) = W2(J1, 1)
         W4(L, 2) = W2(J1, 2)
         IW2(L) = IW1(J1)
         J1 = J1 - 1
  180 CONTINUE
      K = J * (J + 1) / 2
      J2 = NS1
      DO 190 L = K, 1, -1
         W3(L) = W1(J2)
         J2 = J2 - 1
  190 CONTINUE
      NC = IW3(NS3, 1)
      N4 = N - NS3 - NC
      IF (NC .EQ. IW3(NS3, 2)) THEN
         NS1 = J2
         NS2 = J1
         NS3 = NS3 - 1
         IF (IORD .EQ. 1) THEN
            GO TO 170
         END IF
      ELSE
         IW3(NS3, 1) = IW3(NS3, 1) + NI
      END IF
      GO TO 70
  200 IF (NB .NE. NF) THEN
         NB = NB + NI
         GO TO 40
      END IF
      RETURN
      END
      SUBROUTINE OUTR(IMV, RSS, R, NR, C, NC)
C
C        ALGORITHM AS 268.2  APPL.STATIST. (1991), VOL.40, NO.3
C
C        Version A: For all possible subsets regression
C
      INTEGER NC, NR, IMV(NC)
      REAL R(NR), C(NC), RSS
      WRITE (6, '(A,10(1X,I3))') ' MODEL VARIABLES ', IMV
      WRITE (6, '(A,F20.8)') ' RESIDUAL SUM OF SQUARES ', RSS
      WRITE (6, '(A, I3)') ' R MATRIX IN PACKED VECTOR FORM LENGTH ', NR
      WRITE (6, * ) R
      WRITE (6, '(A, I3)') ' C VECTOR LENGTH ', NC
      WRITE (6, * ) C
      RETURN
      END
      SUBROUTINE OUTR(IMV, RSS, R, NR, C, NC)
C
C        ALGORITHM AS 268.3  APPL.STATIST. (1991), VOL.40, NO.3
C
C        Version B:  Best subsets
C
      INTEGER NR, NC, IMV(NC)
      REAL R(NR), C(NC), RSS
      INTEGER IMVB
      REAL BRSS
      COMMON / BEST / BRSS(20), IMVB(210)
      INTEGER I, J
C
      J = NC * (NC - 1) / 2
      DO 10 I = 1, NC
         IMVB(J + I) = IMV(I)
   10 CONTINUE
      BRSS(NC) = RSS
      RETURN
      END
