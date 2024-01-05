C UKC NETLIB DISTRIBUTION COPYRIGHT 1990 RSS
C
        SUBROUTINE SINV(A, K, L, LO, PVT, IFAULT)
C
C         ALGORITHM AS 141 APPL. STATIST. (1979) VOL.28, NO.2
C
C         CALCULATE THE INVERSE OF A SYMMETRIC MATRIX
C         IGNORING A SPECIFIED ROW/COLUMN IF LO .NE. 0,
C         USING EITHER THE ORIGINAL MATRIX OR A COMPLETE INVERSE OF IT.
C
C  Correction mentioned in vol 28 no 3 1979 p 336 applied
C       
        DIMENSION A(L, L)
        REAL A, AA, AIP, BIG, EP, PVT, SMALL, T
        INTEGER P, PM, PP
C
	PARAMETER (BIG=1.0E38, SMALL=1.0E-7)
C
C         PARAMETER CHECKS
C
        IFAULT = 3
        IF (ABS(LO) .GT. K .OR. K .LT. 1 .OR. K .GT. L) RETURN
C
C         INITIAL VALUES
C
        IFAULT = 0
        IF (LO .GE. 0) GOTO 1
        EP = 1.0
        P = -LO
        PVT = ABS(A(P, P))
        T = PVT
        GOTO 3
1       EP = -1.0
        P = 1
        PVT = BIG
C
C         PIVOT BY PIVOT INVERSION
C
2       IF (P .EQ. LO) GOTO 12
        T = ABS(A(P, P))
        IF (T .LT. PVT) PVT = T
3       IF (T .LT. SMALL) IFAULT = 1
        IF (T .EQ. 0.0) GOTO 15
        PM = P - 1
        PP = P + 1
        AA = 1.0 / A(P, P)
        A(P, P) = -AA
        IF (P .EQ. 1) GOTO 8
        DO 7 I = 1, PM
          AIP = A(I, P) * AA
          DO 4 J = I, PM
4         A(I, J) = A(I, J) - AIP * A(J, P)
          IF (P .EQ. K) GOTO 6
          DO 5 J = PP, K
5         A(I, J) = A(I, J) - AIP * A(P, J)
6         A(I, P) = AIP * EP
7       CONTINUE
8       IF (P .EQ. K) GOTO 11
        DO 10 I = PP, K
          AIP = A(P, I) * AA
          DO 9 J = I, K
9         A(I, J) = A(I, J) - AIP * A(P, J)
          A(P, I) = AIP * EP
10      CONTINUE
11      IF (EP .GT. 0.0) RETURN
12      P = P + 1
        IF (P .LE. K) GOTO 2
C
C         SIGN CORRECTION
C
        DO 14 I = 1, K
          DO 13 J = I, K
13        A(I, J) = -A(I, J)
14      CONTINUE
        RETURN
C
C         NIL PIVOT EXIT
C
15      IFAULT = 2
        RETURN
        END
