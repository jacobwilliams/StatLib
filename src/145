C UKC NETLIB DISTRIBUTION COPYRIGHT 1990 RSS
C
        SUBROUTINE MULMAX(N, K, NMIN, NN, NNK, PROB, A, FAC, IFAULT)
C
C         ALGORITHM AS 145 APPL. STATIST. (1979) VOL.28, NO.3
C
C         N IS THE NUMBER OF BALLS
C         K IS THE NUMBER OF BOXES
C         PROB(I) IS THE PROBABILITY THAT THE MAXIMUM
C         NUMBER OF BALLS IN A BOX IS I + NMIN - 1
C
        REAL PROB(NN), FAC(NNK), NUM, NUP
        INTEGER A(K), SUM
C
	PARAMETER (ZERO=0.0,TEN=10.0,EXPLIM=75.0)
C
        IFAULT = 1
        NUP = FLOAT(N) * ALOG10(FLOAT(K))
        IF (NUP .GT. EXPLIM) RETURN
        NUP = TEN ** NUP
        IFAULT = 2
        SUM = (N + K -1) / K
        IF (NMIN .GT. N .OR. NMIN .LT. SUM) RETURN
        IFAULT = 3
        IF (N .LE. 1) RETURN
        IFAULT = 0
        DO 1 I = 1, NN
1       PROB(I) = ZERO
C
C         THIS GENERATES THE FIRST PARTITION
C
        M = 1
        A(1) = N
C
C         A NEW PARTITION IS STORED IN A(1) TO A(M)
C
2       NUM = FAC(N)
C
C         THIS COMPUTES THE LOGARITHM OF THE NUMBER
C         OF PARTITIONS OF BALLS
C
        DO 3 I = 1, M
          NA = A(I)
          NUM = NUM - FAC(NA)
3       CONTINUE
C
C         THIS COMPUTES THE LOGARITHM OF THE NUMBER
C         OF PERMUTATIONS OF BOXES
C
        NUM = NUM + FAC(K)
        KM = K - M
        IF (M .LT. K) NUM = NUM - FAC(KM)
        IF (M .EQ. 1) GOTO 6
        I = 1
        II = 2
4       IF (A(II) .EQ. A(I)) GOTO 5
        IJ = II - I
        NUM = NUM - FAC(IJ)
        I = II
5       II = II + 1
        IF (II .LE. M) GOTO 4
        IJ = II - I
        IF (IJ .GT. 1) NUM = NUM - FAC(IJ)
C
C         THE NUMBER OF ALLOCATIONS GIVING RISE TO THIS PARTITION
C         IS CUMULATED ACCORDING TO THE MAXIMUM FREQUENCY
C
6       KX = A(1) - NMIN + 1
        PROB(KX) = PROB(KX) + EXP(NUM)
C
C         THE MAIN ALGORITHM FOR GENERATING PARTITIONS STARTS HERE
C
7       IF (A(M) .EQ. 1) GOTO 11
        A(M) = A(M) -1
        IF (A(1) .LT. NMIN) GOTO 13
8       M = M + 1
        IF (M .GT. K) GOTO 12
        SUM = A(1)
        L = 2
9       IF (L .GE. M) GOTO 10
        SUM = SUM + A(L)
        L = L + 1
        GOTO 9
10      A(M) = N - SUM
        IF (A(M) .LE. A(M -1)) GOTO 2
        IF (M .GE. K) GOTO 12
        A(M) = A(M -1)
        GOTO 8
11      M = M -1
        GOTO 7
12      M = M -2
        GOTO 7
C
C         ALL PARTITIONS HAVE NOW BEEN GENERATED
C
13      DO 14 I = 1, NN
14      PROB(I) = PROB(I) / NUP
        RETURN
        END
