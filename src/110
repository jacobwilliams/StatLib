                                                                     
                                             
      SUBROUTINE LPEST(N, P, X, Y, MAXIT, A, B, SD, R, RATE, IT, NPO,
     * IFAULT)
C
C FORMAL PARAMETERS
C    N     INTEGER        input : the number of points
C    P     REAL           input : p in the Lp norm
C    X     REAL ARRAY (N) input : the observed values x(i)
C    Y     REAL ARRAY (N) input : the observed values y(i)
C    MAXIT INTEGER        input : maximum allowable number of iterations
C    A     REAL           output : the estimate of alfa
C    B     REAL           output : the estimate of beta
C    SD    REAL           output : the Lp norm
C    R     REAL ARRAY (N) output : the signed residuals
C    RATE  REAL           output : abs(S(k+1)-S(k))/S(k+1), where S(k) is
C                                  the Lp norm of the fit on the kth iteration
C    IT    INTEGER        output : the number of iterations
C    NPO   INTEGER        output : the number of the points on the line
C   IFAULT INTEGER        output :
C                             IFAULT = 0 if the routine converged.
C                                    = 1 if return was due to an increase
C                                        in the norm.
C                                    = 2 if the maximum iterations
C                                        specified was less than 2.
C                                    = 3 if the weighted sample variance of
C                                        x is 0.
C                                    = 4 if N given less than 2.
C                                    = 5 if convergence has not be achieved within
C                                        the maximum number of iterations specified.
C
C  ALGORITHM AS 110 APPL. STATIST. (1977), VOL.26, NO.1
C
C  LP-NORM FIT OF STRAIGHT LINE BY EXTENSION OF SCHLOSSMACHER
C
      implicit INTEGER (h-n)
      implicit REAL (a-g)
      implicit REAL (o-z)
      REAL X(N), Y(N), R(N)
      DATA EPS /1.0E-6/
C
      IF(MAXIT .LT. 2) GOTO 9
      IF(N .LT. 2) GOTO 10
      IFAULT = 0
      WP = P - 2.0
      EPS2 = 2.0 * EPS
      SD = 0.0
      DO 1 I = 1,N
1     R(I) = 1.0
      DO 6 IT = 1, MAXIT
        NPO = 0
C
C          CALCULATE A AND B BY LEAST SQUARES ON WEIGHTED DATA,
C          USING THE HERRAMAN ALGORITHM.
C          OMIT OBSERVATIONS WITH SMALL RESIDUALS.
C
        SW = 0.0
        XMEAN = 0.0
        YMEAN = 0.0
        SSX = 0.0
        SPXY = 0.0
        DO 3 I = 1, N
          ABSRI = DABS(R(I))
          IF(ABSRI .LE. EPS) GOTO 2
          W = ABSRI ** WP
          SW = SW + W
          DIV = W / SW
          XI = X(I) - XMEAN
          YI = Y(I) - YMEAN
          XIW = XI * W
          DX = XI * XIW
          DXY = YI * XIW
          SSX = SSX + DX - DX * DIV
          SPXY = SPXY + DXY - DXY * DIV
          XMEAN = XMEAN + XI * DIV
          YMEAN = YMEAN + YI * DIV
          GOTO 3
2         NPO = NPO + 1
3       CONTINUE
        IF(SSX .LT. EPS) GOTO 11
        B = SPXY / SSX
        A = YMEAN - B * XMEAN
C
C          FORM RESIDUALS AND TEST CONVERGENCE
C
        SD2 = 0.0
        ISW = 0
        DO 4 I = 1, N
          RES = Y(I) - A - B * X(I)
          ABSRI = DABS(RES)
          IF(DABS(ABSRI - DABS(R(I))) .GT. EPS2) ISW = 1
          SD2 = SD2 + ABSRI ** P
          R(I) = RES
4       CONTINUE
        RATE = DABS(SD2 - SD) / SD2
        IF(ISW .EQ. 0) RETURN
        IF(IT .EQ. 1) GOTO 5
C
C        TEST FOR INCREASE IN NORM
C
        IF(SD2 .GT. SD) GOTO 7
C
C ***  CODE ADDED UNDER REMARK AS R29 APPL. STATIST. (1979) VOL.28, NO.1
C
        SD3 = SD2-SD 
        IF (NPO.EQ.0) GOTO 12
        IF (SD3.GT.EPS.AND.SD4.GT.EPS) GOTO 7
        SD4 = SD3 
        GOTO 5 
12      IF (SD3.GT.EPS) GOTO 7 
        SD4= -10.0
C
C ***  END OF CODE ADDITION
C
5       SD = SD2
        A2 = A
        B2 = B
6     CONTINUE
C     FAILED TO CONVERGE IN MAXIT ITERATIONS
C
      IFAULT = 5
      RETURN
C
C      NORM INCREASED, RESTORE A, B, AND R, THEN RETURN
C
7     IFAULT = 1
      A = A2
      B = B2
      DO 8 I = 1, N
8     R(I) = Y(I) - A - B * X(I)
      RETURN
C
C      MAXIT SPECIFIED LESS THAN 2
C
9     IFAULT = 2
      RETURN
C
C      N LESS THAN 2
C
10     IFAULT = 4
       RETURN
C
C      VARIANCE OF WEIGHTED X IS ZERO
C
11    IFAULT = 3
      RETURN
      END

