C This file contains:
C 1. Single-precision version of AS 218 as supplied by the RSS Algorithms
C    Editor
C 2. Double-precision version supplied by Luis Escobar (author)
C 3. A demo program (double precision) supplied by the author
C
C
      SUBROUTINE SEVINF(ITYPE, ZL, ZR, F11, F12, F22, IFAULT)
C
C        ALGORITHM AS218  APPL. STATIST. (1986) VOL. 35, NO. 1
C
C        FISHER INFORMATION MATRIX ELEMENTS FOR TIME (TYPE I) OR
C        FAILURE (TYPE II) CENSORED SMALLEST EXTREME VALUE DATA
C
      REAL BG, DEN, ETA, F11, F12, F22, G, ONE, THET0L, THET0R, THET1L,
     *  THET1R, THET2L, THET2R, X, ZERO, ZL, ZMAX, ZMIN, ZR, ZU
C
      DATA ONE /1.0E0/, ZERO /0.0E0/, ZMAX /3.6E0/, ZMIN /-34.0E0/
C
C        FUNCTION STATEMENTS - PDF (G(X)) AND CDF (BG(X)) FOR A
C        STANDARDIZED SMALLEST EXTREME VALUE RANDOM VARIABLE
C
      G(X) = EXP(X - EXP(X))
      BG(X) = ONE - EXP(-EXP(X))
C
C        CHECK FOR ILLEGAL ITYPE
C
      IF (ITYPE .LT. 1 .OR. ITYPE .GT. 4) GOTO 6
C
C        IDENTIFY THE CENSORING TYPE
C        ITYPE = 1, 2, 3, 4 CORRESPOND TO COMPLETE, RIGHT,
C        LEFT, AND LEFT AND RIGHT CENSORED DATA, RESPECTIVELY
C
      GOTO (1, 2, 3, 4), ITYPE
C
C        COMPLETE DATA - FISHER MATRIX ELEMENTS
C
    1 CALL INTEGR(ZMAX, F11, F12, F22, IFAULT)
      RETURN
C
C        RIGHT CENSORED DATA - FISHER MATRIX ELEMENTS
C
    2 CALL INTEGR(ZR, F11, F12, F22, IFAULT)
      RETURN
C
C        LEFT CENSORED DATA
C        SET THE RIGHT CENSORING BOUND (ZU) EQUAL TO ZMAX
C
    3 ZU = ZMAX
      GOTO 5
C
C        LEFT AND RIGHT CENSORED DATA - CHECK THAT THE LEFT
C        CENSORING BOUND IS NOT GREATER THAN THE RIGHT BOUND
C
    4 IF (ZL .GT. ZR) GOTO 7
      ZU = ZR
C
C        LEFT OR LEFT AND RIGHT CENSORED DATA - FISHER MATRIX ELEMENTS
C        ETA IS THE SECOND TERM IN THE EQUATION FOR F11
C
    5 ETA = ZERO
      DEN = ZERO
      IF (ZL .GT. ZMIN .AND. ZL .LT. ZMAX) DEN = BG(ZL) * EXP(-EXP(ZL))
      IF (DEN .GT. ZERO) ETA = G(ZL) ** 2 / DEN
      CALL INTEGR(ZL, THET0L, THET1L, THET2L, IFAULT)
      IF (IFAULT .EQ. 1) RETURN
      CALL INTEGR(ZU, THET0R, THET1R, THET2R, IFAULT)
      IF (IFAULT .EQ. 1) RETURN
      F11 = THET0R - THET0L + ETA
      F12 = THET1R - THET1L + ZL * ETA
      F22 = THET2R - THET2L + ZL * ZL * ETA
      RETURN
C
C        SET IFAULT
C
    6 IFAULT = 2
      RETURN
    7 IFAULT = 3
      RETURN
      END
C
      SUBROUTINE INTEGR(Z, THETA0, THETA1, THETA2, IFAULT)
C
C        COMPUTES THE INTEGRALS FROM (-INFINITY) TO Z OF THE FOLLOWING
C        FUNCTIONS: G(X)=EXP(X-EXP(X)), G1(X)=(1+X)*G(X), AND
C        G2(X)=(1+X)*G1(X).  G(X) HAS A CLOSED FORM INTEGRAL.
C        G1 AND G2 ARE INTEGRATED BY USING POWER SERIES EXPANSIONS
C        WHEN Z.LE.1, AND BY POWER SERIES EXPANSIONS THROUGH Z=1
C        PLUS A GAUSSIAN QUADRATURE FROM 1 TO Z WHEN Z.GT.1
C
      REAL ASYMP1, ASYMP2, EZ, FACTOR, G1, G1XL, G1XR, G2XL, G2XR, HALF,
     *  ONE, S1, S2, THETA0, THETA1, THETA2, THET11, THET21, TOL, TWO,
     *  X, X1, X2, X3, X4, X5, Z, ZERO, ZMAX, ZMIN
C
      REAL P(8), W(8)
C
C        P AND W ARE CONSTANTS USED IN THE GAUSSIAN QUADRATURES
C
      DATA P(1) /0.4947004674958250E0/, P(2) /0.4722875115366163E0/,
     *  P(3) /0.4328156011939159E0/, P(4) /0.3777022041775015E0/, P(5)
     *  /0.3089381222013219E0/, P(6) /0.2290083888286137E0/, P(7)
     *  /0.1408017753896295E0/, P(8) /0.4750625491881872E-1/
      DATA W(1) /0.1357622970587705E-1/, W(2) /0.3112676196932395E-1/,
     *  W(3) /0.4757925584124639E-1/, W(4) /0.6231448562776694E-1/,
     *  W(5) /0.7479799440828837E-1/, W(6) /0.8457825969750127E-1/, W(7)
     *  /0.9130170752246179E-1/, W(8) /0.9472530522753425E-1/
C
C        ASYMP1, ASYMP2, THET11, AND THET21 ARE THE INTEGRALS
C        FOR G1 AND G2 WHEN Z=+INFINITY AND Z=1 ,RESPECTIVELY
C
      DATA ASYMP1 /0.4227843350984671E0/, ASYMP2 /0.1823680660852879E1/,
     *  HALF /0.5E0/, JMAX /25/, ONE /1.0E0/, TOL /1.0E-11/, THET11 /
     *  0.2720757938345342E0/, THET21 /0.1475933122158450E1/, TWO /
     *  2.0E0/, ZERO /0.0E0/, ZMAX /3.6E0/, ZMIN /-34.0E0/
C
C        FUNCTION REQUIRED IN THE GAUSSIAN QUADRATURE
C
      G1(X) = (ONE + X) * EXP(X - EXP(X))
C
C        CHECK FOR EXTREME VALUES
C
      IFAULT = 0
      IF (Z .GE. ZMAX) GOTO 5
      IF (Z .LE. ZMIN) GOTO 6
C
C        COMPUTATION OF THETA0 - INTEGRAL OF G(X)=EXP(X-EXP(X))
C
      EZ = EXP(Z)
      THETA0 = ONE - EXP(-EZ)
C
C        SELECT INTEGRATION BY POWER SERIES EXPANSIONS OR
C        BY POWER SERIES EXPANSIONS AND GAUSSIAN QUADRATURE
C
      IF (Z .GE. ONE) GOTO 3
C
C        Z IS LESS THAN 1 - INTEGRATION BY POWER SERIES EXPANSIONS
C
C        COMPUTATION OF THETA1 AND THETA2
C
      FACTOR = -ONE
      S1 = ZERO
      S2 = ZERO
      DO 1 J = 1, JMAX
C
C        TERMS TO EVALUATE THETA1(Z) - INTEGRAL OF G1(X)=(1+X)G(X)
C
      X5 = FLOAT(J)
      FACTOR = -FACTOR * EZ / X5
      X1 = Z - ONE / X5
      X3 = FACTOR * X1
      S1 = S1 + X3
C
C        TERMS TO EVALUATE THETA2(Z) - INTEGRAL OF G2(X)=(1+X)G1(X)
C
      X2 = X1 * X1 + ONE / (X5 * X5)
      X4 = FACTOR * X2
      S2 = S2 + X4
C
C        TESTS FOR CONVERGENCE OF THE SERIES -  THE SUMMATIONS STOP WHEN
C        THE ABSOLUTE VALUE OF THE LAST ADDED TERM IS SMALLER THAN TOL
C        FOR BOTH SERIES - A FAULT IS DECLARED IF CONVERGENCE IS NOT
C        REACHED IN A MAXIMUM OF JMAX TERMS
C
      X5 = AMAX1(ABS(X3), ABS(X4))
      IF (X5 .LT. TOL) GOTO 2
    1 CONTINUE
      IFAULT = 1
      RETURN
C
C        ADD TERMS TO OBTAIN INTEGRALS
C
    2 THETA1 = THETA0 + S1
      THETA2 = TWO * THETA1 - THETA0 + S2
      RETURN
C
C        Z IS BETWEEN 1 AND 3.6.
C        THET11 AND THET12 CONTAIN THE INTEGRALS OF G1(X) AND G2(X)
C        OVER (-INFINITY,1) AND THEY ARE DEFINED IN THE DATA STATEMENTS
C
C        GAUSSIAN QUADRATURE TO INTEGRATE G1(X) AND G2(X) OVER (1,Z)
C
    3 S1 = ZERO
      S2 = ZERO
      X1 = HALF * (Z + ONE)
      X2 = Z - ONE
      DO 4 IQUAD = 1, 8
      X3 = X1 - P(IQUAD) * X2
      X4 = X1 + P(IQUAD) * X2
      G1XL = G1(X3)
      G1XR = G1(X4)
      G2XL = (ONE + X3) * G1XL
      G2XR = (ONE + X4) * G1XR
      S1 = S1 + W(IQUAD) * (G1XL + G1XR)
      S2 = S2 + W(IQUAD) * (G2XL + G2XR)
    4 CONTINUE
      S1 = S1 * X2
      S2 = S2 * X2
C
C        ADD TERMS TO OBTAIN INTEGRALS
C
      THETA1 = S1 + THET11
      THETA2 = S2 + THET21
      RETURN
C
C        ASYMPTOTIC VALUES.  Z IS GREATER THAN OR EQUAL TO 3.6
C
    5 THETA0 = ONE
      THETA1 = ASYMP1
      THETA2 = ASYMP2
      RETURN
C
C        ASYMPTOTIC VALUES.  Z IS SMALLER THAN OR EQUAL TO -34
C
    6 THETA0 = ZERO
      THETA1 = ZERO
      THETA2 = ZERO
      RETURN
      END
C
C------------------------- Double precision version -------------------
C
      SUBROUTINE SEVINF(ITYPE,ZL,ZR,F11,F12,F22,IFAULT)
C                                                                        
C        ALGORITHM AS218  APPL. STATIST. (1986) VOL. 35, NO. 1           
C 
C       FISHER INFORMATION MATRIX ELEMENTS FOR TIME (TYPE I) OR          
C       FAILURE (TYPE II) CENSORED SMALLEST EXTREME VALUE DATA           
C
      DOUBLE PRECISION BG, DEN, ETA, F11, F12, F22, G, ONE, THET0L,      
     * THET0R, THET1L, THET1R, THET2L, THET2R, X, ZERO,                  
     * ZL, ZMAX, ZMIN, ZR, ZU
C
      DATA  ONE/1.0D0/, ZERO/0.0D0/, ZMAX/3.6D0/, ZMIN/-34.0D0/          
C                                                                        
C       FUNCTION STATEMENTS - PDF (G(X)) AND CDF (BG(X)) FOR A           
C       STANDARDIZED SMALLEST EXTREME VALUE RANDOM VARIABLE              
C
      G(X) = EXP(X-EXP(X))                                               
      BG(X) = ONE-EXP(-EXP(X))                                           
C
C       CHECK FOR ILLEGAL ITYPE                                          
C
      IF (ITYPE.LT.1 .OR. ITYPE.GT.4) GOTO 6                               
C
C       IDENTIFY THE CENSORING TYPE                                      
C       ITYPE = 1, 2, 3, 4 CORRESPOND TO COMPLETE, RIGHT,                
C       LEFT, AND LEFT AND RIGHT CENSORED DATA, RESPECTIVELY             
C
      GOTO (1,2,3,4), ITYPE                                              
C
C       COMPLETE DATA - FISHER MATRIX ELEMENTS                           
C
    1 CALL INTEGR(ZMAX,F11,F12,F22,IFAULT)                               
      RETURN                                                             
C
C       RIGHT CENSORED DATA - FISHER MATRIX ELEMENTS                     
C
    2 CALL INTEGR(ZR,F11,F12,F22,IFAULT)                                 
      RETURN                                                             
C
C       LEFT CENSORED DATA                                               
C       SET THE RIGHT CENSORING BOUND (ZU) EQUAL TO ZMAX                 
C
    3 ZU = ZMAX                                                            
      GOTO 5                                                             
C
C       LEFT AND RIGHT CENSORED DATA - CHECK THAT THE LEFT               
C       CENSORING BOUND IS NOT GREATER THAN THE RIGHT BOUND              
C
    4 IF (ZL.GT.ZR) GOTO 7                                               
      ZU = ZR                                                              
C
C       LEFT OR LEFT AND RIGHT CENSORED DATA - FISHER MATRIX ELEMENTS    
C       ETA IS THE SECOND TERM IN THE EQUATION FOR F11                   
C                                                                        
    5 ETA = ZERO                                                           
      DEN = ZERO                                                           
      IF (ZL.GT.ZMIN.AND.ZL.LT.ZMAX) DEN = BG(ZL)*DEXP(-DEXP(ZL))          
      IF (DEN.GT.ZERO) ETA = G(ZL)**2/DEN                                  
      CALL INTEGR(ZL,THET0L,THET1L,THET2L,IFAULT)                        
      IF (IFAULT.EQ.1) RETURN                                            
      CALL INTEGR(ZU,THET0R,THET1R,THET2R,IFAULT)                        
      IF (IFAULT.EQ.1) RETURN                                            
      F11 = THET0R-THET0L+ETA                                              
      F12 = THET1R-THET1L+ZL*ETA                                           
      F22 = THET2R-THET2L+ZL*ZL*ETA                                        
      RETURN                                                             
C                                                                        
C       SET IFAULT                                                       
C                                                                        
    6 IFAULT = 2                                                           
      RETURN                                                             
    7 IFAULT = 3                                                           
      RETURN                                                             
      END                                                                

      SUBROUTINE INTEGR(Z,THETA0,THETA1,THETA2,IFAULT)                   
C                                                                        
C       COMPUTES THE INTEGRALS FROM (-INFINITY) TO Z OF THE FOLLOWING    
C       FUNCTIONS: G(X) = EXP(X-EXP(X)), G1(X) = (1+X)*G(X), AND             
C       G2(X) = (1+X)*G1(X).  G(X) HAS A CLOSED FORM INTEGRAL.             
C       G1 AND G2 ARE INTEGRATED BY USING POWER SERIES EXPANSIONS        
C       WHEN Z.LE.1, AND BY POWER SERIES EXPANSIONS THROUGH Z = 1          
C       PLUS A GAUSSIAN QUADRATURE FROM 1 TO Z WHEN Z.GT.1               
C                                                                        
      DOUBLE PRECISION ASYMP1, ASYMP2, EZ, FACTOR, G1, G1XL, G1XR,       
     * G2XL, G2XR, HALF, ONE, S1, S2, THETA0, THETA1,                    
     * THETA2, THET11, THET21, TOL, TWO, X, X1, X2, X3,                  
     * X4, X5, Z, ZERO, ZMAX, ZMIN                                       
C                                                                        
      DOUBLE PRECISION P(8), W(8)                                        
C                                                                        
C       P AND W ARE CONSTANTS USED IN THE GAUSSIAN QUADRATURES           
C                                                                        
      DATA P(1)/0.4947004674958250D0 /, P(2)/0.4722875115366163D0 /,     
     *     P(3)/0.4328156011939159D0 /, P(4)/0.3777022041775015D0 /,     
     *     P(5)/0.3089381222013219D0 /, P(6)/0.2290083888286137D0 /,     
     *     P(7)/0.1408017753896295D0 /, P(8)/0.4750625491881872D-1/      
      DATA W(1)/0.1357622970587705D-1/, W(2)/0.3112676196932395D-1/,     
     *     W(3)/0.4757925584124639D-1/, W(4)/0.6231448562776694D-1/,     
     *     W(5)/0.7479799440828837D-1/, W(6)/0.8457825969750127D-1/,     
     *     W(7)/0.9130170752246179D-1/, W(8)/0.9472530522753425D-1/      
C                                                                        
C       ASYMP1, ASYMP2, THET11, AND THET21 ARE THE INTEGRALS             
C       FOR G1 AND G2 WHEN Z = +INFINITY AND Z = 1 ,RESPECTIVELY             
C                                                                        
      DATA ASYMP1/0.4227843350984671D0/, ASYMP2/0.1823680660852879D1/,   
     *     HALF/0.5D0/, JMAX/25/, ONE/1.0D0/, TOL/1.0D-11/,              
     *     THET11/0.2720757938345342D0/, THET21/0.1475933122158450D1/,   
     *     TWO/2.0D0/, ZERO/0.0D0/, ZMAX/3.6D0/, ZMIN/-34.0D0/           
C                                                                        
C       FUNCTION REQUIRED IN THE GAUSSIAN QUADRATURE                     
C                                                                        
      G1(X) = (ONE+X)*EXP(X-EXP(X))                                      
C                                                                        
C       CHECK FOR EXTREME VALUES                                         
C                                                                        
      IFAULT = 0                                                           
      IF(Z.GE.ZMAX) GOTO 5                                               
      IF(Z.LE.ZMIN) GOTO 6                                               
C                                                                        
C       COMPUTATION OF THETA0 - INTEGRAL OF G(X) = EXP(X-EXP(X))           
C                                                                        
      EZ = DEXP(Z)                                                         
      THETA0 = ONE-DEXP(-EZ)                                               
C                                                                        
C       SELECT INTEGRATION BY POWER SERIES EXPANSIONS OR                 
C       BY POWER SERIES EXPANSIONS AND GAUSSIAN QUADRATURE               
C                                                                        
      IF(Z.GE.ONE) GOTO 3                                                
C                                                                        
C       Z IS LESS THAN 1 - INTEGRATION BY POWER SERIES EXPANSIONS        
C                                                                        
C       COMPUTATION OF THETA1 AND THETA2                                 
C                                                                        
      FACTOR = -ONE                                                        
      S1 = ZERO                                                            
      S2 = ZERO                                                            
      DO 1 J = 1,JMAX                                                      
C                                                                        
C       TERMS TO EVALUATE THETA1(Z) - INTEGRAL OF G1(X) = (1+X)G(X)        
C                                                                        
      X5 = FLOAT(J)                                                       
      FACTOR = -FACTOR*EZ/X5                                               
      X1 = Z-ONE/X5                                                        
      X3 = FACTOR*X1                                                       
      S1 = S1+X3                                                           
C                                                                        
C       TERMS TO EVALUATE THETA2(Z) - INTEGRAL OF G2(X) = (1+X)G1(X)       
C                                                                        
      X2 = X1*X1+ONE/(X5*X5)                                               
      X4 = FACTOR*X2                                                       
      S2 = S2+X4                                                           
C                                                                        
C       TESTS FOR CONVERGENCE OF THE SERIES -  THE SUMMATIONS STOP WHEN  
C       THE ABSOLUTE VALUE OF THE LAST ADDED TERM IS SMALLER THAN TOL    
C       FOR BOTH SERIES - A FAULT IS DECLARED IF CONVERGENCE IS NOT      
C       REACHED IN A MAXIMUM OF JMAX TERMS                               
C                                                                        
      X5 = MAX1(ABS(X3),ABS(X4))                                        
      IF (X5.LT.TOL) GOTO 2                                              
    1 CONTINUE                                                           
      IFAULT = 1                                                           
      RETURN                                                             
C                                                                        
C       ADD TERMS TO OBTAIN INTEGRALS                                    
C                                                                        
    2 THETA1 = THETA0+S1                                                   
      THETA2 = TWO*THETA1-THETA0+S2                                        
      RETURN                                                             
C                                                                        
C       Z IS BETWEEN 1 AND 3.6.                                          
C       THET11 AND THET12 CONTAIN THE INTEGRALS OF G1(X) AND G2(X)       
C       OVER (-INFINITY,1) AND THEY ARE DEFINED IN THE DATA STATEMENTS   
C                                                                        
C       GAUSSIAN QUADRATURE TO INTEGRATE G1(X) AND G2(X) OVER (1,Z)      
C                                                                        
    3 S1 = ZERO                                                            
      S2 = ZERO                                                            
      X1 = HALF*(Z+ONE)                                                    
      X2 = Z-ONE                                                           
      DO 4 IQUAD = 1,8                                                     
      X3 = X1-P(IQUAD)*X2                                                  
      X4 = X1+P(IQUAD)*X2                                                  
      G1XL = G1(X3)                                                        
      G1XR = G1(X4)                                                        
      G2XL = (ONE+X3)*G1XL                                                 
      G2XR = (ONE+X4)*G1XR                                                 
      S1 = S1+W(IQUAD)*(G1XL+G1XR)                                         
      S2 = S2+W(IQUAD)*(G2XL+G2XR)                                         
    4 CONTINUE                                                           
      S1 = S1*X2                                                           
      S2 = S2*X2                                                           
C                                                                        
C       ADD TERMS TO OBTAIN INTEGRALS                                    
C                                                                        
      THETA1 = S1+THET11                                                   
      THETA2 = S2+THET21                                                   
      RETURN                                                             
C                                                                        
C       ASYMPTOTIC VALUES.  Z IS GREATER THAN OR EQUAL TO 3.6            
C                                                                        
    5 THETA0 = ONE                                                         
      THETA1 = ASYMP1                                                      
      THETA2 = ASYMP2                                                      
      RETURN                                                             
C                                                                        
C       ASYMPTOTIC VALUES.  Z IS SMALLER THAN OR EQUAL TO -34            
C                                                                        
    6 THETA0 = ZERO                                                        
      THETA1 = ZERO                                                        
      THETA2 = ZERO                                                        
      RETURN                                                             
      END                                                                
C
C------------------------ Demo program (double precision) -------------
C
C
C         DRIVER PROGRAM FOR SEVINF
C
      PROGRAM DEM218
      DOUBLE PRECISION
     *     COV, DET, F11, F12, F22, ONE, Q1, Q2,
     *     SEVEN, TENTH, VAR1, VAR2, ZERO, ZL, ZR
C
      DATA ONE/1.0D0/, SEVEN/7.0D0/, TENTH/0.1D0/, ZERO/0.0D0/
      DATA IUNIT0/10/, IUNIT1/15/
C
      OPEN(IUNIT0, FILE='HARTER.OUT', STATUS='NEW')
      OPEN(IUNIT1, FILE='MEEKER.OUT', STATUS='NEW')
C
1     FORMAT(1H1, 2X, 'Q1', 4X, 'Q2', 8X, 'VAR(U)', 8X, 'VAR(S)', 8X,
     * 'COV(U,S)', 2X, 'IFAULT')
2     FORMAT(1H , 2F6.2, 3(3X, F11.6), I3)
3     FORMAT(1H1, ' ZETA', 4X, 'VAR(U)', 9X, 'VAR(S)', 9X,
     * 'COV(U,S)', 2X, 'IFAULT')
4     FORMAT(1H , F5.2, 3(3X, F11.5), I3)
C
C     HERE WE REPRODUCE THE FIRST FIVE COLUMNS OF TABLE 1
C     IN HARTER AND MOORE JASA(1968), P. 889-901
C
      WRITE(IUNIT0, 1)
      DO 20 I = 1, 10
	Q1 = TENTH*FLOAT(I-1)
	IF (Q1.GT.ZERO) ZL = LOG(-LOG(ONE-Q1))
	JMAX = 11 - I
	DO 10 J = 1, JMAX
	  Q2 = TENTH*FLOAT(J-1)
	  ITYPE = 4
	  IF (Q1.GT.ZERO.AND.Q2.EQ.ZERO) ITYPE = 3
	  IF (Q1.EQ.ZERO.AND.Q2.GT.ZERO) ITYPE = 2
	  IF (Q1.EQ.ZERO.AND.Q2.EQ.ZERO) ITYPE = 1
	  IF (Q2.GT.ZERO) ZR = LOG(-LOG(Q2))
	  CALL SEVINF(ITYPE, ZL, ZR, F11, F12, F22, IFAULT)
	  DET = F11*F22 - F12*F12
	  VAR1 = F22/DET
	  VAR2 = F11/DET
	  COV = -F12/DET
	  WRITE(IUNIT0, 2) Q1, Q2, VAR1, VAR2, COV, IFAULT
10      CONTINUE
20    CONTINUE
C
C     HERE WE PARTIALLY REPRODUCE TABLE 1 IN MEEKER AND
C     NELSON, TECHNOMETRICS 1977, VOL. 19, NO. 4, 473-476
C
      WRITE(IUNIT1, 3)
      ITYPE = 2
      DO 30 J = 1, 91
	ZR = -SEVEN + TENTH*FLOAT(J-1)
	CALL SEVINF(ITYPE, ZL, ZR, F11, F12, F22, IFAULT)
	DET = F11*F22 - F12*F12
	VAR1 = F22/DET
	VAR2 = F11/DET
	COV = -F12/DET
	WRITE(IUNIT1, 4) ZR, VAR1, VAR2, COV, IFAULT
30    CONTINUE
      STOP
      END
