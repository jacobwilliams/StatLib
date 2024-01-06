      SUBROUTINE RESSMO(X, Y, NIN, CUTOFF, XLO, XHI, BANDW, SMOOTH,             
     +            NOUT, NEWCAL, IWK, W, WA, WB, WC, WK, IFAULT)                 
C                                                                               
C      ALGORITHM AS222  APPL. STATIST. (1987) VOL. 36, NO. 1                    
C      Modified by incorporating AS R73 (Appl. Statist. (1988) vol. 37          
C      (2), p.316)                                                              
C                                                                               
C      Auxiliary routine required: FASTF = algorithm AS 83 as amended           
C      in Griffiths & Hill, 1985 (i.e. with extra final argument)               
C                                                                               
      INTEGER   NIN, NOUT, IWK(NIN), NEWCAL, IFAULT                             
      REAL      X(NIN), Y(NIN), CUTOFF, XLO, XHI, BANDW, SMOOTH(NOUT),          
     +          W(NIN), WA(NOUT), WB(NOUT), WC(NOUT), WK(6, NOUT)               
C                                                                               
C      Local variables                                                          
C                                                                               
      INTEGER NPOW, M, K, N2, JHI, JMAX, I, J, NC, ICAL                         
      REAL    ZERO, ONE, TWO, PI, BIG, RANGE, XSTEP, A, B, XLO1, TEMP,          
     +        RES                                                               
      DATA ZERO, ONE, TWO, PI, BIG /0.0, 1.0, 2.0, 3.14159, 30.0/               
C                                                                               
      RANGE = XHI - XLO                                                         
C                                                                               
C      Check cutoff parameter                                                   
C                                                                               
      IF (CUTOFF .LE. ZERO) THEN                                                
        IFAULT = 1                                                              
        RETURN                                                                  
      END IF                                                                    
C                                                                               
C      Check if NOUT is a power of 2                                            
C                                                                               
      NPOW = 8                                                                  
      M = 16                                                                    
      DO 2 K = 3, 11                                                            
        IF (NPOW .EQ. NOUT) GO TO 3                                             
        M = M * 2                                                               
        NPOW = NPOW * 2                                                         
    2 CONTINUE                                                                  
      IFAULT = 2                                                                
      RETURN                                                                    
C                                                                               
C      Check range                                                              
C                                                                               
    3 IF (RANGE .LE. ZERO) THEN                                                 
        IFAULT = 3                                                              
        RETURN                                                                  
      END IF                                                                    
C                                                                               
C      Check bandwidth                                                          
C                                                                               
      IF (BANDW .LE. ZERO) THEN                                                 
        IFAULT = 4                                                              
        RETURN                                                                  
      END IF                                                                    
C                                                                               
C      Check NEWCAL                                                             
C                                                                               
      IF (NEWCAL .LT. 0 .OR. NEWCAL .GT. 2) THEN                                
        IFAULT = 5                                                              
        RETURN                                                                  
      END IF                                                                    
C                                                                               
      XSTEP = RANGE / FLOAT(NOUT)                                               
      A = FLOAT(NOUT) / (RANGE * FLOAT(NIN))                                    
      N2 = NOUT / 2                                                             
      B = TWO * (PI * BANDW / RANGE)**2                                         
      IF (NEWCAL .EQ. 2) GO TO 8                                                
C                                                                               
      JHI = SQRT(BIG / B)                                                       
      JMAX = MIN(N2-1, JHI)                                                     
      XLO1 = XLO - XSTEP                                                        
      IF (NEWCAL .EQ. 0) THEN                                                   
        DO 5 I = 1, NIN                                                         
          IWK(I) = (X(I) - XLO1) / XSTEP                                        
          W(I) = ONE                                                            
    5   CONTINUE                                                                
      END IF                                                                    
C                                                                               
C      Refresh former results                                                   
C                                                                               
      IF (NEWCAL .EQ. 1) THEN                                                   
        DO 52 J = 1, NOUT                                                       
          WA(J) = WK(5, J)                                                      
          WB(J) = WK(6, J)                                                      
   52   CONTINUE                                                                
      END IF                                                                    
C                                                                               
C      Find numerator of Nadaraya-Watson estimate                               
C                                                                               
      CALL LINSMO(Y, IWK, NIN, SMOOTH, NOUT, NEWCAL, WA, WB, WC, N2,            
     +            JHI, JMAX, A, B, M)                                           
C                                                                               
C     Transfer results to work area                                             
C                                                                               
      DO 54 J = 1, NOUT                                                         
        WK(3, J) = WA(J)                                                        
        WK(4, J) = WB(J)                                                        
        WK(1, J) = SMOOTH(J)                                                    
   54 CONTINUE                                                                  
C                                                                               
C     Refresh former results                                                    
C                                                                               
      IF (NEWCAL .EQ. 1) THEN                                                   
        DO 55 J = 1, NOUT                                                       
          WA(J) = WK(5, J)                                                      
          WB(J) = WK(6, J)                                                      
   55   CONTINUE                                                                
      END IF                                                                    
C                                                                               
C      Find density estimate of marginal distribution of X.                     
C                                                                               
      CALL LINSMO(W, IWK, NIN, SMOOTH, NOUT, NEWCAL, WA, WB, WC, N2,            
     +            JHI, JMAX, A, B, M)                                           
      DO 65 J = 1, NOUT                                                         
        WK(5, J) = WA(J)                                                        
        WK(6, J) = WB(J)                                                        
   65 CONTINUE                                                                  
C                                                                               
C      Compute Nadaraya-Watson estimate                                         
C                                                                               
      DO 7 J = 1, NOUT                                                          
        TEMP = ZERO                                                             
        IF (SMOOTH(J) .GT. ZERO) TEMP = WK(1, J) / SMOOTH(J)                    
        WK(1, J) = TEMP                                                         
    7 CONTINUE                                                                  
C                                                                               
C      Compute Huber's PSI from the residuals                                   
C                                                                               
    8 NC = 0                                                                    
      DO 11 I = 1, NIN                                                          
        RES = Y(I) - WK(1, IWK(I))                                              
        IF (RES .GT. -CUTOFF) GO TO 9                                           
        NC = NC + 1                                                             
        W(I) = -CUTOFF                                                          
        GO TO 11                                                                
    9   IF (RES .LT. CUTOFF) GO TO 10                                           
        NC = NC + 1                                                             
        W(I) = CUTOFF                                                           
        GO TO 11                                                                
   10   W(I) = RES                                                              
   11 CONTINUE                                                                  
C                                                                               
C      Compute non-linear correction                                            
C                                                                               
      IFAULT = -NC                                                              
      ICAL = 0                                                                  
      CALL LINSMO(W, IWK, NIN, SMOOTH, NOUT, ICAL, WA, WB, WC, N2, JHI,         
     +            JMAX, A, B, M)                                                
C                                                                               
C      Store the results                                                        
C                                                                               
      DO 111 J = 1, NOUT                                                        
  111 WK(2, J) = SMOOTH(J)                                                      
C                                                                               
C      Derivative of Huber's PSI from residuals                                 
C                                                                               
      DO 112 I = 1, NIN                                                         
        TEMP = ZERO                                                             
        IF (W(I) .LT. CUTOFF .AND. W(I) .GT. -CUTOFF) TEMP = ONE                
        W(I) = TEMP                                                             
  112 CONTINUE                                                                  
C                                                                               
C      Compute denominatot of non-linear correction                             
C                                                                               
      CALL LINSMO(W, IWK, NIN, SMOOTH, NOUT, ICAL, WA, WB, WC, N2, JHI,         
     +            JMAX, A, B, M)                                                
C                                                                               
C      Compute the full estimator                                               
C                                                                               
      DO 13 J = 1, NOUT                                                         
        TEMP = WK(1, J)                                                         
        IF (SMOOTH(J) .GT. ZERO) TEMP = TEMP + WK(2, J) / SMOOTH(J)             
        SMOOTH(J) = TEMP                                                        
   13 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LINSMO(Y, IWK, NIN, SMOOTH, NOUT, NEWCAL, WA, WB, WC,          
     +            N2, JHI, JMAX, A, B, M)                                       
C                                                                               
C      ALGORITHM AS222   APPL. STATIST. (1987) VOL. 36, NO. 1                   
C                                                                               
      INTEGER   NIN, IWK(NIN), NOUT, NEWCAL, N2, JHI, JMAX, M                   
      REAL      Y(NIN), SMOOTH(NOUT), WA(NOUT), WB(NOUT),                       
     +            WC(NOUT), A, B                                                
C                                                                               
C      Local variables                                                          
C                                                                               
      INTEGER M1, J, I, JI, ITYPE, J1, J2, JHI2, IFAULT                         
      REAL ZERO, TEMP, C                                                        
      DATA ZERO/0.0/                                                            
C                                                                               
      M1 = M / 2                                                                
      IF (NEWCAL .NE. 0) GO TO 30                                               
C                                                                               
C      Transform                                                                
C                                                                               
      DO 10 J = 1, NOUT                                                         
        WA(J) = ZERO                                                            
        WB(J) = ZERO                                                            
   10 CONTINUE                                                                  
C                                                                               
      DO 20 I = 1, NIN                                                          
        JI = IWK(I)                                                             
        IF (JI .LT. 1 .OR. JI .GT. NOUT) GO TO 20                               
        WA(JI) = WA(JI) + A*Y(I)                                                
   20 CONTINUE                                                                  
      ITYPE = 1                                                                 
C                                                                               
      CALL FASTF(WA, WB, M1, ITYPE, IFAULT)                                     
C                                                                               
C      Filter Fourier transform                                                 
C                                                                               
   30 SMOOTH(1) = WA(1)                                                         
      WC(1) = WB(1)                                                             
      DO 40 J = 1, JMAX                                                         
        C = EXP(-B * FLOAT(J*J))                                                
        J1 = J + 1                                                              
        J2 = NOUT - J + 1                                                       
        SMOOTH(J1) = C * WA(J1)                                                 
        WC(J1) = C * WB(J1)                                                     
        SMOOTH(J2) = C * WA(J2)                                                 
        WC(J2) = C * WB(J2)                                                     
   40 CONTINUE                                                                  
C                                                                               
C      Underflow correction                                                     
C                                                                               
      IF (JHI + 1 - N2 .GT. 0) THEN                                             
        TEMP = EXP(-B * FLOAT(N2*N2))                                           
        SMOOTH(N2 + 1) = TEMP * WA(N2 + 1)                                      
        WC(N2 + 1) = TEMP * WB(N2 + 1)                                          
      ELSE IF (JHI + 1 - N2 .LT. 0) THEN                                        
        JHI2 = JHI + 2                                                          
        DO 61 J1 = JHI2, N2                                                     
          J2 = NOUT - J1 + 2                                                    
          SMOOTH(J1) = ZERO                                                     
          WC(J1) = ZERO                                                         
          SMOOTH(J2) = ZERO                                                     
          WC(J2) = ZERO                                                         
   61   CONTINUE                                                                
        SMOOTH(N2 + 1) = ZERO                                                   
        WC(N2 + 1) = ZERO                                                       
      ELSE                                                                      
        SMOOTH(N2 + 1) = ZERO                                                   
        WC(N2 + 1) = ZERO                                                       
      END IF                                                                    
C                                                                               
      ITYPE = -1                                                                
      CALL FASTF(SMOOTH, WC, M1, ITYPE, IFAULT)                                 
C                                                                               
      RETURN                                                                    
      END                                                                       

