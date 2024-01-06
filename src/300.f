      SUBROUTINE SPP (DI, NIN, SPEC, LARGE, NOUT, XSTORE, XOUT, ISEED,
     *                RANNUM, IFAULT)                                                 
C                                                                       
C        ALGORITHM AS 300.1 APPL. STATIST. (1996), VOL.45, NO.1
C                                                                       
C        Efficient algorithm for determining a simulated percentile
C        point with a fixed degree of accuracy
C                                                                       
      INTEGER DI, NIN(4), NOUT(3), ISEED, IFAULT                        
      REAL SPEC(4), XSTORE(DI), XOUT(3), RANNUM,                          
      LOGICAL LARGE                                                     
C
      INTEGER HIGH, LOW, N, N1, N2, NBND, NCHK, NDEX1, NDEX2, NSTORE    
      REAL AN, ANP, C1, C2, D, ERR, HALF, MAXERR, ONE, P, PQ, TWO,      
     *     THREE, XNEW, ZERO                                               
C     
      PARAMETER (ONE = 1.0E0, TWO = 2.0E0, THREE = 3.0E0, HALF = 0.5,   
     *           ZERO = 0.0E0)                                                   
C                                                                       
      EXTERNAL INSERT, RANNUM                                       
C
C        Initialise local variables and check for input failures
C                                                                       
      IFAULT = 0                                                        
      N1 = NIN(1)                                                       
      IF (N1 .LE. 2) IFAULT = IFAULT + 1                                
      N2 = NIN(2)                                                       
      IF (N2 .LE. N1) IFAULT = IFAULT + 2                               
      NCHK = NIN(3)                                                     
      IF (NCHK .LE. N2) IFAULT = IFAULT + 4                             
      NBND = NIN(4)                                                     
      IF (NBND .LE. NCHK) IFAULT = IFAULT + 8                           
      C1 = SPEC(1)                                                      
      IF (C1 .LE. ZERO) IFAULT = IFAULT + 16                            
      C2 = SPEC(2)                                                      
      IF (C2 .LE. C1) IFAULT = IFAULT + 32                              
      MAXERR = SPEC(3)                                                  
      IF (MAXERR .LE. ZERO) IFAULT = IFAULT + 64                        
      P = SPEC(4)                                                       
      IF (P .LE. ZERO .OR. P .GT. HALF) IFAULT = IFAULT + 128           
      IF (DI .LT. TWO * C2 * SQRT(NBND * P * (ONE - P)) + THREE)        
     *   IFAULT = IFAULT + 256                                           
      IF (IFAULT .NE. 0) RETURN                                         
      PQ = P * (ONE - P)                                                
      NDEX1= 1                                                          
      NDEX2 = MAX0(2, N1)                                               
C                                                                       
C        Start iterating simulated values
C                                                                       
      DO 50 N = 1, NBND                                                 
         XNEW = RANNUM(ISEED)                                              
         IF (N .EQ. 1) THEN                                                
            XSTORE(1) = XNEW                                                
         ELSE IF (N .EQ. 2) THEN                                           
            XSTORE(2) = XNEW                                                
            NSTORE = 2                                                      
         IF (XNEW .GT. XSTORE(1) .EQV. LARGE) THEN                       
            XSTORE(2) = XSTORE(1)                                         
            XSTORE(1) = XNEW                                              
         ENDIF                                                           
C                                                                       
C        Stage 1: sort and store all of the first N1 simulated values
C                                                                       
         ELSE IF (N .LE. N1) THEN                                          
            CALL INSERT (NSTORE, DI, XNEW, XSTORE, 0, LARGE)                
C                                                                       
C        Stage 2: sort and store only N1 of the first N2 simulated
C                 values
C                                                                       
         ELSE IF (N .LE. N2) THEN                                          
            IF (XNEW .GT. XSTORE(NDEX2) .EQV. LARGE) THEN                   
            CALL INSERT (NSTORE, DI, XNEW, XSTORE, 2, LARGE)              
         ENDIF                                                           
C                                                                       
C        Stage 3: sort and store only a probabilistically determined
C                 number of simulated values
C            3.1: do not store simulated value but change index
C                 numbers
C                                                                       
         ELSE IF (XNEW .GT. XSTORE(1) .EQV. LARGE) THEN                    
            NDEX1 = NDEX1 + 1                                               
            NDEX2 = NDEX2 + 1                                               
C                                                                       
C            3.2: sort and store the simulated value
C                                                                       
         ELSE IF (XNEW .GT. XSTORE(NSTORE) .EQV. LARGE) THEN               
            AN = N                                                          
            ANP = AN * P                                                    
            D = C2 * SQRT(AN * PQ)                                          
            LOW = ANP - D                                                   
C
C          3.2.1: delete the simulated value that had been in the first
C                 storage location
C                                                                       
            IF (NDEX1 .LT. LOW) THEN                                        
               CALL INSERT (NSTORE, DI, XNEW, XSTORE, 1, LARGE)              
               NDEX1 = NDEX1 + 1                                             
               NDEX2 = NDEX2 + 1                                             
C                                                                       
C          3.2.2: no previously stored simulated value is to be
C                 deleted
C                                                                       
            ELSE                                                            
               HIGH = ANP + D + 1                                            
               IF (NDEX2 .LT. HIGH) THEN                                     
                  NDEX2 = NDEX2 + 1                                           
                  CALL INSERT (NSTORE, DI, XNEW, XSTORE, 0, LARGE)
C                                                                       
C          3.2.3: delete the simulated value that had been in the last
C                 storage location
C                                                                       
            ELSE                                                          
               CALL INSERT (NSTORE, DI, XNEW, XSTORE, 2, LARGE)            
            ENDIF                                                         
         ENDIF                                                           
C                                                                       
C            3.3: do not store simulated value and do not change index
C                 numbers
C                                                                       
      ENDIF                                                             
C                                                                       
C        Stage 4: check to see if the desired precision has been 
C                 achieved at the end of each NCHK iterations
C                                                                       
      IF (N .EQ. (N / NCHK) * NCHK .OR. N .EQ. NBND) THEN               
         AN = N                                                          
         ANP = AN * P                                                    
         D = C1* SQRT(AN * PQ)                                           
         LOW = ANP - D + HALF                                            
         HIGH = ANP + D + ONE + HALF                                     
C
C            4.1: the desired order statistics are not in storage                                                        
C                 (abnormal termination)
C                                                                       
         IF (NDEX1 .GT. LOW .OR. NDEX2 .LT. HIGH) THEN                   
            IFAULT = 999                                                  
            GOTO 100                                                      
         ENDIF                                                           
C
C            4.2: the desired order statistics are in storage                                                        
C                                                                       
         LOW = LOW - NDEX1 + 1                                           
         HIGH = HIGH - NDEX1 + 1                                         
         IF (LARGE) THEN                                                 
            XOUT(3) = XSTORE(LOW)                                         
            XOUT(2) = XSTORE(HIGH)                                        
         ELSE                                                            
            XOUT(3) = XSTORE(HIGH)                                        
            XOUT(2) = XSTORE(LOW)                                         
         ENDIF                                                           
            ERR = (XOUT(3) - XOUT(2)) / TWO                                 
            LOW = ANP - NDEX1 + 1                                           
            XOUT(1) = XSTORE(LOW)                                           
C
C            4.3: the desired precision has been achieved when ERR is                                                       
C                 less than MAXERR
C                                                                       
            IF (ERR .LT. MAXERR) GOTO 100                                   
         ENDIF                                                             
   50 CONTINUE                                                          
C
C            4.4: the desired precision is not achieved within the                                                       
C                 fixed number of iterations
C                                                                       
      N = N - 1                                                         
  100 NOUT(1) = N                                                       
      NOUT(2) = NSTORE                                                  
      NOUT(3) = NDEX1                                                   
C      
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INSERT (NSTORE, DI, XNEW, XSTORE, CONTRL, LARGE)       
C                                                                       
C        ALGORITHM AS 300.2 APPL. STATIST. (1996), VOL.45, NO.1
C                                                                       
C        Divide and conquer algorithm to sort and store a new value
C        (XNEW) among a set of previously sorted and stored values
C                                                                       
      INTEGER NSTORE, DI, CONTRL
      REAL XNEW, XSTORE(DI)                                             
      LOGICAL LARGE                                                     
C                                                                       
      INTEGER I, IBND1, IBND2, ITRY                                     
      REAL XHOLD                                                        
C                                                                       
C        Initialise local variables
C                                                                       
      IBND1 = 0                                                         
      IBND2 = NSTORE + 1                                                
C                                                                       
C        Find the appropriate position to insert XNEW into the
C        sorted set
C                                                                       
   50 ITRY = (IBND1 + IBND2) / 2                                        
      IF (ITRY.LT.0 .OR. ITRY.GT.NSTORE) RETURN                         
      IF (ITRY.EQ.IBND1) GOTO 100                                       
      IF (XSTORE(ITRY) .GT. XNEW .EQV. LARGE) THEN                      
         IBND1 = ITRY                                                    
      ELSE                                                              
         IBND2 = ITRY                                                    
      ENDIF                                                             
      GOTO 50                                                           
  100 IF (CONTRL .NE. 1) THEN                                           
C                                                                       
C        If CONTRL equals 0, add XNEW to the list
C                                                                       
         IF (CONTRL .EQ. 0) NSTORE = NSTORE + 1                          
C                                                                       
C        If CONTRL equals 2, add XNEW to the list and delete the value
C        in the last storage location
C                                                                       
         DO 150 I = IBND2, NSTORE                                        
            XHOLD = XSTORE(I)                                               
            XSTORE(I) = XNEW                                                
            XNEW = XHOLD                                                    
  150    CONTINUE                                                        
C                                                                       
C        If CONTRL equals 1, add XNEW to the list and delete the
C        first storage location
C                                                                       
      ELSE                                                              
         DO 200 I = IBND1, 1, -1                                         
            XHOLD = XSTORE(I)                                               
            XSTORE(I) = XNEW                                                
            XNEW = XHOLD                                                    
  200    CONTINUE                                                        
      ENDIF
C
      RETURN                                                            
      END
