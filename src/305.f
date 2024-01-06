      SUBROUTINE ARLM(K,H,U,DELTA,ARL,IFAULT)  
C
C        ALGORITHM AS 305.1 APPL. STATIST. (1996), VOL.45, NO.4
C
C        Compute the ARL of a CUSUM mean chart with linear drift
C
      DOUBLE PRECISION K,H,U,DELTA,ARL
      INTEGER IFAULT
C
      DOUBLE PRECISION SIGMA
      INTEGER NN, NSIDE, NTYPE
C
      EXTERNAL COMARL
C
      NTYPE=1
      NSIDE=1 
      CALL COMARL(NTYPE,NSIDE,NN,K,H,U,DELTA,SIGMA,ARL,IFAULT)
C
      RETURN
      END
C 
      SUBROUTINE ARLV(NSIDE,NN,K,H,U,SIGMA,ARL,IFAULT)
C
C        ALGORITHM AS 305.2 APPL. STATIST. (1996), VOL.45, NO.4
C
C        Compute the ARL of a CUSUM variance chart
C
      INTEGER NSIDE, NN, IFAULT
      DOUBLE PRECISION K,H,U,SIGMA,ARL
C
      INTEGER NTYPE
      DOUBLE PRECISION DELTA
C
      EXTERNAL COMARL
C
      NTYPE=2
      CALL COMARL(NTYPE,NSIDE,NN,K,H,U,DELTA,SIGMA,ARL,IFAULT)
C
      RETURN
      END

C
      SUBROUTINE COMARL(NTYPE,NSIDE,NN,K,H,U,DELTA,SIGMA,ARL,IFAULT)
C
C        ALGORITHM AS 305.3 APPL. STATIST. (1996), VOL.45, NO.4
C
C        Compute the ARL of a CUSUM chart or a CUSUM variance chart
C
      INTEGER NTYPE, NSIDE, NN, IFAULT
      DOUBLE PRECISION K,H,U,DELTA,SIGMA,ARL
C
      INTEGER I, IA, IB, IP1, IFLAG, IPIVOT(25), J, JJ, JP1, JP2,
     *        MMM, MP1, ND2, NPTS, NQM1, NQUAD, NSIDE2
      DOUBLE PRECISION A(25, 25), ALPHA, B(25), BETA, DMEAN(3000),
     *                 P(25), SIDE, SK, SS, W(25), WK(25), X(25),
     *                 XOLD(25)
C
      DOUBLE PRECISION CDFLG, PDFLG, PDFN
      REAL ALNORM, DNORM
C
      EXTERNAL FACTOR, SUBST, SETUP
C
C        Gaussian quadrature abscissas
C
      P(1)=0.9951872199970213D0
      P(2)=0.9747285559713095D0
      P(3)=0.9382745520027327D0
      P(4)=0.8864155270044010D0
      P(5)=0.8200019859739029D0
      P(6)=0.7401241915785543D0
      P(7)=0.6480936519369755D0
      P(8)=0.5454214713888395D0
      P(9)=0.4337935076260451D0
      P(10)=0.3150426796961634D0
      P(11)=0.1911188674736163D0
      P(12)=0.0640568928626056D0
C
C        Gaussian quadrature weights
C
      W(1)=0.0123412297999872D0
      W(2)=0.0285313886289337D0 
      W(3)=0.0442774388174198D0
      W(4)=0.0592985849154368D0
      W(5)=0.0733464814110803D0
      W(6)=0.0861901615319533D0
      W(7)=0.0976186521041139D0
      W(8)=0.1074442701159656D0
      W(9)=0.1155056680537256D0
      W(10)=0.1216704729278034D0
      W(11)=0.1258374563468283D0
      W(12)=0.1279381953467521D0
C
C        Initialise variables
C
      NPTS=25
      NQUAD=NPTS-1                                                      
      ND2=NQUAD/2                                                       
      SIDE=NSIDE       
      SS=NN 
      ALPHA=(SS-1.0D0)/2.0D0
      BETA=(SIGMA*SIGMA)/ALPHA
      SK=SIDE*K  
C
C        Set IFAULT to be 0, 1, 2, 3, 4, 5, 6
C 
      IFAULT=0
      IF(NTYPE.EQ.1.AND.(DELTA.GE.0.0D0.AND.
     *DELTA.LT.0.0001D0)) IFAULT=1
      IF(NTYPE.EQ.1.AND.DELTA.LT.0.0D0) IFAULT=2
      NSIDE2=NSIDE*NSIDE
      IF(NTYPE.EQ.2.AND.NSIDE2.NE.1) IFAULT=3
      IF(NTYPE.EQ.2.AND.SIGMA.LE.0.0D0) IFAULT=4
      IF(NTYPE.EQ.1.AND.(U.LT.0.0D0.OR.U.GT.H)) IFAULT=5
      IF(NTYPE.EQ.2.AND.NSIDE.EQ.1.AND.
     *(U.LT.0.0D0.OR.U.GT.H)) IFAULT=5
      IF(NTYPE.EQ.2.AND.NSIDE.EQ.-1.AND.
     *(U.LT.-H.OR.U.GT.0.0D0)) IFAULT=5
      IF(NTYPE.EQ.2.AND.NN.LE.1) IFAULT=6
      IF(IFAULT.GE.2) RETURN
C
C        Obtain the full Gaussian abscissas and weights  
C
      DO 2 I=1,ND2                                                      
        P(NPTS-I)=-P(I)                                                
        W(NPTS-I)=W(I)                                                 
    2 CONTINUE                                                          
C
C        Transform the Gaussian abscissas and weights
C        from the interval (-1,1) to (0,h) or (-h,0)
C
      DO 3 I=1,NQUAD                                                    
        W(I)=H*W(I)/2.0D0
        P(I)=H*P(I)/2.0D0 + SIDE*H/2.0D0 
    3 CONTINUE                                                          
C
C        Compute the ARL of a CUSUM variance chart under step shift
C
      IF(NTYPE.EQ.2) THEN 
C
C        Replace the integral equation by a system of linear
C        equations AX=B
C
        B(NPTS)=-1.D0
        A(1,1)=-1.0D0+CDFLG(SK,ALPHA,BETA)
        IF(NSIDE.EQ.-1) A(1,1)=-1.0D0-A(1,1)
        DO 4 I=1,NQUAD
          A(1,I+1)=W(I)*PDFLG(SK+P(I),ALPHA,BETA)
          B(I)=-1.D0
    4   CONTINUE
        DO 5 I=1,NQUAD
          IP1=I+1
          A(IP1,1)=CDFLG(SK-P(I),ALPHA,BETA)
          IF(NSIDE.EQ.-1) A(IP1,1)=1.0D0-A(IP1,1)
    5   CONTINUE
        DO 7 I=1,NQUAD
          DO 6 J=1,NQUAD
            IP1=I+1
            JP1=J+1
            IF(I.EQ.J) THEN
              A(IP1,JP1)=W(I)*PDFLG(SK,ALPHA,BETA)-1.0D0
            ELSE
              A(IP1,JP1)=W(J)*PDFLG(P(J)-P(I)+SK,
     *        ALPHA,BETA)
            ENDIF
    6     CONTINUE
    7   CONTINUE
C
C        Call the subroutines FACTOR and SUBST to solve the system
C        of equations for X. Set IFAULT to be 7 if matrix A is 
C        singular.
C
        CALL FACTOR(A,NPTS,WK,IPIVOT,IFLAG)
        IF(IFLAG.EQ.0) THEN
          IFAULT=7   
          RETURN
        ENDIF
        CALL SUBST(A,IPIVOT,B,NPTS,X)
C
C        Compute L(u) of a CUSUM variance chart
C
        ARL=CDFLG(SK-U,ALPHA,BETA)
        IF(NSIDE.EQ.-1) ARL=1.0D0-ARL
        ARL=1.D0+X(1)*ARL+W(1)*X(2)*PDFLG(P(1)+SK-U,
     *  ALPHA,BETA)
        NQM1=NQUAD-1 
        DO 8 J=1,NQM1
          JP1=J+1   
          JP2=J+2
          ARL=ARL+W(JP1)*X(JP2)*PDFLG(P(JP1)+SK-U,
     *    ALPHA,BETA)
    8   CONTINUE
C
C        Compute the ARL of a CUSUM mean chart under linear drift
C
      ELSE
C
C        Compute the process mean mu_0, mu_1, ..., mu_m
C
        CALL SETUP(DELTA,DMEAN,MP1)
C
C        Replace the integral equation by a system of linear
C        equations AX=B
C
        DO 20 I=1,NPTS  
          B(I)=-1.0D0      
          DO 10 J=1,NPTS                       
            IF(I.LE.NQUAD.AND.J.EQ.I) THEN                              
              A(I,J)=W(J)*PDFN(P(J)+K-P(I),
     *        DMEAN(MP1),1.0D0)-1.0D0
            ELSEIF(I.LE.NQUAD.AND.J.LE.NQUAD) THEN                      
              A(I,J)=W(J)*PDFN(P(J)+K-P(I),
     *        DMEAN(MP1),1.0D0)
            ELSEIF(I.LE.NQUAD.AND.J.EQ.NPTS) THEN    
              DNORM=K-P(I)-DMEAN(MP1)
              A(I,J)=ALNORM(DNORM,.FALSE.)     
            ELSEIF(I.EQ.NPTS.AND.J.LT.I) THEN                           
              A(I,J)=W(J)*PDFN(P(J)+K, 
     *        DMEAN(MP1),1.0D0) 
            ELSE
              DNORM=K-DMEAN(MP1) 
              A(I,J)=ALNORM(DNORM,.FALSE.)-1.0D0   
            ENDIF                                                       
   10    CONTINUE
   20  CONTINUE                                                 
C
C        Call the subroutines FACTOR and SUBST to solve the system
C        of equations for X.  Set IFAULT to be 7 if matrix A is
C        singular.
C
        CALL FACTOR(A,NPTS,WK,IPIVOT,IFLAG)                    
        IF(IFLAG.EQ.0) THEN                                     
          IFAULT=7
          RETURN                                                 
        ENDIF                                                           
        CALL SUBST(A,IPIVOT,B,NPTS,X)                            
C
C        Compute intermediate ARL functions 
C
        MMM=MP1-1
        DO 80 J=2,MMM 
          JJ=MP1-J+1
          DO 60 I=1,NPTS 
            XOLD(I)=X(I)
            X(I)=1.0D0 
   60     CONTINUE 
          DO 75 IA=1,NPTS 
            IF(IA.EQ.NPTS) THEN                                         
              DNORM = K-DMEAN(JJ)
            ELSE                                                     
              DNORM = K-P(IA)-DMEAN(JJ)
            ENDIF                                                       
            X(IA) = X(IA)+ALNORM(DNORM,.FALSE.)*XOLD(NPTS)
            DO 70 IB=1,NQUAD                       
              IF(IA.EQ.NPTS) THEN                  
                X(IA) = X(IA)+W(IB)*XOLD(IB)*      
     *          PDFN(P(IB)+K,DMEAN(JJ),1.0D0)
              ELSE                              
                X(IA) = X(IA)+W(IB)*XOLD(IB)*                         
     *          PDFN(P(IB)+K-P(IA),DMEAN(JJ),1.0D0)
              ENDIF                                                    
   70       CONTINUE                                                    
   75     CONTINUE        
   80   CONTINUE                                                   
C
C        Compute L_0(u,0)                 
C
        ARL=1.0D0                                                 
        DNORM = K-U-DMEAN(1)
        ARL = ARL+ALNORM(DNORM,.FALSE.)*X(NPTS) 
        DO 90 J=1,NQUAD                                            
          ARL=ARL+W(J)*PDFN(P(J)+K-U,DMEAN(1),1.0D0)*X(J)
   90   CONTINUE    
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE SETUP(DELTA,DMEAN,MP1)
C
C        ALGORITHM AS 305.4 APPL. STATIST. (1996), VOL.45, NO.4
C
C
C        Determine MP1 and DMEAN where 
C        DMEAN(j) = process mean at (j-1)th sample;
C                   for j=1, 2, ...,MP1
C
      DOUBLE PRECISION DELTA, DMEAN(3000)
      INTEGER MP1
C
      INTEGER J, JM1
C
      IF(DELTA.LT.0.0001) THEN                                      
        MP1=1                                                      
      ELSEIF(DELTA.LE.0.001) THEN
        MP1=3000                                                   
      ELSEIF(DELTA.LE.0.002) THEN                                   
        MP1=2000                                                   
      ELSEIF(DELTA.LE.0.005) THEN                                   
        MP1=1000                                                   
      ELSEIF(DELTA.LE.0.010) THEN                                   
        MP1=500                                                    
      ELSEIF(DELTA.LE.0.015) THEN                                   
        MP1=300                                                    
      ELSEIF(DELTA.LE.0.020) THEN                                   
        MP1=200                                                    
      ELSEIF(DELTA.LE.0.025) THEN                                   
        MP1=150                                                    
      ELSEIF(DELTA.LE.0.040) THEN                                   
        MP1=100                                                    
      ELSEIF(DELTA.LE.0.050) THEN                                   
        MP1=70                                                     
      ELSEIF(DELTA.LE.0.080) THEN                                   
        MP1=60                                                     
      ELSEIF(DELTA.LE.0.150) THEN                                   
        MP1=40                                                     
      ELSEIF(DELTA.LE.0.200) THEN                                   
        MP1=30
      ELSEIF(DELTA.LE.0.400) THEN                                   
        MP1=20                                                     
      ELSEIF(DELTA.LE.1.000) THEN                                   
        MP1=15                                                     
      ELSE                                                        
        MP1=10                                                     
      ENDIF                                                       
      DMEAN(1)=0.0D0                                         
      IF(MP1.GE.2) THEN
        DO 10 J=2,MP1 
          JM1=J-1
          DMEAN(J)=DMEAN(JM1)+DELTA 
   10   CONTINUE
      ENDIF
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PDFN(X,XMU,SIGMA)
C
C        ALGORITHM AS 305.5 APPL. STATIST. (1996), VOL.45, NO.4
C
C        Compute the normal density function at X where XMU is the
C        mean and SIGMA is the sandard deviation
C
      DOUBLE PRECISION X, XMU, SIGMA
C
      DOUBLE PRECISION ARG
C
      ARG=-0.5D0*(X-XMU)*(X-XMU)/(SIGMA*SIGMA)                       
      PDFN=(3.989422804014327D-1/SIGMA)*DEXP(ARG) 
C      
      RETURN                                                            
      END
