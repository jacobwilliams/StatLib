      SUBROUTINE ULTRA (BW,M,MORD,MUE,N,NUE,OPTIO1,OPTIO2,OPTIO3,
     *                  XIN,XOU,YIN,IWKAR,WKAR,IFAULT,W2,YOU)
C
C        ALGORITHM AS297.1 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Nonparametric regression using ultraspherical
C        polynomials
C
      INTEGER I0, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12,
     *        I13, II1, II2, IFLG, IFLG1, IL, ILB, IR, IRB, IZ, M,
     *        IWKAR(2*M), MORD, N, IFAULT, NUE, OPTIO1, OPTIO2, OPTIO3
      DOUBLE PRECISION BW, PHIJU, SUMW2, XA, XAMB, XAMBN, XAPB, XB, XC,
     *                 XX, XIN(N), XOU(M), W2(M),YIN(N), YOU(M), YOUT,
     *                 YOUT1, WKAR(60*N+2000)
C
      EXTERNAL BOUND, ERROR, INIT, PHIJ, SAVE, WEIGHT
C
C        Error checks
C
      CALL ERROR (BW,M,MORD,MUE,N,NUE,OPTIO1,OPTIO2,XIN,XOU,IFAULT)
      IF (IFAULT.EQ.0.OR.IFAULT.EQ.7) GO TO 5
      GO TO 200
C
C        Convert XIN(N) and XOU(M) to interval [-1,1]
C
    5 XA = XIN(1)-(XIN(2)-XIN(1))/2.0D0
      XB = XIN(N)+(XIN(N)-XIN(N-1))/2.0D0
      XAMB = XB-XA
      XAPB = XA+XB 
      DO 10 I=1, N
         XIN(I) = 2.0D0*XIN(I)/XAMB-XAPB/XAMB
   10 CONTINUE
      DO 15 I=1, M
         XOU(I) = 2.0D0*XOU(I)/XAMB-XAPB/XAMB
   15 CONTINUE   
C
C        Store boundaries of windows and corresponding indices
C   
      I0 = 1
      I1 = I0+N
      I2 = I1+2*M
      CALL BOUND (BW,M,N,OPTIO1,OPTIO2,XIN,XOU,IWKAR(I0),
     *            ILB,IRB,WKAR(I1),WKAR(I2))    
C
C        Calculate the initial values
C
      I3 = I2+N+1
      I4 = I3+MUE+1
      I5 = I4+(MUE+1)*(MUE+1)
      I6 = I5+(MORD+1)*(MORD+1)
      I7 = I6+MORD+1
      CALL INIT(WKAR(I3),BW,MORD,MUE,WKAR(I4),WKAR(I5),WKAR(I6),
     *          WKAR(I7))
      I8 = I7+MUE+1
      I9 = I8+N+1
      I10 = I9+MORD+1
      I11 = I10+(N+1)*(MORD+2*MUE+1)
      I12 = I11+MORD+1
      I13 = I12+N*(MORD+1)
C
C        Calculations for the:
C           1 - left boundary (III = 1)
C           2 - right boundary (III = 2)
C           3 - interior (III = 3)
C 
      DO  120 III=1, 3
         IF (III.EQ.1.AND.ILB.LE.0) GO TO 120
         IF (III.EQ.2.AND.IRB.LE.0) GO TO 120
         IF (III.EQ.3.AND.ILB.GT.IRB-2) GO TO 120
         IF (III.EQ.3) GO TO 30
         IL = IWKAR(1+(M-1)*(III-1))
         IR = IWKAR(M+1+(M-1)*(III-1))
         IZ = IR-IL
         IF (IZ.LE.0) GO TO 120  
         DO 20 I=IL, IR
            WKAR(I8-IL+I) = WKAR(I2+I-1)
   20    CONTINUE
         WKAR(I8) = WKAR(N+1+(M-1)*(III-1))
         WKAR(I8+IZ) = WKAR(N+M+1+(M-1)*(III-1))
         IF (III.EQ.2) GO TO 25
         II1 = 1
         II2 = ILB
         IFLG = 0
         IFLG1 = 0
         XC = WKAR(N+1)+BW
         GO TO 40
   25    IF (ILB.EQ.IRB) IRB = IRB +1
         II1 = IRB
         II2 = M
         IFLG = 0
         IFLG1 = 0
         XC = WKAR(N+2*M)-BW
         GO TO 40
   30    IF ((ILB+1).EQ.IRB) GO TO 120
         II1 = ILB+1
         II2 = IRB -1 
   40    DO 110 II=II1, II2
            XX = XOU(II)
            IF (III.NE.3) GO TO 55
            XC = XX
            IL = IWKAR(II)
            IR = IWKAR(M+II)
            IZ = IR-IL
            IF (IZ.LE.0) GO TO 110
            DO 50 I=IL, IR
               WKAR(I8-IL+I) = WKAR(I2+I-1)
   50       CONTINUE
            WKAR(I8) = WKAR(N+II)
            WKAR(I8+IZ) = WKAR(N+M+II)
            GO TO 60
   55       IF (IFLG.EQ.1) GO TO 65 
   60       CALL SAVE (BW,WKAR(I5),WKAR(I6),IZ,MORD,MUE,WKAR(I7),
     *                 N,WKAR(I8),XC,WKAR(I10),WKAR(I12))
            IFLG = 1
C
C        Calculate array PHIAR(MORD+1)
C
   65       DO 70 JJ=0, MORD
               CALL PHIJ (BW,WKAR(I5),JJ,WKAR(I6),MORD,NUE,XC,XX,
     *                    WKAR(I9),PHIJU)
               WKAR(I11+JJ) = PHIJU
   70       CONTINUE
            YOUT = 0.0D0
            IF (III.EQ.3.AND.OPTIO1.EQ.1) GO TO 130
            IF (OPTIO3.EQ.0) GO TO 80
C
C        Sum of squared weights is requested (OPTIO3 = 1)
C
            CALL WEIGHT (IZ,MORD,N,NUE,WKAR(I11),WKAR(I12),SUMW2,
     *                   WKAR(I0))
            DO 75 I=1, IZ
               YOUT = YOUT+WKAR(I0+I-1)*YIN(IL+I-1)
   75       CONTINUE
            YOU(II) = YOUT
            W2(II) = SUMW2
            GO TO 110
C
C        Sum of squared weights not requested (OPTIO3 = 0)
C
   80       IF (IFLG1.EQ.1.AND.III.EQ.1) GO TO 100
            IF (IFLG1.EQ.1.AND.III.EQ.2) GO TO 100
            DO 95 J=NUE, MORD   
               YOUT1 = 0.0D0
               DO 90 I=1, IZ
                  YOUT1 = YOUT1+YIN(IL+I-1)*WKAR(I12+I-1+N*J)
   90          CONTINUE
               WKAR(I13+J) = YOUT1
   95       CONTINUE
            IFLG1 = 1
  100       DO 105 J=NUE, MORD
               YOUT = YOUT+WKAR(I13+J)*WKAR(I11+J)
  105       CONTINUE
            YOU(II) = YOUT
  110    CONTINUE 
  120 CONTINUE
      GO TO 160  
C
C        Calculation for the interior in the equidistant case
C
  130 CALL WEIGHT (IZ,MORD,N,NUE,WKAR(I11),WKAR(I12),SUMW2,
     *             WKAR(I0))
      DO 150 II=II1, II2
         YOUT = 0.0D0
         IL = IWKAR(II)
         DO 140 I=1, IZ  
            YOUT = YOUT+WKAR(I0+I-1)*YIN(IL+I-1)
  140    CONTINUE
         YOU(II) = YOUT
         W2(II) = SUMW2
  150 CONTINUE
C  
C        Convert XIN(N), XOU(M) and YOU(M) to the original scales
C
  160 CONTINUE 
      DO 170 I=1, N
         XIN(I) = (XAMB*XIN(I)+XAPB)/2.0D0
  170 CONTINUE
      XAMBN = 1.0D0
      IF (NUE.EQ.0) GO TO 185
      DO 180 I=1, NUE
         XAMBN = (2.0D0/XAMB)*XAMBN
  180 CONTINUE
  185 DO 190 I=1, M
         XOU(I) = (XAMB*XOU(I)+XAPB)/2.0D0
         YOU(I) = XAMBN*YOU(I)
  190 CONTINUE
  200 CONTINUE
       RETURN
       END       
C
      SUBROUTINE ERROR (BW,M,MORD,MUE,N,NUE,OPTIO1,OPTIO2,XIN,
     *                  XOU,IFAULT)
C
C        ALGORITHM AS297.2 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Error checking
C
      INTEGER M, MO, MORD, MUE, N, IFAULT, NU, NUE, OPTIO1, OPTIO2
      DOUBLE PRECISION BW, ERR, X1, X2, XIN(N), XOU(M)
C
      IFAULT = 0
C
C        Error code messages
C
      IF (N.GE.4.AND.M.GE.1) GO TO 10
      IFAULT = 1
   10 IF (MORD.LE.30.AND.MUE.LE.10.AND.NUE.LE.3.AND.NUE.LE.MORD)
     *   GO TO 20
      IFAULT = 2
   20 DO 30 I=1, N-1
         IF (XIN(I).LE.XIN(I+1)) GO TO 30
         IFAULT = 3
   30 CONTINUE
      DO 40 I=1, M-1
         IF (XOU(I).LE.XOU(I+1)) GO TO 40
         IFAULT = 4
   40 CONTINUE    
      IF (BW.LE.0.AND.OPTIO2.NE.1)  IFAULT=5
      IF (OPTIO1.NE.1) GO TO 70
      ERR = 0.00001D0
      DO 50 I=1, N
         IF (DABS(XOU(I)-XIN(I)).LT.ERR) GO TO 50
         IFAULT = 6
   50 CONTINUE
      DO 60 I=1, N-2        
         X1 = XIN(I+1)-XIN(I)
         X2 = XIN(I+2)-XIN(I+1)
         IF (DABS(X1-X2).LT.ERR) GO TO 60
         IFAULT = 6
   60 CONTINUE
   70 IF (IFAULT.NE.0) GO TO 80
C
C        Warning message (error code 7)
C
      MO = (MORD/2)*2
      NU = (NUE/2)*2
      IF (MO.EQ.MORD.AND.NU.EQ.NUE.OR.MO.NE.MORD.AND.NU.NE.NUE)
     *   GO TO 80
      IFAULT = 7
   80 CONTINUE
      RETURN
      END
C
      SUBROUTINE BOUND (BW,M,N,OPTIO1,OPTIO2,XIN,XOU,IBND,ILB,
     *                  IRB,BND,S)
C
C        ALGORITHM AS297.3 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Calculate left and right boundaries and indices of windows
C        midpoints
C
      INTEGER IL, IR, IZ, ILB, IRB, M,
     *        IBND(M,2), N, OPTIO1, OPTIO2
      DOUBLE PRECISION B1, B11, B2, B22, BND(M,2), BW, S(N+1),
     *                 XIN(N), XOU(M)
C
      S(1) = XIN(1)-(XIN(2)-XIN(1))/2.0D0
      S(N+1) = XIN(N)+(XIN(N)-XIN(N-1))/2.0D0
      DO 10 I=1, N-1
         S(I+1) = (XIN(I+1)+XIN(I))/2.0D0
   10 CONTINUE  
      IF (OPTIO2.EQ.1) BW = (S(N+1)-S(1))/2.0D0
      B1 = S(1)+BW
      B11 = B1+BW
      B2 = S(N+1)-BW
      B22 = B2-BW 
C
C        Find the indices of the boundaries for the output points
C
       IF (OPTIO1.EQ.0) GO TO 15
C
C        Equidistant case
C        
      IZ = BW/(XIN(2)-XIN(1))
      ILB = 1+IZ
      IRB = N-IZ
      GO TO 25
C
C        Non-equidistant case
C          
   15 DO 20 I=1, M
         IF (XOU(I).LE.B1) ILB = I
         IF (XOU(M-I+1).GE.B2) IRB = M-I+1
   20 CONTINUE
C
C        Find the indices of the boundaries for the input points
C  
       IF (OPTIO1.EQ.0) GO TO 30
C
C        Equidistant case
C
   25 IZ = INT(2*BW/(XIN(2)-XIN(1)))
      IL = 1+IZ
      IR = N+1-IZ
      GO TO 40
C  
C        Non-equidistant case
C 
   30 DO 35 I=1, N
         IF (S(I).LT.B11)    IL = I
         IF (S(N-I+2).GT.B22) IR = N-I+2
   35 CONTINUE
      IL = IL+1
      IF (IR.GT.1) IR = IR-1
      IZ = IR-IL
C
C        Record the boundaries and boundary indices
C
   40 IF (XOU(1).GT.B1) GO TO 45
      IBND(1,1) = 1
      IBND(1,2) = IL
      BND(1,1) = S(1)
      BND(1,2) = B11
   45 IF (XOU(M).LT.B2) GO TO 50
      IBND(M,1) = IR
      IBND(M,2) = N+1
      BND(1,1) = S(1)
      BND(1,2) = B11
      BND(M,1) = B22
      BND(M,2) = S(N+1) 
   50 IF (OPTIO2.EQ.1) GO TO 110
      IF (B1.EQ.B2) GO TO 100
      IF (XOU(M).LT.B2) IRB = M+1
      IF (ILB.GT.IRB-2) GOTO 110
      DO 90 I=ILB+1, IRB-1
         XL = XOU(I)-BW
         XR = XOU(I)+BW
C
C        Computation of relevant indices
C
         IF (OPTIO1.EQ.0) GO TO 55
C
C        Equidistant case
C 
         IZ = INT(BW/(XIN(2)-XIN(1)))
         IL = I-IZ
         IR = I+IZ+1
         GO TO 85
C
C        Non-equidistant case
C
   55    IL = 1
   60    IR = IL
         IF (S(IL).GT.XL) GO TO 65
         IL = IL+1
         GO TO 60
   65    IR = IL-1
   70    IR = IR+1
         IF (IR.EQ.N+1) GO TO 80
         IF (S(IR).LT.XR) GO TO 70
   80    IF (IL.GT.1)  IL = IL-1
         IZ = IR-IL
C
C        Record the indices
C    
   85    IBND(I,1) = IL
         IBND(I,2) = IR
         BND(I,1) = XL
         BND(I,2) = XR
   90 CONTINUE 
      GO TO 110 
  100 OPTIO2 = 1
  110 RETURN
      END
C       
      SUBROUTINE INIT (BIN,BW,MORD,MUE,WKM1,COEF,HJAR,M1BIN)
C
C        ALGORITHM AS297.4 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Computes necessary constants
C
      INTEGER JM, MORD, MUE
      DOUBLE PRECISION AA, BB, BCD, BIN(MUE+1), BW, CC,
     *                 COEF(MORD+1,MORD+1), DD, EE, FF, HJAR(MORD+1),
     *                 HJMU, M1BIN(MUE+1), MUP1, POW2, RR, 
     *                 WKM1(MUE+1,MUE+1)
C
C        Calculate multiple of binomial coefficients
C
      EXTERNAL BINOM, HJ, RGAMMA
      INTRINSIC DSQRT
C
      CALL BINOM (MUE,WKM1,BIN)
      DO 5 K=0, MUE
         M1BIN(K+1) = ((-1)**K)*BIN(K+1)
    5 CONTINUE
C
C        Calculate HJAR(MORD+1)
C  
      MUP1 = (1.0D0)*MUE+1.0D0
      POW2 = 1.0D0
      DO 10 I=1,2*MUE+1
         POW2 = 2.0D0*POW2
   10 CONTINUE
      POW2 = DSQRT(POW2)
      CALL RGAMMA (2*MUP1-1,MUP1,RR)
      DO 15 JJ=0, MORD
         CALL HJ (JJ,MUE,HJMU)
         HJAR(JJ+1) = HJMU*RR/POW2
   15 CONTINUE
C
C        Calculate COEF(MORD+1,MORD+1)
C
      DO 80 JJ=0, MORD
         DO 70 MM=0, JJ/2
            JM = JJ-2*MM
            IF (JJ.GE.3) GO TO 85
            BB = 1.0D0
            IF (MM.EQ.0) GO TO 30
            DO 20 I=1, MM
               BB = BB*(2.0D0*I-1)/(2.0D0*I)
   20       CONTINUE
   30       CC = 1.0D0
            DO 40 UU=1, MUE
               CC = CC*(2*MM+2*UU-1.0D0)/(2*UU-1.0D0)
   40       CONTINUE
            DD = 1.0D0
            IF (JM.EQ.0) GO TO 60
            DO 50 K=1, JM
               DD = DD*(2*MM+2*MUE+2*K-1.0D0)/FLOAT(K)
   50       CONTINUE
   60       DD = ((-1.0D0)**MM)*DD
            BCD = BB*CC*DD
            AA = 1.0D0
            COEF(JJ+1,MM+1) = AA*BCD*DSQRT(BW)
   70    CONTINUE
   80 CONTINUE
   85 IF (MORD.LT.3) GO TO 130
      DO 120 JJ=3, MORD
         DO 110 MM=0, JJ/2
            JM = JJ-2*MM
            EE = COEF(JJ,MM+1)
            IF (MM.EQ.0) GO TO 90
            FF = COEF(JJ-1,MM)
            GO TO 100
   90       FF = 0.0D0
  100       BCD = 2*((JJ+MUE-0.5D0)/FLOAT(JJ))*EE
     *            -((JJ+2*MUE-1.0D0)/FLOAT(JJ))*FF
            COEF (JJ+1,MM+1) = BCD
  110    CONTINUE
  120 CONTINUE
  130 RETURN
      END
C
      SUBROUTINE SAVE (BW,COEF,HJAR,IZ,MORD,MUE,M1BIN,N,S1,
     *                 XC,WKM2,RINT) 
C  
C        ALGORITHM AS297.5 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Saves integrals for the boundaries given in array S1(IZ)
C
      INTEGER  IZ, JMK1, N, MORD, MO2, MUE
      DOUBLE PRECISION AA, BB, BW, CC, COEF(MORD+1,MORD+1),
     *                 HJAR(MORD+1), M1BIN(MUE+1), RINT(N,MORD+1),
     *                 S1(N+1), SUM, XC, WKM2(N+1,MORD+2*MUE+1)
C
C        Calculate the powers
C   
      MO2 = MORD+2*MUE+1
      DO 30 I=1, IZ+1
         WKM2(I,1) = (S1(I)-XC)/BW
         IF (WKM2(I,1).GT.1.0D0) WKM2(I,1) = 1.0D0
         IF (WKM2(I,1).LT.-1.0D0) WKM2(I,1) = -1.0D0     
         DO 20 J=2, MO2
            WKM2(I,J) = WKM2(I,1)*WKM2(I,J-1)
   20    CONTINUE
   30 CONTINUE
C
C        Calculate the integral
C
      DO 70 II=1, IZ
         DO 60 JJ=0, MORD
            AA = HJAR(JJ+1)
            SUM = 0.0D0
            DO 50 MM=0, JJ/2
               BB = COEF(JJ+1,MM+1)
               CC = 0.0D0
               DO 40 K=0, MUE
                  JMK1 = JJ-2*MM+2*K+1
                  CC = CC+(M1BIN(K+1)*(WKM2(II+1,JMK1)
     *                     -WKM2(II,JMK1)))/FLOAT(JMK1)
   40          CONTINUE
               SUM = SUM +BB*CC
   50       CONTINUE
            RINT(II,JJ+1) = AA*SUM
   60    CONTINUE
   70 CONTINUE
      RETURN
      END
C
      SUBROUTINE WEIGHT (IZ,MORD,N,NUE,PHIAR,RINT,SUMW2,W)
C
C        ALGORITHM AS297.6 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Calculate weights W(N) corresponding to estimating NUE-th
C        derivative at point XX using ultraspherical polynomial of
C        order MORD on the interval [XX-BW,XX+BW]
C
      INTEGER IZ, MORD, N, NUE
      DOUBLE PRECISION PHIAR(MORD+1), RINT(N,MORD+1), SUM,
     *                 SUMW2, W(N)
C
      SUMW2 = 0
      DO 10 I=1, N
         W(I) = 0
   10 CONTINUE
      DO 30 I=1, IZ
         SUM = 0.0D0
         DO 20 J=NUE, MORD
            SUM = SUM +PHIAR(J+1)*RINT(I,J+1)
   20    CONTINUE
         W(I) = SUM
         SUMW2 = SUMW2+SUM*SUM
   30 CONTINUE
      RETURN
      END
C
      SUBROUTINE PHIJ (BW,COEF,JJ,HJAR,MORD,NUE,XC,XX,WKAR1,PHIJU)
C
C        ALGORITHM AS297.7 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Calculate NUE-th derivative of the JJ-th ultraspherical
C        polynomial with degree of smoothness MUE, on the interval
C        [XC-BW,XC+BW], evaluated at point XX
C
      INTEGER JJ, JM, MORD, MM, NUE
      DOUBLE PRECISION AA, BB, BW, BWN, CC, CENX, COEF(MORD+1,MORD+1),
     *                 HJAR(MORD+1), HJMU, PHIJU, SUM, XC, XX,
     *                 WKAR1(MORD+1)
C
      CENX = (XX-XC)/BW
      WKAR1(1) = 1.0D0
      DO 10 I=2, JJ+1
         WKAR1(I) = CENX*WKAR1(I-1)
   10 CONTINUE
      BWN = 1.0D0
      DO 20 I=0, NUE
         BWN = BW*BWN
   20 CONTINUE
      SUM = 0.0D0
      DO 50 MM=0, (JJ-NUE)/2
         JM = JJ-2*MM
         BB = COEF(JJ+1,MM+1)
         IF (NUE.EQ.0) GO TO 40
         DO 30 KK=1, NUE
            BB = BB*(JM-KK+1)
   30    CONTINUE 
   40    CC = WKAR1(JM-NUE+1)
         AA = BB*CC
         SUM = SUM+AA
   50 CONTINUE
      HJMU = HJAR(JJ+1)
      PHIJU = HJMU*SUM/BWN
      RETURN
      END
C
       SUBROUTINE HJ(JJ,MUE,HJMU)
C
C        ALGORITHM AS297.8 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Calculates part of the normalizing constant for
C        ultraspherical polynomials of order JJ, constant for
C        JJ-th polynomial
C
      INTEGER JJ, MUE
      DOUBLE PRECISION AA, BB, JJP1, JP1MU, HJMU
C
      EXTERNAL RGAMMA
      INTRINSIC DSQRT
C
      JJP1 = FLOAT(JJ+1)
      JP1MU = JJP1+2*MUE
      CALL RGAMMA (JJP1,JP1MU,BB)
      AA = JJ+JP1MU
      HJMU = DSQRT(AA*BB)
      RETURN
      END
C
      SUBROUTINE GAMMA (NN,GM)
C
C        ALGORITHM AS297.9 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Calculate gamma function at argument which is an
C        integer or integer + 0.5
C
      INTEGER N
      DOUBLE PRECISION ABSN, GM, PI, NN
C
      INTRINSIC DABS, DASIN, DSQRT
C
      PI = 2.0D0*DASIN(1.0D0)
      GM = 1.0D0
      IF (NN.EQ.0.5D0) GO TO 20
      N = NINT(NN+0.00001D0)
      DO 10 I=1, N-1
         GM = GM*(NN-I)
   10 CONTINUE
      ABSN = DABS(N-NN)
      IF (ABSN.GT.0.1D0) GO TO 20
      GO TO 30
   20 GM = GM*DSQRT(PI)
   30 RETURN
      END
C
      SUBROUTINE RGAMMA (MM,NN,RGM)
C
C        ALGORITHM AS297.10 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Calculate the ratio of two gamma functions
C
      DOUBLE PRECISION MM, NN, RGM
C
      RGM = 1.0D0
      IF (MM.EQ.NN) GO TO 40
      IF (MM.LT.NN) GO TO 20
      DO 10 I=1, NINT(MM-NN)
         RGM = RGM*(MM-I)
   10 CONTINUE
      GO TO 40
   20 DO 30 I=1, NINT(NN-MM)
         RGM = RGM*(NN-I)
   30 CONTINUE
      RGM = 1.0D0/RGM
   40 RETURN
      END
C 
      SUBROUTINE BINOM (P,WKM1,BIN)
C
C        ALGORITHM AS297.11 APPL. STATIST. (1995) VOL.44, NO.2
C
C        Binomial coefficients for power P using the method of
C        Pascal triangle
C     
      INTEGER P, P1
      DOUBLE PRECISION BIN(P+1), WKM1(P+1,P+1)
C
      P1 = P+1
      DO 20 I=1, P1
         DO 10 J=1, P1
            WKM1(I,J) = 1.0D0
   10    CONTINUE
   20 CONTINUE
      IF (P.LT.2) GO TO 50
      DO 40 I=3, P1
         DO 30 J=2, I-1
            WKM1(I,J) = WKM1(I-1,J-1)+WKM1(I-1,J)
   30    CONTINUE
   40 CONTINUE
   50 DO 60 J=1, P1
         BIN(J) = WKM1(P1,J)
   60 CONTINUE
      RETURN
      END
