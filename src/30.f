      SUBROUTINE HNPLOT (NTITLE,OBS,N)
C 
C        Algorithm as 30 j.r.statist.soc. c. (1970) vol.19. no.2.
C 
C        Half-normal plotting
C 
      DIMENSION OBS(*), A(6), B(4), P(100), NTITLE(20)
      DIMENSION OUT(101), XPR(11), D(100), YPR(6)
      CHARACTER*1 OUT,DOT,BLANK,PLUS
 10   FORMAT (1H1,20A4,///)
 20   FORMAT (1H ,F11.4,101A1)
 30   FORMAT (1H ,11X,1H.)
 40   FORMAT (1H ,F11.4,101H............................................
     1.........................................................)
 50   FORMAT (1H ,11X,101H.      .       .       .          .         . 
     1         .      .                    .                .)
 60   FORMAT (1H ,8X,F6.2,F7.2,F8.2,F8.2,F11.2,F10.2,F11.2,F7.2,16X,
     1 F5.2,12X,F5.2)
 70   FORMAT (1H ,11X,101A1)
 80   FORMAT (1H ,F11.4,1H.)
C 
C        Set iwrite equal to the output device number (system dependent)
C 
      DATA IWRITE/6/
C 
C     set constants for probability scale, and for printing
C 
      DATA A(1),A(2),A(3),A(4),A(5),A(6)/1.048,-0.8559,0.363,0.108392,
     1 0.328117,1.253314/
      DATA B(1),B(2),B(3),B(4)/-0.001416,-0.039811,-0.6256,0.401703/
      DATA XPR(1),XPR(2),XPR(3),XPR(4),XPR(5),XPR(6),XPR(7),XPR(8),
     1XPR(9),XPR(10)/50.0,60.0,70.0,80.0,90.0,95.0,98.0,99.0,99.9,99.99/
      DATA DOT,BLANK,PLUS/'.',' ','+'/
C 
      OUT(1)=DOT
C 
C        Sort data into ascending order
C 
      DO 100 I=1,N
        DO 100 J=I,N
        IF (OBS(J)-OBS(I)) 90,100,100
 90     TEMP=OBS(I)
        OBS(I)=OBS(J)
        OBS(J)=TEMP
 100  CONTINUE
C 
C     Calculate scales for the axes (the x-axis is assumed to be
C     100 units long, and the y-axis 50 units long
C 
      XSCALE=0.03719
      OBSMAX=OBS(N)
      KOUNT=0
 110  IF (OBSMAX.GT.100.0) GO TO 120
      IF (OBSMAX.GE.10.0) GO TO 130
      OBSMAX=OBSMAX*10.0
      KOUNT=KOUNT-1
      GO TO 110
 120  OBSMAX=OBSMAX/10.0
      KOUNT=KOUNT+1
      GO TO 110
 130  IF (AMOD(OBSMAX,5.0).EQ.0) GO TO 140
      OBSMAX=(AINT(OBSMAX/5.0)+1.0)*5.0
 140  OBSMAX=OBSMAX*(10.0**KOUNT)
      YSCALE=OBSMAX/50.0
      YPR(1)=OBSMAX
      YPR(6)=0.0
      OBSMAX=OBSMAX/5.0
      DO 150 I=2,5
 150    YPR(I)=YPR(I-1)-OBSMAX
C 
C        Calculate positioning of points on probability scale
C 
      X=50.0/(2.0*FLOAT(N))
      P(1)=50.0+X
      DO 160 K=2,N
 160    P(K)=P(K-1)+2.0*X
      DO 200 K=1,N
        W=0.02*P(K)-1.0
        IF (P(K).GT.95.5) GO TO 180
        W2=W*W
        D(K)=A(1)
        DO 170 I=2,6
 170      D(K)=D(K)*W2+A(I)
        D(K)=D(K)*W
        GO TO 200
 180    Z=LOG(1.0-W)
        D(K)=B(1)
        DO 190 I=2,4
 190      D(K)=D(K)*Z+B(I)
 200  CONTINUE
      WRITE (IWRITE,10) (NTITLE(I),I=1,20)
      L=N
      QR=YSCALE/1000.0
C 
C        Print output line by line, except for the last one, with the
C        scale value printed on every tenth line
C 
      DO 280 I=1,50
C 
C        Calculate position on vertical scale, and which points - of
C        those not yet printed - exceed this value
C 
        PR=YPR(1)-FLOAT(I-1)*YSCALE-QR
        IF (OBS(L).LT.PR) GO TO 260
C 
C        Set output array to blanks
C 
        DO 210 IX=2,101
 210      OUT(IX)=BLANK
C 
C        Calculate position of point on probability scale and put print
C        symbol in appropriate array element.  Repeat for further points
C        at this value, and print line of output
C 
 220    JP=(D(L)/XSCALE)+1.0
        OUT(JP)=PLUS
        L=L-1
        IF (L.EQ.0) GO TO 230
        IF (OBS(L).GE.PR) GO TO 220
 230    IF (MOD(I-1,10).NE.0) GO TO 240
        I1=((I-1)/10)+1
        WRITE (IWRITE,20) YPR(I1),(OUT(I2),I2=1,101)
        GO TO 250
 240    WRITE (IWRITE,70) (OUT(I2),I2=1,101)
C 
C        If all points have been printed, or points remain after all
C        lines have been printed (except the last) leave the main loop
C        and complete the printing elsewhere
C 
 250    IF (L.GT.0) GO TO 280
        IF (I.GE.50) GO TO 320
        IL=I+1
        GO TO 290
 260    IF (MOD(I-1,10).NE.0) GO TO 270
        I1=((I-1)/10)+1
        WRITE (IWRITE,80) YPR(I1)
        GO TO 280
 270    WRITE (IWRITE,30)
 280  CONTINUE
      GO TO 320
C 
C        Complete printing when all points have already been output
C 
 290  DO 310 I=IL,50
        IF (MOD(I-1,10).NE.0) GO TO 300
        I1=((I-1)/10)+1
        WRITE (IWRITE,80) YPR(I1)
        GO TO 310
 300    WRITE (IWRITE,30)
 310  CONTINUE
 320  IF (L.LE.0) GO TO 350
C 
C        Print last line, including the probability axis and any
C        remaining points
C 
      DO 330 IX=2,101
 330    OUT(IX)=OUT(1)
 340  JP=(D(L)/XSCALE)+1.0
      OUT(JP)=PLUS
      L=L-1
      IF (L.GT.0) GO TO 340
      WRITE (IWRITE,20) YPR(6),(OUT(I2),I2=1,101)
      GO TO 360
 350  WRITE (IWRITE,40) YPR(6)
C 
C        Print probability scale
C 
 360  WRITE (IWRITE,50)
      WRITE (IWRITE,60) (XPR(IP),IP=1,10)
      RETURN
      END
