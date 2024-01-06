      SUBROUTINE INVMOD(MAT,IMAT,RMOD,CMOD,IWK,NROW,IFAULT)
C
C        ALGORITHM AS 314.1 APPL. STATIST. (1997), VOL.46, NO.2
C
C        Inverts matrix with contents subject to modulo arithmetic
C
      INTEGER   IFAULT,NROW
      INTEGER   CMOD(NROW),IMAT(NROW*NROW),IWK(2*NROW),
     *          MAT(NROW*NROW),RMOD(NROW)
C
      INTEGER   I,IR,J,K,KIR,KJR,N
C
      EXTERNAL  MSORT,MUSORT
C
      INTRINSIC MOD
C
C        Check elements in 'mixed-moduli' positions are all zero
C
      N=0
      DO 26 I=1,NROW
        DO 24 J=1,NROW
          N=N+1
          IF ((RMOD(I).NE.CMOD(J)).AND.(MAT(N).GT.0)) GO TO 1002
          IF ((MAT(N).GT.RMOD(I)).OR.(MAT(N).LT.0)) GO TO 1001
          IMAT(N)=0
   24   CONTINUE
   26 CONTINUE
C
C        Sort rows and columns into ascending order of moduli
C
      CALL MSORT(MAT,IMAT,RMOD,CMOD,IWK,IWK(NROW+1),NROW)
C
C        Complete initialisation of inverse matrix 
C
      DO 28 N=1,NROW*NROW,NROW+1
        IMAT(N)=1
   28 CONTINUE
C
C        Invert the matrix
C
      DO 190 IR=1,NROW
        KIR=(IR-1)*NROW
        IF (MAT(KIR+IR).EQ.0) THEN
C
C        Find a row JR below IR such that K(JR,IR)>0
C
          DO 112 KJR=KIR+NROW+IR,NROW*NROW,NROW
            IF (MAT(KJR).GT.0) GO TO 115
  112     CONTINUE
C
C        Column IR contains all zeros in rows IR or below:
C        look for a row above with zeros to left of column IR 
C        and K(JR,IR)>0
C
          DO 114 KJR=IR,KIR,NROW
            IF (MAT(KJR).GT.0) THEN
              DO 113 I=KJR-IR+1,KJR-1
                IF (MAT(I).GT.0) GO TO 1003
  113         CONTINUE
              GO TO 115
            ENDIF
  114     CONTINUE
C
C        Column IR contains all zeros
C
          GO TO 190
C
C        Switch row JR with row IR
C
  115     KJR=KJR-IR
          DO 116 I=1,NROW
            K=MAT(KIR+I)
            MAT(KIR+I)=MAT(KJR+I)
            MAT(KJR+I)=K
            K=IMAT(KIR+I)
            IMAT(KIR+I)=IMAT(KJR+I)
            IMAT(KJR+I)=K
  116     CONTINUE
        END IF
C
C        Find multiplier N such that N*MAT(IR,IR)=1 mod(P{IR})
C
        K=MAT(KIR+IR)
        DO 122 N=1,RMOD(IR)-1
          IF (MOD(N*K,RMOD(IR)).EQ.1) GO TO 125
  122   CONTINUE
C
C        Multiply row IR by N
C
  125   IF (N.GT.1) THEN
          DO 126 I=KIR+1,IR*NROW
            MAT(I)=MAT(I)*N
            IMAT(I)=IMAT(I)*N
  126     CONTINUE
        END IF
C
C        Subtract MAT(JR,IR) * row IR from each row JR
C
        DO 136 KJR=0,NROW*NROW-1,NROW
          N=RMOD(IR)-MAT(KJR+IR)
          IF ((KJR.NE.KIR).AND.(N.NE.0)) THEN
            DO 132 I=1,NROW
              MAT(KJR+I)=MOD(MAT(KJR+I)+N*MAT(KIR+I),CMOD(I))
              IMAT(KJR+I)=MOD(IMAT(KJR+I)+N*IMAT(KIR+I),CMOD(I))
  132       CONTINUE
          END IF
  136   CONTINUE
  190 CONTINUE
C
C        Check inversion was possible - that result has
C        non-zero elements only on diagonal
C
      IFAULT=0
      DO 202 N=1,NROW*NROW,NROW+1
C
C        Zero diagonal element => left inverse formed
C
        IF (MAT(N).EQ.0) IFAULT=-1
        MAT(N)=-MAT(N)
  202 CONTINUE
      DO 204 N=1,NROW*NROW
        IF (MAT(N).GT.0) GO TO 1003
  204 CONTINUE
      DO 206 N=1,NROW*NROW,NROW+1
        MAT(N)=-MAT(N)
  206 CONTINUE
C
C        Sort rows and columns back into original order
C
      CALL MUSORT(MAT,IMAT,RMOD,CMOD,IWK,IWK(NROW+1),NROW)
C
      RETURN
C
C        Invalid entry (< 0 or > modulus)
C
 1001 IFAULT=1
      RETURN
C
C        Non-zero entry in position where row mod does not 
C        equal column mod
C
 1002 IFAULT=2
      RETURN
C
C        Matrix cannot be inverted
C
 1003 IFAULT=3
C      
      RETURN
      END
C
      SUBROUTINE MSORT(MAT,IMAT,RMOD,CMOD,RSORT,CSORT,NROW)
C
C        ALGORITHM AS 314.2 APPL. STATIST. (1997), VOL.46, NO.2
C
C        Sorts the matrix and associated information to put 
C        rows and columns in ascending order of moduli
C
      INTEGER   NROW
      INTEGER   CMOD(NROW),CSORT(NROW),IMAT(*),MAT(*),RMOD(NROW),
     *          RSORT(NROW)
C
      INTEGER   I,IRC,J,JRC,KIRC,KJRC,P
C
C        Initialise row and column addresses
C
      DO 22 I=1,NROW
        RSORT(I)=I
        CSORT(I)=I
   22 CONTINUE
C
C        Sort rows
C
      DO 48 IRC=1,NROW
C
C        Find next row
C
        JRC=IRC
        P=RMOD(IRC)
        DO 42 I=IRC+1,NROW
          IF (RMOD(I).LT.P) THEN
            P=RMOD(I)
            JRC=I
          END IF
   42   CONTINUE
        IF (IRC.NE.JRC) THEN
          I=RMOD(IRC)
          RMOD(IRC)=RMOD(JRC)
          RMOD(JRC)=I
          I=RSORT(IRC)
          RSORT(IRC)=RSORT(JRC)
          RSORT(JRC)=I
C
C        Switch rows
C
          KIRC=(IRC-1)*NROW
          KJRC=(JRC-1)*NROW
          DO 44 J=1,NROW
            I=MAT(KIRC+J)
            MAT(KIRC+J)=MAT(KJRC+J)
            MAT(KJRC+J)=I
   44     CONTINUE
        END IF
   48 CONTINUE
C
C        Sort columns
C
      DO 68 IRC=1,NROW
C
C        Find next column
C
        JRC=IRC
        P=CMOD(IRC)
        DO 62 I=IRC+1,NROW
          IF (CMOD(I).LT.P) THEN
            P=CMOD(I)
            JRC=I
          END IF
   62   CONTINUE
        IF (IRC.NE.JRC) THEN
          I=CMOD(IRC)
          CMOD(IRC)=CMOD(JRC)
          CMOD(JRC)=I
          I=CSORT(IRC)
          CSORT(IRC)=CSORT(JRC)
          CSORT(JRC)=I
C
C        Switch columns
C
          DO 64 J=0,NROW*NROW-1,NROW
            I=MAT(IRC+J)
            MAT(IRC+J)=MAT(JRC+J)
            MAT(JRC+J)=I
   64     CONTINUE
        END IF
   68 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MUSORT(MAT,IMAT,RMOD,CMOD,RSORT,CSORT,NROW)
C
C        ALGORITHM AS 314.3 APPL. STATIST. (1997), VOL.46, NO.2
C
C        Sorts the inverse matrix and associated information into the 
C        original order of rows and columns
C
      INTEGER   NROW
      INTEGER   CMOD(NROW),CSORT(NROW),IMAT(NROW*NROW),
     *          MAT(NROW*NROW),RMOD(NROW),RSORT(NROW)
C
      INTEGER   I,IRC,J,JRC,KIRC,KJRC
C
C        Sort rows of inverse (=columns of original)
C
      DO 48 IRC=1,NROW
C
C        Find next row
C
        IF (CSORT(IRC).NE.IRC) THEN
          DO 42 JRC=IRC+1,NROW
            IF (CSORT(JRC).EQ.IRC) GO TO 43
   42     CONTINUE
   43     I=CMOD(IRC)
          CMOD(IRC)=CMOD(JRC)
          CMOD(JRC)=I
          I=CSORT(IRC)
          CSORT(IRC)=CSORT(JRC)
          CSORT(JRC)=I
C
C        Switch rows
C
          KIRC=(IRC-1)*NROW
          KJRC=(JRC-1)*NROW
          DO 44 J=1,NROW
            I=IMAT(KIRC+J)
            IMAT(KIRC+J)=IMAT(KJRC+J)
            IMAT(KJRC+J)=I
   44     CONTINUE
        END IF
   48 CONTINUE
C
C        Sort columns of inverse (= rows of original)
C
      DO 68 IRC=1,NROW
C
C        Find next column
C
        IF (RSORT(IRC).NE.IRC) THEN
          DO 62 JRC=IRC+1,NROW
            IF (RSORT(JRC).EQ.IRC) GO TO 63
   62     CONTINUE
   63     I=RMOD(IRC)
          RMOD(IRC)=RMOD(JRC)
          RMOD(JRC)=I
          I=RSORT(IRC)
          RSORT(IRC)=RSORT(JRC)
          RSORT(JRC)=I
C
C        Switch columns
C
          DO 64 J=0,NROW*NROW-1,NROW
            I=IMAT(IRC+J)
            IMAT(IRC+J)=IMAT(JRC+J)
            IMAT(JRC+J)=I
   64     CONTINUE
C
C        Switch diagonal elements of MAT (others are zero)
C
          KIRC=(IRC-1)*NROW+IRC
          KJRC=(JRC-1)*NROW+JRC
          I=MAT(KIRC)
          MAT(KIRC)=MAT(KJRC)
          MAT(KJRC)=I
        END IF
   68 CONTINUE
C
      RETURN
      END
