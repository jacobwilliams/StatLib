      SUBROUTINE BIB(AMAT,A,B,BMAT,C,D,NEWMAT,FLAG)
C      
C        ALGORITHM AS 306 APPL.STATIST. (1996), VOL.45, NO.4
C
C        Calculation of a BIB product operation on any two matrices
C
      INTEGER A, B, C, D
      REAL AMAT(A,B), BMAT(C,D), NEWMAT(A,B*D)
      LOGICAL FLAG
C
      INTEGER CT, I, J, K, PT
C
C     This block of code sets the logical variable FLAG to F and exits
C     the subroutine if the number of non-zero elements in each column
C     of A does not equal the number of rows of B
C      
      DO 20 I = 1, B
         CT=0
         DO 10 J=1, A
            IF (AMAT(J,I).NE.0.0) CT = CT + 1
   10    CONTINUE
         FLAG=CT.EQ.C
         IF(.NOT.FLAG) GOTO 70
   20 CONTINUE
C
C        This block of code computes the BIB product
C
      DO 60 I = 1, B
         PT = 1
         DO 50 J = 1, A
            IF(AMAT(J,I).EQ.0.0) THEN
               DO 30 K = 1, D
                 NEWMAT(J,(I*D)-D+K)=0.0
   30          CONTINUE
            ELSE
               DO 40 K = 1, D
                  NEWMAT(J,(I*D)-D+K)=AMAT(J,I)*BMAT(PT,K)
   40          CONTINUE
               PT=PT+1
            END IF
   50    CONTINUE
   60 CONTINUE
   70 RETURN
      END
