      SUBROUTINE NEXTGJ(A, IA, K, IWK)
C
C     A FORTRAN translation of Garside's algorithm AS38, but without
C     checks and with only the upper triangle of the SSP-matrix stored.
C     Array IWK must contain (K-1) zeroes initially.
C
      DOUBLE PRECISION A(IA)
      INTEGER IWK(K)
C
C     Local variables
C
      DOUBLE PRECISION ONE, AA, AIP
      DATA ONE/1.D0/
C
      IP = K-1
C
C     IREM = 1 if a variable is being removed from the model, 
C          = 0 if it is being added.
C
   10 IREM = IWK(IP)
C
C     Pivot variable no. IP in or out of the model.
C
      IPP = (IP-1)*(K+K+2-IP)/2 + 1
      IJ = IPP+K+1-IP
      AA = -A(IPP)
      A(IPP) = AA
      AA = ONE/AA
      IP1 = IP+1
      DO 30 I = IP1, K
	IPP = IPP+1
	AIP = A(IPP)*AA
	IPJ = IPP
	DO 20 J = I, K
	  A(IJ) = A(IJ) + AIP*A(IPJ)
	  IJ = IJ+1
	  IPJ = IPJ+1
   20   CONTINUE
   30 CONTINUE
C
C     Change the indicator in array IWK for variable IP.
C     If variable was deleted, move onto the next lower numbered
C     variable, otherwise return.
C
      IWK(IP) = 1-IREM
      IF(IREM.EQ.0) RETURN
      IP = IP-1
      IF(IP.GT.0) GO TO 10
C
      RETURN
      END
