      SUBROUTINE SIMDO(QIND, QFOR, IPROD, KDIM, JSUB, IVEC, IFAULT)
C
C     ALGORITHM AS 172  APPL. STATIST. (1982) VOL.31, NO.1
C
C     Generates the subscript vector (index subscript) given the index
C     subscript (subscript vector) for a sequence of Fortran DO-loops
C     nested to depth KDIM
C
      LOGICAL QIND, QFOR
      INTEGER IPROD(KDIM), IVEC(KDIM)
C
      IFAULT = 0
      IF (.NOT. QIND) GO TO 12
C
C     Index subscript to subscript vector conversion
C
      IF (JSUB .LE. IPROD(KDIM)) GO TO 5
      IFAULT = 1
      RETURN
    5 ITEMPV = JSUB - 1
      IJ = KDIM - 1
      DO 10 I = 1, IJ
	IK = KDIM - I
	IVEC(I) = ITEMPV / IPROD(IK)
	ITEMPV = ITEMPV - IPROD(IK) * IVEC(I)
	IVEC(I) = IVEC(I) + 1
   10 CONTINUE
      IVEC(KDIM) = ITEMPV + 1
      IF (QFOR) CALL REVERS(IVEC, KDIM)
      RETURN
C
C     Subscript vector to index subscript conversion
C
   12 IF (IVEC(1) .LE. IPROD(1)) GO TO 14
      IFAULT = 2
      RETURN
   14 DO 15 I = 2, KDIM
	IF (IVEC(I) .LE. IPROD(I)/IPROD(I-1)) GO TO 15
	IFAULT = 2
	RETURN
   15 CONTINUE
      IF (.NOT. QFOR) CALL REVERS(IVEC, KDIM)
      JSUB = IVEC(1)
      DO 20 I = 2, KDIM
	JSUB = JSUB + (IVEC(I) - 1) * IPROD(I-1)
   20 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE REVERS(ITAB, IDIM)
      INTEGER ITAB(IDIM)
C
C     ALGORITHM AS 172.1  APPL. STATIST. (1982) VOL.31, NO. 1
C
C     Reorders subscript vector, if required
C
      ITER = IDIM / 2
      K = IDIM + 1
      DO 10 I = 1, ITER
	ITEMP = ITAB(I)
	IK = K - I
	ITAB(I) = ITAB(IK)
	ITAB(IK) = ITEMP
   10 CONTINUE
C
      RETURN
      END
