      SUBROUTINE YATES (Y, LEVS, KODE, NFAC, X, EXTR)
C
C     ALGORITHM  AS 22 APPL.STATIST. (1969) VOL.18, NO.3
C
C     Yates algorithm for general complete factorial experiments.
C     Computes  Y = (IL..X..IR).Y
C
C     where   . denotes  ordinary matrix product
C            ..          direct matrix product
C            IL          unit matrix of order NLFT
C            IR          unit matrix of order NRGT
C           Y( )         multi-way table stored as a vector
C           X( )         a square matrix stored as a vector
C
      INTEGER LEVS(*), KODE, NFAC
      REAL Y(*), X(*), EXTR(*)
C
C     Local variables
C
      INTEGER NLFT, NRGT, ILEV, NLEV, JUMP, ILFT, IRGT, JUMPHO, IEL,
     *        JEL, JUMPER
C
      NLFT = 1
      NRGT = 1
      DO 3 ILEV = 1, NFAC
	IF (ILEV - KODE .LT. 0) THEN
	  NLFT = NLFT * LEVS(ILEV)
	ELSE IF (ILEV - KODE .GT. 0) THEN
	  NRGT = NRGT * LEVS(ILEV)
	END IF
    3 CONTINUE
      NLEV = LEVS(KODE)
      JUMP = 0
C
C     Loop over left unit matrix
C
      DO 7 ILFT = 1, NLFT
	JUMP = JUMP + 1
C
C     Loop over right unit matrix
C
	DO 6 IRGT = 1, NRGT
	  JUMPHO = JUMP
          JUMP = JUMP - NRGT
C
C     Extract vector elements for linear combination
C
	  DO 4 ILEV = 1, NLEV
	    JUMP = JUMP + NRGT
	    EXTR(ILEV) = Y(JUMP)
    4     CONTINUE
	  JUMP = JUMPHO - NRGT
	  ILEV = 0
C
C     Begin matrix by vector product.
C     Loop over the various contrasts.
C
	  DO 5 IEL = 1, NLEV
	    JUMP = JUMP + NRGT
	    Y(JUMP) = 0.0
C
C     Form linear combination.
C
	    DO 5 JEL = 1, NLEV
	      ILEV = ILEV + 1
	      Y(JUMP) = Y(JUMP) + X(ILEV) * EXTR(JEL)
    5     CONTINUE
          JUMPER = JUMP
	  JUMP = JUMPHO + 1
    6   CONTINUE
	JUMP = JUMPER
    7 CONTINUE
C
      RETURN
      END

