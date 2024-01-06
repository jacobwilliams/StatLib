!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Algorithm AS 319 and test program
! Converted to Fortran 90 free-format style by Alan Miller
! e-mail: Alan.Miller @ vic.cmis.csiro.au
! URL: www.ozemail.com.au/~milleraj

MODULE as319
IMPLICIT NONE

! COMMON /funerr/ler
! COMMON /test/ig,ifn
LOGICAL, SAVE      :: ler
INTEGER, SAVE      :: ig, ifn

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 50)

END MODULE as319



!-----------------------------------------------------------------------

SUBROUTINE varme(fun, npar, b, f0, nsig, maxfn, iout, ier)
!-----------------------------------------------------------------------

!     CALLING SUBROUTINE FOR SUBROUTINE VARMET

!     ALLOWS FOR SETUP OF DEFAULT PARAMETERS
!     AND EFFICIENT USE OF STORAGE
!     AS WELL AS WRITING OF ERROR MESSAGES

!    VERSION 0.1
!    CODED BY JOHN J. KOVAL
!    MARCH 1986

!    VERSION 0.2
!    CODED BY JOHN J. KOVAL
!    JULY 1988

!    VERSION 0.26
!    CODED BY MURRAY ALEXANDER FOR JOHN J. KOVAL
!    JULY 1989

!    VERSION 0.27
!    MODIFIED BY NAZIH HASSAN, JULY 1993

!    VERSION 0.28
!    MODIFIED BY JOHN KOVAL, JUNE 1996
!    BECAUSE OF COMMENTS FROM REVIEWER FOR APPLIED STATISTICS
!    CHANGES TO ORDER OF PARAMETERS IN GRAD

!    PARAMETERS              MEANING                       DEFAULT
!    ----------              -------                       -------

!     FUN         NAME OF FUNCTION TO BE MINIMIZED

!     NPAR        ORDER OF PARAMETER VECTOR
!                 (NUMBER OF UNKNOWNS)

!     B           ARRAY CONTAINING INITIAL ESTIMATES
!                 ON OUTPUT CONTAINING FINAL ESTIMATES

!     F0          VALUE OF FUNCTION AT MINIMUM

!     NSIG        MACHINE ACCURACY AS NEGATIVE POWER       10 OR 5
!                 OF TEN

!     MAXFN       MAXIMUM NUMBER OF FUNCTION EVALUATIONS     1000
!                 (DOES INCLUDE EVALUATIONS BY SUBROUTINE
!                 GRAD WHICH CALCULATES APPROXIMATE GRADIENT)

!     IOUT        OUTPUT CHANNEL FOR ERROR MESSAGES            0
!                 (IF 0, THEN MESSAGES NOT WRITTEN)

!      IER        ERROR INDICATOR                              0
!                 INTEGER

!-----------------------------------------------------------------------

USE as319
IMPLICIT NONE
INTEGER, INTENT(IN)               :: npar
REAL (dp), INTENT(IN OUT)         :: b(:)
REAL (dp), INTENT(IN OUT)         :: f0
INTEGER, INTENT(IN OUT)           :: nsig
INTEGER, INTENT(IN OUT)           :: maxfn
INTEGER, INTENT(IN OUT)           :: iout
INTEGER, INTENT(OUT)              :: ier

INTERFACE
  SUBROUTINE fun(nord, bp, q)
    USE as319
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: nord
    REAL (dp), INTENT(IN)         :: bp(:)
    REAL (dp), INTENT(OUT)        :: q
  END SUBROUTINE fun
END INTERFACE

REAL (dp) :: gradtl

INTEGER, PARAMETER :: maxf = 1000, msig = 10

!     INITIALIZE

ier=0

IF(nsig == 0) nsig = msig
gradtl = 1.0 / (10.0**(nsig))

IF(gradtl < 0.0) THEN
  IF(iout > 0) WRITE(iout, 300) nsig, gradtl
  300    FORMAT(' NSIG VALUE OF ', i3, ' CREATES NEGATIVE VALUE OF',  &
      ' GRADTL, NAMELY, ', g12.5)
  gradtl = 1.0/(10**(msig))
  IF(iout > 0) WRITE(iout, 310) msig, gradtl
  310    FORMAT(' PROGRAM SUBSTITUTES NSIG VALUE OF ', i3, ' WHICH',  &
      ' GIVES GRADTL VALUE OF ', g12.5)
END IF

IF(maxfn == 0) maxfn = maxf

!      NOW WE ARE READY TO CALL THE MINIMIZATION SUBROUTINE

CALL varmet(fun, npar, b, f0, gradtl, maxfn, ier)

IF(ier > 0.AND.iout > 0)THEN
  WRITE(iout, 30) ier
  30   FORMAT(/' SUBROUTINE VARMET ERROR NUMBER ', i3)
  IF(ier == 1) THEN
    WRITE(iout, 40)
  ELSE IF(ier == 2) THEN
    WRITE(iout, 60)
  ELSE IF(ier == 3) THEN
    WRITE(iout, 70)
  ELSE IF(ier == 4) THEN
    WRITE(iout, 80)
  END IF
  40   FORMAT(' FUNCTION UNDEFINED AT INITIAL VALUE         ')
  60   FORMAT(' GRADIENT UNDEFINED IN TOO MANY DIMENSIONS   ')
  70   FORMAT(' FUNCTON NOT MINIMIZED BUT'  &
      /' UNABLE TO FIND MINIMUM IN DIRECTION OF SEARCH')
  80   FORMAT(' TOO MANY FUNCTION EVALUATIONS REQUIRED      ')
  
END IF

RETURN

CONTAINS


!----------------------------------------------------------------------

SUBROUTINE varmet(fun, npar, b, f0, gradtl, maxfn, ifault)

!       ALGORITHM AS 319 APPL.STATIST. (1997), VOL.46, NO.4

!               VARIABLE METRIC FUNCTION MINIMISATION

INTEGER, INTENT(IN)               :: npar
REAL (dp), INTENT(IN OUT)         :: b(:)
REAL (dp), INTENT(OUT)            :: f0
REAL (dp), INTENT(OUT)            :: gradtl
INTEGER, INTENT(OUT)              :: maxfn
INTEGER, INTENT(OUT)              :: ifault

INTERFACE
  SUBROUTINE fun(nord, bp, q)
    USE as319
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: nord
    REAL (dp), INTENT(IN)         :: bp(:)
    REAL (dp), INTENT(OUT)        :: q
  END SUBROUTINE fun
END INTERFACE

REAL (dp)            :: d1, s, ck, f1, d2
INTEGER              :: i, ic, icount, ifn, ig, ilast, j, k, np
INTEGER, PARAMETER   :: icmax=20
REAL (dp), PARAMETER :: toler=0.00001, w=0.2
REAL (dp)            :: g(npar), h(npar,npar), c(npar), d(npar), t(2*npar)

ig = 0
ifn = 0
ler = .false.
ifault = 0
np = npar + 1

IF (maxfn == 0) maxfn = 1000
IF (gradtl == 0.0) gradtl = 1.0D-10

CALL fun(npar, b, f0)
IF(ler) THEN
  ifault = 1
  RETURN
END IF
ifn = ifn + 1

CALL grad(fun, npar, b, f0, g, t(np:), gradtl, ifault)
IF(ifault > 0) RETURN

ig = ig + 1
ifn = ifn + npar
IF(ifn > maxfn) THEN
  ifault = 4
  RETURN
END IF

10 DO k = 1, npar
  h(k,1:npar) = 0.0
  h(k,k) = 1.00
END DO
ilast = ig

40 DO i = 1, npar
  d(i) = b(i)
  c(i) = g(i)
END DO

d1 = 0.0
DO i = 1, npar
  s = - DOT_PRODUCT( h(i,1:npar), g(1:npar) )
  t(i) = s
  d1 = d1 - s*g(i)
END DO

IF(d1 <= 0.0) THEN
  IF(ilast == ig) THEN
    RETURN
  END IF
  GO TO 10
ELSE
  ck = 1.0
  ic = 0
  90    icount = 0
  DO i = 1,npar
    b(i) = d(i) + ck*t(i)
    IF(b(i) == d(i)) THEN
      icount = icount + 1
    END IF
  END DO
  
  IF(icount >= npar) THEN
    IF(ilast == ig) THEN
      RETURN
    END IF
    GO TO 10
  ELSE
    CALL fun(npar, b, f1)
    
    ifn = ifn + 1
    IF(ifn > maxfn) THEN
      ifault = 4
      RETURN
    ELSE IF(ler) THEN
      ck = w * ck
      ic = ic+1
      IF(ic > icmax) THEN
        ifault = 3
        RETURN
      END IF
      GO TO 90
      
    ELSE IF(f1 >= f0 - d1*ck*toler) THEN
      ck = w * ck
      GO TO 90
    ELSE
      f0 = f1
      CALL grad(fun, npar, b, f0, g, t(np:), gradtl, ifault)
      IF(ifault > 0) THEN
        RETURN
      END IF
      ig = ig + 1
      ifn = ifn + npar
      IF(ifn > maxfn) THEN
        ifault = 4
        RETURN
      END IF
      
      d1 = 0.0
      DO i = 1, npar
        t(i) = ck*t(i)
        c(i) = g(i) - c(i)
        d1 = d1 + t(i)*c(i)
      END DO
      
      IF(d1 <= 0.0) THEN
        GO TO 10
      END IF
      
      d2 = 0.0
      DO i = 1, npar
        s = 0.0
        DO j = 1, npar
          s = s + h(i,j)*c(j)
        END DO
        d(i) = s
        d2 = d2 + s*c(i)
      END DO
      d2 = 1.0 + d2/d1
      
      DO i = 1, npar
        DO j = 1, npar
          h(i,j) = h(i,j) - (t(i)*d(j) + d(i)*t(j) - d2*t(i)*t(j))/d1
        END DO
      END DO
    END IF
  END IF
END IF
GO TO 40
END SUBROUTINE varmet


SUBROUTINE grad(f, npar, b, f0, g, sa, er, ifault)

!     CALCULATE APPROXIMATE GRADIENT

INTEGER, INTENT(IN)           :: npar
REAL (dp), INTENT(IN OUT)     :: b(:)
REAL (dp), INTENT(IN)         :: f0
REAL (dp), INTENT(OUT)        :: g(:)
REAL (dp), INTENT(OUT)        :: sa(:)
REAL (dp), INTENT(IN)         :: er
INTEGER, INTENT(OUT)          :: ifault

INTERFACE
  SUBROUTINE f(nord, bp, q)
    USE as319
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: nord
    REAL (dp), INTENT(IN)         :: bp(:)
    REAL (dp), INTENT(OUT)        :: q
  END SUBROUTINE f
END INTERFACE

REAL (dp) :: h, f1
INTEGER   :: i, jc, jcmax

jcmax=npar - 2
jc = 0

DO i = 1, npar
  h =(ABS(b(i)) + SQRT(er)) * SQRT(er)
  sa(i) = b(i)
  b(i) = b(i) + h
  CALL f(npar, b, f1)
  b(i) = sa(i)
  
  IF(ler) THEN
    f1 = f0 + h
    jc = jc + 1
  END IF
  
  g(i) = (f1 - f0)/h
END DO

IF(jc > jcmax) ifault = 2
RETURN
END SUBROUTINE grad

END SUBROUTINE varme


!----------------------------------------------------------------------
PROGRAM  var
!----------------------------------------------------------------------
!       A PROGRAM TO IMPLEMENT A QUASI-NEWTON METHOD.
!       USING NEW ALGORITHM VARMET   JULY 1994
!----------------------------------------------------------------------

USE as319
IMPLICIT NONE
INTEGER, PARAMETER :: n=2
INTEGER            :: ier
INTEGER, SAVE      :: gradtl = 12, maxfn = 1000, mess = 6
REAL (dp)          :: x(n), xtmp(n), fp

INTERFACE
  SUBROUTINE fun(nord, bp, q)
    USE as319
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: nord
    REAL (dp), INTENT(IN)         :: bp(:)
    REAL (dp), INTENT(OUT)        :: q
  END SUBROUTINE fun

  SUBROUTINE varme(fun, npar, b, f0, nsig, maxfn, iout, ier)
    USE as319
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: npar
    REAL (dp), INTENT(IN OUT)         :: b(:)
    REAL (dp), INTENT(IN OUT)         :: f0
    INTEGER, INTENT(IN OUT)           :: nsig
    INTEGER, INTENT(IN OUT)           :: maxfn
    INTEGER, INTENT(IN OUT)           :: iout
    INTEGER, INTENT(OUT)              :: ier
    INTERFACE
      SUBROUTINE fun(nord, bp, q)
        USE as319
        IMPLICIT NONE
        INTEGER, INTENT(IN)           :: nord
        REAL (dp), INTENT(IN)         :: bp(:)
        REAL (dp), INTENT(OUT)        :: q
      END SUBROUTINE fun
    END INTERFACE
  END SUBROUTINE varme
END INTERFACE

WRITE(*,*)'  '
WRITE(*,*) 'INPUT YOUR STARTING GUESS: '
READ(*,*) xtmp(1:n)
WRITE(*,*)' '
WRITE(*,*)'INITIALIZATION COMPLETE.'
WRITE(*,*)'***************************************'

x(1:n)=xtmp(1:n)
ifn = 0
CALL varme(fun, n, x, fp, gradtl, maxfn, mess, ier)

WRITE(*,*)' '
IF (ier /= 0) WRITE(*, *) '** IER =', ier, ' **'

WRITE(*,*)' THE NUMBER OF FUNCTION EVALUATIONS IS ', ifn
WRITE(*,*)' '
WRITE(*,*)'THE MINIMUM FOUND IS ', x(1:n)
WRITE(*,*)' '
CALL fun(n, x, fp)
WRITE(*, *)'THE FUNCTION VALUE IS: ', fp
STOP
END PROGRAM  var
!----------------------------------------------------------------------

SUBROUTINE fun(nord, bp, q)
!----------------------------------------------------------------------

USE as319
IMPLICIT NONE
INTEGER, INTENT(IN)           :: nord
REAL (dp), INTENT(IN)         :: bp(:)
REAL (dp), INTENT(OUT)        :: q

q=100.*(bp(2) - bp(1)**2)**2 + (bp(1) - 1.)**2
ifn = ifn + 1
ler = .false.
IF (nord < 1) WRITE(*, *)'** NORD must be > 0, actual value =', nord
RETURN
END SUBROUTINE fun

