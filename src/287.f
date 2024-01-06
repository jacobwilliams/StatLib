C        .. ULCC Toolpack/1 2.1
      SUBROUTINE initial(ns, m, emax, x, hx, hpx, lb, xlb, ub, xub,
     *                   ifault, iwv, rwv)
C
C        ALGORITHM AS 287.1 APPL.STATIST. (1993), VOL.42, NO.4
C
C        This subroutine takes as input the number of starting values M
C        and the starting points x(i), hx(i), hpx(i)  i=1, m. As output
C        we have pointer ipt along with ilow and ihigh and the lower
C        and upper hulls defined  by z, hz, scum, cu, hulb, huub stored
C        in working vectors iwv and rwv
C
      INTEGER ns, m, ifault, iwv( * )
      LOGICAL lb, ub
      REAL emax, x( * ), hx( * ), hpx( * ), xlb, xub, rwv( * )
C
      INTEGER i, ihigh, ihpx, ihuz, ihx, iipt, ilow, iscum, ix, iz, nn
      LOGICAL HORIZ
      REAL a, b, alcu, cu, eps, expon, hulb, huub, huzmax, zlog, zmax
C
      zlog(a) = alog(a)
      zmax(a, b) = amax1(a, b)
      eps = expon(-emax, emax)
C
      ifault = 0
      ilow = 1
      ihigh = 1
      nn = ns + 1
C
C        at least one starting point
C
      IF (m .LT. 1) ifault = 1
C
      huzmax = hx(1)
C
      IF ( .NOT. ub) xub = 0.0
      IF ( .NOT. lb) xlb = 0.0
C
      hulb = (xlb - x(1)) * hpx(1) + hx(1)
      huub = (xub - x(1)) * hpx(1) + hx(1)
C
C        if bounded on both sides
C
      IF ((ub) .AND. (lb)) THEN
         huzmax = zmax(huub, hulb)
         horiz = (abs(hpx(1)) .LT. eps)
         IF (horiz) THEN
            cu = expon((huub + hulb) * 0.5 - huzmax, emax) * (xub - xlb)
         ELSE
            cu = expon(huub - huzmax, emax) *
     *           (1 - expon(hulb - huub, emax)) / hpx(1)
         END IF
      ELSE IF ((ub) .AND. ( .NOT. lb)) THEN
C
C        if bounded on the right and unbounded on the left
C
         huzmax = huub
         cu = 1.0 / hpx(1)
      ELSE IF (( .NOT. ub) .AND. (lb)) THEN
C
C        if bounded on the left and unbounded on the right
C
         huzmax = hulb
         cu = -1.0 / hpx(1)
C
C        if unbounded at least 2 starting points
C
      ELSE
         cu = 0.0
         IF (m .LT. 2) ifault = 1
      END IF
C
      IF (cu .GT. 0.0) alcu = zlog(cu)
C
C        set pointers
C
      iipt = 6
      iz = 9
      ihuz = nn + iz
      iscum = nn + ihuz
      ix = nn + iscum
      ihx = nn + ix
      ihpx = nn + ihx
C
C        store values in working vectors
C
      iwv(1) = ilow
      iwv(2) = ihigh
      iwv(3) = ns
      iwv(4) = 1
C
      IF (lb) THEN
         iwv(5) = 1
      ELSE
         iwv(5) = 0
      END IF
C
      IF (ub) THEN
         iwv(6) = 1
      ELSE
         iwv(6) = 0
      END IF
C
      IF (ns .LT. m) ifault = 2
C
      iwv(iipt + 1) = 0
      rwv(1) = hulb
      rwv(2) = huub
      rwv(3) = emax
      rwv(4) = eps
      rwv(5) = cu
      rwv(6) = alcu
      rwv(7) = huzmax
      rwv(8) = xlb
      rwv(9) = xub
      rwv(iscum + 1) = 1.0
      DO 10 i = 1, m
         rwv(ix + i) = x(i)
         rwv(ihx + i) = hx(i)
         rwv(ihpx + i) = hpx(i)
   10 CONTINUE
C
C        create lower and upper hulls
C
      i = 1
   20 IF (i .LT. m) THEN
         CALL update(iwv(4), iwv(1), iwv(2), iwv(iipt + 1),
     *               rwv(iscum + 1), rwv(5), rwv(ix + 1), rwv(ihx + 1),
     *               rwv(ihpx + 1), rwv(iz + 1), rwv(ihuz + 1), rwv(7),
     *               rwv(3), lb, rwv(8), rwv(1), ub, rwv(9), rwv(2),
     *               ifault, rwv(4), rwv(6))
         i = iwv(4)
         IF (ifault .NE. 0) RETURN
         GO TO 20
      END IF
C
C        test for wrong starting points
C
      IF (( .NOT. lb) .AND. (hpx(iwv(1)) .LT. eps)) ifault = 3
      IF (( .NOT. ub) .AND. (hpx(iwv(2)) .GT. -eps)) ifault = 4
C
      RETURN
      END
      SUBROUTINE sample(iwv, rwv, eval, beta, ifault)
C
C        ALGORITHM AS 287.2 APPL.STATIST. (1993), VOL.42, NO.4
C
      INTEGER iwv( * ), ifault
      REAL rwv( * ), beta
C
      INTEGER ihpx, ihuz, ihx, iipt, iscum, ix, iz, nn, ns
      LOGICAL lb, ub
      EXTERNAL EVAL
C
C        set pointers
C
      iipt = 6
      iz = 9
      ns = iwv(3)
      nn = ns + 1
      ihuz = nn + iz
      iscum = nn + ihuz
      ix = nn + iscum
      ihx = nn + ix
      ihpx = nn + ihx
      lb = .FALSE.
      ub = .FALSE.
C
      IF (iwv(5) .EQ. 1) lb = .TRUE.
      IF (iwv(6) .EQ. 1) ub = .TRUE.
C
C        call sampling subroutine
C
      CALL spl1(ns, iwv(4), iwv(1), iwv(2), iwv(iipt + 1),
     *          rwv(iscum + 1), rwv(5), rwv(ix + 1), rwv(ihx + 1),
     *          rwv(ihpx + 1), rwv(iz + 1), rwv(ihuz + 1), rwv(7), lb,
     *          rwv(8), rwv(1), ub, rwv(9), rwv(2), eval, beta, ifault,
     *          rwv(3), rwv(4), rwv(6))
C
      RETURN
      END
      SUBROUTINE spl1(ns, n, ilow, ihigh, ipt, scum, cu, x, hx, hpx, z,
     *                huz, huzmax, lb, xlb, hulb, ub, xub, huub, eval,
     *                beta, ifault, emax, eps, alcu)
C
C        ALGORITHM AS 287.3 APPL.STATIST. (1993), VOL.42, NO.4
C
C        this subroutine performs the adaptive rejection sampling, it
C        calls subroutine splhull to sample from the upper hull, if the
C        sampling involves a function evaluation it calls the updating
C        subroutine
C
      INTEGER ns, n, ilow, ihigh, ipt( * ), ifault
      LOGICAL lb, ub
      REAL scum( * ), cu, x( * ), hx( * ), hpx( * ), z( * ), huz( * ),
     *     huzmax, xlb, hulb, xub, huub, beta, emax, eps, alcu
C
      INTEGER i, j, l, n1
      LOGICAL sampld
      REAL a, alhl, alhu, alu1, fx, random, u1, u2, zlog
      EXTERNAL eval
C
      zlog(a) = alog(a)
C
      sampld = .FALSE.
   10 IF ( .NOT. sampld) THEN
         u2 = random(L)
C
C        test for zero random number
C
         IF (u2 .EQ. 0.0) THEN
            ifault = 6
            RETURN
         END IF
         CALL splhull(u2, ipt, ilow, lb, xlb, hulb, huzmax, alcu, x, hx,
     *                hpx, z, huz, scum, eps, emax, beta, i, j)
C
C        sample u1 to compute rejection
C
         u1 = random(l)
         IF (u1 .EQ. 0.0) ifault = 6
         alu1 = zlog(u1)
C
C        compute alhu: upper hull at point u1
C
         alhu = hpx(i) * (beta - x(i)) + hx(i) - huzmax
         IF ((beta .GT. x(ilow)) .AND. (beta .LT. x(ihigh))) THEN
C
C        compute alhl: value of the lower hull at point u1
C
            IF (beta .GT. x(i)) THEN
               j = i
               i = ipt(i)
            END IF
            alhl = hx(i) + (beta - x(i)) * (hx(i) - hx(j)) /
     *             (x(i) - x(j)) - huzmax
C
C        squeezing test
C
            IF ((alhl - alhu) .GT. alu1) THEN
               sampld = .TRUE.
            END IF
         END IF
C
C        if not sampled evaluate the function, do the rejection test
C        and update
C
         IF ( .NOT. sampld) THEN
            n1 = n + 1
            x(n1) = beta
            CALL eval(x(n1), hx(n1), hpx(n1))
            fx = hx(n1) - huzmax
            IF (alu1 .LT. (fx - alhu)) sampld = .TRUE.
C
C        update while the number of points defining the hulls is lower
C        than ns
C
            IF (n .LT. ns) THEN
               CALL update(n, ilow, ihigh, ipt, scum, cu, x, hx, hpx, z,
     *                     huz, huzmax, emax, lb, xlb, hulb, ub, xub,
     *                     huub, ifault, eps, alcu)
            END IF
            IF (ifault .NE. 0) RETURN
         END IF
         GO TO 10
      END IF
C
      RETURN
      END
      SUBROUTINE splhull(u2, ipt, ilow, lb, xlb, hulb, huzmax, alcu, x,
     *                   hx, hpx, z, huz, scum, eps, emax, beta, i, j)
C
C        ALGORITHM AS 287.4 APPL.STATIST. (1993), VOL.42, NO.4
C
C        this subroutine samples beta from the normalised upper hull
C
      REAL u2, xlb, hulb, huzmax, alcu, x( * ), hx( * ), hpx( * ),
     *     z( * ), huz( * ), scum( * ), eps, emax, beta
      INTEGER ipt( * ), ilow, i, j
      LOGICAL lb
C
      REAL a, eh, expon, logdu, logtg, sign, zlog
      LOGICAL horiz
C
      zlog(a) = alog(a)
      i = ilow
C
C        find from which exponential piece you sample
C
   10 IF (u2 .GT. scum(i)) THEN
         j = i
         i = ipt(i)
         GO TO 10
      END IF
C
      IF (i .EQ. ilow) THEN
C
C        sample below z(ilow), depending on the existence of a lower
C        bound
C
         IF (lb) THEN
            eh = hulb - huzmax - alcu
            horiz = (abs(hpx(ilow)) .LT. eps)
            IF (horiz) THEN
               beta = xlb + u2 * expon(-eh, emax)
            ELSE
               sign = abs(hpx(i)) / hpx(i)
               logtg = zlog(abs(hpx(i)))
               logdu = zlog(u2)
               eh = logdu + logtg - eh
               IF (eh .LT. emax) THEN
                  beta = xlb + zlog(1.0 + sign * expon(eh, emax)) /
     *                   hpx(i)
               ELSE
                  beta = xlb + eh / hpx(i)
               END IF
            END IF
         ELSE
C
C        hpx(i) must be positive, x(ilow) is left of the mode
C
            beta = (zlog(hpx(i) * u2) + alcu - hx(i) + x(i) * hpx(i) +
     *             huzmax) / hpx(i)
         END IF
      ELSE
C
C        sample above (j)
C
         eh = huz(j) - huzmax - alcu
         horiz = (abs(hpx(i)) .LT. eps)
         IF (horiz) THEN
            beta = z(j) + (u2 - scum(j)) * expon(-eh, emax)
         ELSE
            sign = abs(hpx(i)) / hpx(i)
            logtg = zlog(abs(hpx(i)))
            logdu = zlog(u2 - scum(j))
            eh = logdu + logtg - eh
            IF (eh .LT. emax) THEN
               beta = z(j) + (zlog(1.0 + sign * expon(eh, emax))) /
     *                hpx(i)
            ELSE
               beta = z(j) + eh / hpx(i)
            END IF
         END IF
      END IF
C
      RETURN
      END
      SUBROUTINE intersection(x1, y1, yp1, x2, y2, yp2, z1, hz1, eps,
     *                        ifault)
C
C        ALGORITHM AS 287.5 APPL.STATIST. (1993), VOL.42, NO.4
C
C        computes the intersection (z1,hz1) between 2 tangents defined
C        by x1, y1, yp1 and x2, y2, yp2
C
      REAL x1, y1, yp1, x2, y2, yp2, z1, hz1, eps
      INTEGER ifault
C
      REAL y12, dh
C
C        first test for non-concavity
C
      y12 = y1 + yp1 * (x2 - x1)
      y21 = y2 + yp2 * (x1 - x2)
      IF ((y21 .LT. y1) .OR. (y12 .LT. y2)) THEN
         ifault = 5
         RETURN
      END IF
C
      dh = yp2 - yp1
C
C        IF the lines are nearly parallel,
C        the intersection is taken at the mid-point
C
      IF (abs(dh) .LE. eps) THEN
         z1 = 0.5 * (x1 + x2)
         hz1 = 0.5 * (y1 + y2)
C
C        Else compute from the left or the right for greater
C        numerical precision
C
      ELSE IF (abs(yp1) .LT. abs(yp2)) THEN
         z1 = x2 + (y1 - y2 + yp1 * (x2 - x1)) / dh
         hz1 = yp1 * (z1 - x1) + y1
      ELSE
         z1 = x1 + (y1 - y2 + yp2 * (x2 - x1)) / dh
         hz1 = yp2 * (z1 - x2) + y2
      END IF
C
C        test for misbehaviour due to numerical imprecision
C
      IF ((z1 .LT. x1) .OR. (z1 .GT. x2)) ifault = 7
C
      RETURN
      END
      SUBROUTINE update(n, ilow, ihigh, ipt, scum, cu, x, hx, hpx, z,
     *                  huz, HUZMAX, EMAX, lb, xlb, hulb, ub, xub, huub,
     *                  ifault, eps, alcu)
C
C        ALGORITHM AS 287.6 APPL.STATIST. (1993), VOL.42, NO.4
C
C        this subroutine increments n and updates all the parameters
C        which define the lower and the upper hull
C
      INTEGER n, ilow, ihigh, ipt( * ), ifault
      LOGICAL lb, ub
      REAL scum( * ), cu, x( * ), hx( * ), hpx( * ), z( * ), huz( * ),
     *     huzmax, emax, xlb, hulb, xub, huub, eps, alcu
C
      INTEGER i, j
      LOGICAL horiz
      REAL a, b, dh, expon, u, zlog, zmax
      EXTERNAL eval
C
      zlog(a) = alog(a)
      zmax(a, b) = amax1(a, b)
C
      n = n + 1
C
C        update z, huz and ipt
C
      IF (x(n) .LT. x(ilow)) THEN
C
C        insert x(n) below x(ilow), test for non-concavity
C
         IF (hpx(ilow) .GT. hpx(n)) ifault = 5
C
         ipt(n) = ilow
         CALL intersection(x(n), hx(n), hpx(n), x(ilow), hx(ilow),
     *                     hpx(ilow), z(n), huz(n), eps, ifault)
C
         IF (ifault .NE. 0) RETURN
         IF (lb) hulb = hpx(n) * (xlb - x(n)) + hx(n)
C
         ilow = n
      ELSE
         i = ilow
         j = i
C
C        find where to insert x(n)
C
   10    IF ((x(n) .GE. x(i)) .AND. (ipt(i) .NE. 0)) THEN
            j = i
            i = ipt(i)
            GO TO 10
         END IF
C
         IF (x(n) .GE. x(i)) THEN
C
C        insert above x(ihigh), test for non-concavity
C
            IF (hpx(i) .LT. hpx(n)) ifault = 5
C
            ihigh = n
            ipt(i) = n
            ipt(n) = 0
            CALL intersection(x(i), hx(i), hpx(i), x(n), hx(n), hpx(n),
     *                        z(i), huz(i), eps, ifault)
C
            IF (ifault .NE. 0) RETURN
C
            huub = hpx(n) * (xub - x(n)) + hx(n)
            z(n) = 0.0
            huz(n) = 0.0
         ELSE
C
C        insert x(n) between x(j) and x(i), test for non-concavity
C
            IF ((hpx(j) .LT. hpx(n)) .OR.
     *          (hpx(i) .GT. hpx(n))) ifault = 5
C
            ipt(j) = n
            ipt(n) = i
C
C        insert z(j) between x(j) and x(n)
C
            CALL intersection(x(j), hx(j), hpx(j), x(n), hx(n), hpx(n),
     *                        z(j), huz(j), eps, ifault)
C
            IF (ifault .NE. 0) RETURN
C
C        insert z(n) between x(n) and x(i)
C
            CALL intersection(x(n), hx(n), hpx(n), x(i), hx(i), hpx(i),
     *                        z(n), huz(n), eps, ifault)
C
            IF (ifault .NE. 0) RETURN
C
         END IF
      END IF
C
C        update huzmax
C
      j = ilow
      i = ipt(j)
      huzmax = huz(j)
C
   20 IF ((huz(j) .LT. huz(i)) .AND. (ipt(i) .NE. 0)) THEN
         j = i
         i = ipt(i)
         huzmax = zmax(huzmax, huz(j))
         GO TO 20
      END IF
C
      IF (lb) huzmax = zmax(huzmax, hulb)
      IF (ub) huzmax = zmax(huzmax, huub)
C
C        update cu, scum receives area below exponentiated upper hull
C        left of z(i)
C
      i = ilow
      horiz = (abs(hpx(ilow)) .LT. eps)
      IF (( .NOT. lb) .AND. ( .NOT. horiz)) THEN
         cu = expon(huz(i) - huzmax, emax) / hpx(i)
      ELSE IF (lb .AND. horiz) THEN
         cu = (z(ilow) - xlb) * expon(hulb - huzmax, emax)
      ELSE IF (lb .AND. ( .NOT. horiz)) THEN
         dh = hulb - huz(i)
         IF (dh .GT. emax) THEN
            cu = -expon(hulb - huzmax, emax) / hpx(i)
         ELSE
            cu = expon(huz(i) - huzmax, emax) * (1 - expon(dh, emax)) /
     *           hpx(i)
         END IF
      ELSE
         cu = 0
      END IF
      scum(i) = cu
      j = i
      i = ipt(i)
C
   30 IF (ipt(i) .NE. 0) THEN
         dh = huz(j) - huz(i)
         horiz = (abs(hpx(i)) .LT. eps)
         IF (horiz) THEN
            cu = cu + (z(i) - z(j)) * expon((huz(i) + huz(j)) * 0.5 -
     *           huzmax, emax)
         ELSE
            IF (dh .LT. emax) THEN
               cu = cu + expon(huz(i) - huzmax, emax) *
     *              (1 - expon(dh, emax)) / hpx(i)
            ELSE
               cu = cu - expon(huz(j) - huzmax, emax) / hpx(i)
            END IF
         END IF
         j = i
         i = ipt(i)
         scum(j) = cu
         GO TO 30
      END IF
C
      horiz = (abs(hpx(i)) .LT. eps)
C
C        if the derivative is very small the tangent is nearly
C        horizontal
C
      IF ( .NOT. (ub .OR. horiz)) THEN
         cu = cu - expon(huz(j) - huzmax, emax) / hpx(i)
      ELSE IF (ub .AND. horiz) THEN
         cu = cu + (xub - x(i)) * expon((huub + hx(i)) * 0.5 - huzmax,
     *        emax)
      ELSE IF (ub .AND. ( .NOT. horiz)) THEN
         dh = huz(j) - huub
         IF (dh .GT. emax) THEN
            cu = cu - expon(huz(j) - huzmax, emax) / hpx(i)
         ELSE
            cu = cu + expon(huub - huzmax, emax) *
     *           (1 - expon(dh, emax)) / hpx(i)
         END IF
      END IF
C
      scum(i) = cu
      IF (cu .GT. 0) alcu = zlog(cu)
C
C        normalize scum to obtain a cumulative probability while
C        excluding unnecessary points
C
      i = ilow
      u = (cu - scum(i)) / cu
C
      IF ((u .EQ. 1.0) .AND. (hpx(ipt(i)) .GT. 0)) THEN
         ilow = ipt(i)
         scum(i) = 0.0
      ELSE
         scum(i) = 1.0 - u
      END IF
C
      j = i
      i = ipt(i)
C
   40 IF (ipt(i) .NE. 0) THEN
         j = i
         i = ipt(i)
         u = (cu - scum(j)) / cu
         IF ((u .EQ. 1.0) .AND. (hpx(i) .GT. 0)) THEN
            ilow = i
         ELSE
            scum(j) = 1.0 - u
         END IF
         GO TO 40
      END IF
C
      scum(i) = 1.0
      IF (ub) huub = hpx(ihigh) * (xub - x(ihigh)) + hx(ihigh)
      IF (lb) hulb = hpx(ilow) * (xlb - x(ilow)) + hx(ilow)
C
      RETURN
      END
      REAL FUNCTION expon(x, emax)
C
C        ALGORITHM AS 287.7 APPL.STATIST. (1993), VOL.42, NO.4
C
C        performs an exponential without underflow
C
      REAL x, emax
C
      REAL ZEXP
C
      zexp(x) = exp(x)
C
      IF (x .LT. -emax) THEN
         expon = 0.0
      ELSE
         expon = zexp(x)
      END IF
C
      RETURN
      END

