c This file contains two versions of the algorithm.   The first is that
c published in the journal, including several auxiliary routines which
c are not in the journal, but converted to double precision.
c The second version is one which calls IMSL routines.
c
      subroutine mulnor(a, b, sig, eps, n, inf, prob, bound, ifault)
c
c        Algorithm AS 195  appl. statist. (1984) vol.33, no.1
c
c        Computes multivariate normal distribution function and computes
c        the probability that a multivariate normal vector falls in a
c        rectangle in n-space with error less than eps.
c
c     Auxiliary routines required: BIVNOR (CACM algorithm 462) which
c         needs DERF, PPND = AS111 (or PPND7 or PPND16 from AS241), and
c         SYMINV = AS7 which calls CHOL = AS6.   BIVNOR and DERF are
c         included in this file.
c
      double precision
     *     a(*), b(*), sig(*), c(7), d(7), co(25), sd(2), coef(5, 3),
     *  binc(7, 5), bl(7, 5), br(7, 35, 5), r(36, 5), s(5, 27),
     *  xss(6, 6), pr2(6), prep(6), preb(6), cv(5, 15), sinv(28),
     *  condl(5), xm(5), cond(5), beta(5, 5), bh(5), fe(5), sigma(28),
     *  ep(5), del(3, 5), bou4(5), bou5(5), ans(6), fact(5), prod(5),
     *  bint(6), bcs(5), bcn(7), eo(13), do(25)
      integer inf(7), intvl(5), ind(5), ksa(5), num(5), itype(5)
      logical simps, ipr
      double precision eps, bound, prob, cons2, cons3, cons4, epsmin,
     +    epssim, zero, p05, p15, half, one, two, three, six, eight,
     +    ten, twelve, fiften, forty5, ninety, nine45, bou1, bou2, bou3,
     +    cup, det, eplos, epsi, ept, fac, fsa, rho, sigc, t, tem, temb,
     +    temp, wb, wl, wt, wu, x1, x2, xs, y1, y2, z
      double precision alnorm, bivnor, f2, f3, f4, f5, f6, ppnd, phi, x,
     +    y
      data coef
     *  /0.311111111111111D0, 1.422222222222222D0, 0.533333333333333D0,
     *  1.422222222222222D0, 0.311111111111111D0, 0.333333333333333D0,
     *  0.0D0, 1.333333333333333D0, 0.0D0, 0.333333333333333D0, 0.5D0,
     *  0.0D0, 0.0D0, 0.0D0, 0.5D0/
      data bcs /1.0D0, 4.0D0, 6.0D0, 4.0D0, 1.D0/
      data bcn /1.0D0, 6.0D0, 15.0D0, 20.0D0, 15.0D0, 6.0D0, 1.D0/
      data do /-3.75043972D0,-3.3242574D0, -2.85697D0, -2.36675941D0,
     *  -2.3344142D0, -1.8891759D0, -1.7320508D0, -1.3556262D0,
     *  -1.15440539D0, -1.D0, -0.74196378D0, -0.61670659D0,0.0D0,
     *  0.61670659D0, 0.74196378D0, 1.0D0, 1.15440539D0,1.3556262D0,
     *  1.7320508D0, 1.8891759D0, 2.3344142D0, 2.36675941D0, 2.85697D0,
     *  3.3242574D0,3.75043972D0/
      data eo / -2.85697D0, -2.3344142D0, -1.7320508D0, -1.3556262D0,
     *  -1.0D0, -0.74196378D0, 0.0D0, 0.74196378D0, 1.0D0, 1.3556262D0,
     *  1.7320508D0, 2.3344142D0, 2.85697D0/
      data cons2 /6.879833D0/, cons3 /4.517004D0/, cons4 /76.371214D0/,
     *  epsmin /1.0d-8/, epssim /6.0d-5/
      data zero, p05, p15, half, one, two, three, six, eight, ten,
     *  twelve, fiften, forty5, ninety, nine45 /0.0D0, 0.05D0, 0.15D0,
     *  0.5D0, 1.0D0, 2.0D0, 3.0D0, 6.0D0, 8.0D0, 10.0D0, 12.0D0,
     *  15.0D0, 45.0D0, 90.D0, 945.D0/
c
c        Functions needed to calculate the first six derivatives of
c        the normal density
c
      f2(x, y) = abs((x * x - one) / (y * y))
      f3(x, y) = abs(-x * (x * x - three) / (y * y * y))
      f4(x, y) = abs((three + x * x * (-six + x * x)) / (y ** 4))
      f5(x, y) = abs(x * (-fiften + x * x * (ten - x * x))) / (y ** 5)
      f6(x, y) = abs(-fiften + x * x * (forty5 - x * x * (-fiften +
     *  x * x))) / (y ** 6)
c        Checking for faulty data
c
      ifault = 0
      if (eps .le. epsmin) ifault = 2
      if (n .le. 0 .or. n .gt. 7) ifault = 3
      if (ifault .ne. 0) return
      do 1 i = 1, n
    1 if (inf(i) .eq. 2 .and. a(i) .lt. b(i)) ifault = 100 + i
      if (ifault .ne. 0) return
      bound = eps
c        Finding z such that p(n(0,1).gt.z) .lt. 0.15*eps/n
c        co will contain the roots of the first 5 or 7 Hermite
c        polynomials
      ept = eps * p15 / float(n)
      z = -ppnd(ept, ifault) + epsmin
      if (ifault .ne. 0) return
      cup = alnorm(z, .true.)
c
c        inverting sig and the n-2 lower right hand principal minors
c
      ik = 0
      ij = 0
      do 3 i = 1, n
        do 3 j = 1, i
          ik = ik + 1
          if (i .eq. j) goto 2
          ij = ij + 1
          sigma(ik) = sig(ij)
          goto 3
    2     sigma(ik) = one
    3 continue
      if (n .le. 2) goto 4
      call invert(sigma, sinv, cv, n, det, ifault)
      if (ifault .ne. 0) return
      simps = .true.
      if (det .lt. p05 .or. eps .le. epssim) simps = .false.
      prob = zero
      det = sinv(1) * sinv(3) - sinv(2) * sinv(2)
      sd(1) = sqrt(sinv(3) / det)
      sd(2) = sqrt(sinv(1) / det)
      rho = -sinv(2) / (sd(1) * sd(2) * det)
      if (abs(rho) .gt. one) goto 400
    4 nm2 = n - 2
      nm1 = n - 1
c
c        checking whether upper and lower endpoints are too big
c
      eplos = zero
      do 6 l = 1, n
        c(l) = max(b(l), -z)
        d(l) = min(a(l), z)
        if (inf(l) .eq. 0) d(l) = z
        if (inf(l) .eq. 1) c(l) = -z
        if (a(l) .gt. z .or. inf(l) .eq. 0) eplos = eplos + cup
        if (b(l) .lt. -z .or. inf(l) .eq. 1) eplos = eplos + cup
        if (c(l) .ge. d(l)) return
    6 continue
      if (n .eq. 1) goto 350
      fac = one
      ipr = .false.
      if (inf(1) .ne. 1 .or. inf(2) .ne. 1) goto 7
      ipr = .true.
      eplos = eplos - two * cup
      goto 8
    7 if (inf(1) .ne. 0 .or. inf(2) .ne. 0) goto 8
      fac = -one
      ipr = .true.
      d(1) = c(1)
      d(2) = c(2)
      eplos = eplos - two * cup
    8 if (n .eq. 2) goto 360
      ifault = 5
      epsi = (eps - eplos) / float(nm2)
c        Finding regression coefficients (beta,bh) and bounds on the
c        conditional integrals (binc)
c        cond(l)=conditional variance of variable n-l+1 given later
c        variables
c
      do 15 l = 1, nm2
        cond(l) = one / sqrt(cv(l, 1))
        condl(l) = log(cond(l))
        do 10 i = 1, l
        beta(l, i) = zero
        do 10 j = 1, l
          jk = (l - i + 1) * (l - i) / 2 + j
          if (j .gt. l - i + 1) jk = j * (j - 1) / 2 + l - i + 1
          jn = (n - l + j) * (n - l + j - 1) / 2 + n - l
          beta(l, i) = beta(l, i) + sigma(jn) * cv(l, jk)
   10   continue
        k = n - l - 1
        bh(k + 1) = beta(l, l)
        do 11 i = 1, k
          bh(i) = zero
          do 11 j = 1, l
            jn = (j + n - l) * (j + n - l - 1) / 2 + n - l
            ijk = 1 + j * (j - 1) / 2
            bh(i) = bh(i) + sigma(jn) * cv(l, ijk)
   11   continue
        k = 0
        sigc = zero
        do 12 j = 1, l
          do 12 i = 1, j
            k = k + 1
            sigc = sigc + bh(i) * bh(j) * sinv(k)
   12   continue
        binc(1, l) = one
        binc(2, l) = sqrt(sigc)
        binc(3, l) = two * sigc
        binc(4, l) = cons2 * sigc * binc(2, l)
        binc(5, l) = twelve * sigc * sigc
        if (simps) goto 13
        binc(6, l) = sigc * sigc * binc(1, l) * cons3
        binc(7, l) = (sigc ** 3) * cons4
   13   if (l .lt. nm2) goto 15
        do 14 i = 1, nm2
          bh(i) = zero
          do 14 j = 1, nm2
            jk = (l - i + 1) * (l - i) / 2 + j
            if (j .gt. l - i + 1) jk = j * (j - 1) / 2 + l - i + 1
            jn = (2 + j) * (1 + j) / 2 + 1
            bh(i) = bh(i) + sigma(jn) * cv(l, jk)
   14   continue
   15 continue
      l = 1
c
c        co will contain the roots of the first 5 or 7 Hermite
c        polynomials
c
      if (simps) goto 50
      do 40 i = 1, 25
   40 co(i) = do(i)
      iend = 25
      ien = 7
      goto 60
   50 do 55 i = 1, 13
   55 co(i) = eo(i)
      iend = 13
      ien = 5
c
c        Initialising values.   xss contains partial sums used for
c        calculating conditional means.
c
   60 do 70 i = 1, nm1
   70 xss(i, 1) = zero
      xm(1) = zero
      prod(1) = one
      pr2(1) = one
      do 80 i = 1, nm2
        ni = n - i + 1
        pr2(i + 1) = pr2(i) * (d(ni) - c(ni))
   80 continue
c
c        bint(l) is a bound on the error accumulated at levels l and
c                deeper.
c        ans(l)  is the accumulated integral at level l.
c        prep(l) contains the integrand at level l.
c        preb(l) bounds the error accumulated at levels deeper than l.
c
      bint(nm1) = zero
   90 intvl(l) = 2
      ans(l) = zero
      bou4(l) = zero
      bint(l) = zero
      prep(l) = zero
      preb(l) = zero
      k = 1
c
c        Finding which of the co are in the current interval.
c        s(l,.) are the endpoints of intervals on which the integrand
c        at level l and its derivatives are monotone.
c        num(l) is the number of such intervals.
c
      nl = n - l + 1
      s(l, 1) = c(nl) - xm(l)
      s(l, iend + 2) = d(nl) - xm(l)
      num(l) = iend + 2
      do 91 i = 1, iend
        njs = i
        if (s(l, 1) .lt. co(i) * cond(l)) goto 92
        num(l) = num(l) - 1
   91 continue
   92 if (num(l) .eq. 2) goto 99
      do 94 i = njs, iend
        mjs = iend - i + njs
        if (s(l, iend + 2) .ge. co(mjs) * cond(l)) goto 96
        num(l) = num(l) - 1
   94 continue
   96 if (num(l) .eq. 2) goto 99
      do 98 i = njs, mjs
        inj = i - njs + 2
        s(l, inj) = co(i) * cond(l)
   98 continue
   99 numl = num(l)
      s(l, numl) = s(l, iend + 2)
c
c        ep(l)  is an upper limit on the allowable error at level l.
c        r(k,l) is the right end-point of the current sub-interval.
c        fe(l)  is the left end-point of the current sub-interval.
c        ind(l)-1  indicates which point of the Newton-Cotes formula
c        we are dealing with.
c
      ep(l) = epsi / pr2(l + 1)
      r(1, l) = s(l, 2)
      ind(l) = 6
c
c        Bounding derivatives at left end-point of current interval.
c        bl(i,l) is a bound on the ith derivative of the normal
c        density at level l at the left end-point.
c
      fe(l) = s(l, 1)
      t = fe(l) / cond(l)
      bl(1, l) = phi(t, condl(l))
      bl(2, l) = bl(1, l) * abs(t / cond(l))
      bl(3, l) = bl(1, l) * f2(t, cond(l))
      bl(4, l) = bl(1, l) * f3(t, cond(l))
      bl(5, l) = bl(1, l) * f4(t, cond(l))
      if (simps) goto 100
      bl(6, l) = bl(1, l) * f5(t, cond(l))
      bl(7, l) = bl(1, l) * f6(t, cond(l))
c        Bounding derivatives at right end-point of sub-interval.
c        br(i,l) is a bound on the ith derivative of the normal
c        density at level l at the right end-point.
c
  100 t = r(k, l) / cond(l)
      br(1, k, l) = phi(t, condl(l))
      br(2, k, l) = br(1, k, l) * abs(t / cond(l))
      br(3, k, l) = br(1, k, l) * f2(t, cond(l))
      br(4, k, l) = br(1, k, l) * f3(t, cond(l))
      br(5, k, l) = br(1, k, l) * f4(t, cond(l))
      if (simps) goto 104
      br(6, k, l) = br(1, k, l) * f5(t, cond(l))
      br(7, k, l) = br(1, k, l) * f6(t, cond(l))
  104 r(k + 1, l) = (fe(l) + r(k, l)) * half
      bou5(l) = ep(l) * (r(k, l) - s(l, 1))
      del(2, l) = r(k + 1, l) - fe(l)
c
c        Checking the bound for the trapezoidal rule
c
      del(3, l) = two * del(2, l)
      bou1 = max(br(1, k, l), bl(1, l)) * binc(3, l) +
     *  two * max(br(2, k, l), bl(2, l)) * binc(2, l) +
     *  max(br(3, k, l), bl(3, l)) * binc(1, l)
      bou3 = bou4(l) + bou1 * (del(3, l) ** 3) * prod(l) / twelve
      itype(l) = 3
      if (bou3 .le. bou5(l)) goto 200
c
c        Checking the bound for Simpsons rule.
c
      bou1 = zero
      do 110 ij = 1, 5
        jk = 6 - ij
        bou2 = max(br(ij, k, l), bl(ij, l))
        bou1 = bou1 + bou2 * binc(jk, l) * bcs(ij)
  110 continue
      bou3 = bou4(l) + bou1 * (del(2, l) ** 5) * prod(l) / ninety
      itype(l) = 2
      if (bou3 .le. bou5(l)) goto 200
      if (simps) goto 130
c
c        Checking the bound for boules rule, if necessary.
c
      del(1, l) = half * del(2, l)
      bou1 = zero
      itype(l) = 1
      do 120 ij = 1, 7
        jk = 8 - ij
        bou2 = max(br(ij, k, l), bl(ij, l))
        bou1 = bou1 + bou2 * binc(jk, l) * bcn(ij)
  120 continue
      bou3 = bou4(l) + bou1 * (del(1, l) ** 7) * prod(l) * eight /
     *  nine45
      if (bou3 .le. bou5(l)) goto 200
c
c        Sub-dividing further at level l when the bound is too big.
c
  130 k = k + 1
      if (k .gt. 35) return
      goto 100
  200 bint(l) = bint(l) + bou3 - bou4(l)
      bou4(l) = bou3
      ksa(l) = k
      if (ind(l) .eq. 6) goto 202
      if (itype(l) - 2) 205, 206, 210
c
c        The next 30 lines condition on the value xs and go to level l+1
c
  202 ind(l) = 5
      xs = fe(l)
      fact(l) = bl(1, l)
  203 xss(nm1, l + 1) = xss(nm1, l) + bh(l) * (xs + xm(l))
      do 204 ll = l, nm2
  204 xss(ll, l + 1) = xss(ll, l) + beta(ll, l) * (xs + xm(l))
      if (l .eq. nm2) goto 300
c
c        xm is the mean of the next variable given those fixed so far.
c
      xm(l + 1) = xss(l, l + 1)
      prod(l + 1) = prod(l) * fact(l)
      l = l + 1
      goto 90
  205 ind(l) = 4
      k = ksa(l)
      xs = half * (fe(l) + r(k + 1, l))
      goto 207
  206 ind(l) = 3
      k = ksa(l)
      xs = r(k + 1, l)
  207 t = xs / cond(l)
      fact(l) = phi(t,condl(l))
      goto 203
  208 ind(l) = 2
      k = ksa(l)
      xs = half * (r(k, l) + r(k + 1, l))
      goto 207
  210 ind(l) = 1
      k = ksa(l)
      xs = r(k, l)
      fact(l) = br(1, k, l)
      goto 203
c
c        Evaluate conditional bivariate probabilities at deepest level
c
  300 x1 = fac * (xss(nm1, nm1) - d(1)) / sd(1)
      x2 = fac * (xss(nm2, nm1) - d(2)) / sd(2)
      l = nm1
      ans(l) = bivnor(x1, x2, rho)
      if (ipr) goto 310
      y1 = (xss(nm1, nm1) - c(1)) / sd(1)
      y2 = (xss(nm2, nm1) - c(2)) / sd(2)
      wu = bivnor(y1, y2, rho)
      wt = bivnor(x1, y2, rho)
      wb = bivnor(y1, x2, rho)
      ans(l) = ans(l) + wu - wt - wb
310   if (l .eq. 1) goto 340
      l = l - 1
c
c        Advancing the integration at the current level
c
      indl = ind(l)
      numl = num(l)
      ity = itype(l)
      temp = fact(l) * ans(l + 1)
      temb = bint(l + 1)
      fsa = one
      if (indl .ne. 1) goto 315
      tem = temp
      temp = temp + prep(l)
      temb = temb + preb(l)
      prep(l) = tem
      preb(l) = bint(l + 1)
      fsa = two
  315 ans(l) = ans(l) + coef(indl, ity) * temp * del(ity, l)
      bint(l) = bint(l) + coef(indl, ity) * temb * del(ity, l)
c
c        Making use of the error which did not accumulate at level l+1
c
      ep(l) = ep(l) + del(ity, l) * coef(indl, ity) * (fsa * float(nm2 -
     *  l) * epsi / pr2(l + 1) - temb) / (s(l, numl) - s(l, 1))
      if (indl .eq. 1) goto 320
      igo = indl - (1 + (itype(l) * (itype(l) - 1)) / 2)
      goto (210, 208, 206, 205), igo
c
c        Un-subdividing at level l.
  320 k = ksa(l)
      do 322 i = 1, ien
  322 bl(i, l) = br(i, k, l)
      ind(l) = 5
      fe(l) = r(k, l)
      if (k .eq. 1) goto 326
      k = k - 1
      goto 104
  326 if (intvl(l) .eq. num(l)) goto 310
      intvl(l) = intvl(l) + 1
      intl = intvl(l)
      r(1, l) = s(l, intl)
      goto 100
c
c        Completion of integration and bounding.
c
 340  ifault = 0
      prob = ans(1)
      bound = bint(1) + eplos
      return
c
c        special cases --
c        label 350 - n=1
c        label 360 - n=2
c
 350  prob = alnorm(d(1), .false.) - alnorm(c(1), .false.)
      return
 360  rho = sigma(2)
      if (abs(rho) .gt. one) goto 400
      y1 = -d(1) * fac
      y2 = -d(2) * fac
      prob = bivnor(y1, y2, rho)
      if (ipr) return
      x1 = -c(1)
      x2 = -c(2)
      wl = bivnor(x1, x2, rho)
      wt = bivnor(x1, y2, rho)
      wb = bivnor(y1, x2, rho)
      prob = wl - wt - wb + prob
      return
c
c        Error return for covariance not positive definite.
c
  400 ifault = 4
      return
      end
c
c
c
      double precision function phi(x, y)
c
c        algorithm as 195.1  appl. statist. (1984) vol.33, no.1
c
c        computes univariate normal density
c        xlow=log(smallest floating point number)
c        sq2p=log(sqrt(two*pi))
c
      double precision arg, half, sq2p, x, xlow, y, zero
c
      data xlow /-87.0D0/, sq2p /0.91893853320467274D0/, zero /0.D0/,
     *  half /0.5D0/
      phi = zero
      arg = -half * x * x - sq2p - y
      if (arg .gt. xlow) phi = exp(arg)
      return
      end
c
c
c
      DOUBLE PRECISION FUNCTION BIVNOR(AH, AK, R)
C
C     BIVNOR is a controlled precision Fortran function to calculate
C     the bivariate normal upper right area, viz. the probability for
C     two normal variates X and Y whose correlation is R, that X > AH
C     and Y > AK.
C     The accuracy is specified as the number of decimal digits, IDIG.
C
C     Reference:
C        Donnelly, T.G. (1973) Algorithm 462, Bivariate normal
C        distribution, Comm. A.C.M., vol.16, p. 638.
C
C     Corrected for the case AH = 0, AH non-zero by reversing AH & AK
C     in such cases.
C
      DOUBLE PRECISION AH, AK, R
C
C     Local variables
C
      DOUBLE PRECISION TWOPI, B, XAH, XAK, GH, GK, RR, GAUSS, DERF, H2,
     +   A2, H4, EX, W2, AP, S2, SP, S1, SN, SGN, SQR, CON, WH, WK, GW,
     +   T, G2, CONEX, CN, TWO, ZERO, ONE, FOUR, QUART, HALF, EXPLIM
      INTEGER IDIG, IS
      DATA TWO/2.D0/, ZERO/0.D0/, ONE/1.D0/, FOUR/4.D0/, QUART/0.25D0/,
     +   HALF/0.5D0/
      GAUSS(T) = (ONE + DERF(T/SQRT(TWO)))/TWO
C
C     GAUSS is a univariate lower normal tail area calculated here from
C     the central error function, DERF.
C     It may be replaced by the algorithm in Hill, I.D. and Joyce, S.A.
C     Algorithm 304, Normal curve integral (S15), Comm. A.C.M. (10),
C     (June 1967), p.374 or with Applied Statistics algorithm AS66.
C     Source: Owen, D.B. Ann. Math. Statist., vol.27 (1956), p.1075.
C
      DATA TWOPI/6.2831 85307 17958 7D0/, IDIG/15/, EXPLIM/80.D0/
C
      B = ZERO
      IF (AH .EQ. ZERO) THEN
        XAH = AK
        XAK = AH
      ELSE
        XAH = AH
        XAK = AK
      END IF

      GH = GAUSS(-XAH) / TWO
      GK = GAUSS(-XAK) / TWO
      IF (R .EQ. ZERO) THEN
        B = FOUR * GH * GK
        GO TO 350
      END IF

      RR = ONE - R*R
      IF (RR .LT. ZERO) THEN
        WRITE(*, *)'Error in BIVNOR, R = ', R
        GO TO 390
      END IF
      IF (RR .GT. ZERO) GO TO 100
C
C     R^2 = 1.0
C
      IF (R .GE. ZERO) GO TO 70
      IF (XAH + XAK .GE. ZERO) GO TO 350
      B = TWO * (GH + GK) - ONE
      GO TO 350
   70 IF (XAH - XAK .LT. ZERO) THEN
        B = TWO * GK
      ELSE
        B = TWO * GH
      END IF
      GO TO 350
C
C     Regular case, R^2 < 1
C
  100 SQR = SQRT(RR)
      IF (IDIG .EQ. 15) THEN
        CON = TWOPI * 1.D-15 / TWO
      ELSE
        CON = TWOPI / TWO / 10**IDIG
      END IF

      IF (XAH .NE. ZERO) GO TO 170
      IF (XAK .NE. ZERO) GO TO 190
      B = ATAN(R/SQR) / TWOPI + QUART
      GO TO 350
  170 B = GH

      IF (XAH*XAK) 180, 200, 190
  180 B = B - HALF
  190 B = B + GK
  200 WH = -XAH
      WK = (XAK/XAH - R)/SQR
      GW = TWO * GH
      IS = -1
  210 SGN = -ONE
      T = ZERO
      IF (WK .EQ. ZERO) GO TO 320
      IF (ABS(WK) - ONE) 270, 230, 240
  230 T = WK * GW * (ONE - GW) / TWO
      GO TO 310
  240 SGN = -SGN
      WH = WH * WK
      G2 = GAUSS(WH)
      WK = ONE / WK
      IF (WK .LT. ZERO) B = B + HALF
      B = B - (GW + G2)/TWO + GW*G2
  270 H2 = WH * WH
      A2 = WK * WK
      H4 = H2 / TWO
      IF (H4 .LT. EXPLIM) THEN
        EX = EXP(-H4)
      ELSE
        EX = ZERO
      END IF
      W2 = H4 * EX
      AP = ONE
      S2 = AP - EX
      SP = AP
      S1 = ZERO
      SN = S1
      CONEX = ABS(CON / WK)
      GO TO 290

  280 SN = SP
      SP = SP + ONE
      S2 = S2 - W2
      W2 = W2 * H4 / SP
      AP = -AP*A2
  290 CN = AP * S2 / (SN + SP)
      S1 = S1 + CN
      IF (ABS(CN) - CONEX .GT. ZERO) GO TO 280

      T = (ATAN(WK) - WK*S1) / TWOPI
  310 B = B + SGN*T
  320 IF (IS .GE. 0) GO TO 350
      IF (XAK .NE. ZERO) THEN
        WH = -XAK
        WK = (XAH/XAK - R) / SQR
        GW = TWO * GK
        IS = 1
        GO TO 210
      END IF

  350 IF (B .LT. ZERO) B = ZERO
      IF (B .GT. ONE) B = ONE
  390 BIVNOR = B
C
      RETURN
      END

C
C
C
      SUBROUTINE INVERT(A, AINV, C, N, DET, IFAULT)
C
C     Invert the NxN symmetric matrix, A, of which only the lower
C     triangle is stored by rows.   The inverse is returned in AINV.
C     The inverse of the last I rows and columns, for I = 1, 2, ...,
C     N-2, is returned in C(I,*).   DET = the determinant of A.
C     IFAULT = 4 if A is not positive definite.
C
      DOUBLE PRECISION A(*), AINV(*), C(5,*), DET
      INTEGER N, IFAULT
C
C     Local variables
C
      DOUBLE PRECISION ONE, WK(15), W(7)
      INTEGER NULLTY, NN, I, ROW1, WKPOS, ROW, APOS, COL
      DATA ONE/1.D0/
C
      DO 50 I = 1, N-2
C
C     Copy the last I rows and columns of A into WK.
C
        NN = I * (I + 1) / 2
        ROW1 = N + 1 - I
        WKPOS = 1
        DO 20 ROW = ROW1, N
          APOS = ROW * (ROW - 1) / 2 + ROW1
          DO 10 COL = ROW1, ROW
            WK(WKPOS) = A(APOS)
            WKPOS = WKPOS + 1
            APOS = APOS + 1
   10     CONTINUE
   20   CONTINUE
C
C     Call SYMINV to invert WK in situ.
C
        CALL SYMINV(WK, I, NN, WK, W, NULLTY, IFAULT)
        IF (IFAULT .NE. 0) GO TO 70
C
C     Copy the inverse into the I-th row of C.
C
        DO 40 APOS = 1, NN
   40   C(I, APOS) = WK(APOS)
   50 CONTINUE
C
C     Use AS6 (CHOL) to calculate the determinant.
C     This is inefficient as AS7 will call AS6 again to repeat the same
C     calculations.   AS6/7 can easily be converted to calculate the
C     determinant.
C
      NN = N * (N + 1) / 2
      CALL CHOL(A, N, NN, AINV, NULLTY, IFAULT)
      IF (IFAULT .NE. 0) GO TO 70
      DET = ONE
      APOS = 1
      DO 60 ROW = 1, N
        DET = DET * AINV(APOS)
        APOS = APOS + ROW + 1
   60 CONTINUE
      DET = DET * DET
C
C     Invert the full matrix.
C
      CALL SYMINV(A, N, NN, AINV, W, NULLTY, IFAULT)
      IF (IFAULT .EQ. 0) RETURN
C
   70 IFAULT = 4
      RETURN
      END
C
C
C
      SUBROUTINE CALERF(ARG,RESULT,JINT)
C------------------------------------------------------------------
C
C THIS PACKET COMPUTES THE ERROR AND COMPLEMENTARY ERROR FUNCTIONS
C   FOR REAL ARGUMENTS  ARG.  IT CONTAINS TWO FUNCTION TYPE
C   SUBPROGRAMS,  ERF  AND  ERFC  (OR  DERF  AND  DERFC),  AND ONE
C   SUBROUTINE TYPE SUBPROGRAM,  CALERF.  THE CALLING STATEMENTS
C   FOR THE PRIMARY ENTRIES ARE
C
C                   Y=ERF(X)     (OR   Y=DERF(X) )
C   AND
C                   Y=ERFC(X)    (OR   Y=DERFC(X) ).
C
C   THE ROUTINE  CALERF  IS INTENDED FOR INTERNAL PACKET USE ONLY,
C   ALL COMPUTATIONS WITHIN THE PACKET BEING CONCENTRATED IN THIS
C   ROUTINE.  THE FUNCTION SUBPROGRAMS INVOKE  CALERF  WITH THE
C   STATEMENT
C          CALL CALERF(ARG,RESULT,JINT)
C   WHERE THE PARAMETER USAGE IS AS FOLLOWS
C
C      FUNCTION                     PARAMETERS FOR CALERF
C       CALL              ARG                  RESULT          JINT
C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
C     ERFC(ARG)     ABS(ARG) .LT. XMAX        ERFC(ARG)         1
C
C   THE MAIN COMPUTATION EVALUATES NEAR MINIMAX APPROXIMATIONS
C   FROM "RATIONAL CHEBYSHEV APPROXIMATIONS FOR THE ERROR FUNCTION"
C   BY W. J. CODY, MATH. COMP., 1969, PP. 631-638.  THIS
C   TRANSPORTABLE PROGRAM USES RATIONAL FUNCTIONS THAT THEORETICALLY
C   APPROXIMATE  ERF(X)  AND  ERFC(X)  TO AT LEAST 18 SIGNIFICANT
C   DECIMAL DIGITS.  THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC
C   SYSTEM, THE COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER
C   SELECTION OF THE MACHINE-DEPENDENT CONSTANTS.
C
C*******************************************************************
C
C EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
C
C   XSMALL = ARGUMENT BELOW WHICH ERF(X) MAY BE REPRESENTED
C            BY   2*X/SQRT(PI)  AND ABOVE WHICH  X*X  WILL
C            NOT UNDERFLOW.  A CONSERVATIVE VALUE IS THE
C            LARGEST X SUCH THAT   1.0 + X = 1.0   TO MACHINE
C            PRECISION.
C   XMAX   = LARGEST ARGUMENT ACCEPTABLE TO  ERFC;  SOLUTION TO
C            EQUATION:  W(X) * (1-0.5/X**2) = XMIN,  WHERE
C            W(X) = EXP(-X*X)/(X*SQRT(PI)),  AND XMIN IS THE
C            SMALLEST POSITIVE MACHINE NUMBER (SEE TABLE BELOW).
C
C     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
C
C                          XSMALL     XMAX     XMIN
C
C    IBM 195     (D.P.)   1.39D-17   13.306   5.40D-79
C    CDC 7600    (S.P.)   7.11E-15   25.922   3.13E-294
C    CRAY-1      (S.P.)   7.11E-15   75.326   4.58E-2467
C    UNIVAC 1108 (D.P.)   1.73D-18   26.582   2.78D-309
C    VAX 11/780  (S.P.)   5.96E-8     9.269   2.94E-39
C    VAX 11/780  (D.P.)   1.39D-17    9.269   2.94D-39
C    IBM PC      (S.P.)   5.96E-8     9.194   1.18E-38
C    IBM PC      (D.P.)   1.11D-16   26.543   2.23D-308
C
C*******************************************************************
C
C ERROR RETURNS
C
C  THE PROGRAM RETURNS  ERFC = 0  FOR  ARG .GT. XMAX.
C
C
C OTHER SUBPROGRAMS REQUIRED (SINGLE PRECISION VERSION)
C
C     ABS, EXP
C
C OTHER SUBPROGRAMS REQUIRED (DOUBLE PRECISION VERSION)
C
C     DABS, DEXP
C
C
C  AUTHOR: W. J. CODY
C          MATHEMATICS AND COMPUTER SCIENCE DIVISION
C          ARGONNE NATIONAL LABORATORY
C          ARGONNE, IL 60439
C
C  LATEST MODIFICATION: JANUARY 8, 1985
C
C------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL             A,ARG,B,C,D,FOUR,HALF,P,ONE,Q,RESULT,SQRPI,
CS   1               TWO,THRESH,X,XMAX,XDEN,XNUM,XSMALL,Y,YSQ,ZERO
      DOUBLE PRECISION A,ARG,B,C,D,FOUR,HALF,P,ONE,Q,RESULT,SQRPI,
     1               TWO,THRESH,X,XMAX,XDEN,XNUM,XSMALL,Y,YSQ,ZERO
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
C------------------------------------------------------------------
C  MATHEMATICAL CONSTANTS
C------------------------------------------------------------------
CS    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/
CS    DATA SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/
      DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/
      DATA SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/
C------------------------------------------------------------------
C  MACHINE-DEPENDENT PARAMETERS
C------------------------------------------------------------------
CS    DATA XSMALL/4.2E-16/, XMAX/9.269E0/
      DATA XSMALL/4.2D-16/, XMAX/9.269D0/
C------------------------------------------------------------------
C  COEFFICIENTS FOR APPROXIMATION TO DERF IN FIRST INTERVAL
C------------------------------------------------------------------
CS    DATA A/3.16112374387056560E00,1.13864154151050156E02,
CS   1       3.77485237685302021E02,3.20937758913846947E03,
CS   2       1.85777706184603153E-1/
CS    DATA B/2.36012909523441209E01,2.44024637934444173E02,
CS   1       1.28261652607737228E03,2.84423683343917062E03/
      DATA A/3.16112374387056560D00,1.13864154151050156D02,
     1       3.77485237685302021D02,3.20937758913846947D03,
     2       1.85777706184603153D-1/
      DATA B/2.36012909523441209D01,2.44024637934444173D02,
     1       1.28261652607737228D03,2.84423683343917062D03/
C------------------------------------------------------------------
C  COEFFICIENTS FOR APPROXIMATION TO DERFC IN SECOND INTERVAL
C------------------------------------------------------------------
CS    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
CS   1       6.61191906371416295E01,2.98635138197400131E02,
CS   2       8.81952221241769090E02,1.71204761263407058E03,
CS   3       2.05107837782607147E03,1.23033935479799725E03,
CS   4       2.15311535474403846E-8/
CS    DATA D/1.57449261107098347E01,1.17693950891312499E02,
CS   1       5.37181101862009858E02,1.62138957456669019E03,
CS   2       3.29079923573345963E03,4.36261909014324716E03,
CS   3       3.43936767414372164E03,1.23033935480374942E03/
      DATA C/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/
      DATA D/1.57449261107098347D01,1.17693950891312499D02,
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/
C------------------------------------------------------------------
C  COEFFICIENTS FOR APPROXIMATION TO DERFC IN THIRD INTERVAL
C------------------------------------------------------------------
CS    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
CS   1       1.25781726111229246E-1,1.60837851487422766E-2,
CS   2       6.58749161529837803E-4,1.63153871373020978E-2/
CS    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
CS   1       5.27905102951428412E-1,6.05183413124413191E-2,
CS   2       2.33520497626869185E-3/
      DATA P/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/
      DATA Q/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/
C------------------------------------------------------------------
      X = ARG
CS    Y = ABS(X)
      Y = DABS(X)
      IF (Y .GT. FOUR) GO TO 200
      IF (Y .GT. THRESH) GO TO 100
C------------------------------------------------------------------
C  EVALUATE ERF FOR ABS(X) .LE. 0.46875
C------------------------------------------------------------------
      YSQ = ZERO
      IF (Y .GT. XSMALL) YSQ = Y * Y
      XNUM = A(5)*YSQ
      XDEN = YSQ
      DO 20 I = 1, 3
         XNUM = (XNUM + A(I)) * YSQ
         XDEN = (XDEN + B(I)) * YSQ
   20 CONTINUE
      RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
      IF (JINT .NE. 0) RESULT = ONE - RESULT
      GO TO 800
C------------------------------------------------------------------
C  EVALUATE ERFC FOR 0.46875 .LT. ABS(X) .LE. 4.0
C------------------------------------------------------------------
  100 YSQ = Y * Y
      XNUM = C(9)*Y
      XDEN = Y
      DO 120 I = 1, 7
         XNUM = (XNUM + C(I)) * Y
         XDEN = (XDEN + D(I)) * Y
  120 CONTINUE
CS    RESULT = EXP(-YSQ) * (XNUM + C(8)) / (XDEN + D(8))
      RESULT = DEXP(-YSQ) * (XNUM + C(8)) / (XDEN + D(8))
      GO TO 300
C------------------------------------------------------------------
C  EVALUATE ERFC FOR ABS(X) .GT. 4.0
C------------------------------------------------------------------
  200 RESULT = ZERO
      IF (Y .GE. XMAX) GO TO 300
  220 YSQ = ONE / (Y * Y)
      XNUM = P(6)*YSQ
      XDEN = YSQ
      DO 240 I = 1, 4
         XNUM = (XNUM + P(I)) * YSQ
         XDEN = (XDEN + Q(I)) * YSQ
  240 CONTINUE
      RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
CS    RESULT = (EXP(-Y*Y) / Y) * (SQRPI - RESULT)
      RESULT = (DEXP(-Y*Y) / Y) * (SQRPI - RESULT)
C------------------------------------------------------------------
C  FIX UP FOR NEG. ARG., ERF, ETC.
C------------------------------------------------------------------
  300 IF (JINT .EQ. 0) GO TO 350
      IF (X .LT. ZERO) RESULT = TWO - RESULT
      GO TO 800
  350 RESULT = (HALF - RESULT) + HALF
      IF (X .LT. ZERO) RESULT = -RESULT
C------------------------------------------------------------------
  800 RETURN
C---------- LAST CARD OF CALERF ----------
      END




CS    REAL FUNCTION ERF(X)
      DOUBLE PRECISION FUNCTION DERF(X)
C------------------------------------------------------------------
C
C  PROGRAM TO COMPUTE THE ERROR FUNCTION
C
C   AUTHOR - W. J. CODY
C
C   DATE - JANUARY 8, 1985
C
C------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 0
      CALL CALERF(X,RESULT,JINT)
CS    ERF = RESULT
      DERF = RESULT
      RETURN
C---------- LAST CARD OF DERF ----------
      END




CS    REAL FUNCTION ERFC(X)
      DOUBLE PRECISION FUNCTION DERFC(X)
C------------------------------------------------------------------
C
C  PROGRAM TO COMPUTE THE COMPLEMENTARY ERROR FUNCTION
C
C   AUTHOR - W. J. CODY
C
C   DATE - JANUARY 8, 1985
C
C------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 1
      CALL CALERF(X,RESULT,JINT)
CS    ERFC = RESULT
      DERFC = RESULT
      RETURN
C---------- LAST CARD OF DERFC ----------
      END
c
c-------------------------------------------------------------------------
c
      subroutine mulnor(a, b, sig, eps, n, inf, prob, bound, ifault)
c
c        algorithm as 195  appl. statist. (1984) vol.33, no.1
c
c        computes multivariate normal distribution function and computes
c        the probability that a multivariate normal vector falls in a
c        rectangle in n-space with error less than eps.
c
c     Auxiliary functions required: ANORDF, ANORIN, BNRDF, LFDRG, LFTRG
c     and LINRG from the IMSL Stat/Library.
c
        dimension
     *     a(*), b(*), sig(*), c(7), d(7), co(25), sd(2), coef(5, 3),
     *  binc(7, 5), bl(7, 5), br(7, 35, 5), r(36, 5), s(5, 27),
     *  xss(6, 6), pr2(6), prep(6), preb(6), cv(5, 15), sinv(28),
     *  condl(5), xm(5), cond(5), beta(5, 5), bh(5), fe(5), sigma(28),
     *  ep(5), del(3, 5), bou4(5), bou5(5), ans(6), fact(5), prod(5),
     *  bint(6), bcs(5), bcn(7), eo(13), do(25)
      dimension inf(7), intvl(5), ind(5), ksa(5), num(5), itype(5)
      logical simps, ipr
      data coef(1, 1), coef(2, 1), coef(3, 1), coef(4, 1), coef(5, 1),
     *  coef(1, 2), coef(2, 2), coef(3, 2), coef(4, 2), coef(5, 2),
     *  coef(1, 3), coef(2, 3), coef(3, 3), coef(4, 3), coef(5, 3)
     *  /0.311111111111111, 1.422222222222222, 0.533333333333333,
     *  1.422222222222222, 0.311111111111111, 0.333333333333333, 0.0,
     *  1.333333333333333, 0.0, 0.333333333333333, 0.5, 0.0, 0.0, 0.0,
     *  0.5/
      data bcs(1), bcs(2), bcs(3), bcs(4), bcs(5) /1.0, 4.0, 6.0, 4.0,
     *  1.0/
      data bcn(1), bcn(2), bcn(3), bcn(4), bcn(5), bcn(6), bcn(7)
     *  /1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0/
      data do(1), do(2), do(3), do(4), do(5), do(6), do(7), do(8),
     *  do(9), do(10), do(11), do(12), do(13), do(14), do(15), do(16),
     *  do(17), do(18), do(19), do(20), do(21), do(22), do(23), do(24),
     *  do(25) / - 3.75043972,-3.3242574, -2.85697, -2.36675941,
     *  -2.3344142, -1.8891759, -1.7320508, -1.3556262,-1.15440539,-1.,
     *  -0.74196378, -0.61670659,0.0,
     *  0.61670659, 0.74196378, 1.0, 1.15440539,1.3556262,  1.7320508,
     *  1.8891759,2.3344142, 2.36675941,2.85697, 3.3242574,3.75043972/
      data eo(1), eo(2), eo(3), eo(4), eo(5), eo(6), eo(7), eo(8),
     *  eo(9), eo(10), eo(11), eo(12), eo(13) / - 2.85697, -2.3344142,
     *  -1.7320508, -1.3556262, -1.0, -0.74196378, 0.0, 0.74196378, 1.0,
     *  1.3556262, 1.7320508, 2.3344142, 2.85697/
      data cons2 /6.879833/, cons3 /4.517004/, cons4 /76.371214/,
     *  epsmin /1.0e-8/, epssim /6.0e-5/
      data zero, p05, p15, half, one, two, three, six, eight, ten,
     *  twelve, fiften, forty5, ninety, nine45 /0.0, 0.05, 0.15, 0.5,
     *  1.0, 2.0, 3.0, 6.0, 8.0, 10.0, 12.0, 15.0, 45.0, 90, 945/
c
c        functions needed to calculate the first six derivatives of
c        the normal density
c
      f2(x, y) = abs((x * x - one) / (y * y))
      f3(x, y) = abs(-x * (x * x - three) / (y * y * y))
      f4(x, y) = abs((three + x * x * (-six + x * x)) / (y ** 4))
      f5(x, y) = abs(x * (-fiften + x * x * (ten - x * x))) / (y ** 5)
      f6(x, y) = abs(-fiften + x * x * (forty5 - x * x * (-fiften +
     *  x * x))) / (y ** 6)
c        checking for faulty data
c
      ifault = 0
      if (eps .le. epsmin) ifault = 2
      if (n .le. 0 .or. n .gt. 7) ifault = 3
      if (ifault .ne. 0) return
      do 1 i = 1, n
    1 if (inf(i) .eq. 2 .and. a(i) .lt. b(i)) ifault = 100 + i
      if (ifault .ne. 0) return
      bound = eps
c        finding z such that p(n(0,1).gt.z).lt. 0.15*eps/n
c        co will contain the roots of the first 5 or 7 hermite
c        polynomials
      ifault=0.
      ept = eps * p15 / float(n)
      z = -anorin(ept) + epsmin
      cup=1-anordf(z)
c
c        inverting sig and the n-2 lower right hand principal minors
c
      ik = 0
      ij = 0
      do 3 i = 1, n
        do 3 j = 1, i
          ik = ik + 1
          if (i .eq. j) goto 2
          ij = ij + 1
          sigma(ik) = sig(ij)
          goto 3
    2     sigma(ik) = one
    3 continue
      if (n .le. 2) goto 4
      call invert(sigma, sinv, cv, n, det, ifault)
      if (ifault .ne. 0) return
      simps = .true.
      if (det .lt. p05 .or. eps .le. epssim) simps = .false.
      prob = zero
      det = sinv(1) * sinv(3) - sinv(2) * sinv(2)
      sd(1) = sqrt(sinv(3) / det)
      sd(2) = sqrt(sinv(1) / det)
      rho = -sinv(2) / (sd(1) * sd(2) * det)
      if (abs(rho) .gt. one) goto 400
    4 nm2 = n - 2
      nm1 = n - 1
c
c        checking whether upper and lower endpoints are too big
c
      eplos = zero
      do 6 l = 1, n
        c(l) = max(b(l), -z)
        d(l) = min(a(l), z)
        if (inf(l) .eq. 0) d(l) = z
        if (inf(l) .eq. 1) c(l) = -z
        if (a(l) .gt. z .or. inf(l) .eq. 0) eplos = eplos + cup
        if (b(l) .lt. -z .or. inf(l) .eq. 1) eplos = eplos + cup
        if (c(l) .ge. d(l)) return
    6 continue
      if (n .eq. 1) goto 350
      fac = one
      ipr = .false.
      if (inf(1) .ne. 1 .or. inf(2) .ne. 1) goto 7
      ipr = .true.
      eplos = eplos - two * cup
      goto 8
    7 if (inf(1) .ne. 0 .or. inf(2) .ne. 0) goto 8
      fac = -one
      ipr = .true.
      d(1) = c(1)
      d(2) = c(2)
      eplos = eplos - two * cup
    8 if (n .eq. 2) goto 360
      ifault = 5
      epsi = (eps - eplos) / float(nm2)
c        finding regression coefficients (beta,bh) and bounds on the
c        conditional integrals (binc)
c        cond(l)=conditional variance of variable n-l+1 given later
c        variables
c
      do 15 l = 1, nm2
        cond(l) = one / sqrt(cv(l, 1))
        condl(l) = log(cond(l))
        do 10 i = 1, l
        beta(l, i) = zero
        do 10 j = 1, l
          jk = (l - i + 1) * (l - i) / 2 + j
          if (j .gt. l - i + 1) jk = j * (j - 1) / 2 + l - i + 1
          jn = (n - l + j) * (n - l + j - 1) / 2 + n - l
          beta(l, i) = beta(l, i) + sigma(jn) * cv(l, jk)
   10   continue
        k = n - l - 1
        bh(k + 1) = beta(l, l)
        do 11 i = 1, k
          bh(i) = zero
          do 11 j = 1, l
            jn = (j + n - l) * (j + n - l - 1) / 2 + n - l
            ijk = 1 + j * (j - 1) / 2
            bh(i) = bh(i) + sigma(jn) * cv(l, ijk)
   11   continue
        k = 0
        sigc = zero
        do 12 j = 1, l
          do 12 i = 1, j
            k = k + 1
            sigc = sigc + bh(i) * bh(j) * sinv(k)
   12   continue
        binc(1, l) = one
        binc(2, l) = sqrt(sigc)
        binc(3, l) = two * sigc
        binc(4, l) = cons2 * sigc * binc(2, l)
        binc(5, l) = twelve * sigc * sigc
        if (simps) goto 13
        binc(6, l) = sigc * sigc * binc(1, l) * cons3
        binc(7, l) = (sigc ** 3) * cons4
   13   if (l .lt. nm2) goto 15
        do 14 i = 1, nm2
          bh(i) = zero
          do 14 j = 1, nm2
            jk = (l - i + 1) * (l - i) / 2 + j
            if (j .gt. l - i + 1) jk = j * (j - 1) / 2 + l - i + 1
            jn = (2 + j) * (1 + j) / 2 + 1
            bh(i) = bh(i) + sigma(jn) * cv(l, jk)
   14   continue
   15 continue
      l = 1
c
c        co will contain the roots of the first 5 or 7 hermite
c        polynomials
c
      if (simps) goto 50
      do 40 i = 1, 25
   40 co(i) = do(i)
      iend = 25
      ien = 7
      goto 60
   50 do 55 i = 1, 13
   55 co(i) = eo(i)
      iend = 13
      ien = 5
c
c        initialising values.   xss contains partial sums used for
c        calculating conditional means.
c
   60 do 70 i = 1, nm1
   70 xss(i, 1) = zero
      xm(1) = zero
      prod(1) = one
      pr2(1) = one
      do 80 i = 1, nm2
        ni = n - i + 1
        pr2(i + 1) = pr2(i) * (d(ni) - c(ni))
   80 continue
c
c        bint(l) is a bound on the error accumulated at levels l and
c                deeper.
c        ans(l)  is the accumulated integral at level l.
c        prep(l) contains the integrand at level l.
c        preb(l) bounds the error accumulated at levels deeper than l.
c
      bint(nm1) = zero
   90 intvl(l) = 2
      ans(l) = zero
      bou4(l) = zero
      bint(l) = zero
      prep(l) = zero
      preb(l) = zero
      k = 1
c
c        finding which of the co are in the current interval.
c        s(l,.) are the endpoints of intervals on which the integrand
c        at level l and its derivatives are monotone.
c        num(l) is the number of such intervals.
c
      nl = n - l + 1
      s(l, 1) = c(nl) - xm(l)
      s(l, iend + 2) = d(nl) - xm(l)
      num(l) = iend + 2
      do 91 i = 1, iend
        njs = i
        if (s(l, 1) .lt. co(i) * cond(l)) goto 92
        num(l) = num(l) - 1
   91 continue
   92 if (num(l) .eq. 2) goto 99
      do 94 i = njs, iend
        mjs = iend - i + njs
        if (s(l, iend + 2) .ge. co(mjs) * cond(l)) goto 96
        num(l) = num(l) - 1
   94 continue
   96 if (num(l) .eq. 2) goto 99
      do 98 i = njs, mjs
        inj = i - njs + 2
        s(l, inj) = co(i) * cond(l)
   98 continue
   99 numl = num(l)
      s(l, numl) = s(l, iend + 2)
c
c        ep(l)  is an upper limit on the allowable error at level l.
c        r(k,l) is the right end-point of the current sub-interval.
c        fe(l)  is the left end-point of the current sub-interval.
c        ind(l)-1  indicates which point of the newton-cotes formula
c        we are dealing with.
c
      ep(l) = epsi / pr2(l + 1)
      r(1, l) = s(l, 2)
      ind(l) = 6
c
c        bounding derivatives at left end-point of current interval.
c        bl(i,l) is a bound on the ith derivative of the normal
c        density at level l at the left end-point.
c
      fe(l) = s(l, 1)
      t = fe(l) / cond(l)
      bl(1, l) = phi(t, condl(l))
      bl(2, l) = bl(1, l) * abs(t / cond(l))
      bl(3, l) = bl(1, l) * f2(t, cond(l))
      bl(4, l) = bl(1, l) * f3(t, cond(l))
      bl(5, l) = bl(1, l) * f4(t, cond(l))
      if (simps) goto 100
      bl(6, l) = bl(1, l) * f5(t, cond(l))
      bl(7, l) = bl(1, l) * f6(t, cond(l))
c        bounding derivatives at right end-point of sub-interval.
c        br(i,l) is a bound on the ith derivative of the normal
c        density at level l at the right end-point.
c
  100 t = r(k, l) / cond(l)
      br(1, k, l) = phi(t, condl(l))
      br(2, k, l) = br(1, k, l) * abs(t / cond(l))
      br(3, k, l) = br(1, k, l) * f2(t, cond(l))
      br(4, k, l) = br(1, k, l) * f3(t, cond(l))
      br(5, k, l) = br(1, k, l) * f4(t, cond(l))
      if (simps) goto 104
      br(6, k, l) = br(1, k, l) * f5(t, cond(l))
      br(7, k, l) = br(1, k, l) * f6(t, cond(l))
  104 r(k + 1, l) = (fe(l) + r(k, l)) * half
      bou5(l) = ep(l) * (r(k, l) - s(l, 1))
      del(2, l) = r(k + 1, l) - fe(l)
c
c        checking the bound for the trapezoidal rule
c
      del(3, l) = two * del(2, l)
      bou1 = max(br(1, k, l), bl(1, l)) * binc(3, l) +
     *  two * max(br(2, k, l), bl(2, l)) * binc(2, l) +
     *  max(br(3, k, l), bl(3, l)) * binc(1, l)
      bou3 = bou4(l) + bou1 * (del(3, l) ** 3) * prod(l) / twelve
      itype(l) = 3
      if (bou3 .le. bou5(l)) goto 200
c
c        checking the bound for simpsons rule.
c
      bou1 = zero
      do 110 ij = 1, 5
        jk = 6 - ij
        bou2 = max(br(ij, k, l), bl(ij, l))
        bou1 = bou1 + bou2 * binc(jk, l) * bcs(ij)
  110 continue
      bou3 = bou4(l) + bou1 * (del(2, l) ** 5) * prod(l) / ninety
      itype(l) = 2
      if (bou3 .le. bou5(l)) goto 200
      if (simps) goto 130
c
c        checking the bound for boules rule, if necessary.
c
      del(1, l) = half * del(2, l)
      bou1 = zero
      itype(l) = 1
      do 120 ij = 1, 7
        jk = 8 - ij
        bou2 = max(br(ij, k, l), bl(ij, l))
        bou1 = bou1 + bou2 * binc(jk, l) * bcn(ij)
  120 continue
      bou3 = bou4(l) + bou1 * (del(1, l) ** 7) * prod(l) * eight /
     *  nine45
      if (bou3 .le. bou5(l)) goto 200
c
c        sub-dividing further at level l when the bound is too big.
c
  130 k = k + 1
      if (k .gt. 35) return
      goto 100
  200 bint(l) = bint(l) + bou3 - bou4(l)
      bou4(l) = bou3
      ksa(l) = k
      if (ind(l) .eq. 6) goto 202
      if (itype(l) - 2) 205, 206, 210
c
c        the next 30 lines condition on the value xs and go to level l+1
c
  202 ind(l) = 5
      xs = fe(l)
      fact(l) = bl(1, l)
  203 xss(nm1, l + 1) = xss(nm1, l) + bh(l) * (xs + xm(l))
      do 204 ll = l, nm2
  204 xss(ll, l + 1) = xss(ll, l) + beta(ll, l) * (xs + xm(l))
      if (l .eq. nm2) goto 300
c
c        xm is the mean of the next variable given those fixed so far.
c
      xm(l + 1) = xss(l, l + 1)
      prod(l + 1) = prod(l) * fact(l)
      l = l + 1
      goto 90
  205 ind(l) = 4
      k = ksa(l)
      xs = half * (fe(l) + r(k + 1, l))
      goto 207
  206 ind(l) = 3
      k = ksa(l)
      xs = r(k + 1, l)
  207 t = xs / cond(l)
      fact(l) = phi(t,condl(l))
      goto 203
  208 ind(l) = 2
      k = ksa(l)
      xs = half * (r(k, l) + r(k + 1, l))
      goto 207
  210 ind(l) = 1
      k = ksa(l)
      xs = r(k, l)
      fact(l) = br(1, k, l)
      goto 203
c
c        evaluate conditional bivariate probabilities at deepest level
c
  300 x1 = fac * (xss(nm1, nm1) - d(1)) / sd(1)
      x2 = fac * (xss(nm2, nm1) - d(2)) / sd(2)
      l = nm1
      ans(l) = bivnor(x1, x2, rho)
      if (ipr) goto 310
      y1 = (xss(nm1, nm1) - c(1)) / sd(1)
      y2 = (xss(nm2, nm1) - c(2)) / sd(2)
      wu = bivnor(y1, y2, rho)
      wt = bivnor(x1, y2, rho)
      wb = bivnor(y1, x2, rho)
      ans(l) = ans(l) + wu - wt - wb
310   if (l .eq. 1) goto 340
      l = l - 1
c
c        advancing the integration at the current level
c
      indl = ind(l)
      numl = num(l)
      ity = itype(l)
      temp = fact(l) * ans(l + 1)
      temb = bint(l + 1)
      fsa = one
      if (indl .ne. 1) goto 315
      tem = temp
      temp = temp + prep(l)
      temb = temb + preb(l)
      prep(l) = tem
      preb(l) = bint(l + 1)
      fsa = two
  315 ans(l) = ans(l) + coef(indl, ity) * temp * del(ity, l)
      bint(l) = bint(l) + coef(indl, ity) * temb * del(ity, l)
c
c        making use of the error which did not accumulate at level l+1
c
      ep(l) = ep(l) + del(ity, l) * coef(indl, ity) * (fsa * float(nm2 -
     *  l) * epsi / pr2(l + 1) - temb) / (s(l, numl) - s(l, 1))
      if (indl .eq. 1) goto 320
      igo = indl - (1 + (itype(l) * (itype(l) - 1)) / 2)
      goto (210, 208, 206, 205), igo
c
c        un-subdividing at level l.
  320 k = ksa(l)
      do 322 i = 1, ien
  322 bl(i, l) = br(i, k, l)
      ind(l) = 5
      fe(l) = r(k, l)
      if (k .eq. 1) goto 326
      k = k - 1
      goto 104
  326 if (intvl(l) .eq. num(l)) goto 310
      intvl(l) = intvl(l) + 1
      intl = intvl(l)
      r(1, l) = s(l, intl)
      goto 100
c
c        completion of integration and bounding.
c
 340  ifault = 0
      prob = ans(1)
      bound = bint(1) + eplos
      return
c
c        special cases --
c        label 350 - n=1
c        label 360 - n=2
c
350   prob=anordf(d(1))-anordf(c(1))
      return
  360 rho = sigma(2)
      if (abs(rho) .gt. one) goto 400
      y1 = -d(1) * fac
      y2 = -d(2) * fac
      prob = bivnor(y1, y2, rho)
      if (ipr) return
      x1 = -c(1)
      x2 = -c(2)
      wl = bivnor(x1, x2, rho)
      wt = bivnor(x1, y2, rho)
      wb = bivnor(y1, x2, rho)
      prob = wl - wt - wb + prob
      return
c
c        error return for covariance not positive definite.
c
  400 ifault = 4
      return
      end
c
c
c
      function phi(x, y)
c
c        algorithm as 195.1  appl. statist. (1984) vol.33, no.1
c
c        computes univariate normal density
c        xlow=log(smallest floating point number)
c        sq2p=log(sqrt(two*pi))
c
      real arg, half, sq2p, x, xlow, y, zero
c
      real exp
c
      data xlow /-87.0/, sq2p /0.91893853320467274/, zero /0.0/,
     *  half /0.5/
      phi = zero
      arg = -half * x * x - sq2p - y
      if (arg .gt. xlow) phi = exp(arg)
      return
      end
c
c
c
        real function bivnor(x,y,r)
        data range/25.0/
c
        xx=max(-range,min(range,x))
        yy=max(-range,min(range,y))
        bivnor=1+bnrdf(xx,yy,r)-anordf(xx)-anordf(yy)
        return
        end
c
c
c
        subroutine invert(a,ai,c,n,det,ier)
        dimension a(*),ai(*),c(5,1),s(5,5),ipiv(20),aa(5,5)
c
        l=0
        do 3 i=1,n
          do 3 j=1,i
            l=l+1
            aa(i,j)=a(l)
3       aa(j,i)=aa(i,j)
        call linrg(n,aa,5,s,5)
        l=0
        do 4 i=1,n
          do 4 j=1,i
            l=l+1
4       ai(l)=s(i,j)
        do 1 nm=1,(n-2)
          ii=n+1-nm
          l1=0
          do 2 i=ii,n
            l1=l1+1
            l2=0
            do 2 j=ii,n
              l2=l2+1
2         s(l1,l2)=aa(i,j)
          call linrg(nm,s,5,s,5)
          l=0
          do 1 i=1,l1
            do 1 j=1,i
              l=l+1
1       c(nm,l)=s(i,j)
        call lftrg(n,aa,5,aa,5,ipiv)
        call lfdrg(n,aa,5,ipiv,det1,det2)
        det=det1*10.0**det2
        ier=0
        return
        end
c
c
c
        subroutine matinv(a,nd,n,ifault)
        dimension a(nd,1)
c
        ifault=0
        call linrg(n,a,nd,a,nd)
        return
        end
