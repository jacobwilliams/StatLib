       subroutine monte(est, conf, nsims, simval, ordval, lsb, sbound,
     # lbb, bbound, alow, ahi, blow, bhi, clow, chi, dlow, dhi, 
     # ifault)
       implicit double precision (a-h,o-z)
c
c         Algorithm AS214 Appl. Statistics. (1985) vol. 34, no. 3
c
c         Sets up monte carlo confidence intervals
c         uses function ppnd - algorithm AS 111 (see also AS241)
c
      logical lsb, lbb
      dimension simval(nsims), ordval(nsims)
      data least /5/
      data zero, half, one, two, hun, big /0.0d0,
     #  0.5d0, 1.0d0, 2.0d0, 1.0d2, 1.0d20 /
c
      ifault = 0
      sbnd = sbound
      bbnd = bbound
      if (conf .ge. hun) ifault = 6
      if (conf .le. zero) ifault = 7
      if (ifault .gt. 0) return
c
c        Symmetric mci
c
c        Find mean and variance
c
      v1 = zero
      v2 = zero
      do 10 j = 1, nsims
        fj = float(j)
        v2 = v2 + (fj - one) * (simval(j) - v1) ** 2 / fj
        v1 = (simval(j) + (fj - one) * v1) / fj
 10   continue
      stderr = zero
      if (v2 .gt. zero) stderr = sqrt(v2 / float(nsims - 1))
      alpha = half * (hun - conf) / hun
      z = -ppnd(alpha, ifault)
      if (ifault .eq. 0) goto 20
      ifault = 6
      return
 20   ahi = z * stderr
      alow = est -ahi
      ahi = est + ahi
c
c        Calculate bias adjustment, so that the number of values that
c        must be ordered is known
c
      call biasad(est, nsims, simval, z, limit1, limit2, ifault)
      limitl = alpha * float(nsims + 1) + half
      limitu = (one-alpha) * float(nsims + 1) + half
      if (limitl .lt. least) ifault = 5
      if (limitl .eq. limitu) ifault = 7
      if (ifault .gt. 0) return
      l1 = max0(limit1, 2 * limitl)
      l2 = min0(limit2, 2 * limitu - nsims - 1)
      if (lsb) sbnd = -big
      if (lbb) bbn = big
c
c        Select and order the l1 smallest and the nsims+1-l2 largest
c        values
c
      call order(simval,ordval, nsims, l1, l2)
      if (sbnd .le. ordval(1)) goto 40
      if (lsb) goto 30
      ifault = 1
      return
 30   sbnd = ordval(1) * two
      if (sbnd .gt. zero) sbnd = -sbnd
 40   if (bbnd .ge. ordval(nsims)) goto 60
      if (lbb) goto 50
      ifault = 2
      return
 50   bbnd = ordval(nsims) * two
      if (bbnd .lt. zero) bbnd = -bbnd
c
c        Equal tails mci ( percentile method )
c
 60   blow = ordval(limitl)
      bhi = ordval(limitu)
c
c        Bias-corrected percentile method
c
      clow = sbnd
      if (limit1 .gt. 0) clow = ordval(limit1)
      chi = bbnd
      if (limit2 .le. nsims) chi = ordval(limit2)
c
c        Minimum lenght mci
c
      call lmin(ordval, nsims, sbnd, bbnd, limitl, limitu, dlow, dhi)
      return
      end
c
c
c
c
c
      subroutine biasad(est, nsims, simval, z, limit1, limit2, ifault)
      implicit double precision (a-h,o-z)
c
c        Finds adjustment required to implement bias corrected 
c        percentile method
c
c        Uses functions alnorm and ppnd -
c        algorithms  AS 66 and AS 111
c
      dimension simval(nsims)
      data half, two/0.5d0, 2.0d0/
c
      j = 0
      k = 0
c
c        m = (number of values less than est)+(half the number equal to est),
c        rounded up if not integral
c
      do 30  m = 1, nsims
        if (simval(m) - est) 10, 20, 30
 10     j = j + 1
 20     k = k + 1
 30   continue
      m = (j + k + 1) / 2
      if (m .eq. 0) ifault = 3
      if (m .eq. nsims) ifault = 4
      if (ifault .gt. 0) return
      zed = two * ppnd(float(m) / float(nsims), ifail)
c
c ifail cannot exceed 0 since m .ge. 1 and m .le. (nsims-1)
c
      fn1 = float(nsims + 1)
      limit1 = fn1 * alnorm(zed - z, .false.) + half
      limit2 = fn1 * alnorm(zed + z, .false.) + half
      return
      end
c
c
c
c
c
c
      subroutine order(simval, ordval, nsims, l1, l2)
      implicit double precision (a-h, o-z)
c
c        Orders the smallest and largest values from simval and
c        puts them in ordval
c
      dimension simval(nsims), ordval(nsims)
c
      ordval(1) = simval(1)
      ordval(nsims) = ordval(1)
      lm1 = 1
      lm2 = nsims
      do 100 j = 2, nsims
        if (simval(j) .lt. ordval(lm1)) goto 10
        if (lm1 .eq. l1) goto 50
        lm1 = lm1 +1
        ordval(lm1) = simval(j)
        goto 50
 10     if (lm1 .lt. l1) lm1 =lm1 +1
        if (lm1 .eq. 2) goto 30
        ll =lm1 -2
        do 20 k = 1, ll
          kk = lm1 - k
          ordval(kk + 1) = ordval(kk)
          if (ordval(kk-1) .le. simval(j)) goto 40
 20     continue
 30     ordval(2) = ordval(1)
        kk=1
 40     ordval(kk) = simval(j)
 50     if (simval(j) .gt. ordval(lm2)) goto 60
        if (lm2 .eq. l2) goto 100
        lm2 = lm2 - 1
        ordval(lm2) = simval(j)
        goto 100
 60     if (lm2 .gt. l2) lm2 = lm2 - 1
        if (lm2 .eq. nsims - 1) goto 80
        ll = nsims - 2
        do 70 k = lm2, ll
          kk = k + 1
          ordval(k) = ordval(kk)
          if (ordval(kk + 1) .ge. simval(j)) goto 90
 70     continue
 80     kk = nsims
        ordval(kk -1) = ordval(kk)
 90     ordval(kk) = simval(j)
 100  continue
      return
      end
c
c
c
      subroutine lmin(ordval, nsims, sbnd, bbnd, limitl, limitu,
     #dlow,dhi)
      implicit double precision (a-h, o-z)
c
c        Finds the interval with minimum length
c
      dimension ordval(nsims)
      length = limitu - limitl
      numb = nsims - length
      diff = ordval(nsims) - ordval(1)
      l = 0
      do 10 k = 1, numb
        kk = k + length
        if (diff .le. ordval(kk) - ordval(k)) goto 10
        l = k
        diff = ordval(kk) - ordval(k)
 10   continue
      dlow = ordval(l)
      l = l + length
      dhi = ordval(l)
      l = numb + 1
      if (dhi - dlow .le. bbnd - ordval(l)) goto 20
      dlow = ordval(l)
      dhi = bbnd
  20  l = length
      if (dhi - dlow .le. ordval(l) -sbnd) return
      dlow = sbnd
      dhi = ordval (l)
      return
      end
