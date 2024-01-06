      subroutine spect(a, m, bigt, mu, nusq, lognu, u, iu, v, iv, z,
     +		ifault)
c
c     Algorithm AS 193  Appl. Statist. (1983) vol. 32, no. 3
c
c     A revised algorithm for the spectral test
c
      integer bigt, iu, iv, t, t1, ifault
      double precision a, m, mu(bigt), nusq(bigt), lognu(bigt),
     +		u(iu, bigt), v(iv, bigt), z(bigt)
      double precision h, hprime, mmax, mmax2, msq, p, pi, pprime, q,
     +		vprod
      data zero/0.d0/, one/1.d0/, two/2.d0/, four/4.d0/
c
c     Suitable values for
c     1) IBM REAL*8
c     data mmax/33554432.d0/
c     2) IBM REAL*16
c
      data mmax/9007199254740992.d0/
c
c     3) CDC 7600 Double precision
c     data mmax/35184372088832.d0/
c
c     Test the validity of the input parameters
c
      mmax2 = mmax + mmax
      ifault = 0
      if (bigt .lt. 2) ifault = 1
      if (a .ge. m .or. a .le. zero .or. m .le. zero) ifault = 2
      if (m .gt. mmax) ifault = 3
      if (ifault .gt. 0) return
c
c     Check that A and M are relatively prime.
c     Needs valid A and M.   Uses Euclid's algorithm.
c
      h = a
      hprime = m
   10 r = mod(hprime, h)
      if (r .eq. zero) go to 20
      hprime = h
      h = r
      go to 10
   20 if (h .ne. one) then
	ifault = 4
	return
      end if
      msq = m * m
c
c     All steps refer to those in Knuth's algorithm.
c     Step 1 - initialization.
c
      h = a
      hprime = m
      p = one
      pprime = zero
      r = a
      s = one + a * a
c
c     Step 2 - Euclidean step
c
   30 q = int(hprime / h)
      uc = hprime - q * h
      vc = pprime - q * p
      w = uc * uc + vc * vc
      if (w .ge. s) go to 40
      s = w
      hprime = h
      h = uc
      pprime = p
      p = vc
      go to 30
c
c     Step 3 - compute nu(2)
c
   40 uc = uc - h
      vc = vc - p
      w = uc * uc + vc * vc
      if (w .ge. s) go to 50
      s = w
      hprime = uc
      pprime = vc
   50 nusq(2) = s
c
c     Initialize U and V matrices.
c     N.B. We store by columns whereas Knuth stores by rows.
c
      t = 2
      u(1,1) = -h
      u(1,2) = -hprime
      u(2,1) = p
      u(2,2) = pprime
      sign = one
      if (pprime .gt. zero) sign = -one
      v(1,1) = sign * pprime
      v(1,2) = -sign * p
      v(2,1) = sign * hprime
      v(2,2) = -sign * h
c
c     Step 4 - advance T
c
   60 if (t .eq. bigt) go to 200
      t1 = t
      t = t + 1
      r = mod(a*r, m)
      u(1,t) = -r
      u(t,t) = one
      u(t,1) = zero
      v(1,t) = zero
      v(t,t) = m
      do 70 i = 2, t1
	u(i,t) = zero
	u(t,i) = zero
	v(i,t) = zero
   70 continue
      do 90 i = 1, t1
	qtemp = v(1,i) * r
	q = nint(qtemp / m)
	v(t,i) = qtemp - q * m
	do 80 i2 = 1, t
   80   u(i2, t) = u(i2, t) + q * u(i2, i)
   90 continue
      s = min(s, vprod(u(1,t), u(1,t), t))
      k = t
      j = 1
c
c     Step 5 - transform
c
  100 do 120 i = 1, t
	if (i .eq. j) go to 120
	vij = vprod(v(1,i), v(1,j), t)
	vjj = vprod(v(1,j), v(1,j), t)
	if (two * abs(vij) .le. vjj) go to 120
	q = nint(vij / vjj)
	do 110 i2 = 1, t
	  v(i2, i) = v(i2, i) - q * v(i2, j)
	  u(i2, j) = u(i2, j) + q * u(i2, i)
  110   continue
	k = j
  120 continue
c
c     Step 6 - examine new bound
c
      if (k .eq. j) s = min(s, vprod(u(1,j), u(1,j), t))
c
c     Step 7 - advance J
c
      j = j + 1
      if (j .eq. t + 1) j = 1
      if (j .ne. k) go to 100
c
c     Step 8 - prepare for search
c
c     MU and LOGNU are used to store Knuth's X and Y respectively
c
      do 130 i = 1, t
	mu(i) = zero
	lognu(i) = zero
	qtemp = vprod(v(1,i), v(1,i), t)
	if (qtemp .gt. mmax2) go to 240
	qtemp = qtemp / msq
	z(i) = int(sqrtdble((int(qtemp * s))))
  130 continue
      k = t
c
c     Step 9 - advance XK
c
  140 if (mu(k) .eq. z(k)) go to 190
      mu(k) = mu(k) + one
      do 150 i = 1, t
  150 lognu(i) = lognu(i) + u(i,k)
c
c     Step 10 - advance K
c
  160 k = k + 1
      if (k .gt. t) go to 180
      mu(k) = -z(k)
      do 170 i = 1, t
  170 lognu(i) = lognu(i) - two * z(k) * u(i,k)
      go to 160
  180 s = min(s, vprod(lognu, lognu, t))
c
c     Step 11 - decrease K
c
  190 k = k - 1
      if (k .ge. 1) go to 140
      nusq(t) = s
      go to 60
c
c     Calculate NU and log(NU)
c
  200 do 210 i = 2, bigt
	mu(i) = sqrt(nusq(i))
	lognu(i) = log(mu(i)) / log(two)
  210 continue
c
c     Calculate transformed MU values
c
      pi = four * atan(one)
      q = one
      do 220 t = 2, bigt, 2
	q = q * pi * two / t
	mu(t) = q * mu(t) ** t / m
  220 continue
      if (bigt .eq. 2) return
      q = two
      do 230 t = 3, bigt, 2
	q = q * pi * two / t
	mu(t) = q * mu(t) ** t / m
  230 continue
      return
c
  240 ifault = 5
      return
      end
c
c
      double precision function vprod(u, v, t)
c
c     Algorithm AS 193.1  Appl. Statist. (1983) vol.32, no.3
c
c     Auxiliary function to calculate the inner product of vectors
c     U and V of length T.
c     N.B. Equivalent to DDOT from LINPACK.
c
      integer t
      double precision u(t), v(t), zero
      data zero/0.d0/
c
      vprod = zero
      do 10 i = 1, t
   10 vprod = vprod + u(i) * v(i)
      return
      end
