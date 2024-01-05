      subroutine logccs(ns, nca, nct, nimax, nmax2, nvmax, nv, nv1, z,
     +		ivar, covi, w, ww, dl, iw, is, b, cov, chi2, st, ifault)
c
c     ALGORITHM AS 162  APPL. STATIST. (1981) VOL. 30, NO. 2
c
c     Logistic analysis of case-control studies
c
c     Auxiliary routine required: SYMINV from AS 7.
c     N.B. The description in the journal states that AS 6 is also
c          needed, but where ?
c
      integer	nca(ns), nct(ns), ivar(nv), iw(nmax2), is(ns)
      real	z(nvmax, nimax), covi(nv1), w(nv), ww(nv), dl(nv),
     +		b(nv), cov(nv1)
      real	one, zero, two
      logical	ifg, id
      data maxit /20/, eps /1.e-6/, one/1.0/, zero/0.0/, two/2.0/
c
c     Initial settings
c
      rlikp = one
      ifault = 0
      if (nvmax .lt. nv) go to 27
      do 1 i = 1, ns
	if (nmax2 .lt. nca(i) + nct(i) + 2) go to 27
    1 continue
      is(1) = 0
      IF(NS .EQ. 1) GOTO 28
      do 2 j = 2, ns
	j1 = j - 1
	is(j) = is(j1) + nca(j1) + nct(j1)
    2 continue
   28 if (nimax .lt. is(ns) + nca(ns) + nct(ns)) go to 27
      its = 0
c
c     Start of iteration
c
    3 its = its + 1
      if (its .gt. maxit) go to 24
      rlik = zero
      k = 0
      do 4 j = 1, nv
	dl(j) = zero
	do 4 jj = 1, j
	  k = k + 1
	  covi(k) = zero
    4 continue
c
c     Loop through strata
c
      do 17 i = 1, ns
	if (nca(i) * nct(i) .eq. 0) go to 17
	ifg = .false.
	sx = zero
	k = 0
	do 5 j = 1, nv
	  w(j) = zero
	  ww(j) = zero
	  do 5 jj = 1, j
	    k = k + 1
	    cov(k) = zero
    5	continue
	m = nca(i)
	n = m + nct(i)
	xx = one
	x = zero
	iw(1) = n + 1
	iw(n + 2) = -2
	kk = nct(i) + 1
	do 6 j = 2, kk
    6	iw(j) = 0
	do 7 j = 1, m
	  jj = kk + j
	  iw(jj) = j
    7	continue
c
c     Calculator numerator of terms of likelihood
c
	do 9 k = 1, nv
	  l = ivar(k)
	  bk = b(k)
	  wk = zero
	  do 8 j = 1, m
	    ji = is(i) + j
	    wk = wk + z(l, ji)
	    x = x + bk * z(l, ji)
    8	  continue
	  w(k) = wk
    9	continue
c
c     Go through all possible combinations to calculate denominator
c     terms of likelihood
c
   10   xx = exp(x)
	sx = sx + xx
	if (ifg) go to 12
	rlik = rlik + log(sx)
	do 11 k = 1, nv
   11   dl(k) = dl(k) + w(k)
	ifg = .true.
   12   l = 0
	do 13 k = 1, nv
	  ww(k) = ww(k) + xx * w(k)
	  do 13 kk = 1, k
	    l = l + 1
	    cov(l) = cov(l) + xx * w(k) * w(kk)
   13   continue
	call twidl(ips, im, iz, id, iw, nmax2)
	if (id) go to 15
c
c     Use the special features of TWIDL that only one element is altered
c     at a time, to calculate contribution of succeeding combinations
c     with minimal arithmetic.
c
        ips1i = is(i) + n - ips + 1
	im1i = is(i) + n - im + 1
	do 14 k = 1, nv
	  l = ivar(k)
	  zc = z(l, ips1i) - z(l, im1i)
	  w(k) = w(k) + zc
	  x = x + b(k) * zc
   14   continue
	go to 10
c
   15	rlik = rlik - log(sx)
	l = 0
        do 16 j = 1, nv
	  dl(j) = dl(j) - ww(j) / sx
	  do 16 k = 1, j
	    l = l + 1
	    covi(l) = covi(l) + (sx * cov(l) - ww(j) * ww(k)) / sx**2
   16   continue
   17 continue
c
      if (its .eq. 1) rliks = rlik
      call syminv(covi, nv, cov, w, nullty, ifault, nv1)
      if (ifault .ne. 0) go to 25
c
c     Calculate new parameter estimates.
c
      do 20 i = 1, nv
	w(i) = zero
	i2 = i * (i - 1) / 2
	do 18 j = 1, i
	  k = i2 + j
	  w(i) = w(i) + dl(j) * cov(k)
   18   continue
	i1 = i + 1
	if (i1 .gt. nv) go to 20
	do 19 k = i1, nv
	  j = k * (k - 1) / 2 + i
	  w(i) = w(i) + dl(k) * cov(j)
   19   continue
   20 continue
      do 21 i = 1, nv
   21 b(i) = b(i) + w(i)
      if (its .ne. 1) go to 23
c
c     Calculate score test
c
      st = zero
      do 22 i = 1, nv
   22 st = st + w(i) * dl(i)
c
c     Test for convergence
c
   23 rlik = rlik - rliks
      if (abs(rlikp - rlik) .le. eps) go to 26
      rlikp = rlik
      go to 3
   24 ifault = 1
      return
   25 ifault = 2
      return
   26 chi2 = two * rlik
      return
   27 ifault = 3
      return
      end
c
c
c
      subroutine twidl(x, y, z, done, p, n2)
      integer x, y, z, n2, p(n2)
      logical done
c
c        Algorithm AS 162.1  Appl. Statist. (1981) Vol.30, No. 2
c
c        This subroutine is a fortran version of 
c        CACM Algorithm 382 for generating all combinations
c        of M out of N objects.  All subscripts in the array P have
c        been increased by unity to avoid reference to P(0).
c
c        Ref:- P.J. Chase (1970), Comm. ACM, Vol.13, No.6, p368.
c
c        Parameters:-
c       
c        done         logical     input: Initially set to .false.
c                                 output: If all combinations have been
c                                         found then done is set to .true.
c                                         Otherwise done=.false.
c 
c        p            integer array  input: Initially, p(1)=N+1, p(N+2)=-2
c                     of length             p(2,...,N-M+1)=0,
c                     N+2                   p(N-M+2,...,N+1)=1,...,M
c                                           If M=0, set p(2)=1.
c
c       n2            integer       input: n2=n+2
c
c       x, y, z       integers      output: See below.
c
c       Initially, let A(1,...,N)=1,...,N
c                  let C(1,...,M)=A(N-M+1),...,A(N) = N-M+1,...,N
c                      be the initial combination.
c       Then the next combination is given by setting C(Z)=A(X)
c
c       Alternatively, initially let B(1,...,N-M)=0 and B(N-M+1,...,N)=1
c                                be the initial sequence of zeros and ones.
c       Then the next sequence of zeros and ones is obtained by setting
c       B(X)=1 and B(Y)=0.
c
				     
      j = 0
    1 j = j + 1
      if (p(j + 1) .le. 0) goto 1
      if (p(j) .ne. 0) goto 4
      if (j .lt. 3) goto 3
      do 2 i = 3, j
    2 p(i) = -1
    3 p(j + 1) = 0
      p(2) = 1
      x = 1
      z = 1
      y = j
      goto 10
    4 if (j .gt. 1) p(j) = 0
    5 j = j + 1
      j1 = j + 1
      if (p(j1) .gt. 0) goto 5
      i = j - 1
      k = i
    6 i = i + 1
      i1 = i + 1
      if (p(i1) .ne.0) goto 7
      p(i1) = -1
      goto 6
    7 if (p(i1) .ne. -1) goto 8
      z = p(k + 1)
      p(i1) = z
      x = i
      y = k
      p(k + 1) = -1
      goto 10
    8 if (i .ne. p(1)) goto 9
      done = .true.
      goto 10
    9 z = p(i1)
      p(j1) = z
      p(i1) = 0
      x = j
      y = i
   10 return
      end
c
c
c*******************************************************************
c
c  The routine below initializes TWIDL.   It is NOT part of AS 162.
c  It enables TWIDL to be used for other applications.
c
      subroutine initp(n, m, p, n2, done)
      integer p(n2)
      logical done
c
      done = .false.
      n1 = n + 1
      nm1 = n - m + 1
      p(1) = n1
      p(n+2) = -2
      do 1 i = 2, nm1
    1 p(i) = 0
      do 2 i = 1,m
      i1 = nm1 + i
    2 p(i1) = i
      return
      end
