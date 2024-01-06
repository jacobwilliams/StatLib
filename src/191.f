      subroutine sarmas(z, nz, n, beta, nbeta, ip, iq, ips, iqs, isea,
     +		iqap, maxit, a, s, sm, w, nw, ifault)
c
c     ALGORITHM AS 191  APPL. STATIST. (1983) VOL. 32, NO. 2
c
      double precision z(nz), a(nz), w(nw), beta(nbeta), s, sm
c
c     Local variables
c
      logical switch
      double precision zero, one, oneneg, etol, sprev, detm, detms,
     +		relerr, temp
c
c     Initialize numerical constants
c
      data zero/0.d0/, one/1.d0/, oneneg/-1.d0/
c
c     ETOL = error tolerance in convergence criterion
c
      data etol/1.d-10/
c
      iter = 0
      switch = .false.
      sprev = zero
      ipq = ip + iq
      ipqps = ipq + ips
      ipsts = ips + isea
      ipsts1 = ipsts + 1
      iqsts = iqs + isea
      iqsts1 = iqsts + 1
      ipsqs = ips + iqs
      iqap2 = iqap
      if (ip .eq. 0 .and. ips .eq. 0) iqap2 = min(iqap, iq + iqsts)
      maxit2 = maxit
      if (iq .eq. 0 .and. iqs .eq. 0) maxit2 = 0
      if (maxit2 .eq. 0) switch = .true.
      nby2 = n / 2
c
c     Input validation
c
      ifault = 0
      ir = max(iq + iqsts, ip + ipsts)
      if (ir .ge. n) ifault = 5
      if (maxit .gt. 0 .and. ir .gt. iqap) ifault = 7
      if (ipq + ipsqs .gt. nbeta) ifault = 4
      if (n + iqap2 .gt. nz) ifault = 4
      if (nw .lt. max(nz, 1 + ipq * (ipq + 3)/2, 
     +		1 + ipsqs * (ipsqs + 3)/2)) ifault = 4
      if (min(ip, iq, ips, iqs, isea, iqap, maxit) .lt. 0) ifault = 6
      if (ifault .ge. 1) return
c
c     Obtain necessary determinants.
c     Check for stationarity/invertibility
c
      detm = one
      detms = one
      ier = 0
      if (ipq .ne. 0) call dtarma(beta, ipq, ip, iq, w, nw, detm, ier)
      if (ier .gt. 0) go to 340
      if (ipsqs .eq. 0) go to 20
      ii = ipq
      do 10 i = 1, ipsqs
	ii = ii + 1
	a(i) = beta(ii)
   10 continue
      call dtarma(a, ipsqs, ips, iqs, w, nw, detms, ier)
      if (ier .gt. 0) go to 340
c
c     If IQAP2 = 0, use conditional sum of squares method
c
   20 if (iqap2 .eq. 0) go to 200
c
c     If no seasonal component and no moving-average component, proceed
c     directly to back-forecasting step (Y and E-series not needed)
c
      if (ipsqs .eq. 0 .and. iq .eq. 0) go to 110
c
c     Calculate Y-series, use W-vector
c
      do 60 i = 1, n
	w(i) = zero
	if (i .le. ip) go to 60
	w(i) = z(i)
	if (ip .eq. 0) go to 40
	do 30 j = 1, ip
	  iii = i - j
	  w(i) = w(i) - beta(j) * z(iii)
   30   continue
   40   l = min(iq, i - 1)
	if (l .eq. 0) go to 60
	do 50 j = 1, l
	  jj = ip + j
	  iii = i - j
	  w(i) = w(i) + beta(jj) * w(iii)
   50   continue
   60 continue
c
c     Calculate E-series, use A-vector
c
      lqs = iqs
      do 100 i = 1, n
	a(i) = zero
	if (i .le. ipsts) go to 100
	a(i) = w(i)
	if (ips .eq. 0) go to 80
	iii = i
	jj = ipq
	do 70 j = 1, ips
	  iii = iii - isea
	  jj = jj + 1
	  a(i) = a(i) - beta(jj) * w(iii)
   70   continue
   80   if (iqs .eq. 0) go to 100
	if (i .le. iqsts1) lqs = (i - 1) / isea
	if (lqs .eq. 0) go to 100
	iii = i
	do 90 j = 1, lqs
	  iii = iii - isea
          jj = ipqps + j
	  a(i) = a(i) + beta(jj) * a(iii)
   90   continue
  100 continue
c
c     Back-forecast Y-series, use w(n+1), w(n+2), ...
c
  110 do 150 i = 1, iqap2
	npi = n + i
	w(npi) = zero
	a(npi) = zero
	if (i .gt. iqsts) go to 130
	iii = npi
	do 120 j = 1, iqs
	  iii = iii - isea
	  jj = ipqps + j
	  w(npi) = w(npi) - beta(jj) * a(iii)
  120   continue
  130   if (ips .eq. 0) go to 150
	iii = npi
	do 140 j = 1, ips
	  iii = iii - isea
	  jj = ipq + j
	  w(npi) = w(npi) + beta(jj) * w(iii)
  140   continue
  150 continue
c
c     Back-forecast Z-series, use z(n+1), z(n+2), ...
c
      do 190 i = 1, iqap2
	npi = n + i
	z(npi) = w(npi)
	if (iq .eq. 0) go to 170
	do 160 j = 1, iq
	  npimj = npi - j
	  jj = ip + j
	  z(npi) = z(npi) - beta(jj) * w(npimj)
  160   continue
  170   if (ip .eq. 0) go to 190
	do 180 j = 1, ip
	  npij = npi - j
	  z(npi) = z(npi) + beta(j) * z(npij)
  180   continue
  190 continue
c
c     Calculate X-series, use W-vector
c
  200 npqap = n + iqap2
      ii = npqap + 1
      do 240 i = 1, npqap
	ii = ii - 1
	w(ii) = z(ii)
	im1 = i - 1
	l = min(im1, ip)
	if (l .eq. 0) go to 220
	iii = ii
	do 210 j = 1, l
	  iii = iii + 1
	  w(ii) = w(ii) - beta(j) * z(iii)
  210   continue
  220   l = min(im1, iq)
	if (l .eq. 0) go to 240
	iii = ii
	do 230 j = 1, l
	  iii = iii + 1
	  jj = ip + j
	  w(ii) = w(ii) + beta(jj) * w(iii)
  230   continue
  240 continue
c
c     Calculate A-series, use A-vector
c
      ii = npqap + 1
      do 280 i = 1, npqap
	ii = ii - 1
	a(ii) = w(ii)
	if (isea .eq. 0) go to 280
	if (i .le. ipsts1) lps = (i - 1) / isea
	if (lps .eq. 0) go to 260
	iii = ii
	do 250 j = 1, lps
	  iii = iii + isea
	  jj = ipq + j
	  a(ii) = a(ii) - beta(jj) * w(iii)
  250   continue
  260   if (i .le. iqsts1) lqs = (i - 1) / isea
	if (lqs .eq. 0) go to 280
	iii = ii
	do 270 j = 1, lqs
	  iii = iii + isea
	  jj = ipqps + j
	  a(ii) = a(ii) + beta(jj) * a(iii)
  270   continue
  280 continue
c
c     Calculate the sum of squares
c
      s = zero
      do 300 i = 1, npqap
  300 s = s + a(i) * a(i)
c
c     Test for convergence
c
      if (iqap2 .eq. 0) go to 330
      if (switch) go to 310
      ifault = 0
      relerr = (s - sprev) / s
      if (abs(relerr) .le. etol) go to 330
c
c     Convergence not obtained
c
  310 ifault = 1
      if (iter .ge. maxit2) go to 330
c
c     Reverse the series and proceed to the forecasting step.
c
      sprev = s
      ii = n
      do 320 i = 1, nby2
	temp = w(ii)
	w(ii) = w(i)
	w(i) = temp
	temp = a(ii)
	a(ii) = a(i)
	a(i) = temp
	temp = z(ii)
	z(ii) = z(i)
	z(i) = temp
	ii = ii - 1
  320 continue
      if (switch) iter = iter + 1
      switch = .not. switch
      go to 110
c
c     Modified sum of squares
c
  330 temp = oneneg / float(n)
      sm = s * detm**temp * detms**(float(isea) * temp)
      if (maxit2 .eq. 0) ifault = 0
      return
c
c     Model is nonstationary or noninvertible
c
  340 ifault = ier + 1
      return
      end
c
c
c
      subroutine dtarma(beta, nbeta, ip, iq, ws, nws, detm, ifault)
c
c     ALGORITHM AS 191.1  APPL. STATIST. (1983) VOL. 32, NO. 2
c
      double precision beta(nbeta), ws(nws), detm
c
c     Local variables
c
      double precision zero, one, oneneg, det, det1
      data zero/0.d0/, one/1.d0/, oneneg/-1.d0/
c
      ifault = 0
      if (nbeta .lt. ip + iq) go to 140
      nwchek = 1 + nbeta * (nbeta + 3)/2
      if (nwchek .gt. nws) go to 140
      det = one
      det1 = one
      ir = ip + iq
      irs = nws - ir - 1
      irsp1 = irs + 1
      ws(irsp1) = oneneg
      isw = 0
      if (ip .eq. 0) isw = 1
      iloop = ip
   10 if (isw .eq. 1) iloop = iq
      if (iloop .eq. 0) go to 120
      if (isw .eq. 2) go to 30
      do 20 i = 1, iloop
	irspi = irs + i + 1
	ippi = isw * ip + i
	ws(irspi) = beta(ippi)
   20 continue
      go to 60
   30 if (ip .eq. 0) go to 120
      iloop = ir
c
c     Multiply the autoregressive and moving-average operators to obtain
c     coefficients in the left-adjoint AR(ip+iq) model
c
      do 50 i = 1, ir
	ii = irs + 1 + i
	ws(ii) = zero
	imiq = i - iq
	j1 = max(0, imiq) + 1
	j2 = min(i, ip) + 1
	do 40 j = j1, j2
	  jm1 = j - 1
	  if (j .eq. 1) then
	    bj = oneneg
	  else
	    bj = beta(jm1)
	  end if
	  imj = i - j + 1
	  ippimj = ip + imj
	  if (imj .eq. 0) then
	    bi = oneneg
	  else
	    bi = beta(ippimj)
	  end if
	  ws(ii) = ws(ii) - bi * bj
   40   continue
   50 continue
c
c     Form the Schur matrix
c
   60 m = 0
      iend = iloop + 1
      do 90 i = 1, iloop
	do 80 j = 1, i
	  m = m + 1
	  ws(m) = zero
	  l = min(i, j)
	  do 70 k = 1, l
	    irsi = irs + i - k + 1
	    irsj = irs + j - k + 1
	    irspi = irs + iend - i + k
	    irspj = irs + iend - j + k
	    ws(m) = ws(m) + ws(irsi) * ws(irsj)
	    ws(m) = ws(m) - ws(irspi) * ws(irspj)
   70     continue
   80   continue
   90 continue
c
c     Calculate the determinant using the modified Cholesky
c     decomposition
c
      call mchol(ws, nws, iloop, det, ifault)
      if (ifault .gt. 0) go to 130
c
      if (isw .ge. 1) go to 110
      isw = 1
      det1 = det * det
      go to 10
c
  110 if (isw .eq. 2) go to 120
      isw = 2
      det1 = det1 * det * det
      go to 10
c
  120 detm = det1 / det
      return
c
  130 ifault = isw + 1
      return
  140 ifault = 3
      return
      end
c
c
c
      subroutine mchol(a, na, n, det, ifault)
c
c        Algorithm AS 191.2  Appl. Statist. (1983) Vol. 32, No. 2
c
c     This algorithm performs the modified Cholesky decomposition
c     A=LDL' where L is lower triangular and D is diagonal.  This
c     form of decomposition avoids the square-root computation of
c     the standard decomposition.
c
c     Parameters:-
c      
c     a     Double precision array (na)
c                               input:  Positive definite matrix stored in
c                                       symmetric mode as a11, a21, a22,...
c                                       amm.
c                               output: Modified Cholesky decomposition
c                                       stored as a one-dimensional array
c                                       as d11,l21,d22,l31,l32,d33,....,
c                                       lm1,lm2,...,lmm-1,dmm.
c     na     Integer            input:  m*(m+1)/2
c     n      Integer            input:  m, order of input matrix
c     det    Double precision   output: Determinant of a
c     ifault Integer            output: Fault Indicator
c                                       =1 if na or n is invalid
c                                       =2 if matrix a is not positive-definite
c                                       =0 otherwise
c
      double precision a(na), eta, one, det, w, t, tt
      data one /1.0d0/
c
c        eta  - largest number such that 1.0 + eta = 1.0
c        (depends upon machine precision)
c
       data eta /1.0d-14/
c
      ifault = 1
      det = one
      if (n .le. 0) goto 70
      if (na .lt. n*(n+1)/2) goto 70
      ifault = 2
      j = 1
      k = 0
      do 60 irow = 1,n
        l = 0
        do 20 icol = 1,irow
          k = k+1
          w = a(k)
          m = j
          if (irow .eq. icol) goto 30
          do 10 i = 1,icol
            l = l+1
            if (i .eq. icol) goto 20
            w = w-a(l)*a(m)
            m = m+1
   10     continue
   20     a(k) = w
c
   30   ii = 0
        do 40 i = 1,icol
          if(i .eq. icol) goto 50
          ii = ii+i
          t = a(m)
          tt = a(m) / a(ii)
          w = w - t * tt
          a(m) = tt
          m  =  m+1
   40   continue
   50   if (w .lt. eta*abs(a(k))) goto 70
        a(k) = w
        det = det*w
        j = j+irow
   60 continue
      ifault = 0
   70 return
      end
