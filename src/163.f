	subroutine gtrans(r, nr, m, n, ifault)
c
c	Algorithm AS163 Applied Statistics (1981) vol.30, no. 2
c
c	Transposes rows m and m+1 in the standardized lower-triangular
c	working matrix r and restores the triangular structure by plane
c	rotations.   r is of length nr = n(n+1)/2 and is stored row by
c	row, shortest row first.
c
	double precision r(nr), a, b, c, d, zero, one
c
	data zero, one/0.d0, 1.d0/
c
	ifault = 1
	if (n .lt. 2 .or. m .lt. 1 .or. m .ge. n) return
	ifault = 0
	ma = m * (m + 1) / 2
	mb = ma + m
c
c	Exchange columns
c
	i = mb
   20	i = i - 1
	if (i .eq. ma) go to 30
	j = i - m
	a = r(i)
	r(i) = r(j)
	r(j) = a
	go to 20
   30	mc = mb + 1
c
c	a will now be set equal to element (m,m)
c	b = element (m,m+1), and c = element (m+1,m+1).
c	If a and c are zero, there is nothing to do.
c
	a = r(ma)
	c = r(mc)
	if (a .le. zero .and. c .le. zero) return
c
c	Other special cases where one or more of a, b or c is zero.
c
	b = r(mb)
	j = m + 2
	mi = mb + j
	if (b .ne. zero) go to 50
	r(mc) = a
	r(ma) = c
	if (j .gt. n) return
	do 40 i = j, n
	  b = r(mi-1)
	  r(mi-1) = r(mi)
	  r(mi) = b
	  mi = mi + i
   40	continue
	return

   50	if (c .gt. zero) go to 70
	r(ma) = a * b * b
	r(mb) = one / b
	if (j .gt. n) return
	do 60 i = j, n
	  r(mi-1) = r(mi-1) / b
	  mi = mi + i
   60	continue
	return
c
c	General case - a, b, c all non-zero
c
   70	d = a * b * b + c
	r(ma) = d
	c = c / d
	r(mc) = a * c
	a = a * b / d
	r(mb) = a
	if (j .gt. n) return
	do 80 i = j, n
	  d = r(mi-1)
	  r(mi-1) = a * d + c * r(mi)
	  r(mi) = d - b * r(mi)
	  mi = mi + i
   80	continue
	return
	end
c
c
	subroutine update(r, x, w, n, nr)
c
c	Algorithm AS163.1 Applied Statistics (1981) vol.30, no.2
c
c	Uses plane rotations to update the standardized square root
c	(Cholesky factor) of a matrix of sums of squares and products
c	when another data vector x is to be included.   r is of length
c	nr = n(n+1)/2 and is stored row by row, shortest row first.
c
	double precision x(n), r(nr), w, y, d, c, s, eps
c
	data eps/0.d0/
c
	l = 0
	do 20 i = 1, n
	  if (w .le. eps) return
	  l = l + i
	  y = x(i)
	  if (abs(y) .le. eps) go to 20
	  d = r(l) + w * y * y
	  if (abs(d) .le. eps) go to 20
	  c = r(l) / d
	  r(l) = d
	  if (i .eq. n) return
	  s = w * y / d
	  w = c * w
	  k = l
	  m = i + 1
	  do 10 j = m, n
	    k = k + j - 1
	    d = x(j)
	    r(k) = c * r(k) + s * d
   10	  continue
   20	continue
	return
	end
