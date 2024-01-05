      subroutine median (m, n, x, norm, c, f, med, ifault, 
     *                   maxm, maxn, md, g, q, l)
c
c        Algorithm AS 143  Appl. Statist. (1979) Vol. 28, No. 3
c
c        The Mediancentre
c        At least 10 times faster than AS 78.
c
c        Parameters:
c   
c        m       integer     input: number of dimensions
c        n       integer     input: number of sample points
c        x       real array  input: data, n points, m dimensions
c        norm    real        output: norm of gradient
c        c       real        output: multiplicity of m in S ; 0 if m not in S
c        med     real array  output: coordinates of the mediancentre
c        ifault  integer     output: -1 if m or n out of range
c                                     0 if matrix of second derivatives
c                                       is positive definite
c                                     1 if matrix of second derivatives
c                                       is not positive definite
c
      implicit double precision (a-h, o-z)
      dimension x(maxn,maxm), med(maxm), md(maxm), g(maxm),
     *          q(maxm,maxm), l(maxn)
      double precision norm, med, md
c
      data zero, one /0.0d0, 1.0d0/
      data eps /0.0001d0/
c
      ifault = -1
      if (m .lt. 1 .or. m .gt. maxm .or. n .lt. 1 .or. n .gt. maxn)
     *    return
      ifault = 0
      do 1 i = 1, m
    1 med(i) = zero
      nn = n
      do 2 i = 1, n
    2 l(i) = i
c
c        computation of (residual) mean
c
    3 f = one / float(nn)
      do 5 i = 1, m
        s = zero
        g(i) = zero
        do 4 j = 1, nn
          k = l(j)
          s = s + x(k, i)
    4   continue
        md(i) = s * f
        med(i) = med(i) + md(i)
    5 continue
c
c        Computation of function, gradient, and norm of gradient
c
      c = zero
      f = zero
      s = zero
      do 9 i = 1, n
        t = zero
        do 6 j = 1, n
          x(i, j) = x(i, j) - md(j)
          t = t + x(i, j) * x(i, j)
    6   continue
        if (t .gt. zero) goto 7
        c = c + one
        goto 9
    7   t = sqrt(t)
        r = one / t
        s = s + r
        f = f + t
        do 8 j = 1, m
          g(j) = g(j) - r * x(i, j)
    8   continue
    9 continue
      norm = zero
      do 10 i = 1, m
        norm = norm + g(i) * g(i)
   10 continue
      norm = sqrt(norm)
c
c        Check for extremal points
c
      if (norm .le. c + eps) return
      if (nn .eq. 1) goto 13
c
c        Simplex
c
      nt = nn
      nn = 0
      do 12 i = 1, nt
        k = l(i)
        ang = zero
        do 11 j = 1, m
          ang = ang + g(j) * x(k, j)
   11   continue
        if (ang .gt. zero) goto 12
        nn = nn + 1
        l(nn) = k
   12 continue
      if (nn .gt. 0) goto 3
c
c        Starting value
c
   13 r = (c / norm - one) / s
      do 14 i = 1, m
        md(i) = r * g(i)
   14 continue
c
c        Newton - Raphson procedure
c
   15 c = zero
      f = zero
      s = zero
      do 16 i = 1, m
        g(i) = zero
        med(i) = med(i) + md(i)
        do 16 j = i, m
          q(i, j) = zero
   16 continue
      do 19 i = 1, n
        t = zero
        do 17 j = 1, m
          x(i, j) = x(i, j) - md(j)
          t = t + x(i, j) * x(i, j)
   17   continue
        if (t .eq. zero) goto 19
        t = sqrt(t)
        r = one / t
        s = s + r
        rr = r * r * r
        f = f + t
        do 18 j = 1, m
          g(j) = g(j) - r * x(i, j)
          do 18 k = j, m
            q(j, k) = q(k, j) - rr * x(i, j) * x(i, k)
   18   continue
   19 continue
      norm = zero
      do 20 j = 1, m
        norm = norm + g(j) * g(j)
   20 continue
      norm = sqrt(norm)
      if (norm .le. eps) return
c
c        Cholesky and solution of equation
c
      ifault = 1
      if (q(1,1) + s .le. zero) return
      q(1,1) = sqrt(q(1,1) + s) 
      do 21 i = 2, m
        q(i, 1) = q(1, i) / q(1, 1)
   21 continue
      do 25 j = 2, m
        tmp = s + q(j, j)
        j1 = j - 1
        do 22 i = 1, j1
          tmp = tmp - q(j, i) * q(j, i)
   22   continue
        if (tmp .le. zero) return
        q(j, j) = sqrt(tmp)
        j2 = j + 1
        if (j2 .gt. m) goto 25
        do 24 k = j2, m
          tmp = q(j, k)
          do 23 i = 1, j1
            tmp = tmp - q(j, i) * q(k, i)
   23     continue
          q(k, j) = tmp / q(j, j)
   24   continue
   25 continue
      md(1) = -g(1) / q(1, 1)
      do 27 j = 2, m
	j1 = j - 1
        tmp = -g(j)
        do 26 i = 1, j1
          tmp = tmp - q(j, i) * md(i)
   26   continue
        md(j) = tmp / q(j, j)
   27 continue
      md(m) = md(m) / q(m, m)
      j1 = m - 1
      do 29 j2 = 1, j1
        i = m - j2
        tmp = md(i)
        k = i + 1
        do 28 j = k, m
          tmp = tmp - q(j, i) * md(j)
   28   continue
        md(i) = tmp / q(i, i)
   29 continue
      ifault = 0
      goto 15
      end
