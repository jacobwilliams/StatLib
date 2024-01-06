      subroutine wext(x, n, ssq, a, n2, eps, w, pw, ifault)
c
c        Algorithm AS 181   Appl. Statist. (1982) Vol. 31, No. 2
c
c        Calculates Shapiro and Wilk's W statistic and its sig. level
c
c     Auxiliary routines required: ALNORM = algorithm AS 66 and NSCOR2
c     from AS 177.
c
      real x(n), a(n2), lamda, wa(3), wb(4), wc(4), wd(6), we(6), wf(7),
     *  c1(5, 3), c2(5, 3), c(5), unl(3), unh(3)
      integer nc1(3), nc2(3)
      logical upper
      data wa(1), wa(2), wa(3)
     *  /0.118898, 0.133414, 0.327907/,
     *     wb(1), wb(2), wb(3), wb(4)
     *  /-0.37542, -0.492145, -1.124332, -0.199422/,
     *     wc(1), wc(2), wc(3), wc(4)
     *  /-3.15805, 0.729399, 3.01855, 1.558776/,
     *     wd(1), wd(2), wd(3), wd(4), wd(5), wd(6)
     *  /0.480385, 0.318828, 0.0, -0.0241665, 0.00879701, 0.002989646/,
     *     we(1), we(2), we(3), we(4), we(5), we(6)
     *  /-1.91487, -1.37888, -0.04183209, 0.1066339, -0.03513666,
     *   -0.01504614/,
     *     wf(1), wf(2), wf(3), wf(4), wf(5), wf(6), wf(7)
     *  /-3.73538, -1.015807, -0.331885, 0.1773538, -0.01638782,
     *   -0.03215018, 0.003852646/
      data c1(1,1), c1(2,1), c1(3,1), c1(4,1), c1(5,1),
     *     c1(1,2), c1(2,2), c1(3,2), c1(4,2), c1(5,2),
     *     c1(1,3), c1(2,3), c1(3,3), c1(4,3), c1(5,3) /
     *     -1.26233, 1.87969, 0.0649583, -0.0475604, -0.0139682,
     *     -2.28135, 2.26186, 0.0, 0.0, -0.00865763,
     *     -3.30623, 2.76287, -0.83484, 1.20857, -0.507590/
      data c2(1,1), c2(2,1), c2(3,1), c2(4,1), c2(5,1),
     *     c2(1,2), c2(2,2), c2(3,2), c2(4,2), c2(5,2),
     *     c2(1,3), c2(2,3), c2(3,3), c2(4,3), c2(5,3) /
     *     -0.287696, 1.78953, -0.180114, 0.0, 0.0,
     *     -1.63638, 5.60924, -3.63738, 1.08439, 0.0,
     *     -5.991908, 21.04575, -24.58061, 13.78661, -2.835295/
      data unl(1), unl(2), unl(3) /-3.8,-3.0, -1.0/,
     *     unh(1), unh(2), unh(3) / 8.6, 5.8,  5.4/
      data nc1(1), nc1(2), nc1(3) /5, 5, 5/,
     *     nc2(1), nc2(2), nc2(3) /3, 4, 5/
      data pi6 /1.90985932/, stqr /1.04719755/, upper /.true./,
     *  zero/0.0/, tqr /0.75/, one /1.0/, onept4 /1.4/, three/3.0/,
     *  five/5.0/
c
      ifault = 1
      pw = one
      w = one
      if(n.le.2) return
      ifault = 3
      if(n/2 .ne. n2) return
      ifault = 2
      if(n.gt.2000) return
c
c        Calculate W
c
      ifault = 0
      w = zero
      an = n
      i = n
      do 10 j = 1,n2
        w = w+a(j) * (x(i)-x(j))
        i = i-1
   10 continue
      w = w * w / ssq
      if (w .lt. one) goto 20
      w = one
      return
c
c        Get significance level of W
c
   20 if (n .le. 6) goto 100
c
c        N between 7 and 2000 ... Transform W to Y, get mean and sd,
c        standardize and get significance level
c
      if(n.gt.20) goto 30
      al = log(an) - three
      lamda  =  poly(wa,3,al)
      ybar = exp(poly(wb,4,al))
      sdy = exp(poly(wc,4,al))
      goto 40
   30 al = log(an) - five
      lamda = poly(wd,6,al)
      ybar = exp(poly(we,6,al))
      sdy = exp(poly(wf,7,al))
   40 y  =  (one-w) ** lamda
      z  =  (y - ybar) / sdy
      pw  =  alnorm(z,upper)
      return
c
c        Deal with N less than 7 (Exact significance level for N = 3).
c
  100 if(w .le. eps) goto 160
      ww = w
      if(n.eq.3) goto 150
      un = log((w-eps)/(one-w))
      n3 = n-3
      if(un .lt. unl(n3)) goto 160
      if (un .ge. onept4) goto 120
      nc = nc1(n3)
      do 110 i = 1,nc
  110 c(i) = c1(i,n3)
      eu3 = exp(poly(c,nc,un))
      goto 140
  120 if (un .gt. unh(n3)) return
      nc = nc2(n3)
      do 130 i = 1,nc
  130 c(i) = c2(i,n3)
      un = log(un)
      eu3 = exp(exp(poly(c,nc,un)))
  140 ww = (eu3+tqr)/(one+eu3)
  150 pw = pi6 *(atan(sqrt(ww/(one-ww)))-stqr)
      return
  160 pw = zero
      return
      end
c
c
      subroutine wcoef(a, n, n2, eps, ifault)
c
c        Algorithm AS 181.1   Appl. Statist.  (1982) Vol. 31, No. 2
c
c        Obtain array A of weights for calculating W
c
      real a(n2), c4(2), c5(2), c6(3)
      data c4(1), c4(2) /0.6869, 0.1678/, c5(1), c5(2) /0.6647, 0.2412/,
     *  c6(1), c6(2), c6(3) /0.6431, 0.2806, 0.0875/
      data rsqrt2 /0.70710678/, zero/0.0/, half/0.5/, one/1.0/,
     *  two/2.0/, six/6.0/, seven/7.0/, eight/8.0/, thirt /13.0/
      ifault = 1
      if(n.le.2) return
      ifault = 3
      if(n/2 .ne. n2) return
      ifault = 2
      if(n.gt.2000) return
      ifault = 0
      if(n .le. 6) goto 30
c
c        N .GT. 6   Calculate rankits using approximate routine nscor2
c                   (AS177)
c
      call nscor2(a, n, n2, ifault)
      sastar = zero
      do 10 j = 2,n2
   10 sastar = sastar+a(j)*a(j)
      sastar = sastar*eight
      nn = n
      if (n .le. 20) nn = nn-1
      an = nn
      a1sq = exp(log(six*an+seven)-log(six*an+thirt)
     *  +half*(one+(an-two)*log(an+one)-(an-one)
     *  *log(an+two)))
      a1star = sastar/(one/a1sq-two)
      sastar = sqrt(sastar+two*a1star)
      a(1) = sqrt(a1star)/sastar
      do 20 j = 2,n2
   20 a(j) = two*a(j)/sastar
      goto 70
c
c        n.le.6   Use exact values for weiughts
c
   30 a(1) = rsqrt2
      if(n.eq.3) goto 70
      n3 = n-3
      goto (40,50,60), n3
   40 do 45 j = 1,2
   45 a(j) = c4(j)
      goto 70
   50 do 55 j = 1,2
   55 a(j) = c5(j)
      goto 70
   60 do 65 j = 1,3
   65 a(j) = c6(j)
c
c        Calculate the minimum possible value of W
c
   70 eps = a(1)*a(1)/(one-one/float(n))
      return
      end
c
c
      function poly(c, nord, x)
c
c
c        Algorithm AS 181.2   Appl. Statist.  (1982) Vol. 31, No. 2
c
c        Calculates the algebraic polynomial of order nored-1 with
c        array of coefficients c.  Zero order coefficient is c(1)
c
      real c(nord)
      poly = c(1)
      if(nord.eq.1) return
      p = x*c(nord)
      if(nord.eq.2) goto 20
      n2 = nord-2
      j = n2+1
      do 10 i = 1,n2
      p = (p+c(j))*x
      j = j-1
   10 continue
   20 poly = poly+p
      return
      end
c
c
c
      subroutine wgp(x, n, ssq, gp, h, a, n2, eps, w, u, p, ifault)
c
c        AS R63 Appl. Statist. (1986) Vol. 35, No.2
c
c        A remark on AS 181
c
c        Calculates Sheppard corrected version of W test.
c
c     Auxiliary functions required: ALNORM = algorithm AS 66, and
c     PPND = algorithm AS 111 (or PPND7 from AS 241).
c
      real x(n), a(n2)
      data twelve /1.2e01/, ten/1.0e1/, five/5.0e0/, one/1.0e0/,
     *     zero/0.0e0/
c
      zbar = zero
      zsd = one
      ifault = 1
      if (n.lt.7) return
c
c        No correction applied if gp=0.
c
      if (gp .le. zero) goto 10
      an1 = float(n-1)
c
c        correct ssq and find standardized grouping interval (h)
c
      ssq = ssq - an1 * gp * gp / twelve
      h = gp / sqrt(ssq / an1)
      ifault = 4
      if (h .gt. 1.5) return
   10 call wext(x, n, ssq, a, n2, eps, w, p, ifault)
      if (ifault .ne. 0) return
      if (p .gt. zero .and. p .lt. one) goto 20
      u = five - ten * p
      return
   20 if (gp .le. zero) goto 30
c
c        correct u for grouping interval (n .le. 100 and n .gt. 100
c        separately)
c
      hh = sqrt(h)
      if (n .gt. 100) goto 25
      zbar = -h*(1.07457 + hh*(-2.8185 + hh*1.8898))
      zsd = one + h*(0.50933 + hh*(-0.98305 + hh*0.7408))
      goto 30
   25 zbar = -h*(0.96436 + hh*(-2.1300 + hh*1.3196))
      zsd = one + h*(0.2579 + h*0.15225)
c
c        ppnd is AS 111 (Beasley and Springer, 1977)
c
   30 u = (-ppnd(p, ifppnd) - zbar) / zsd
c
c       alnorm is AS 66 (Hill, 1973)
c
      p = alnorm(u, .true.)
      return
      end
