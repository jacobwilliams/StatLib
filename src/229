      subroutine rq(m,n,m5,n2,a,b,t,toler,ift,x,e,s,wa,wb,
     1  nsol,sol,lsol)
      implicit double precision (a-h,o-z)
c     
c     Algorithm AS229 Appl. Statist. (1987) Vol. 36, No. 3
c     
      integer i,j,k,kl,kount,kr,l,lsol,m,m1,m2,m3,m4,m5,
     1  ift
      integer n,n1,n2,nsol,out,s(m)
      logical stage,test,init,iend
      double precision a1,aux,b1,big,d,dif,pivot,smax,t,t0,t1,tnt
      double precision min,max,toler,zero,half,one,two,three,four,five
      double precision b(m),sol(n2,nsol),a(m,n),x(n),wa(m5,n2),wb(m)
      double precision sum,e(m)
      double precision xx,y
      data big /1.0d37/
      data zero /0.0d0/
      data half /0.5d0/
      data one /1.0d0/
      data two /2.0d0/
      data three /3.0d0/
      data four /4.0d0/
      data five /5.0d0/
c
c     check dimension parameters
c
      ift = 0
      if (m5.ne.m + 5) ift = 3
      if (n2.ne.n + 2) ift = 4
      if (m.le.zero.or.n.le.zero) ift = 5
      if (ift.gt.two) return
c
c     initialization
c 
      m1 = m + 1
      n1 = n + 1
      m2 = m + 2
      m3 = m + 3
      m4 = m + 4
      do 2 i = 1,m
        wb(i) = b(i)
        do 1 j = 1,n
          wa(i,j) = a(i,j)
   1    continue
   2  continue
      wa(m2,n1) = zero
      dif = zero
      iend = .true.
      if (t.ge.zero.and.t.le.one) goto 3
      t0 = one / float(m) - toler
      t1 = one - t0
      t = t0
      iend = .false.
   3  continue
      init = .false.
      lsol = 1
      kount = 0
      do 9 k = 1,n
        wa(m5,k) = zero
        do 8 i = 1,m
          wa(m5,k) = wa(m5,k) + wa(i,k)
   8    continue
        wa(m5,k) = wa(m5,k) / float(m)
   9  continue
      do 10 j = 1,n
        wa(m4,j) = j
        x(j) = zero
  10  continue
      do 40 i = 1,m
        wa(i,n2) = n + i
        wa(i,n1) = wb(i)
        if (wb(i).ge.zero) goto 30
        do 20 j = 1,n2
          wa(i,j) = -wa(i,j)
  20    continue
  30    e(i) = zero
  40  continue
      do 42 j = 1,n
        wa(m2,j) = zero
        wa(m3,j) = zero
        do 41 i = 1,m
          aux = sign(one,wa(m4,j)) * wa(i,j)
          wa(m2,j) = wa(m2,j) + aux * (one - sign(one,wa(i,n2)))
          wa(m3,j) = wa(m3,j) + aux * sign(one,wa(i,n2))
  41    continue
        wa(m3,j) = two * wa(m3,j)
  42  continue
      goto 48
  43  continue
      lsol = lsol + 1
      do 44 i = 1,m
        s(i) = zero
  44  continue
      do 45 j = 1,n
        x(j) = zero
  45  continue
c
c     compute next t
c
      smax = two
      do 47 j = 1,n
        b1 = wa(m3,j)
        a1 = (-two - wa(m2,j)) / b1
        b1 = -wa(m2,j) / b1
        if (a1.lt.t) goto 46
        if (a1.ge.smax) goto 46
        smax = a1
        dif = (b1 - a1) / two
  46    if (b1.le.t) goto 47
        if (b1.ge.smax) goto 47
        smax = b1
        dif = (b1 - a1) / two
  47  continue
      tnt = smax + toler * (one + dabs(dif))
      if (tnt.ge.t1 + toler) iend = .true.
      t = tnt
      if (iend) t = t1
  48  continue
c
c     compute new marginal costs
c
      do 49 j = 1,n
        wa(m1,j) = wa(m2,j) + wa(m3,j) * t
  49  continue
      if (init) goto 265
c
c     stage 1
c
c      determine the vector to enter the basis
      stage = .true.
      kr = 1
      kl = 1
  70  max = -one
      do 80 j = kr,n
        if (abs(wa(m4,j)).gt.n) goto 80
        d = abs(wa(m1,j))
        if (d.le.max) goto 80
        max = d
        in = j
  80  continue
      if (wa(m1,in).ge.zero) goto 100
      do 90 i = 1,m4
        wa(i,in) = -wa(i,in)
  90  continue
c
c      determine the vector to leave the basis
c
 100  k = 0
      do 110 i = kl,m
        d = wa(i,in)
        if (d.le.toler) goto 110
        k = k + 1
        wb(k) = wa(i,n1) / d
        s(k) = i
        test = .true.
 110  continue
 120  if (k.gt.0) goto 130
      test = .false.
      goto 150
 130  min = big
      do 140 i = 1,k
        if (wb(i).ge.min) goto 140
        j = i
        min = wb(i)
        out = s(i)
 140  continue
      wb(j) = wb(k)
      s(j) = s(k)
      k = k - 1
c
c     check for linear dependence in stage 1
c
 150  if (test.or. .not.stage) goto 170
      do 160 i = 1,m4
        d = wa(i,kr)
        wa(i,kr) = wa(i,in)
        wa(i,in) = d
 160  continue
      kr = kr + 1
      goto 260
 170  if (test) goto 180
      wa(m2,n1) = two
      goto 390
 180  pivot = wa(out,in)
      if (wa(m1,in) - pivot - pivot .le. toler) goto 200
      do 190 j = kr,n1
        d = wa(out,j)
        wa(m1,j) = wa(m1,j) - d - d
        wa(m2,j) = wa(m2,j) - d - d
        wa(out,j) = -d
 190  continue
      wa(out,n2) = -wa(out,n2)
      goto 120
c
c      pivot on wa(out,in)
c
 200  do 210 j = kr,n1
        if (j.eq.in) goto 210
        wa(out,j) = wa(out,j) / pivot
 210  continue
      do 230 i = 1,m3
        if (i.eq.out) goto 230
        d = wa(i,in)
        do 220 j = kr,n1
          if (j.eq.in) goto 220
          wa(i,j) = wa(i,j) - d * wa(out,j)
 220    continue
 230  continue
      do 240 i = 1,m3
        if (i.eq.out) goto 240
        wa(i,in) = -wa(i,in) / pivot
 240  continue
      wa(out,in) = one / pivot
      d = wa(out,n2)
      wa(out,n2) = wa(m4,in)
      wa(m4,in) = d
      kount = kount + 1
      if (.not.stage) goto 270
c
c       interchange rows in stage 1
c
      kl = kl + 1
      do 250 j = kr,n2
        d = wa(out,j)
        wa(out,j) = wa(kount,j)
        wa(kount,j) = d
 250  continue
 260  if (kount + kr.ne.n1) goto 70
c
c       stage 2
c
 265  stage = .false.
c
c      determine the vector to enter the basis
c
 270  max = -big
      do 290 j = kr,n
        d = wa(m1,j)
        if (d.ge.zero) goto 280
        if (d.gt. (-two)) goto 290
        d = -d - two
 280  if (d.le.max) goto 290
      max = d
      in = j
 290  continue
      if (max.le.toler) goto 310
      if (wa(m1,in) .gt.zero) goto 100
      do 300 i = 1,m4
        wa(i,in) = -wa(i,in)
 300  continue
      wa(m1,in) = wa(m1,in) - two
      wa(m2,in) = wa(m2,in) - two
      goto 100
c
c       compute quantiles
 310  continue
      do 320 i = 1,kl - 1
        k = wa(i,n2) * sign(one,wa(i,n2))
        x(k) = wa(i,n1) * sign(one,wa(i,n2))
 320  continue
      sum = zero
      do 330 i = 1,n
        sum = sum + x(i) * wa(m5,i)
        sol(i + 2,lsol) = x(i)
 330  continue
      sol(1,lsol) = t
      sol(2,lsol) = sum
      if (iend) goto 340
      init = .true.
      goto 43
 340  continue
      if (lsol.le.2) goto 355
      do 350 i = 2,lsol
        sol(1,i-1) = sol(1,i)
 350  continue
      lsol = lsol - 1
      sol(1,lsol) = one
 355  continue
      l = kl - 1
      do 370 i = 1,l
        if (wa(i,n1).ge.zero) goto 370
        do 360 j = kr,n2
          wa(i,j) = -wa(i,j)
 360    continue
 370  continue
      wa(m2,n1) = zero
      if (kr.ne.1) goto 390
      do 380 j = 1,n
        d = abs(wa(m1,j))
        if (d.le.toler .or. two - d .le. toler) goto 390
 380  continue
      wa(m2,n1) = one
 390  do 400 i = kl,m
        k = wa(i,n2) * sign(one,wa(i,n2))
        d = wa(i,n1) * sign(one,wa(i,n2))
        k = k - n
        e(k) = d
 400  continue
      wa(m2,n2) = kount
      wa(m1,n2) = n1 - kr
      sum = zero
      do 410 i = 1,m
        sum = sum + e(i) * (half + sign(one,e(i)) * (t - half))
 410  continue
      wa(m1,n1) = sum
      return
      end
