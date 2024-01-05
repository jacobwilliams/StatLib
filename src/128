      subroutine covmat(v,n,mdim,v11,ex1,ex2,summ2,ifault)
c
c      Appl. Statist. algorithm as 128 (1978), vol. 27
c      Davis C.S. and Stephens M.A.
c
c      Computes and normalises the David-Johnson approximation
c      for the covariance matrix of normal order statistics.
c
c      Auxiliary function equired: PPND = algorithm AS111
c      N.B. This file also includes ASR72 to calculate v11
c
      integer n,mdim,ifault
      double precision v(mdim,n),v11,ex1,ex2,summ2
c
c      local integer variables
c
      integer i,j,k,ni,nj
c
c      local real variables
c
      double precision cnst,dxr,d2xr,d3xr,d4xr,d5xr,dxs,d2xs,
     1       d3xs,d4xs,d5xs,pr,ps,qr,rn,rn1,rn2,rn22,rn23,sum,
     2       two,xr,xs,zero
c
      common /cons/ rn2,rn22,rn23
c
c     initialise constants
c
      data zero /0.0d0/,half /0.5d0/, one /1.0d0/ ,two/2.0d0/
c
      ifault = 1
      if(n .gt. mdim .or. n. lt. 2) return
      ifault = 0
      rn=n
      rn1=rn+one
      rn2=rn+two
      rn22=rn2*rn2
      rn23=rn22*rn2
      nhalf1=(n+1) / 2
c
c     the elements of the upper triangle
c     are first computed
c
      ni=n
      do 50 i=1,nhalf1
         pr=float(i)/rn1
         qr=one-pr
c
c        in function ppnd, xr is computed to satisfy
c        prob(z.lt.pr) = pr, where z is n(0,1)
c
         xr = ppnd(pr, ifault)
         call der(xr,dxr,d2xr,d3xr,d4xr,d5xr)
         do 40 j=i,ni
            if (i .ne. j) goto 30
c
c     if i is equal to j, var(xr) is calculated
c
            v(i,j)=var(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr)
            go to 40
c
c     if i is not equal to j, cov(xr,xs) is calculated
c
   30       ps=float(j)/rn1
            xs=ppnd(ps,ifault)
            call der(xs,dxs,d2xs,d3xs,d4xs,d5xs)
            v(i,j)=cov(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr,
     1                    dxs,d2xs,d3xs,d4xs,d5xs,ps)
            v(j,i)=v(j,i)
   40    continue
         ni=ni-1
   50 continue
c
c        By symmetry the other elements of v will now be filled
c
      nj = n
      do 70 i = 2,n
      njm1 = nj - 1
      im1 = i - 1
      do 60 j = nj, n
      v(i,j) = v(im1, njm1)
      im1 = im1 - 1
   60 continue
      nj = nj - 1
   70 continue
c
c     insert exact values of v(1,1) and v(1,2)
c
      v(1,1)=v11
      v(n,n)=v(1,1)
      v(1,2)=v(1,1)+ex1*(ex1-ex2)-one
      v(2,1)=v(1,2)
      nsub1=n-1
      v(n,nsub1)=v(1,2)
      v(nsub1,n)=v(1,2)
c
c     normalise the first row of v, leaving
c     v(1,1) and v(1,2) fixed
c
      if(n.eq.2) return
      sum=zero
      do 80 j=3,n
         sum=sum+v(1,j)
   80 continue
      cnst=(one-v(1,1)-v(1,2))/sum
      nj=n-2
      do 90 j=3,n
         v(1,j)=v(1,j)*cnst
         v(j,1)=v(1,j)
	 v(n,nj)=v(1,j)
	 v(nj,n)=v(1,j)
         nj=nj-1
   90 continue
c
c     normalise rows 2 through n-1 of v
c
      call rwnorm(v,n,mdim,0)
c
c     modify v(2,2) and its equal
c     v(n-1,n-1) so the trace identity is satisfied
c
      sum=zero
      do 100 k=1,n
         if(k.eq.2.or.k.eq.nsub1) go to 100
         sum=sum+v(k,k)
  100 continue
      v(2,2)=half*(float(n)-summ2-sum)
      v(nsub1,nsub1)=v(2,2)
c
c     renormalise rows 2 through n-1 of v,
c     leaving diagonal elements fixed,
c
      call rwnorm(v,n,mdim,1)
      return
      end
c
      subroutine rwnorm(v,n,mdim,id)
c
c      Appl. Statist algorithm as 128.4 (1978), vol. 27
c      Davis C.S. and Stephens M.A.
c
c      Normalises rows of covariance matrix of normal order statistics
c      so that sum of row elements equals one.
c
c      arguments :      v - array (mdim,mdim) containing the covariance
c                           matrix approximation.
c                       n - sample size.
c                    mdim - row dimension of v in the calling program.
c                      id - 
c
      integer n,mdim,id
      double precision v(mdim,n)
c
c      local integer variables
c
      integer i,j,k,l,m,ni,nj
c
c      local real variables
c
      double precision cnst,one,small,sum,term,zero
c
      data zero /0.0d0/, small /1.0d-12/, one /1.0d0/
c
      nhalf1 = (n+1)/2
      ni=n-1
      do 75 i=2,nhalf1
c
c     find sums of computed terms in each row
c
         sum=zero
         do 55 j=i,ni
            sum=sum+v(i,j)
   55    continue
         if(id .ne.0) sum = sum - v(i,i)
         if(abs(sum).lt.small) go to 75
c
c     normalise rows leaving appropriate elements fixed
c
         k=i-1
         if(id .ne. 0) k = i
         term=zero
         do 60 j=1,k
   60       term=term+v(i,j)
         l=ni+1
         do 65 j=l,n
   65       term=term+v(i,j)
         cnst=(one-term)/sum
         m=i
         if(id.ne.0) m=i+1
         nj=n-m+1
         do 70 j=m,ni
            v(i,j)=v(i,j)*cnst
            v(j,i)=v(i,j)
            v(ni,nj)=v(i,j)
            v(nj,ni)=v(i,j)
            nj=nj-1
   70    continue
         ni=ni-1
   75 continue
      return
      end
c
      subroutine der(x,dx,d2x,d3x,d4x,d5x)
c
c      Appl. Statist. algorithm as 128.1 (1978), vol 27.
c      Davis C.S. and Stephens M.A.
c
c      Computes derivatives for the david-johnson approximation to the
c      variances and covariances of normal order statistics.
c
c      arguments : x - real number at which derivative is calculated.
c                 dx - first derivative of normal probability integral
c                      evaluated at x.
c                  :                 :               :
c                d5x - fifth derivative of normal probability integral
c                      evaluated at x.
c
      double precision x,dx,d2x,d3x,d4x,d5x
c
c      local real variables
c
      double precision forty6,one,onept5,rad2pi,seven,six,
     1       term,twent4,two,twopi,x2
c
c     initialise constants
c
      data one /1.0d0/, onept5 /1.5d0/, two /2.0d0/,
     1     six /6.0d0/,seven /7.0d0/, twent4 /24.0d0/,
     2     forty6 /46.0d0/, rad2pi /2.506628274631d0/,
     3     twopi/6.2831853071796d0/
      x2=x*x
      dx=rad2pi*exp(x2/two)
      d2x=twopi*x*exp(x2)
      d3x=twopi*rad2pi*(two*x2+one)*exp(onept5*x2)
      term=twopi*twopi*exp(two*x2)
      d4x=term*x*(six*x2+seven)
      d5x=term*dx*(x2*(twent4*x2+forty6)+seven)
      return
      end
c
      double precision function var(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr)
c
c      Appl. Statist. algorithm 128.2 (1978), vol 27
c
c      Computes David-Johnson approximation for the variance of
c      the rth largest order statistic from the normal dist. for
c      a sample size n.
c
c      arguments : dxr - first derivative of normal probability integral
c                        evaluated at xr.
c                   :                   :                 :
c                 d5xr - fifth derivative of normal probability integral
c                        evaluated at xr.
c                   pr - expected value of rth largest order statistic
c                        from uniform dist. ( = r/(n+1) r=1,...,n).
c                   qr - 1-pr
c                        n.b. xr is the inverse normal probability integral
c                        of pr.
c
      double precision dxr,d2xr,d3xr,d4xr,d5xr,pr,qr
c
c      local real variables
c
      double precision d2xr2,dxr2,fiveth,fourth,half,onept5,prqr,
     1       qrmpr,rn2,rn22,rn23,two,three
c
      common /cons/ rn2,rn22,rn23
c
c      initialise constants
c
      data fourth /0.25d0/, half /0.5d0/, onept5 /1.5d0/,
     1     fiveth /1.6666666667d0/, two /2.0d0/ , three /3.0d0/
      dxr2=dxr*dxr
      prqr=pr*qr
      var=prqr*dxr2/rn2
c
c     to order (n+2)**(-2)
c
      qrmpr=qr-pr
      d2xr2=d2xr*d2xr
      var = var+prqr/rn22*(two*qrmpr*dxr*d2xr+prqr*
     1    (dxr*d3xr+half*d2xr2))
c
c     to order (n+2)**(-3)
c
      var = var+prqr/rn23*(-two*qrmpr*dxr*d2xr+
     1    (qrmpr*qrmpr-prqr)*(two*dxr*d3xr+onept5*d2xr2)
     2    +prqr*qrmpr*(fiveth*dxr*d4xr+three*d2xr*d3xr)
     3    +fourth*prqr*prqr*(dxr*d5xr+two*d2xr*d4xr+
     4    fiveth*d3xr*d3xr))
      return
      end
c
      double precision function cov(dxr,d2xr,d3xr,d4xr,d5xr
     1                     ,pr,qr,dxs,d2xs,d3xs,d4xs,d5xs,ps)
c
c      Appl. statist. algorithm as 128.3 (1978), vol. 27
c      Davis C.S. and Stephens M.A.
c
c      Computes David-Johnson approximation for covariance between rth
c      and sth order statistics from the normal dist. for a sample size n.
c
c      arguments : dxr - first derivative of normal probability integral
c                        evaluated at xr.
c                   :                     :               :
c                 d5xr - fifth derivative of normal probability integral
c                        evaluated at xr.
c                   pr - expected value of rth order statistic from
c                        uniform dist. ( = r/(n+1) r=1,...n ).
c                   qr - 1-pr.
c                  dxs - first derivative of normal probability integral
c                        evaluated at xs.
c                   :                     :             :
c                 d5xs - fifth derivative of normal probability integral
c                        evaluated at xs.
c                   ps - expected value of sth order statistic from
c                        uniform distribution. ( = s/(n+1) s=1,...,n).
c                        n.b. xr is the inverse normal probability
c                        integral of pr etc.
c
      double precision dxr,d2xr,d3xr,d4xr,d5xr,pr,qr,dxs,d2xs,
     1         d3xs,d4xs,d5xs,ps
c
c      local real variables
c
      double precision eigth,five6,fourth,half,one,onept5,pr2,
     1       prqr,prqs,ps2,psqr,psqs,qr2,qrmpr,qs,qs2,qsmps,
     2       rn2,rn22,rn23,term1,term2,term3,term4,term5,
     3       three,twelth,two
c
      common /cons/ rn2,rn22,rn23
c
c      initialise constants
c
      data twelth /0.0833333333333d0/, eigth /0.125d0/,
     1     fourth /0.25d0/, half /0.5d0/,
     2     five6 /0.833333333333d0/, one /1.0d0/,
     3     onept5 /1.5d0/, two /2.0d0/, three /3.0d0/
c
      qs=one-ps
      prqs=pr*qs
      cov=prqs*dxr*dxs/rn2
c
c     to order (n+2)**(-2)
c
      qrmpr=qr-pr
      qsmps=qs-ps
      prqr=pr*qr
      psqs=ps*qs
      cov = cov+prqs/rn22*(qrmpr*d2xr*dxs+
     1    qsmps*dxr*d2xs+half*prqr*d3xr*dxs+
     2    half*psqs*dxr*d3xs+half*prqs*d2xr*d2xs)
c
c     to order (n+2)**(-3)
c
      pr2=pr*pr
      qr2=qr*qr
      ps2=ps*ps
      qs2=qs*qs
      psqr=ps*qr
      term1=-d2xr*dxs*qrmpr-qsmps*dxr*d2xs+
     1      (qrmpr*qrmpr-prqr)*d3xr*dxs
      term2=(qsmps*qsmps-psqs)*dxr*d3xs+(onept5*qrmpr*
     1       qsmps+half*psqr-two*prqs)*d2xr*d2xs
      term3=five6*(prqr*qrmpr*d4xr*dxs+psqs*qsmps*dxr*
     1      d4xs)+(prqs*qrmpr+half*prqr*qsmps)*d3xr*d2xs
      term4=(prqs*qsmps+half*psqs*qrmpr)*d2xr*d3xs+
     1      eigth*(pr2*qr2*d5xr*dxs+ps2*qs2*dxr*d5xs)
      term5=fourth*(pr2*qr*qs*d4xr*d2xs+pr*ps*qs2*
     1      d2xr*d4xs)+twelth*(two*pr2*qs2+three*pr*qr*
     2      ps*qs)*d3xr*d3xs
      cov = cov+prqs/rn23*(term1+term2+term3+term4+term5)
      return
      end
c
	double precision function v11(n, ifault)
c
c	ASR 72 (Remark on AS 128) Applied Stats. (1988) vol. 37 (1)
c
c	Calculates an approximation to the variance of the largest
c	normal order statistic.
c
	integer n, ifault
c
c	Local variables
c
	double precision zero, one, a0, a1, a2, a3, a4, a5, a6, x,
     +  	d0, d1, d2, d3,	d4, d5, d6, pt09, c0, c1, c2, c3,
     +          c4, c5, c6, c7, c8, c9, mpt15, b0, b1, b2, b3, b4
	parameter (mpt15 = -0.15, b0 = -0.934d-4, zero = 0.d0,
     +          one = 1.d0,
     +		b1 = -0.5950321d0, b2 = 0.0165504d0, b3 = 0.0056975d0,
     +		pt09 = 0.091105452691946d0, c0 = 0.7956d-11,
     +		c1 = -0.595628869836878d0, c2 = 0.08967827948053d0,
     +		c3 = -0.007850066416039d0, c4 = -0.296537314353d-3,
     +		c5 = 0.215480033104d-3, c6 = -0.33811291323d-4,
     +		c7 = 0.2738431187d-5, c8 = -0.106432868d-6,
     +		c9 = 0.1100251d-8, a0 = 0.04619831847696d0,
     +		a1 = -0.147930264017706d0, a2 = -0.451288155800301d0,
     +		a3 = 0.010055707621709d0, a4 = 0.007412441980877d0,
     +		a5 = -0.001143407259055d0, a6 = 0.54428754576d-4,
     +		d0 = 0.093256818332708d0, d1 = 1.336952989217635d0,
     +		d2 = -1.783195691545387d0, d3 = 0.488682076188729d0,
     +		d4 = -0.078737246197474d0, d5 = 0.00662561987806d0,
     +		d6 = -0.226486218258d-3, b4 = -0.8531d-3)
c
	v11 = zero
	ifault = 1
	if (n .lt. 1) return
	ifault = 0
	if (n .eq. 1) then
	  v11 = one
	  return
	end if
c
	x = n
	if (n .gt. 370) then
	  x = (x**mpt15 - one) / mpt15
	  v11 = exp(b0 + x*(b1 + x*(b2 + x*(b3 + x*b4))))
	else if (n .le. 100) then
	  x = (x**pt09 - one) / pt09
	  v11 = exp(c0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + x*
     +		(c6 + x*(c7 + x*(c8 + x*c9)))))))))
	else if (n .le. 200) then
	  x = log(a0 + x)
	  v11 = exp(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + x*a6)))))
	else
	  x = log(d0 + x)
	  v11 = exp(d1 + x*(d2 + x*(d3 + x*(d4 + x*(d5 + x*d6)))))
	end if
c
	return
	end
