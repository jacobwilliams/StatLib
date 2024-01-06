      subroutine dssp(x,xmean,xssp,wt,sumwt,nvar,nunit,ifault)
c
c        Algorithm AS 41 j.r.statist.soc.c. (1971) vol. 20 no.2
c
c          This subroutine updates the mean vector xmean (length nvar)
c          and the matrix of corrected sums of squares and products xssp
c          (length nvar(nvar+1)/2, stored by lower triangle), when a
c          data vector x (length nvar) with weight wt is either included
c          (wt.gt.0) or excluded (wt.lt.0).  sumwt is the current sum of
c          weights on entry and the updated sum on exit and nunit is 
c          the current and updated sample size.  ifault=0 indicates normal
c          exit,  ifault=1 indicates zero or negative value of sumwt,
c          ifault=2 indicates zero or negative nunit, ifault=3 indicates
c          nvar.lt.1.  Note that x, xmean, xssp, wt and sumwt are double
c          precision and must be declared as such in the calling program.
c
      dimension x(*), xmean(*), xssp(*)
      double precision x, xmean, xssp, wt, sumwt, b, c, co
      data co/0.0d0/
c
c          Check variates, weights and sample size
c
      ifault = 0
      if(nvar.lt.1) goto 103
      if(wt)107,100,106
  107 nunit = nunit-1
      go to 105
  106 nunit = nunit+1
  105 sumwt = sumwt+wt
      if(sumwt.le.co) go to 101
      k = 0
      b = wt/sumwt
      if(nunit-1)102,120,110
c
c          Update means and ssp for sample size greater than 1
c
  110 c = wt-b*wt
      do 111 i = 1,nvar
        x(i) = x(i)-xmean(i)
        xmean(i) = xmean(i)+b*x(i)
        do 111 j = 1,i
          k = k+1
          xssp(k) = xssp(k)+c*x(i)*x(j)
  111 continue
      return
c
c          Initialise means and ssp for sample size  =  1
c
  120 do 121 i = 1,nvar
        if(wt.lt.co) goto 122
        xmean(i) = x(i)
        goto 123
  122   xmean(i) = xmean(i)+b*(x(i)-xmean(i))
  123   x(i) = co
        do 121 j = 1,i
          k = k+1
          xssp(k) = co
  121 continue
      return
c
c          Set fault indicators
c
  103 ifault = ifault+1
  102 ifault = ifault+1
  101 ifault = ifault+1
  100 return
      end
