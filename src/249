c This file contains algorithm AS 249 for the evaluation of the mean and
c covariance of the truncated multinormal distribution, followed by code
c by the same authors for the numerical estimation of the observed
c information matrix (= maximum likelihood covariance estimates).  As
c the latter is not in the published algorithm, a test program is included
c showing its use.   For further details, contact Phil Leppard, Statistics
c Department, University of Adelaide.   E-mail: PLEPPARD@f.ua.oz
c
c
c
        subroutine mvntrc(n,bounds,u,v,ndim,eps,tmean,tvar,ifault)
c
c	ALGORITHM AS 249  APPL. STATIST. (1989) VOL. 38, NO. 3 (?)
c
c	Auxiliary routines required: MULNOR and MATINV from AS 195.
c       MATINV calls LINRG from the IMSL Stat/Library, but AS 7 can
c       be substituted.
c
        dimension bounds(ndim),u(ndim),v(ndim,ndim),tmean(ndim,1),
     +  tvar(ndim,ndim,1),r(5,5),rinv(5,5),rq(5,5),rqr(5,5),
     +  ind(5),iord(5),sd(5),a(5),aq(5),aqr(5),
     +  phiu(5),phib(5,5),phin1(5),ephin1(5),phin2(5,5),ephin2(5,5),
     +  alpha(3),sum(3),w(10),ww(5,5)
        data zero,one,two/0.0, 1.0, 2.0/
        data  twopi/6.283185308/
        data rtwopi/2.506628275/
        data alphas/1.e-3/
        data ind/5*0/
c
        ifault=0
        if(n.gt.5) go to 100
c
c       Calculate normalised truncation points and correlation
c       matrix R
c
        do 1 i=1,n
          iord(i)=i
          if(v(i,i).le.zero) go to 200
1       sd(i)=sqrt(v(i,i))
        do 2 i=1,n
          a(i)=(bounds(i)-u(i))/sd(i)
          do 2 j=1,n
2       r(i,j)=v(i,j)/(sd(i)*sd(j))
c
c       Order variates by decreasing size of range of integration
c       to increase the speed of MULNOR at highest level
c       (Smallest standardised boundary = largest range, etc)
c
        if(n.gt.2) then
        do 3 i=1,n-1
          do 3 j=i+1,n
            if(a(i).gt.a(j)) then
            aa=a(i)
            a(i)=a(j)
            a(j)=aa
            k=iord(i)
            iord(i)=iord(j)
            iord(j)=k
            do 4 k=1,n
              aa=r(i,k)
              r(i,k)=r(j,k)
4           r(j,k)=aa
            do 5 k=1,n
              aa=r(k,i)
              r(k,i)=r(k,j)
5           r(k,j)=aa
            end if
3       continue
        end if
c
c
c       Evaluate n dimensional integral alpha = Phi(a:R) and error bound
c
c       alpha(1)=alpha+error   alpha(2)=alpha   alpha(3)=alpha-error
c
        l=0
        do 6 i=1,n
          do 6 j=1,i-1
            l=l+1
6       w(l)=r(i,j)
        call mulnor(w,a,w,eps,n,ind,alpha(2),error,i)
        alpha(1)=alpha(2)+error
        alpha(3)=alpha(2)-error
        if(i.ne.0.or.alpha(2).le.alphas) go to 400
c
c
c       Evaluate univariate and bivariate normal densities
c       Initialise multinormals of dimension n-1 and n-2
c
c
        do 7 i=1,n
          phiu(i)=exp(-a(i)*a(i)/two)/rtwopi
          phin1(i)=one
          ephin1(i)=zero
          do 7 j=i+1,n
            phib(i,j)=(a(i)*a(i)-two*r(i,j)*a(i)*a(j)+a(j)*a(j))/
     +      (two*(1-r(i,j)*r(i,j)))
            phib(i,j)=exp(-phib(i,j))/(twopi*sqrt(one-r(i,j)*r(i,j)))
            phib(j,i)=phib(i,j)
            phin2(i,j)=one
            phin2(j,i)=one
            ephin2(i,j)=zero
7       ephin2(j,i)=zero
c
c
c       Calculation of n-1 dimensional integrals Phi(Aq:Rq)
c
c
        if(n.gt.1) then
        do 8 i=1,n
          do 8 j=1,n
8       rinv(i,j)=r(i,j)
        call matinv(rinv,5,n,ifault)
        if(ifault.gt.0) go to 300
        do 9 iq=1,n
          m1=0
c
c       Determine bounds of integration,vector Aq
c
          do 10 i=1,n
            if(i.ne.iq) then
            m1=m1+1
            aq(m1)=(a(i)-r(i,iq)*a(iq))/sqrt(one-r(i,iq)*r(i,iq))
            m2=0
            do 11 j=1,n
              if(j.ne.iq) then
              m2=m2+1
              ww(m1,m2)=rinv(i,j)
              end if
11          continue
            end if
10        continue
c
c       Determine n-1 x n-1 correlation matrix Rq
c
          call matinv(ww,5,n-1,ifault)
          if(ifault.gt.0) go to 300
          do 13 i=1,m1
            do 13 j=1,m1
13          rq(i,j)=ww(i,j)/sqrt(ww(i,i)*ww(j,j))
c
c       Evaluate n-1 dimensional integral  Phi(Aq:Rq) and error bound
c
            l=0
            do 14 i=1,m1
              do 14 j=1,i-1
                l=l+1
14          w(l)=rq(i,j)
            call mulnor(w,aq,w,eps,n-1,ind,phin1(iq),ephin1(iq),i)
9       if(i.ne.0) go to 500
        end if
c
c       Calculation of n-2 dimensional integrals Phi(Aqr:Rqr)
c
        if(n.gt.2) then
        do 15 iq=1,n
          do 15 ir=iq+1,n
c
c       Determine n-2 x n-2 correlation matrix Rqr
c
          m1=0
          do 16 i=1,n
            if(i.ne.iq.and.i.ne.ir) then
            m1=m1+1
            m2=0
            do 17 j=1,n
              if(j.ne.iq.and.j.ne.ir) then
              m2=m2+1
              ww(m1,m2)=rinv(i,j)
              end if
17          continue
            end if
16        continue
          call matinv(ww,5,n-2,ifault)
          if(ifault.gt.0) go to 300
          do 19 i=1,m1
            do 19 j=1,m1
19        rqr(i,j)=ww(i,j)/sqrt(ww(i,i)*ww(j,j))
c
c       Determine bounds of integration,vector Aqr
c
          m2=0
          do 20 i=1,n
            if(i.ne.iq.and.i.ne.ir) then
            m2=m2+1
            bsqr=(r(iq,i)-r(iq,ir)*r(ir,i))/(one-r(iq,ir)*r(iq,ir))
            bsrq=(r(ir,i)-r(iq,ir)*r(iq,i))/(one-r(iq,ir)*r(iq,ir))
            rsrq=(one-r(i,iq)*r(i,iq))*(one-r(ir,iq)*r(ir,iq))
            rsrq=(r(i,ir)-r(i,iq)*r(iq,ir))/sqrt(rsrq)
            aqr(m2)=a(i)-bsqr*a(iq)-bsrq*a(ir)
            aqr(m2)=aqr(m2)/sqrt((one-r(i,iq)*r(i,iq))*(one-rsrq*rsrq))
            end if
20        continue
c
c       Evaluate n-2 dimensional integral Phi(Aqr:Rqr) and error bound
c
          l=0
          do 21 i=1,m2
            do 21 j=1,i-1
              l=l+1
21        w(l)=rqr(i,j)
          call mulnor(w,aqr,w,eps,n-2,ind,phin2(iq,ir),ephin2(iq,ir),i)
          if(i.ne.0) go to 600
          phin2(ir,iq)=phin2(iq,ir)
15      ephin2(ir,iq)=ephin2(iq,ir)
        end if
c
c       Calculation of E(Xi) ,with upper and lower bounds
c       Tallis(1961) , equation (3).
c
        do 22 i=1,n
          do 23 j=1,3
23        sum(j)=zero
          do 24 j=1,n
            aa=r(i,j)*phiu(j)
            do 25 k=1,3
25          w(k)=aa*(phin1(j)+(k-2)*ephin1(j))
            if(w(1).gt.w(3)) then
            aa=w(1)
            w(1)=w(3)
            w(3)=aa
            end if
            do 24 k=1,3
24        sum(k)=sum(k)+w(k)
          k=iord(i)
          do 22 j=1,3
22      tmean(k,j)=sum(j)/alpha(j)
c
c       Calculation of E(Xi,Xj), with upper and lower bounds
c       Tallis(1961), equation (4).
c
        do 26 i=1,n
          do 26 j=i,n
            do 27 k=1,3
27          sum(k)=zero
            do 28 k=1,n
              aa=r(k,i)*r(k,j)*a(k)*phiu(k)
              do 29 m1=1,3
29            w(m1)=aa*(phin1(k)+(m1-2)*ephin1(k))
              if(w(1).gt.w(3)) then
              aa=w(1)
              w(1)=w(3)
              w(3)=aa
              end if
              do 30 m1=1,3
30            sum(m1)=sum(m1)+w(m1)
              do 28 l=1,n
                if(l.ne.k) then
                aa=(r(l,j)-r(k,l)*r(k,j))*r(k,i)*phib(k,l)
                do 31 m1=1,3
31              w(m1)=aa*(phin2(k,l)+(m1-2)*ephin2(k,l))
                if(w(1).gt.w(3)) then
                aa=w(1)
                w(1)=w(3)
                w(3)=aa
                end if
                do 32 m1=1,3
32              sum(m1)=sum(m1)+w(m1)
                end if
28          continue
            m1=iord(i)
            m2=iord(j)
            do 26 m3=1,3
              tvar(m1,m2,m3)=r(i,j)+sum(m3)/alpha(m3)
26      tvar(m2,m1,m3)=tvar(m1,m2,m3)
c
c       Calculate truncated covariance matrix, with upper and
c       lower bounds. Original scale restored.
c
        do 33 i=1,n
          do 33 j=i,n
            tvar(i,j,2)=(tvar(i,j,2)-tmean(i,2)*tmean(j,2))*sd(i)*sd(j)
            l=0
            do 34 m1=1,3,2
              do 34 m2=1,3,2
                do 34 m3=1,3,2
                  l=l+1
34          w(l)=(tvar(i,j,m1)-tmean(i,m2)*tmean(j,m3))*sd(i)*sd(j)
            tvar(i,j,1)=min(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8))
            tvar(i,j,3)=max(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8))
            do 33 k=1,3
33      tvar(j,i,k)=tvar(i,j,k)
c
c       Original location and scale restored to truncated mean
c
        do 35 i=1,n
          do 35 j=1,3
35      tmean(i,j)=u(i)+sd(i)*tmean(i,j)
c
c       Normal termination of subroutine
c
        return
c
c       Error conditions
c
c
100     ifault=100
c       Dimension n greater than 5
        return
c
200     ifault=200+i
c       i-th variance less than or equal to zero
        return
c
300     ifault=300+i
c       Error in inversion of i x i matrix
        return
c
400     ifault=400+i
c       Error in Phi(n)...MULNOR error number i
c               or
c       truncated region, ie alpha, "too small" (then ifault=400)
c
        return
c
500     ifault=500+i
c       Error in Phi(n-1)...MULNOR error number i
        return
c
600     ifault=600+i
c       Error in Phi(n-2)...MULNOR error number i
        return
c
        end
c
c  **   End of AS 249
c
c-------------------------------------------------------------------------
c
        program covtest
        implicit double precision (a-h,o-z)
        dimension cov(5,5),a(5),b(5),x(5)
        common r(1000),nn,y(1000,3)
        external fnorm,fmult
        data a,b/10*777.d0/

C       TEST RUN FROM N(0,4)

        nn=100
        x(1)=0.
        x(2)=0.
        do 1 i=1,nn
          r(i)=dnorin((i-.5)/nn)*2
          x(1)=x(1)+r(i)/nn
1       x(2)=x(2)+r(i)*r(i)
        x(2)=(x(2)-nn*x(1)*x(1))/nn
        print 5,x(1),x(2)
5       format(//5x'Mle = ',5f10.3)
        call mlecov(fnorm,x,2,a,b,cov,5,10,1.d-06,1,ifault)
        print 3,ifault
3       format(//5x'Final estimate   ifault='i4)
        do 4 i=1,2
4       print 2,(cov(i,j),j=1,2)
2       format(5x,3f10.5)
        do 20 i=1,2
20      cov(i,i)=sqrt(cov(i,i))
        cov(1,2)=cov(1,2)/(cov(1,1)*cov(2,2))
        cov(2,1)=cov(1,2)
        print 26
26      format(40x'Standard errors and correlations')
        do 21 i=1,2
21      print 22,(cov(i,j),j=1,2)
22      format(40x,3f10.5)

C       TEST RUN FROM MULTINOMIAL P=(.5,.2,.2,.1)

        nn=4
        r(1)=500
        r(2)=200
        r(3)=200
        r(4)=100
        x(1)=.5
        x(2)=.2
        x(3)=.2
        x(4)=.1
        print 5,(x(i),i=1,3)
        call mlecov(fmult,x,3,a,b,cov,5,20,1.d-10,1,ifault)
        print 3,ifault
        do 8 i=1,3
8       print 2,(cov(i,j),j=1,3)
        do 23 i=1,3
23      cov(i,i)=sqrt(cov(i,i))
        do 25 i=1,3
        do 25 j=1,3
25      if(i.ne.j) cov(i,j)=cov(i,j)/(cov(i,i)*cov(j,j))
        print 26
        do 24 i=1,3
24      print 22,(cov(i,j),j=1,3)
        end
c
c
c
        subroutine matinv(cov,nd,n,ifault)
        implicit double precision (a-h,o-z)
        dimension cov(nd,1)
c
        call dlinrg(n,cov,nd,cov,nd)
        ifault=0
        return
        end
c
c
c
        double precision function fnorm(x)
        implicit double precision (a-h,o-z)
        common r(1000),nn
        dimension x(1)
c
        fnorm=0.
        do 73 i=1,nn
73      fnorm=fnorm-(r(i)-x(1))*(r(i)-x(1))
        fnorm=fnorm/(2*x(2))-nn*log(x(2))/2
        return
        end
c
c
c
        double precision function fmult(x)
        implicit double precision (a-h,o-z)
        dimension x(1),p(5)
        common r(1000),nn
c
        fmult=0.
        p(nn)=1.
        do 1 i=1,nn-1
        p(i)=x(i)
1       p(nn)=p(nn)-p(i)
        do 2 i=1,nn
2       fmult=fmult+r(i)*log(p(i))
        return
        end
c
c
c  **   End of test program and functions
c----------------------------------------
c
        subroutine mlecov(f,theta,n,a,b,cov,nd,maxh,eps,iprint,ifault)
        implicit double precision (a-h,o-z)
        dimension cov(nd,1),theta(*),a(*),b(*),x(5),h(5),store(5,5)
	external f
        data zero,small,two,hun/0.d0, .001d0, 2.d0, 100.d0/
c
        ifault=0
        if(n.gt.5.or.n.gt.nd) go to 100
c
c       Determine initial values for h(i), i=1,..n
c
        do 1 i=1,n
          h(i)=small
          h(i)=max(h(i),abs(theta(i))/hun)
          if(a(i).eq.b(i)) go to 1
          if(theta(i).lt.a(i).or.theta(i).gt.b(i)) go to 200
          h(i)=min(h(i),abs(a(i)-theta(i)),abs(b(i)-theta(i)))
1       continue
        f00=f(theta)
        ncy=0
2       ncy=ncy+1
        if(ncy.gt.maxh) go to 300
c
c       Evaluation of diagonal elements of -I(theta)
c
        do 3 i=1,n
          cov(i,i)=-two*f00
          do 4 j=1,n
4         x(j)=theta(j)
          x(i)=x(i)-h(i)
          do 3 j=1,2
            cov(i,i)=cov(i,i)+f(x)
3       x(i)=x(i)+two*h(i)
c
c       Evaluation of off-diagonal elements of -I(theta)
c
        do 5 i=1,n
          jj=i+1
          if(jj.gt.n) go to 5
          do 7 j=jj,n
            cov(i,j)=cov(i,i)+cov(j,j)+two*f00
            do 6 k=1,n
6           x(k)=theta(k)
            x(i)=x(i)-h(i)
            x(j)=x(j)-h(j)
            do 8 k=1,2
              cov(i,j)=cov(i,j)-f(x)
              x(i)=x(i)+two*h(i)
8           x(j)=x(j)+two*h(j)
            cov(i,j)=cov(i,j)/(two*h(i)*h(j))
7         cov(j,i)=cov(i,j)
5       cov(i,i)=-cov(i,i)/(h(i)*h(i))
c
c       In-place inversion of COV
c
        call matinv(cov,nd,n,ifault)
        if(ifault.ne.0) go to 400
c
c       Check for convergence
c
        s=zero
        if(ncy.eq.1) go to 9
        do 10 i=1,n
          do 10 j=i,n
10      s=max(s,abs(cov(i,j)-store(i,j)))
c
c       Display option if requested
c
9       if(iprint.eq.0) go to 11
        print 12,ncy-1,s
12      format(//5x'Covariance estimate at halving ',i3/
     +  5x'Maximum absolute difference from previous estimate = ',e14.7)
        do 13 i=1,n
13      print 14,(cov(i,j),j=1,n)
14      format(5x,5(e14.7,2x))
c
c       Vector h halved and process repeated
c
11      if(s.le.eps.and.ncy.gt.1) return
        do 15 i=1,n
          h(i)=h(i)/two
          do 15 j=i,n
15      store(i,j)=cov(i,j)
        go to 2
c
c       Error conditions
c
100     ifault=100
c       Dimension or problem size misspecified
        return
200     ifault=200+i
c       Maximum likelihood estimate for i-th parameter outside bounds
        return
300     ifault=300
c       Convergence not achieved in maxh iterations
        return
400     ifault=400+ncy
c       Matrix inversion difficulties at iteration ncy
        return
        end
