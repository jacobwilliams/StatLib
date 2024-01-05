      Program driver
C
C             Program to test and illustrate algorithm BEST.  Program
C             DRIVER is not necessary for using BEST.
C  
C             At the prompt, input values for p, q, and x.  The value of
C             p is the dimension of the covariance matrix and the value
C             of q is the degrees of freedom.  The program outputs an
C             approximation to Pr(GG < x), where GG is the Greenhouse-
C             Geiser epsilon.
C
C             The program is intended to be run interactively.  It expects
C             that unit 5 is input and unit 6 is output.
C
      real x,best,pval
      integer ip,iq
      write(6,10)
 10   format(/,' To exit program, type Control C',/)
 20   write(6,30)
 30   format(/,' Input p, q, and x')
      read(5,*,err=50) ip,iq,x
      pval=best(x,ip,iq,ifault)
      write(6,40) pval,ifault
 40   format(' Pr(GG < x) = ',f7.5,', Ifault = ',i1)
      goto 20
 50   stop
      end



      real function best(x,ip,iq,ifault)
C
C      ALGORITHM AS 284, APPL. STATIST. (1993), VOL. 42, NO. 3
C
C             Computes The Null Distribution Function of The Test
C             Statistic For The Locally Best Invariant Test of
C             Sphericity And Additivity
C
      real x,p,q,a,b,t,s1,s2,s3,s4,s5,zero,one,two,three,four,epsi,
     $ delta(4),fact(0:6),pocham,zreal,alogam,betain,beta,small,lag,
     $ nbit,base
      intrinsic abs,exp,log,sign,sqrt
C
C             Fixed Constants
C
      data zero,one,two,three,four,small,(fact(i),i=0,6) /0.0e0, 
     $ 1.0e0, 2.0e0, 3.0e0, 4.0e0, 1.0e-2, 1.0e0, 1.0e0, 2.0e0,
     $ 6.0e0, 2.4e1, 1.2e2, 7.2e2/
C
C             Machine Dependent Constants
C
      data base,nbit/2.0e0, 23.0e0/
      zreal(i)=real(i)
      epsi=one/base**(nbit-one)
C
C             Check Validity of Input Arguments
C
      best=one
      ifault=4
      if(ip.le.1.or.iq.lt.ip) return
      p=zreal(ip)
      q=zreal(iq)
      ifault=5
      if(x.ge.one-epsi) return
      best=zero
      if(x*p.le.one+epsi) return
C
C             Use Chebyshev's Inequality to Check for Probabilities
C             Near 0 or 1
C
      ifault=3
      t=(one/x-one)/(p-one)
      a=(p+two)*(p*(p*q-q+two*four)-three*four)/(four*p*(q+two))
      b=(q-one)*a*p/(p+two)
      s1=a/(a+b)-t
      s2=(log(a)+log(b)-log(a+b+one))/two-log(a+b)-log(small)
      if(exp(s2).le.abs(s1)) then
         if(s1.lt.zero) then
            return
         else
            best=one
            return
         endif
      endif
      ifault=0
C
C             Use Exact Beta(1,(q-1)/2) Distribution if p = 2
C
      if(ip.eq.2) then
          best=exp((q-one)*log(one-t)/two)
C
C             Check Magnitude of Higher Order Terms.
C
      else
          call moment(p,q,delta)
          s1=abs(delta(1))+abs(delta(2))+abs(delta(3))+abs(delta(4))
          if(s1.ge.epsi**(one/three)) then
C
C             Use Exact Distribution from John (1972) if p = 3
C             and Magnitude of Higher Order Terms is Large
C
              if(ip.eq.3) then
                  a=q/two
                  b=q-one
C
C             Compute log Beta Function.  No Need to Check Error
C             Indicator, it is Known to be 0.
C
                  beta=alogam(a,ier)+alogam(b,ier)-alogam(a+b,ier)
                  s1=sqrt(t)
                  s2=(one+two*s1)/three
                  best=betain(s2,a,b,beta,ier)
     $             -three*exp(a*log(s2)+b*log(one-s2)-beta)
                  if(t.lt.one/four-epsi) then
                      s3=(one-two*s1)/three
C
C             Compute Probability.  No Need to Check Error
C             Indicator, it is Known to be 0.
C
                      best=best-betain(s3,a,b,beta,ier)
     $                 +three*exp(a*log(s3)+b*log(one-s3)-beta)
                  end if
                  best=one-best
C
C             Use Six moment Jacobi Polynomial Expansion if p > 3
C             And Magnitude of Higher Order Terms is Large
C
              else
                  ifault=1
C
C             Compute log Beta Function And Probability.  No Need to
C             Check Error Indicator, it is Known to be 0.
C
                  beta=alogam(a,ier)+alogam(b,ier)-alogam(a+b,ier)
                  best=betain(t,a,b,beta,ier)
                  do 30 n=3,6
                      s1=pocham(a,n-1)+pocham(a+b+zreal(n),n-1)
     $                 -pocham(a+b,n)+pocham(a+b+zreal(n-1),n-1)
     $                 -pocham(b,n)-log(fact(n))
                      s2=zero
                      do 10 i=0,n-3
                          s3=pocham(zreal(n-i+1),i)-log(fact(i))
     $                     +log(a+zreal(n-1))
     $                     -pocham(a+b+zreal(2*n-i-1),i)
     $                     +pocham(a+b+zreal(n-i),i)
     $                     +log(abs(delta(n-i-2)))
     $                     +log(a+b+zreal(2*n-2))
                          s2=s2+exp(s3)*sign(one,delta(n-i-2))
     $                     *(-one)**i
 10                   continue
                      if(abs(s2).le.epsi) goto 30
                      s3=zero
                      do 20 i=1,n
                          s4=pocham(zreal(n-i+1),i)-log(fact(i))
     $                     -pocham(a+b+zreal(2*n-i-1),i)
     $                     +pocham(a+b+zreal(n-i),i)
                          do 20 j=1,i
                              s5=pocham(a+b,n-j)-pocham(a,n-j+1)
     $                         +(a+zreal(n-j))*log(t)+b*log(one-t)
     $                         -beta
                              s3=s3+exp(s5+s4+log(a+b+zreal(2*n-1)))
     $                         *(-one)**i
 20                   continue
                      if(abs(s3).le.epsi) goto 30 
                      s4=exp(s1+log(abs(s2))+log(abs(s3)))
                      best=best+s4*sign(one,s2*s3)
 30               continue
                  best=one-best
              end if
C
C         Use Three moment Laguerre Polynomial Approximation if
C         Magnitude of Higher Order Terms is Small
C
          else
              best=lag(x,ip,iq,ifault)
          end if
      end if
      if(best.lt.zero) best=zero
      if(best.gt.one) best=one
      return
      end


      real function lag(x,ip,iq,ifault)
C
C             Computes Three moment Laguerre Polynomial Approximation
C             to the Distribution Function
C
      real x,p,q,pq,a,b,s1,s2,t,delta,gammad,alogam,zreal,
     $ zero,one,two,three,four,six,cn(9),cd(20)
      intrinsic abs,exp,log,sign
C
C             Fixed Constants
C
      data zero,one,two,three,four,six/ 0.0e0, 1.0e0, 2.0e0, 3.0e0,
     $ 4.0e0, 6.0e0/
      data (cn(i),i=1,9)/ 0.0e0, 1.2e2, 1.8e1, 4.8e1, -1.2e1, -3.0e0,
     $ 4.0e1, 0.0e0, -1.0e0/
      data (cd(i),i=1,20)/ 1.152e3, 3.84e2, -1.84e2, 1.2e1, -1.152e3,
     $ -6.72e2, 6.4e1, 8.0e0, -2.88e2, -1.92e2, 5.4e1, 9.0e0, 2.88e2,
     $ 2.96e2, 4.8e1, 2.0e0, 6.4e1, 8.8e1, 1.8e1, 1.0e0/
      zreal(i)=real(i)
      ifault=2
      p=zreal(ip)
      q=zreal(iq)
      pq=p*q
      a=two*four
      b=-log(a/pq+one)-log((a+two)/pq+one)-three*log(pq)
     $ +log(q-one)+log(q+two)+log(zreal(32))
      b=-exp(b)
C
C             Compute Standardized Third moment
C
      n1=0
      n2=0
      s1=zero
      do 20 j=0,2
          a=zero
          do 10 k=0,2
              n1=n1+1
              a=a/pq+cn(n1)
 10       continue
          s1=s1/p+a
 20   continue
      s2=zero
      do 40 j=0,4
          a=zero
          do 30 k=0,3
              n2=n2+1
              a=a/pq+cd(n2)
 30       continue
          s2=s2/p+a
 40   continue
      delta=b*s1/s2
      a=(p-one)*(p+two)*(p*q+four)*(p*q+six)
     $ /(four*(q-one)*(q+two)*p**2)
      b=(p*q+two)*a/(p+two)
      t=b*(one/x-one)/(p-one)
C
C             Compute Log Gamma Function.  No Need to Check Error
C             Indicator, it is Known to be 0.
C
      s1=a*log(t)-t-alogam(a,ier)+log(abs(delta))-log(six)
      s2=exp(s1+two*log(abs(t-a-two)))-exp(s1+log(abs(a-two)))
C
C             Compute Probability.  No Need to Check Error
C             Indicator, it is Known to be 0.
C
      lag=one-gammad(t,a,ier)+s2*sign(one,delta)
      return
      end

      real function pocham(x,i)
C
C             Compute Natural Log of Pochhammer's Symbol: (x) sub i
C
      real x,zreal,zero
      intrinsic log
      data zero /0.0e0/
      zreal(i)=real(i)
      pocham=zero
      if(i.eq.0) return
      do 10 j=1,i
 10   pocham=pocham+log(x+zreal(i)-zreal(j))
      return
      end

      Subroutine moment(p,q,delta)
      real p,q,pq,num,den,s,zero,one,two,three,four,cn(240),cd(154),
     $ g(4),f(4),delta(4)
      intrinsic exp,log
C
C             Compute Standardized moments of LBI Statistic
C
C             Fixed Constants
C
      data zero,one,two,three,four /0.0e0, 1.0e0, 2.0e0, 3.0e0, 4.0e0/
      data (g(i),i=1,4) /3.2e1, 1.28e2, 6.4e1, 1.28e2/
      data (cn(i),i=1,110)/ 0.0e0, -2.8e1, 0.0e0, 1.0e1, -5.0e0, -1.0e0,
     $ 0.0e0, -1.12896e5, 1.94784e5, -3.3888e4, -3.2448e4, 0.0e0,
     $ 4.032e4, -9.7968e4, 7.6856e4, 2.096e3, -1.732e4, -1.008e3,
     $ 3.552e3, -1.4604e4, 1.002e4, 2.42e3, -3.268e3, -3.28e2, 3.28e2,
     $ -4.34e2, 1.65e2, 3.19e2, -2.45e2, -3.7e1, -8.2e1, 1.16e2, -2.9e1,
     $ 1.0e0, -5.0e0, -1.0e0, 0.0e0, -3.90168576e9, 1.4764032e10,
     $ -1.7053949952e10, 3.842187264e9, 3.411984384e9, -9.30594816e8,
     $ -3.29416704e8, 0.0e0, 1.3934592e9, -5.74580736e9, 1.010539008e10,
     $ -8.198281728e9, 7.0906368e8, 2.382832128e9, -4.64467968e8,
     $ -2.78255616e8, -9.068544e6, -3.787776e7, -5.95800576e8,
     $ 1.800523008e9, -1.554594816e9, -2.6491904e7, 6.44580096e8,
     $ -9.4626048e7, -9.609088e7, -5.374464e6, -1.525248e6, 4.829952e6,
     $ 7.0553536e7, -9.968272e7, -1.3282944e7, 8.1219456e7, -8.673152e6,
     $ -1.7555424e7, -1.35904e6, -6.517632e6, 1.8864096e7, -2.4434256e7,
     $ 6.285264e6, 2.517872e6, 3.552544e6, -3.14064e5, -1.801424e6,
     $ -1.7992e5, 2.71712e5, 1.948336e6, -4.017808e6, 1.359176e6,
     $ 8.19112e5, -2.2064e5, -1.072e4, -9.9608e4, -1.2696e4, -9.7584e4,
     $ 2.53584e5, -2.20896e5, 1.557e4, 7.9326e4, -2.446e4, -2.58e3,
     $ -2.55e3, -4.1e2, 5.056e3, -1.1944e4, 8.55e3, -2.077e3, 8.75e2,
     $ -3.7e2, -6.0e1, -2.5e1, -5.0e0, 0.0e0, -4.12018016256e13/
      data (cn(i),i=111,170) / 2.544631676928e14, -5.8340956962816e14,
     $ 5.799309115392e14, -1.667991011328e14, -1.0232596463616e14,
     $ 5.784648155136e13, 9.4170710016e12, -5.5084843008e12,
     $ -1.03842054144e12, 0.0e0, 1.4714929152e13, -8.270448427008e13,
     $ 2.1422392786944e14, -3.36219910668288e14, 2.93880843804672e14,
     $ -6.9645704552448e13, -7.5702487277568e13, 3.8893641228288e13,
     $ 9.646000324608e12, -4.669805494272e12, -1.170392481792e12,
     $ -2.569273344e10, -4.8813539328e12, 1.241719676928e13,
     $ 4.57903300608e12, -4.737808889856e13, 6.1560777191424e13,
     $ -1.9353568569344e13, -2.04090638848e13, 1.1918790731776e13,
     $ 3.827214399488e12, -1.7948152832e12, -5.66963986432e11,
     $ -2.1138112512e10, 1.2911837184e11, 2.002896101376e12,
     $ -5.52099796992e12, 2.399627067392e12, 4.85868078592e12,
     $ -4.4841321088e12, -2.005616851968e12, 2.151339017216e12,
     $ 7.65947608064e11, -4.0405611776e11, -1.56711717888e11,
     $ -7.879186432e9, -7.0632824832e10, 3.56417660928e11,
     $ -1.269488833536e12, 1.782649662464e12, -3.94358012672e11,
     $ -8.37563249024e11, 1.6834389056e11, 2.62308470784e11,
     $ 7.3064108032e10, -5.86671328e10, -2.7335101312e10,
     $ -1.713773312e9, 3.708598272e10, -5.6949588992e10/
      data (cn(i),i=171,240)/ -1.01034381696e11, 3.1574732672e11,
     $ -1.58358678688e11, -1.11234267456e11, 7.1553216224e10,
     $ 2.4591867072e10, -2.46103264e8, -5.746261696e9, -3.116678112e9,
     $ -2.36430272e8, -2.898390528e9, 3.32693408e9, -1.5851460992e10,
     $ 3.3324926128e10, -1.7291739616e10, -1.069476608e10,
     $ 8.859817856e9, 2.11305632e9, -8.7746096e8, -3.95701248e8,
     $ -2.31185472e8, -2.1039408e7, 4.6900352e8, -3.39609568e8,
     $ -8.90319408e8, 1.59414452e9, -7.83590964e8, -5.89984244e8,
     $ 5.16508044e8, 1.60988444e8, -9.91687e7, -2.0192284e7,
     $ -1.070398e7, -1.17714e6, -7.802352e7, 1.48137968e8, -2.5871288e7,
     $ -7.2898636e7, 2.0969796e7, -5.79076e5, 5.97858e6, 8.300068e6,
     $ -4.904884e6, -7.84316e5, -2.86284e5, -3.8408e4, 4.14784e6,
     $ -1.4920968e7, 2.0030724e7, -1.0498098e7, -8.32893e5, 2.605797e6,
     $ -5.28245e5, 1.13753e5, -9.6255e4, -1.6865e4, -4.115e3, -6.75e2,
     $ -8.3304e4, 2.76992e5, -3.30074e5, 1.73416e5, -5.7129e4, 2.7445e4,
     $ -6.681e3, 1.85e2, -6.75e2, -1.45e2, -2.5e1, -5.0e0/
      data (cd(i),i=1,124)/ 5.76e2, -7.68e2, -1.44e2, 2.56e2, 6.4e1,
     $ -1.92e2, 6.4e1, 4.8e1, 4.8e1, 1.6e1, 1.2e1, 8.0e0, 9.0e0, 2.0e0,
     $ 1.0e0, 1.3824e4, -3.456e4, 1.344e4, 1.632e4, -4.48e3, -3.84e3,
     $ -5.12e2, -1.0368e4, 1.4016e4, 2.528e3, -2.352e3, -2.096e3,
     $ -1.152e3, -1.92e2, 2.208e3, -5.92e2, -4.56e2, -8.6e2, -3.6e2,
     $ -1.08e2, -2.4e1, -1.2e2, -9.2e1, -1.1e2, -3.7e1, -2.1e1, -3.0e0,
     $ -1.0e0, 3.31776e5, -1.327104e6, 1.456128e6, 1.8432e5, -8.0256e5,
     $ -6.144e4, 1.61792e5, 4.9152e4, 4.096e3, -4.42368e5, 1.179648e6,
     $ -5.28384e5, -4.66944e5, 4.7104e4, 1.08544e5, 6.912e4, 2.048e4,
     $ 2.048e3, 1.98144e5, -2.79552e5, -4.6336e4, 1.6896e4, 5.5136e4,
     $ 3.7824e4, 1.2128e4, 3.072e3, 3.84e2, -3.3792e4, 8.192e3, 6.4e3,
     $ 1.6896e4, 8.064e3, 3.712e3, 1.056e3, 1.92e2, 3.2e1, 1.68e3,
     $ 1.408e3, 1.752e3, 7.2e2, 4.41e2, 1.0e2, 3.8e1, 4.0e0, 1.0e0,
     $ 7.962624e6, -4.644864e7, 9.068544e7, -4.902912e7, -3.9020544e7,
     $ 3.236352e7, 1.3006848e7, -5.44768e6, -3.35872e6, -5.7344e5,
     $ -3.2768e4, -1.65888e7, 7.133184e7, -8.626176e7, -8.41728e5,
     $ 3.8708736e7, 5.992192e6, -5.542144e6, -4.4032e6, -1.6896e6,
     $ -3.072e5, -2.048e4, 1.271808e7, -3.621888e7, 1.8376704e7,
     $ 1.2434432e7, 1.097984e6, -3.187328e6, -3.21408e6, -1.3088e6,
     $ -3.4944e5, -6.4e4, -5.12e3, -4.3776e6, 6.517248e6, 9.78944e5/
      data (cd(i),i=125,154)/ 8.0128e4, -1.564096e6, -1.15376e6,
     $ -5.1376e5, -1.768e5, -3.824e4, -6.4e3, -6.4e2, 6.48576e5,
     $ -1.53792e5, -1.22944e5, -3.89792e5, -2.0356e5, -1.1454e5,
     $ -3.78e4, -1.06e4, -2.28e3, -3.0e2, -4.0e1, -3.024e4, -2.7024e4,
     $ -3.4624e4, -1.612e4, -1.041e4, -2.961e3, -1.225e3, -2.1e2,
     $ -6.0e1, -5.0e0, -1.0e0/
      pq=p*q
      s=two*four
      f(1)=log(q-one)+log(q+two)+log(p-two)-log(pq)
     $ -log(p-one)-log(pq+s)-log(pq+s+two)
      do 10 i=2,4
          s=s+four
          f(i)=f(i-1)+three*log(p)+two*log(q)
     $     -log(p-one)-log(pq+s)-log(pq+s+two)
 10   continue
      do 20 i=1,4
 20   f(i)=exp(f(i)+log(g(i)))*(-one)**i
      n1=0
      n2=0
      do 70 ii=1,4
          i=ii+2
          num=zero
          ind=3*(i-3)+1
          do 40 j=0,ind
              s=zero
              do 30 k=0,ind+1
                  n1=n1+1
                  s=s/p+cn(n1)
 30           continue
              num=num/pq+s
 40       continue
          den=zero
          ind=i-1
          do 60 j=0,ind
              s=zero
              do 50 k=0,2*ind
                  n2=n2+1
                  s=s/p+cd(n2)
 50           continue
              den=den/pq+s
 60       continue
          delta(ii)=f(ii)*num/den
 70   continue
      return
      end

C
C
C    *******************************************************************
C    *                                                                 *
C    *    The remaining subprograms are included for referees use      *
C    *    only.  They are not intended for publication.                *
C    *                                                                 *
C    *******************************************************************
C
C
      real function betain(x,p,q,beta,ifault)
C
C         Algorithm AS 63 Appl. Statist. (1973) Vol. 22, p. 409
C
C         Computes Incomplete Beta Function Ratio for Arguments
C         x between zero and one, p and q positive.
C         Log of Complete Beta Function, Beta, is Assumed to
C         be Known.
C
      logical index
      real acu,ai,beta,cx,one,p,pp,psq,q,qq,rx
      real temp,term,x,xx,zero,zabs,zexp,zlog
C
C         Define Accuracy and Initialize
C
      data acu,one,zero/1.0e-7, 1.0e0, 0.0e0/
      zabs(X)=abs(X)
      zexp(X)=exp(X)
      zlog(X)=alog(X)
      betain=x
C
C         Test for Admissibility of Arguments
C
      Ifault=1
      If(p.le.zero.or.q.le.zero) return
      Ifault=2
      If(x.lt.zero.or.x.gt.one) return
      Ifault=0
      if(x.eq.0.or.x.eq.one) return
C
C         Change Tail If Necessary and Determine s
C
      psq=p+q
      cx=one-x
      if(p.ge.psq*x) goto 1
      xx=cx
      cx=x
      pp=q
      qq=p
      index=.true.
      goto 2
 1    xx=x
      pp=p
      qq=q
      index=.false.
 2    term=one
      ai=one
      betain=one
      ns=qq+cx*psq
C
C         Use Reduction Formulae of Soper
C
      rx=xx/cx
 3    temp=qq-ai
      if(ns.eq.0) rx=xx
 4    term=term*temp*rx/(pp+ai)
      betain=betain+term
      temp=zabs(term)
      if(temp.le.acu.and.temp.le.acu*betain) goto 5
      ai=ai+one
      ns=ns-1
      if(ns.ge.0) goto 3
      temp=psq
      psq=psq+one
      goto 4
C
C         Calculate Result
C
 5    betain=betain*zexp(pp*zlog(xx)+(qq-one)*zlog(cx)-beta)/pp
      if(index) betain=one-betain
      return
      end

      real function gammad(x,p,ifault)
C
C          Algorithm AS239 Appl. Statist. (1988) Vol. 37, No. 3
C
C          Computation of the Incomplete Gamma Integral
C
      integer ifault
      real pn1,pn2,pn3,pn4,pn5,pn6,x,tol,oflo,xbig,arg,c,rn,p,a,
     $ b,one,zero,alogam,an,two,elimit,plimit,alnorm,three,nine
      parameter (zero=0.0e0, one=1.0e0, two=2.0e0, oflo=1.3e19,
     $ three=3.0e0, nine=9.0e0, tol=1.0e-7, plimit=1.0e3,
     $ xbig=1.0e6, elimit=-8.8e1)
      intrinsic abs,log,exp,sqrt,min
      external alogam,alnorm
C
      gammad=zero
C
C        Check that we have valid values for X and P
C
      if(p.le.zero.or.x.lt.zero) then
      ifault=1
      return
      endif
      ifault=0
      if(x.eq.zero) then
      gammad=zero
      return
      endif
C
C         Use a Normal Approximation if P .gt. PLIMIT
C
      if(p.gt.plimit) then
      pn1=three*sqrt(p)*((x/p)**(one/three)+one/(nine*p)-one)
      gammad=alnorm(pn1,.false.)
      return
      endif
C
C         If X is extremely large compared to P then set GAMMAD to ONE
C
      if(x.gt.xbig) then
      gammad=one
      return
      endif
C
      if(x.le.one.or.x.lt.p) then
C
C        Use Pearson's series expansion (Note that P is not large enough
C        to force overflow in ALOGAM)  No need to test IFAULT on exit
C        since P .gt. xero
C
      arg=p*log(x)-x-alogam(p+one,ifault)
      c=one
      gammad=one
      a=p
 40   a=a+one
      c=c*x/a
      gammad=gammad+c
      if(c.gt.tol) goto 40
      arg=arg+log(gammad)
      gammad=zero
      if(arg.ge.elimit) gammad=exp(arg)
C
      else
C
C        Use a continued fraction expansion
C
      arg=p*log(x)-x-alogam(p,ifault)
      a=one-p
      b=a+x+one
      c=zero
      pn1=one
      pn2=x
      pn3=x+one
      pn4=x*b
      gammad=pn3/pn4
 60   a=a+one
      b=b+two
      c=c+one
      an=a*c
      pn5=b*pn3-an*pn1
      pn6=b*pn4-an*pn2
      if(abs(pn6).gt.zero) then
      rn=pn5/pn6
      if(abs(gammad-rn).le.min(tol,tol*rn)) goto 80
      gammad=rn
      endif
C
      pn1=pn3
      pn2=pn4
      pn3=pn5
      pn4=pn6
      if(abs(pn5).ge.oflo) then
C
C         Re-scale terms in continued fraction if terms are large
C
      pn1=pn1/oflo
      pn2=pn2/oflo
      pn3=pn3/oflo
      pn4=pn4/oflo
      endif
      goto 60
 80   arg=arg+log(gammad)
      gammad=one
      if(arg.ge.elimit) gammad=one-exp(arg)
      endif
C
      return
      end

      real function alnorm(x,upper)
C
C        Algorithm AS 66 Appl. Statist. (1973) Vol. 22, p. 424
C
C        Evaluates the tail of the standardized normal curve
C        from x to infinity if .upper. is true or from minus
C        infinity to x if upper is .false.
C
      real ltone,utzero,zero,half,one,con,a1,a2,a3,a4,a5,a6,a7,
     $ b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,x,y,z,zexp
      logical upper,up
C
C        Ltone and utzero must be set to suit the particular computer
C        (See introductory text)
C
      data ltone,utzero/5.6e0, 1.2e1/
      data zero,half,one,con/ 0.0e0, 5.0e-1, 1.0e0, 1.28e0/
      data a1,a2,a3,a4,a5,a6,a7/ 3.98942280444e-1, 3.99903438504e-1,
     $ 5.75885480458e0, 2.98213557808e1, 2.62433121679e0,
     $ 4.86959930692e1, 5.92885724438e0/
      data b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12/ 3.98942280385e-1,
     $ 3.8052e-8, 1.00000615302e0, 3.98064794e-4, 1.98615381364e0,
     $ 1.51679116635e-1, 5.29330324926e0, 4.8385912808e0,
     $ 1.51508972451e1, 7.42380924027e-1, 3.0789933034e1,
     $ 3.99019417011e0/
C
      zexp(z)=exp(z)
C
      up=upper
      z=x
      if(z.ge.zero) goto 10
      up=.not. up
      z=-z
 10   if(z.le.ltone.or.up.and.z.le.utzero) goto 20
      alnorm=zero
      goto 40
 20   y=half*z*z
      if(z.gt.con) goto 30
C
      alnorm=half-z*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
      goto 40
C
 30   alnorm=b1*zexp(-y)/(z-b2+b3/(z+b4+b5/(z-b6+b7/
     $ (z+b8-b9/(z+b10+b11/(z+b12))))))
C
 40   if(.not.up) alnorm=one-alnorm
      return
      end
      real function alogam(x,ifault)
C
C            Algorithm ACM 291, Comm. ACM. (1966) Vol. 9, p. 684
C
C            Evaluates Natural Logarithm of Gamma(x)
C
      real a1,a2,a3,a4,a5,f,x,y,z,zlog,half,zero,one,seven
C
C         The Following Constants are Alog(2pi)/2,
C         1/1680, 1/1260, 1/360, and 1/12
C
      data a1,a2,a3,a4,a5/ 9.18938533204673e-1, 5.95238095238e-4,
     $ 7.93650793651e-4, 2.777777777778e-3,
     $ 8.3333333333333e-2/
      data half,zero,one,seven/ 0.5, 0.0, 1.0, 7.0/
C
      zlog(f)=alog(f)
C
      alogam=zero
      ifault=1
      if(x.le.zero) return
      ifault=0
      y=x
      f=zero
      if(y.ge.seven) goto 30
      f=y
 10   y=y+one
      if(y.ge.seven) goto 20
      f=f*y
      goto 10
 20   f=-zlog(f)
 30   z=one/(y*y)
      alogam=f+(y-half)*zlog(y)-y+a1
     $ +(((-a2*z+a3)*z-a4)*z+a5)/y
      return
      end

