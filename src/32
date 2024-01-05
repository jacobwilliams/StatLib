      double precision function gamain (x,p,g,ifault)
      implicit double precision (a-h,o-z)
c
c        Algorithm AS 32 J.R. Statist. Soc. C.  (1970) Vol.19 No. 3
c	 Algorithm AS 239 is recommended as an alternative.
c
c        Computes incomplete gamma ratio for positive values of
c        arguments x and p.  G must be calculated and should be equal
c        to ln(gamma(p)).   Algorithm AS 245 may be used for this.
c        ifault=1 if p.le.0 else 2 if x.lt.0 else 0.
c        Uses series expansion if p.gt.x or x.le.1, otherwise
c        continued fraction approximation.
c
c	 Revised to incorporate the recommendations of Rice & Das,
c	 Appl. Statist., 34, 326, 1985, and of Cran, Appl. Statist.,
c	 38, 423, 1989.
c
      dimension pn(6)
      data zero/0.0d0/, one/1.0d0/, uflo/1.0d-37/, two/2.0d0/
      data acu/1.0d-8/,oflo/1.0d37/
c
      g=alngam(p)
c
c        Define accuracy and initialize
c
      gin=zero
      ifault=0
c
c        Test for admissibility of arguments
c
      if (p.le.zero) ifault=1
      if(x.lt.zero) ifault=2
      if (ifault .gt. 0 .or. x .eq. zero) go to 50
      arg=p*log(x)-x-g
      if(arg.lt.log(uflo)) then
	ifault=3
	go to 50
      end if
      factor=exp(arg)
      if(x.gt.one .and. x.ge.p) goto 30
c
c        Calculation by series expansion
c
      gin=one
      term=one
      rn=p
   20 rn=rn+one
      term=term*x/rn
      gin=gin+term
      if(term.gt.acu) goto 20
      gin=gin*factor/p
      goto 50
c
c        Calculation by continued fraction
c
   30 a=one-p
      b=a+x+one
      term=zero
      pn(1)=one
      pn(2)=x
      pn(3)=x+one
      pn(4)=x*b
      gin=pn(3)/pn(4)
c
   32 a=a+one
      b=b+two
      term=term+one
      an=a*term
      do 33 i=1,2
   33 pn(i+4)=b*pn(i+2)-an*pn(i)
      if(pn(6).eq.zero) goto 35
      rn=pn(5)/pn(6)
      dif=abs(gin-rn)
      if(dif.gt.acu) goto 34
      if(dif.le.acu*rn) goto 42
   34 gin=rn
   35 do 36 i=1,4
   36 pn(i)=pn(i+2)
      if(abs(pn(5)).lt.oflo) goto 32
      do 41 i=1,4
   41 pn(i)=pn(i)/oflo
      goto 32
   42 gin=one-factor*gin
   50 gamain=gin
      return
      end
