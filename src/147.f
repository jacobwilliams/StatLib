      double precision function gammds (y,p,ifault)
c
c        Algorithm AS 147  Appl. Statist. (1980) Vol. 29, No. 1
c
c        Computes the incomplete gamma integral for positive
c        parameters y,p using an infinite series
c
c        Auxiliary function required: ALNGAM = CACM algorithm 291
c
c	 AS239 should be considered as an alternative to AS147
c
      implicit double precision (a-h,o-z)
      data e/1.0d-9/, zero/0.0d0/, one/1.0d0/, uflo/1.0d-37/
c
c        Checks admissibility of arguments and value of f
c
      ifault = 1
      gammds = zero
      if(y.le.zero .or. p.le.zero) return
      ifault = 2
c
c        alngam is natural log of gamma function
c
      arg = p*log(y)-alngam(p+one)-y
      if(arg.lt.log(uflo)) return
      f = exp(arg)
      if(f.eq.zero) return
      ifault = 0
c
c          Series begins
c
      c = one
      gammds = one
      a = p
    1 a = a+one
      c = c*y/a
      gammds = gammds+c
      if (c/gammds.gt.e) goto 1
      gammds = gammds*f
      return
      end
