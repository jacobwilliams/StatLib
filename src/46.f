      subroutine gmsmod (x,n,m,v,c,d,nullty)
c
c        Algorithm AS 46 Applied Statistics (J.R.Statist.Soc C),
c       (1971) Vol.20, No.3
c       
c       Expresses the n*m matrix x as the product of the n*m matrix v
c       and the m*m upper unit triangular matrix c.  v-transpose.v is
c       a diagonal matrix and is returned in d.  nullty is the nullity
c       of d.  v may coincide with x.
c
        double precision x(*),v(*),c(*),d(*)
        double precision w0, w1, w2, w3, eta
        double precision test, zero, one
c
      data eta/1.0d-12/, zero/0.0d0/, one/1.0d0/
c
      nullty=0
      ic=0
      i1=-n
      do 70 icol=1,m
        i1=i1+n
        w0=zero
        i2=i1
        do 10 irow=1,n
          i2=i2+1
          w1=x(i2)
          v(i2)=w1
   10   w0=w0+w1*w1
        i3=0
        i4=0
        do 70 icolw=1,icol
          i2=i1
          w2=zero
          do 20 irow=1,n
            i2=i2+1
            i3=i3+1
            w3=v(i2)
            w1=v(i3)
   20     w2=w2+w1*w3
          if(icolw.lt.icol) goto 30
          test=w2/w0
          if(test.ge.eta) goto 25
          do 21 irow=1,n
            i4=i4+1
   21     v(i4)=zero
          i4=i4-n
          nullty=nullty+1
          w2=zero
   25     d(icol)=w2
          i4=i4+n
          w1=one
          goto 60
   30     w1=zero
          w3=d(icolw)
          if(w3.eq.zero) goto 40
          w1=w2/w3
   40     i2=i1
          do 50 irow=1,n
            i2=i2+1
            i4=i4+1
   50     v(i2)=v(i2)-w1*v(i4)
   60     ic=ic+1
          c(ic)=w1
   70 continue
      return
      end
