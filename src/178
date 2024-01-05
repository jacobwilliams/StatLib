      subroutine gsweep (s, t, k, l, m, n, e, ifault)
c
c        Algorithm AS 178   Appl. Statist.  (1982)  Vol. 31, No. 2
c
c        Performs Gauss-Jordan pivot for row/col K in N by N working
c        array stored lower triangle only row by row, shortest row
c        first, in locations 1 to (N*N+N)/2 of T
c
      double precision a, b, e, s(n), t(m), zero, one
      data zero/0.0d0/, one/1.0d0/
c
      ifault=1
      if(n.lt.1 .or. m.lt.(n*n+n)/2) return
      if(k.lt.1 .or. k.gt.n) return
      if(e.lt.zero) return
      ifault=0
c
c        Parameters in range so test for collinearity
c
      l=k
      kk=(k*k+k)/2
      if(t(kk).lt.zero) goto 20
      if(t(kk).lt.e*s(k)) return
      ii=0
      ik=kk-k
      do 10 l=1,n
        ii=ii+l
        ik=ik+1
        if (l.gt.k) ik=ik+l-2
        if (t(ii).ge.zero) goto 10
        if(one/(t(ik)*t(ik)/t(kk)-t(ii)).lt.e*s(l)) return
   10 continue
c
c        No collinearity so update triangle
c
   20 l=0
      t(kK)=-one/t(kk)
      a=abs(t(kk))
      ik=kk-k
      ij=0
      do 90 i=1,n
        ik=ik+1
        if(i-k) 50,30,40
   30   ij=ij+k
        goto 90
   40   ik=ik+i-2
   50   b=t(ik)
        if(t(kk).lt.zero) b=-b
        t(ik)=a*t(ik)
        jk=kk-k
        do 80 j=1,i
          ij=ij+1
          jk=jk+1
          if(j-k) 70,80,60
   60     jk=jk+j-2
   70     t(ij)=t(ij)+b*t(jk)
   80   continue
   90 continue
      return
      end
