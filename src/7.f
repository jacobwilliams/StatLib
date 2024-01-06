c This file contains AS7 and an enhanced alternative - ASR44.  See also AS6.
c
c
        subroutine syminv(a, n, nn, c, w, nullty, ifault)
c
c       Algorithm AS7, Applied Statistics, vol.17, 1968, p.198.
c
c       Forms in c( ) as lower triangle, a generalised inverse
c       of the positive semi-definite symmetric matrix a( )
c       order n, stored as lower triangle.
c
c       arguments:-
c       a()     = input, the symmetric matrix to be inverted, stored in
c                 lower triangular form
c       n       = input, order of the matrix
c       nn      = input, the size of the a and c arrays     n*(n+1)/2
c       c()     = output, the inverse of a (a generalized inverse if c is
c                 singular), also stored in lower triangular.
c                 c and a may occupy the same locations.
c       w()     = workspace, dimension at least n.
c       nullty  = output, the rank deficiency of a.
c       ifault  = output, error indicator
c                       = 1 if n < 1
c                       = 2 if a is not +ve semi-definite
c                       = 3 if nn < n*(n+1)/2
c                       = 0 otherwise
c
c***************************************************************************
c
        double precision a(nn), c(nn), w(n), x, zero, one
c
        data zero, one /0.0d0, 1.0d0/
c
c       cholesky factorization of a, result in c
c
        call chol(a, n, nn, c, nullty, ifault)
        if(ifault.ne.0) return
c
c       invert c & form the product (cinv)'*cinv, where cinv is the inverse
c       of c, row by row starting with the last row.
c       irow = the row number, ndiag = location of last element in the row.
c
        irow=n
        ndiag=nn
   10   l=ndiag
        if (c(ndiag) .eq. zero) goto 60
        do 20 i=irow,n
          w(i)=c(l)
          l=l+i
   20   continue
        icol=n
        jcol=nn
        mdiag=nn
   30   l=jcol
        x=zero
        if(icol.eq.irow) x=one/w(irow)
        k=n
   40   if(k.eq.irow) go to 50
        x=x-w(k)*c(l)
        k=k-1
        l=l-1
        if(l.gt.mdiag) l=l-k+1
        go to 40
   50   c(l)=x/w(irow)
        if(icol.eq.irow) go to 80
        mdiag=mdiag-icol
        icol=icol-1
        jcol=jcol-1
        go to 30
   60   do 70 j=irow,n
          c(l)=zero
          l=l+j
   70   continue
   80   ndiag=ndiag-irow
        irow=irow-1
        if(irow.ne.0) go to 10
        return
        end
c
c
c
c
        subroutine subinv(a, nm, b, n, c, w, nullty, ifault, ndim, det)
        implicit double precision (a-h,o-z)
c
c       Remark asr 44.1   Appl. Statist. (1982) Vol.31, No.3
c
c       A revised and enhanced version of
c       algorithm as7, applied statistics, vol.17, 1968.
c
c       Forms in c() as lower triangle, a generalised inverse
c       of a sub-matrix of the positive semi-definite symmetric 
c       matrix a() of order n, stored as lower triangle.
c       c() may co-incide with a().  nullty is returned as the nullity
c       of a().  ifault is returned as 1 if n.lt.1, 2 if the
c       submatrix is not positive semi-definite, 3 if the elements
c       of b are inadmissible, otherwise zero.
c       w() is a work array of length at least n that is allocated by
c       the calling routine.
c
c       arguments:-
c       a()     = input, the symmetric matrix to be inverted, stored in
c                 lower triangular form
c       nm      = input, the order of a
c       n       = input, order of the sub-matrix to be inverted
c       b       = input, an array containing the numbers of rows and
c                 columns of a() to be included.
c       c()     = output, the inverse of a (a generalized inverse if c is
c                 singular), also stored in lower triangular.
c                 c and a may occupy the same locations.
c       w()     = workspace, dimension at least n.
c       nullty  = output, the rank deficiency of a.
c       ifault  = output, error indicator
c                       = 1 if n < 1
c                       = 2 if the submatrix of a is not +ve semi-definite
c                       = 3 if elements of b are inadmissible
c                       = 0 otherwise
c
c     Auxiliary routine required: SUBCHL from ASR44 (see AS6)
c
c***************************************************************************
c
        dimension a(ndim),c(ndim),w(n)
        integer b(n)
        data zero/0.0d0/, one/1.0d0/
c
        ifault = 3
        if (n .gt. nm) return
        if (b(1) .lt.1 .or. b(1) .gt. nm-n+1) return
        if (n .eq. 1) goto 19
        do 18 i = 2, n
          if (b(i) .le. b(i-1) .or. b(i) .gt. nm-n+i) return
   18   continue
   19   nrow=n
        ifault=1
        if(nrow.le.0) go to 100
        ifault=0
c
c       Cholesky factorization of a, result in c
c
	call subchl(a,b,nrow,c,nullty,ifault,ndim,det)
        if(ifault.ne.0) go to 100
c
c       invert c & form the product (cinv)'*cinv, where cinv is the inverse
c       of c, row by row starting with the last row.
c       irow = the row number, ndiag = location of last element in the row.
c
        nn=nrow*(nrow+1)/2
        irow=nrow
        ndiag=nn
   10   if(c(ndiag).eq.zero) go to 60
        l=ndiag
        do 20 i=irow,nrow
          w(i)=c(l)
          l=l+i
   20   continue
        icol=nrow
        jcol=nn
        mdiag=nn
   30   l=jcol
        x=zero
        if(icol.eq.irow) x=one/w(irow)
        k=nrow
   40   if(k.eq.irow) go to 50
        x=x-w(k)*c(l)
        k=k-1
        l=l-1
        if(l.gt.mdiag) l=l-k+1
        go to 40
   50   c(l)=x/w(irow)
        if(icol.eq.irow) go to 80
        mdiag=mdiag-icol
        icol=icol-1
        jcol=jcol-1
        go to 30
   60   l=ndiag
        do 70 j=irow,nrow
          c(l)=zero
          l=l+j
   70   continue
   80   ndiag=ndiag-irow
        irow=irow-1
        if(irow.ne.0) go to 10
  100   return
        end
