      subroutine clustr (x,d,dev,b,f,e,i,j,n,nz,k)
      implicit double precision (a-h,o-z)
c
c     Algorithm AS 58  Appl. Statist. (1973) Vol.22, No.1
c
c     Given a matrix of i observations on j variables, the
c     observations are allocated to n clusters in such a way that the
c     within-cluster sum of squares is minimised.
c
      integer b(i),e(k)
      dimension x(i,j),d(k,j),dev(k),f(i)
      data zero/0.0d0/, one/1.0d0/, two/2.0d0/
      data big/1.0d10/
c
      do 10 ia=1,n
        e(ia)=0
   10 continue
c
c     For each observation, calculate the distance from each cluster
c     centre, and assign to the nearest
c
      do 40 ic=1,i
        f(ic)=zero
        da=big
        do 30 id=1,n
          db=zero
          do 20 ie=1,j
            dc=x(ic,ie)-d(id,ie)
            db=db+dc*dc
            if (db.ge.da) go to 30
   20     continue
          da=db
          b(ic)=id
   30   continue
        ig=b(ic)
        e(ig)=e(ig)+1
   40 continue
c
c     Calculate the mean and sum of squares for each cluster
c
      do 50 ix=1,n
        dev(ix)=zero
        do 50 iy=a1,j
          d(ix,iy)=zero
   50 continue
      do 60 ic=1,i
        ig=b(ic)
        do 60 ih=1,j
          d(ig,ih)=d(ig,ih)+x(ic,ih)
   60 continue
      do 80 ij=1,j
        do 70 ii=1,n
          d(ii,ij)=d(ii,ij)/dble(e(ii))
   70   continue
        do 80 ik=1,i
          il=b(ik)
          da=x(ik,ij)-d(il,ij)
          db=da*da
          f(ik)=f(ik)+db
          dev(il)=dev(il)+db
   80 continue
      do 85 ik=1,i
        il=b(ik)
        fl=e(il)
        if (fl.ge.two) f(ik)=f(ik)*fl/(fl-one)
   85 continue
c
c     Examine each observation in turn to see if it should be
c     reassigned to a different cluster
c
      if (nz.le.0) nz=1
   90 iw=0
      do 140 ik=1,i
      il=b(ik)
      ir=il
c
c     If the number of cluster points is less than or equal to the
c     specified minimum, nz, bypass this iteration
c
      if (e(il).le.nz) go to 140
      fl=e(il)
      dc=f(ik)
      do 100 in=1,n
        if (in.eq.il) goto 100
        fm=e(in)
        fm=fm/(fm+one)
        de=zero
        do 95 ip=1,j
          da=x(ik,ip)-d(in,ip)
          de=de + da *da * fm
          if (de.ge.dc) goto 100
   95 continue
      dc=de
      ir=in
  100 continue
      if (ir.eq.il) goto 140
c
c     Reassignment is made here if necessary
c
      fq=e(ir)
      dev(il)=dev(il)-f(ik)
      dev(ir)=dev(ir)+dc
      e(ir)=e(ir)+1
      e(il)=e(il)-1
      do 110 is=1,j
        d(il,is) = (d(il,is) * fl - x(ik,is)) / (fl-one)
        d(ir,is) = (d(ir,is) * fq + x(ik,is)) / (fq + one)
  110 continue
      b(ik) = ir
      do 130 it=1,i
        ij = b(it)
        if (ij.ne.il. and. ij.ne.ir) goto 130
        f(it)=zero
        do 120 iu=1,j
          da = x(it,iu) - d(ij,iu)
          f(it)=f(it) + da * da
  120   continue
        fl=e(ij)
        f(it) = f(it)*fl/(fl-one)
  130 continue
      iw=iw+1
  140 continue
c
c     Return to calling program if no reassignments were made during
c     this iteration
c
      if (iw.eq.o) return
      goto 90
      end
