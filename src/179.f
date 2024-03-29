! ##### User Supplied Comment and Possible Correction #####
!   I think I've found a bug in routine 193 permut. There is an
! array overflow caused statement label "7", which currently reads

!     7 do 8 j = 1,t

! t is equal to n at the end of the algorithm, therefore j 
! equals n.  However, the a(k,n) array is referenced in the following lines
! as a(i,j+1) which is outside the array bounds.

! The correct line should read:

!     7 do 8 j = 1,t-1

! I've tested the change and the program appears to behave the same.  The
! program appeared to work previously but this relied on IFAULT being zero
! and immediately after A in memory.

! cheers,
! Mark
! Mark Lakata <lakata@sseos.lbl.gov>
! ##### 

      subroutine permut(n,k,r,a,ifault)
c     
c     Algorithm AS 179   Appl. Statist. (1982) Vol.31, No.2
c
c     A single call generates all possible permutations of N
c     objects, partitioned into K non-empty subsets, and passed
c     as array R.
c
      integer a(k,n),r(k),t,t2,tj1
c
c     check for fault conditions
c
      ifault = 1
      do 1 i = 1,k
        if (r(i).eq.0) return
    1 continue
      ifault = 2
      isum = 0
      do 2 i = 1,k
    2 isum = isum + r(i)
      if (isum.ne.n) return
      ifault = 0
c
c     step 1.  To initialise array A, fill rows, 1 to k-1, with
c              paired marks.
c
      kount = 1
      l = n
      k1 = k-1
      do 5 i = 1,k1
        ir = r(i)
        do 3 j = 1,ir
    3   a(i,j) = i
        ir1 = r(i) + 1
        do 4 j = ir1,l
    4   a(i,j) = i+1
        l = l - r(i)
    5 continue
c
c     step 2.  Begin generation of permutations by setting looping
c              indices i and t.  i indicates a row in array a, j
c              indicates a column in array a, and t indicates the
c              number of marks in row i.
c              If this is the first permutation, go to step 5.
c
    6 i = k - 1
      t = r(k) + r(i)
      if (kount.eq.1) goto 14
c
c     step 3.  Search for an i mark followed by an i+1 mark in
c              row i and interchange the two marks.  If an 
c              interchange is made, left-shift any marks out of
c              sequence, otherwise go to step 4.
c
    7 do 8 j = 1,t
        j1 = j + 1
        if (a(i,j).ne.i.or.a(i,j1).ne.i+1) goto 8
        limit = j - 2
        a(i,j) = a(i,j1)
        a(i,j1)= a(i,j) - 1
        if (limit) 14, 14, 12
    8 continue
c
c     step 4.  Interchange not made, so return row i to its initial
c              configuration in step 1 and go to step 1.
c
      if (t.eq.1) goto 10
      t2 = t/2
      do 9 j = 1,t2
        tj1 = t - j + 1
        itemp = a(i,j)
        a(i,j) = a(i,tj1)
        a(i,tj1) = itemp
    9 continue
c
c     Reset looping indices i and t.  Return if i = 0, otherwise
c     go to step 3.
c
   10 i = i -1
      if (i.gt.0) goto 11
      return
   11 t = t + r(i)
      goto 7
c
c     Interchange made, so left-shift any marks out of sequence
c     and go to step 5.
c
   12 iflag = 0
      do 13 j = 1,limit
        j1 = j + 1
        if (a(i,j).ne.i + 1.or.a(i,j1).ne.i) goto 13
        a(i,j) = a(i,j1)
        a(i,j1)= a(i,j) + 1
        iflag = 1
   13 continue
      if (iflag.eq.1) goto 12
c
c     step 5.  Fill last row (k) of array a with marks from row 1
c              and generate current permutation in last row (k)
c              of array a.
c
   14 do 15 j = 1,n
   15 a(k,j) = a(1,j)
c
      if (k.eq.2) goto 18
      do 17 l = 2,k1
        m = 1
        do 16 j = 1,n
          if (a(k,j).ne.l) goto 16
          a(k,j) = a(l,m)
          m = m + 1
   16   continue
   17 continue
c
c     step 6.  At this point, call subroutine job to process the
c              current permutation, or execute equivalent
c              statements, and go to step 2.
c
   18 call job(n,k,a,kount)
c
      kount = kount + 1
      goto 6
      end
      
