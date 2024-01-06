      subroutine allnr(n, r, j, ifault)
c
c        Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3
c
c        When called once, generates all possible combinations
c        from a group of N items.  Each combination (represented in j as
c        r ordered integers between 1 and n) is processed within allnr.
c
c        Parameters:-
c       
c        n        integer             input:  The size of the group from which
c                                             the combinations are selected.
c
c        r        integer             input:  The size of each comination.
c
c        j        integer array(r)  workspace: Used by allnr to store
c                                              combinations.
c
c        ifault   integer            output:  Fault indicator, equal to:
c                                             0 if 1 le R le N;
c                                             1 otherwise.
c
      integer r, j(r)
c
      ifault = 1
      if (r .lt.1 .or. r .gt. n) return
      ifault = 0
      kount = 0
      nmr = n - r
c
c        Initialize J(1) to lower limit separately, since lower limit for
c        each index depends on lower limit for previous index
c
      i = 1
      j(1) = 1
c
c        Initialize indices for loops i=1,...,r to lower limits
c
    1 if (i .eq. r) goto 3
      ip1 = i + 1
      do 2 l = ip1, r
    2 j(l) = j(l - 1) + 1
c
c        Update the count (kount) of combinations and process the current
c        combination.  The call to Subroutine job may be replaced by
c        statements to process the current combination.
c
    3 kount = kount + 1
      call job(n, r, j, kount)
c
c        Increment the first possible index (of loop i) among indices of
c        loops R, R-1,...,1
c
      i = r
    4 if (j(i) .lt. nmr + i) goto 5
      i = i - 1
c
c        Return after all indices have achieved their upper limits
c
      if (i .le. 0) return
      goto 4
    5 j(i) = j(i) + 1
      goto 1
      end
c
c
c
      subroutine job (n, r, j, kount)
      integer r, j(r)
      return
      end
