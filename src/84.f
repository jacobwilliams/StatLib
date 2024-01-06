      subroutine skewkt(x, ip, n, b1, b2, work, ssq, nn, irank, ifault)
c
c       Algorithm AS 84  Appl. Statist. (1975) Vol.24, No. 2
c
c       For a given observation matrix X(ip, n), the subroutine
c       calculates measures of multivariate skewness and kurtosis
c       and stores them in B1 and B2
c
c       It is important that the value of irank be tested on exit as
c       the null distribution of the measures under multinormality
c       depends on the rank of the sample covariance matrix
c
c
c       Formal Parameters:-
c
c       X           Real Array (IP, N)     Input: Array containing the data,
c                                                 Each observn occupying 
c                                                 one column
c                                          Output: Mean subtracted from data
c       IP          Integer                Input: Number of variables p
c       N           Integer                Input: Number of observations n
c       B1          Real                   Output: Skewness b1(1,p)
c       B2          Real                   Output: Kurtosis b(2,p)
c       WORK        Real Array (IP)        Workspace: Used by SYMINV
c       SSQ         Real Array (NN)        Output: Contains inverse A of
c                                                  corrected ssp matrix n*S
c                                                  stored in the sequence
c                                                  a11,a21,a22,a31,a32,a33,...
c       NN          Integer                Input: NN=IP*(IP+1)/2
c       IRANK       Integer                Output: Rank r of matrix S
c       IFAULT      Integer                Output: Failure Indicator:-
c                                          =0 None
c                                          =1 S is singular, g-inverse used
c                                          =2 Error returned from Syminv
c                                          =3 IP LE 1
c                                          =4 NN NE IP*(IP+1)/2
c                                          =5 N LT IP
c
c
      double precision x(ip,n), ssq(nn), work(ip), sumwts, xm, b1,
     #                 b2, qq, qr, xki, xkj, ssqm, ss, zero, two
      data zero, two /0.0d0, 2.0d0/
c
      ifault = 0
      if (ip .le. 1) ifault = 3
      if (nn .ne. ip * (ip + 1) / 2) ifault = 4
      if (n .lt. ip) ifault = 5
      if (ifault .ne. 0) return
      sumwts = n
      do 32 j = 1, ip
        xm = zero
        do 25 i = 1, n
   25   xm = xm + x(j, i)
        xm = xm / sumwts
c
c       xm at this stage contains the sample mean of the j-th variable,
c       where j is the control variable of the do 32 loop
c
        do 30 i = 1, n
   30   x(j, i) = x(j, i) - xm
   32 continue
      l = 0
      do 40 j = 1, ip
        do 40 k = 1, j
          l = l + 1
          xm = zero
          do 35 i = 1, n
   35     xm = xm + x(j, i) * x(k, i)
          ssq(l) = xm
   40 continue
c
c       ssq at this stage contains the sample matrix of corrected sums
c       of squares and products, stored as lower triangle
c
      call syminv (ssq, ip, nn, ssq, work, irank, ifault)
      if (ifault .ne. 0) return
      if (irank .ne. 0) ifault = 1
      irank = ip - irank
      b1 = zero
      b2 = zero
      do 100 i = 1, n
        xm = ssq(1) * x(1, i) ** 2
        m = 1
        do 70 k = 2, ip
          qq = zero
          kk = k - 1
          do 50 l = 1, kk
            m = m + 1
            qq = qq + ssq(m) * x(l, i)
   50     continue
          m = m + 1
          xki = x(k, i)
          xm = xm + (ssq(m) * xki + two * qq) * xki
   70   continue
        qq = xm ** 2
        b2 = b2 + qq
        b1 = b1 + qq * xm
  100 continue
      do 150 i = 2, n
        ii = i - 1
        ss = ssq(1) * x(1, i)
        do 150 j = 1, ii
          xm = ss * x(1, j)
          m = 1
          do 120 k = 2, ip
            xki = x(k, i)
            xkj = x(k, j)
            qq = zero
            qr = zero
            kk = k - 1
            do 110 l = 1, kk
              m = m + 1
              ssqm = ssq(m)
              qq = qq + ssqm * x(l, j)
              qr = qr + ssqm * x(l, i)
  110       continue
            qq = qq * xki + qr * xkj
            m = m + 1
            xm = xm + qq + ssq(m) * xki * xkj
  120     continue
          b1 = b1 + two * xm ** 3
  150 continue
      b1 = b1 * sumwts
      b2 = b2 * sumwts
      return
      end
        
