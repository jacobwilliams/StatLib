      double precision function trigam(x, ifault)
      implicit double precision (a-h,o-z) 
c
c        algorithm as121   Appl. Statist. (1978) vol 27, no. 1
c
c        calculates trigamma(x) = d**2(log(gamma(x))) / dx**2
c
      double precision a, b, one, half, b2, b4, b6,b8, x, y, z, zero
      data a, b, one, half /1.0d-4, 5.0d0, 1.0d0, 0.5d0/
      data zero /0.0d0/
c
c        b2, b4, b6 and b8 are Bernoulli numbers
c
      data b2, b4, b6,b8
     */0.1666666667d0, -0.03333333333d0, 0.02380952381, -0.03333333333/
c
c        check for positive value of x
c
      trigam = zero
      ifault = 1
      if (x.le.zero) return
      ifault = 0
      z = x
c
c        use small value approximation if x .le. a
c
      if (z .gt. a) goto 10
      trigam = one / (z * z)
      return
c
c        increase argument to (x+i) .ge. b
c
   10 if (z .ge. b) goto 20
      trigam = trigam + one / (z * z)
      z = z + one
      goto 10
c
c        apply asymptotic formula if argument .ge. b
c
   20 y = one / (z * z)
      trigam = trigam + half * y +
     * (one + y * (b2 + y * (b4 + y * (b6 + y * b8)))) / z
      return
      end
                                   
