      double precision function rngpi(t, n, ifault)
c
c        Algorithm AS 126  Appl. Statist. (1978) Vol. 27, No. 2
c
c        Computes the probability of the normal range given t, the
c        upper limit of integration, and n, the sample size.
c
c     Auxiliary function required: ALNORM = algorithm AS66
c
      implicit double precision (a-h,o-z)
      dimension g(8), h(8)
c
      data g(1), g(2), g(3), g(4), g(5), g(6), g(7), g(8)
     */0.4947004675d0, 0.4722875115d0, 0.4328156012d0, 0.3777022042d0,
     * 0.3089381222d0, 0.2290083888d0, 0.1408017754d0, 0.04750625492d0/
c
      data h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8)
     */0.01357622971d0, 0.03112676197d0, 0.04757925584d0,
     * 0.06231448563d0,
     * 0.07479799441d0, 0.08457825969d0, 0.09130170752d0,
     * 0.09472530523d0/
c
      data zero, half, two, eight /0.0d0, 0.5d0, 2.0d0, 8.0d0/
c
      risf(x) = 0.3989422804d0 * exp(-half*x*x) * (alnorm(x, .false.)
     *  - alnorm(x-t, .false.)) ** (n-1)
c
      ifault = 0
      rngpi = zero
      if (t .le. zero .or. n .le. 1) return

      ifault = 0
      xl = half * t
      a = half * (eight + xl)
      b = eight - xl
      y = zero
      do 10 i = 1, 8
      c = b * g(i)
      y = y + h(i) * (risf(a + c) + risf(a - c))
   10 continue
      rngpi = (two * (alnorm(xl, .false.) - half)) ** n +
     *  two * b * y * n
      return
      end
