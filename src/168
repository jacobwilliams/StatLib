      subroutine scale(fmn, fmx, n, mpv, valmin, step, nvals,
     *  ir, ifault)
c
c        algorithm as 168  appl. statist. (1981) vol.30, no.3
c
c        calculates neat values for lowest print position and
c        step between printed values on a scale to include fmn
c        and fmx
c
      real unit(12)
      real tol, bias
      data nunit /12/
      data unit(1), unit(2), unit(3), unit(4), unit(5), unit(6),
     *   unit(7), unit(8), unit(9), unit(10), unit(11), unit(12)
     *  /12.0, 15.0, 20.0, 25.0, 30.0, 40.0,
     *  50.0, 60.0, 80.0, 100.0, 120.0, 150.0/
      data tol /5.0e-6/, bias /1.0e-5/
      data minn /2/, maxn /10000/, cover /0.7/
c
      fmax = fmx
      fmin = fmn
c
c        test for valid parameters
c
      ifault = 0
      if (fmax .lt. fmin) ifault = ifault + 1
      if (n .lt. minn .or. n .gt. maxn) ifault = ifault + 2
      if (mpv .le. 0 .or. mpv .ge. n) ifault = ifault + 4
      if (ifault .ne. 0) return
      nvals = (n - 1) / mpv + 1
c
c        test for values effectively equal
c
      if (fmax - fmin .gt. tol * amax1(abs(fmax), abs(fmin))) goto 20
      ifault = -1
      if (fmax) 5, 10, 15
    5 fmax = 0.0
      goto 20
   10 fmax = 1.0
      goto 20
   15 fmin = 0.0
c        find neat trial step size
c
   20 finter = float(n) / float(mpv)
      s = (fmax - fmin) * (1.0 + 2.0 * bias) / finter
      ir = 0
   25 if (s .gt. 10.0) goto 30
      s = s * 10.0
      ir = ir + 1
      goto 25
   30 if (s .le. 100.0) goto 35
      s = s / 10.0
      ir = ir - 1
      goto 30
   35 do 40 i = 1, nunit
      if (s .le. unit(i)) goto 45
   40 continue
   45 step = 10.0 ** (-ir) * unit(i)
c
c        find neat trial start value
c
      aj = 0.0
   50 aj = aj + 1.0
      if (unit(i) - 0.1 .gt. aint((unit(i) + 0.1) / aj) * aj) goto 50
      tstep = step / aj
      temp = fmin / tstep + aj * (0.5 / float(mpv) - finter * bias)
      valmin = aint(temp) * tstep
      if (temp .lt. 0.0 .and. temp .ne. aint(temp))
     * valmin = valmin - tstep
c
c        test whether fmax is in scale
c
      if (fmax .lt. valmin + step *
     *  (finter * (1.0 - bias) - 0.5 / float(mpv))) goto 55
c
c        try new step or start value
c
      if (unit(i) / unit(i + 1) * (1.0 - 1.0 / (aj * finter)) .lt.
     *  cover) goto 50
      i = i + 1
      goto 45
c
c        get position of least significant figure on scale
c
   55 do 60 j = 1, 2
        aj = aj * 10.0
        if (unit(i) - 0.1 .lt. aint((unit(i) + 0.1) / aj) * aj)
     *          ir = ir - 1
   60 continue
      return
      end
c
      subroutine axis(valmin, step, nvals, maxpr, ir, irprin, offset,
     *  ifact, vals, iv, ifault)
c
c        algorithm as 168.1 appl. statist. (1981) vol.30, no.3
c
c        sets up values and formats for printing on an axis
c
      real vals(iv)
      data irmax /20/, mprmax /20/
c
c        check for valid parameters
c
      ifault = 0
      if (nvals .lt. 2) ifault = ifault + 1
      fmax = valmin + step * float(nvals - 1)
      if (nvals .ge. 2 .and. fmax .le. valmin) ifault = ifault + 2
      if (maxpr .lt. 2 .or. maxpr .gt. mprmax) ifault = ifault + 4
      if (nvals .gt. iv) ifault = ifault + 8
      if (ir .gt. irmax) ifault = ifault + 16
      if (ifault .gt. 0) return
c        find position of most significant digit(il) and number of
c        significant digits overall(is) and varying(it)
c
      tmax = 10.0 ** maxpr
      fl = abs(fmax)
      fs = abs(valmin)
      il = 0
   10 if (fl .lt. 1.0 .and. fs .lt. 1.0) goto 20
      fl = fl / 10.0
      fs = fs / 10.0
      il = il + 1
      goto 10
c
   20 if (fl .ge. 0.1 .or. fs .ge. 0.1) goto 30
      fl = fl * 10.0
      fs = fs * 10.0
      il = il - 1
      goto 20
c
   30 is = il + ir
      it = is
      if (valmin .le. 0.0 .and. fmax .ge. 0.0) goto 50
   40 fl = mod(fl, 1.0) * 10.0
      fs = mod(fs, 1.0) * 10.0
      if (it .le. 0) goto 1016
      if (int(fl) .ne. int(fs)) goto 50
      it = it - 1
      goto 40
c
c        decide on printing format
c
   50 ifact = 0
      offset = 0.0
      irprin = max(ir, 0)
      ilprin = max(il, 0)
c
      if (irprin + ilprin .le. maxpr) goto 70
      if (is .le. maxpr) goto 60
      irprin = maxpr - 1
      ifact = max(it, maxpr) - 1 - ir
      goto 70
   60 ifact = il - 1
      irprin = is - 1
   70 fs = 10.0 ** (-ifact)
      vstep = step * fs
      vmin = valmin * fs
      if (is .le. maxpr) goto 80
      offset = aint(vmin / 10.0) * 10.0
      vmin = vmin - offset
c
c        write values for axis
c
   80 do 90 i = 1, nvals
        vals(i) = vmin
        vmin = vmin + vstep
   90 continue
c
c        check that all values can be printed
c
      fs = 0.1 ** irprin
      if (abs(vals(1)) * fs + 0.5 .lt. tmax .and.
     *  abs(vals(nvals)) * fs + 0.5 .lt. tmax) return
      il = il + 1
      is = is + 1
      it = it + 1
      goto 50
c
c        error indicator
c
 1016 ifault = 16
      return
      end
