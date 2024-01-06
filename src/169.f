      subroutine scatpl(a, n, m, icy, ncy, icx, ny, nx, scaley,
     *  scalex, istand, ifault, iwrite)
c
c        algorithm as 169  appl. statist. (1981) vol.30, no.3
c
c        produces a scatter plot of one variable against several
c
      character*1 iout,intch,markch,iblank,idot,icolon
      character*1 icomma,iapost,isemi,itwo,idash,form21*14,form22*6
      character iform1*20,iform2*21,form11*8,form12*11
      dimension iout(161), vals(20), a(n, m), icy(ncy), scalex(2),
     *  scaley(2), intch(11), markch(5)
      data maxwid /132/, maxht /50/, maxy /5/
      data mpvx /10/, mpvy /5/
      data intch(1), intch(2), intch(3), intch(4), intch(5), intch(6),
     *  intch(7), intch(8), intch(9), intch(10), intch(11)
     *  /'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '9'/
      data markch(1), markch(2), markch(3), markch(4), markch(5)
     *  /'*', '0', '+', 'x', '='/
      data iblank /' '/, idot /'.'/, icolon /':'/, icomma /','/,
     *  iapost /''''/, isemi /';'/, itwo /'2'/, idash /'-'/
c
c        formats for printing
c
      data iform1 /'(1h ,f8.0,1x,152a1) '/
      data iform2 /'(1h ,5x,16(f8.0,2x)) '/
c
    1 format(' ', 10x, ':', 151a1)
    2 format(' times 10**', i3)
    3 format(' offset', f10.0)
    4 format(' ', 14x, 'times 10**', i3)
    5 format(' ', 14x, 'offset', f10.0)
    6 format(' ', 2x, 16(9x, a1))
c
c        test for valid parameters
c
      ifault = 0
      if (n .lt. 1) ifault = ifault + 1
      if (m .lt. 2) ifault = ifault + 2
      if (icx .lt. 1 .or. icx .gt. m) ifault = ifault + 4
      if (ncy .le. 0 .or. ncy .gt. maxy) ifault = ifault + 8
      if (ifault .gt. 0) return
      do 10 i = 1, ncy
        if (icy(i) .lt. 1 .or. icy(i) .gt. m) goto 1016
   10 continue
c
c        set plot size
c
      nly = maxht - 5
      if (nly .gt. ny) nly = ny
      if (nly .le. mpvy) nly = mpvy + 1
      nlx = maxwid - 11
      if (nlx .gt. nx) nlx = nx
      if (nlx .le. mpvx) nlx = mpvx + 1
c
c        calculate horizontal scale
c
      xmin = scalex(1)
      xmax = scalex(2)
      if (xmax .gt. xmin) goto 30
      xmin = a(1, icx)
      xmax = xmin
      if (n .eq. 1) goto 30
      do 20 i = 2, n
        ai = a(i, icx)
        if (ai .lt. xmin) xmin = ai
        if (ai .gt. xmax) xmax = ai
   20 continue
c
   30 call scale(xmin, xmax, nlx, mpvx, temp, xvstep,
     *  nxvals, irx, ifail)
      if (ifail .gt. 0) goto 1032
      xmin = temp
      xstep = xvstep / float(mpvx)
c
c        calculate vertical scale
c
      ymin = scaley(1)
      ymax = scaley(2)
      if (ymax .gt. ymin) goto 50
      k = icy(1)
      ymin = a(1, k)
      ymax = ymin
      do 40 j = 1, ncy
        k = icy(j)
        do 40 i = 1, n
          ai = a(i, k)
          if (ai .lt. ymin) ymin = ai
          if (ai .gt. ymax) ymax = ai
   40 continue
c
   50 call scale(ymin, ymax, nly, mpvy, temp, yvstep,
     *  nyvals, iry, ifail)
      if (ifail .gt. 0) goto 1064
      ymin = temp
      ystep = yvstep / float(mpvy)
c
c        calculate printed values and set mark for 2 points
c
      call axis(ymin, yvstep, nyvals, 6, iry, irpr, offset,
     *  ifact, vals, 20, ifail)
      if (ifail .gt. 0) goto 1064
      read(iform1,'(a8,1x,a11)')form11,form12
      write(iform1,'(a8,a1,a11)')form11,intch(irpr + 1),form12
      if (ifact .ne. 0) write (iwrite, 2) ifact
      if (offset .ne. 0.0)  write (iwrite, 3) offset
      if (istand .eq. 0) intch(3) = isemi
c
c        for each line of output
c
      iplted = 0
      do 140 i = 1, nly
        iy = nly - i
        do 60 ix = 1, nlx
   60   iout(ix) = iblank
c
c        scan data for points on current line
c
        do 120 l = 1, n
          indx = (a(l, icx) - xmin) / xstep + 1.5
          if (indx .lt. 1 .or. indx .gt. nlx) goto 120
          do 110 j = 1, ncy
            k = icy(j)
            y = (a(l, k) - ymin) / ystep
            indy = y + 0.5
            if (indy .ne. iy) goto 110
            iplted = iplted + 1
            if (iout(indx) .ne. iblank) goto 80
c
c        single point
c
            if (istand .eq. 0) goto 70
            iout(indx) = markch(j)
            goto 110
   70       iout(indx) = icomma
            if (int(y) .eq. iy) iout(indx) = iapost
            goto 110
c
c        multiple points
c
   80       do 90 ic = 3, 10
              if (iout(indx) .eq. intch(ic)) goto 100
   90       continue
            ic = 2
  100       iout(indx) = intch(ic + 1)
  110     continue
  120   continue
c
c        print line
c
        if (mod(iy, mpvy) .eq. 0) goto 130
        write (iwrite, 1) (iout(ix), ix = 1, nlx)
        goto 140
  130   write (iwrite, iform1) vals(nyvals), idash, icolon,
     *       (iout(ix), ix = 1, nlx)
        nyvals = nyvals - 1
  140 continue
c
c        print horizontal axis using variable formats
c
      write (iwrite, 1) (idot, i = 1, nlx)
      call axis(xmin, xvstep, nxvals, 6, irx, irpr, offset,
     *  ifact, vals, 20, ifail)
      intch(3) = itwo
      read(iform2,'(a14,1x,a6)')form21,form22
      write(iform2,'(a14,a1,a6)')form21,intch(irpr + 1),form22
      if (ifail .gt. 0) goto 1032
      write (iwrite, 6) (icolon, i = 1, nxvals)
      write (iwrite, iform2) (vals(i), i = 1, nxvals)
      if (ifact .ne. 0) write (iwrite, 4) ifact
      if (offset .ne. 0.0) write (iwrite, 5) offset
      ifault = iplted - n * ncy
      return
c
c        set error indicator
c
 1064 ifault = ifault + 32
 1032 ifault = ifault + 16
 1016 ifault = ifault + 16
      return
      end
