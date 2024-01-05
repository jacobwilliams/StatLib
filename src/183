      real function random()
c
c     Algorithm AS 183 Appl. Statist. (1982) vol.31, no.2
c
c     Returns a pseudo-random number rectangularly distributed
c     between 0 and 1.   The cycle length is 6.95E+12 (See page 123
c     of Applied Statistics (1984) vol.33), not as claimed in the
c     original article.
c
c     IX, IY and IZ should be set to integer values between 1 and
c     30000 before the first entry.
c
c     Integer arithmetic up to 30323 is required.
c
      integer ix, iy, iz
      common /randc/ ix, iy, iz
c
      ix = 171 * mod(ix, 177) - 2 * (ix / 177)
      iy = 172 * mod(iy, 176) - 35 * (iy / 176)
      iz = 170 * mod(iz, 178) - 63 * (iz / 178)
c
      if (ix .lt. 0) ix = ix + 30269
      if (iy .lt. 0) iy = iy + 30307
      if (iz .lt. 0) iz = iz + 30323
c
c     If integer arithmetic up to 5212632 is available, the preceding
c     6 statements may be replaced by:
c
c     ix = mod(171 * ix, 30269)
c     iy = mod(172 * iy, 30307)
c     iz = mod(170 * iz, 30323)
c
      random = mod(float(ix) / 30269. + float(iy) / 30307. +
     +                        float(iz) / 30323., 1.0)
      return
      end
c
c
c
c
	real function uniform()
c
c	Generate uniformly distributed random numbers using the 32-bit
c	generator from figure 3 of:
c	L'Ecuyer, P. Efficient and portable combined random number
c	generators, C.A.C.M., vol. 31, 742-749 & 774-?, June 1988.
c
c	The cycle length is claimed to be 2.30584E+18
c
c	Seeds can be set by calling the routine set_uniform
c
c	It is assumed that the Fortran compiler supports long variable
c	names, and integer*4.
c
	integer*4 z, k, s1, s2
	common /unif_seeds/ s1, s2
	save /unif_seeds/
c
	k = s1 / 53668
	s1 = 40014 * (s1 - k * 53668) - k * 12211
	if (s1 .lt. 0) s1 = s1 + 2147483563
c
	k = s2 / 52774
	s2 = 40692 * (s2 - k * 52774) - k * 3791
	if (s2 .lt. 0) s2 = s2 + 2147483399
c
	z = s1 - s2
	if (z .lt. 1) z = z + 2147483562
c
	uniform = z / 2147483563.
	return
	end


	subroutine set_uniform(seed1, seed2)
c
c	Set seeds for the uniform random number generator.
c
	integer*4 s1, s2, seed1, seed2
	common /unif_seeds/ s1, s2
	save /unif_seeds/

	s1 = seed1
	s2 = seed2
	return
	end
