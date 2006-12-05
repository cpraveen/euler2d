      subroutine solid_flux(x1, x2, qc, res)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qc(nvar), res(nvar)

      double precision r, u, v, p, q2, fx, fy

C     Left state
      r = qc(1)
      u = qc(2)/r
      v = qc(3)/r
      q2= u**2 + v**2
      p = gamma1*( qc(4) - 0.5d0*r*q2 )

      fx=  p*(x2(2) - x1(2))
      fy= -p*(x2(1) - x1(1))

      res(2) = res(2) + fx
      res(3) = res(3) + fy

      return
      end
