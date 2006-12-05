      subroutine killnormalvel(x1, x2, x3, qv)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), x3(2), qv(nvar)

      double precision dx1, dy1, dx2, dy2, ds1, ds2, nx1, ny1, nx2, ny2,
     +                 nxa, nya, nsa, nx, ny, un

      dx1 = x2(1) - x1(1)
      dy1 = x2(2) - x1(2)
      ds1 = dsqrt(dx1*dx1 + dy1*dy1)
      nx1 =-dy1/ds1
      ny1 = dx1/ds1

      dx2 = x3(1) - x2(1)
      dy2 = x3(2) - x2(2)
      ds2 = dsqrt(dx2*dx2 + dy2*dy2)
      nx2 =-dy2/ds2
      ny2 = dx2/ds2

      nxa = 0.5d0*(nx1 + nx2)
      nya = 0.5d0*(ny1 + ny2)
      nsa = dsqrt(nxa**2 + nya**2)

      nx  = nxa/nsa
      ny  = nya/nsa

      un    = qv(2)*nx + qv(3)*ny
      qv(2) = qv(2) - un*nx
      qv(3) = qv(3) - un*ny

      return
      end
