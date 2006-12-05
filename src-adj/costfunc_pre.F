C.....Calculate lift and drag coefficients
      subroutine costfunc(x1, x2, qv1, qv2, cost)
      implicit none
#if 1
#include<common.h>
#include<inf.h>
#endif
      double precision x1(2), x2(2), qv1(nvar), qv2(nvar), cost

      double precision xf, yf, p, p1, p2, dx, dy, nx, ny, dcl, dcd

      p1 = GAMMA1*( qv1(4) - 0.5d0*(qv1(2)**2 + qv1(3)**2)/qv1(1) )
      p2 = GAMMA1*( qv2(4) - 0.5d0*(qv2(2)**2 + qv2(3)**2)/qv2(1) )
      p  = 0.5d0*(p1 + p2)
      dx = x2(1) - x1(1)
      dy = x2(2) - x1(2)
      nx = dy
      ny =-dx
      xf = p*nx
      yf = p*ny
c     dc = (-dsin(aoa)*xf + dcos(aoa)*yf)/(0.5d0*r_inf*q_inf**2)
      dcl= yf/(0.5d0*r_inf*q_inf**2)
      dcd= xf/(0.5d0*r_inf*q_inf**2)

      cost = cost + dcl

      return
      end
