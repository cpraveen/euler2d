C.....Calculate lift and drag coefficients
      subroutine costfunc(x1, x2, qc, cost)
      implicit none
#if 1
#include<common.h>
#include<inf.h>
#endif
      double precision x1(2), x2(2), qc(nvar), cost

      double precision xf, yf, p, dx, dy, nx, ny, dcl, dcd

      p  = GAMMA1*( qc(4) - 0.5d0*(qc(2)**2 + qc(3)**2)/qc(1) )
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
