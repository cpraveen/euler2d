C.....Cost function is L2 norm of pressure difference
      subroutine costfunc_pre(cp0, qv, cost)
      implicit none
#if 1
#include<common.h>
#include<inf.h>
#endif
      double precision cp0, qv(nvar), cost

      double precision p, cp, dc

      
      p    = GAMMA1*( qv(4) - 0.5d0*(qv(2)**2 + qv(3)**2)/qv(1) )
      cp   = -(p - pre_inf)/(0.5d0*r_inf*q_inf**2)
      dc   = 0.5d0*(cp - cp0)**2
      cost = cost + dc

      return
      end
