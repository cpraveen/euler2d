

C-----------------------------------------------------------------------------
C Flux for a farfield edge
C-----------------------------------------------------------------------------
      subroutine farfield_flux(x1, x2, qc, cl, cd, res)
      implicit none
      include 'common.h'
      include 'inf.h'
      double precision x1(2), x2(2), qc(nvar), cl, cd, res(nvar)

      double precision dx, dy, dr, nx, ny, r, u, v, p, q2, a,
     +                 un, ut, un_inf, ut_inf, l1, l2, R1, R2, S, 
     +                 unf, utf, af, rf, uf, vf, pf, ef
      double precision dref, theta, circ, fact1, fact2, fact3, fact,
     +                 uinf, vinf, qinf, pinf1, pinf2, pinf, rinf, ainf

      if(mach_inf .lt. 1.0d0 .and. vortex .eq. yes)then
         dx    = 0.5d0*( x1(1) + x2(1) ) - xref
         dy    = 0.5d0*( x1(2) + x2(2) ) - yref
         dref  = dsqrt(dx**2 + dy**2)
         theta = datan2(dy, dx)
         circ  = 0.5d0*q_inf*Cl
         fact1 = circ*dsqrt(1.0d0 - mach_inf**2)
         fact2 = 2.0d0*M_PI*dref
         fact3 = 1.0d0 - (mach_inf*dsin(theta - aoa))**2
         fact  = fact1/(fact2*fact3)
         uinf  = u_inf + fact*dsin(theta)
         vinf  = v_inf - fact*dcos(theta)
         qinf  = dsqrt(uinf**2 + vinf**2)
         pinf1 = p_inf**(GAMMA1/GAMMA)
         pinf2 = 0.5d0*(GAMMA1/GAMMA)*r_inf*(q_inf**2 -
     &           qinf**2)/p_inf**(1.0d0/GAMMA)
         pinf  = (pinf1 + pinf2)**(GAMMA/GAMMA1)
         rinf  = r_inf*(pinf/p_inf)**(1.0d0/GAMMA)
         ainf  = dsqrt(GAMMA*pinf/rinf)
      else
         uinf  = u_inf
         vinf  = v_inf
         pinf  = p_inf
         rinf  = r_inf
         ainf  = a_inf
      endif

      dx =  x2(1) - x1(1)
      dy =  x2(2) - x1(2)
      dr =  dsqrt(dx**2 + dy**2)
      nx =  dy/dr
      ny = -dx/dr

      r = qc(1)
      u = qc(2)/r
      v = qc(3)/r
      q2= u**2 + v**2
      p = gamma1*( qc(4) - 0.5d0*r*q2 )
      a = dsqrt(gamma*p/r)

      un     =  u*nx     + v*ny
      ut     = -u*ny     + v*nx
      un_inf =  uinf*nx + vinf*ny
      ut_inf = -uinf*ny + vinf*nx

      l1 = un - a
      l2 = un + a

c First Riemann invariant
      if( l1 .gt. 0.0d0)then
         R1 = 0.5d0*un     - a/gamma1
      else
         R1 = 0.5d0*un_inf - ainf/gamma1
      endif

c Second Riemann invariant
      if( l2 .gt. 0.0d0)then
         R2 = 0.5d0*un     + a/gamma1
      else
         R2 = 0.5d0*un_inf + ainf/gamma1
      endif

c Entropy and tangential velocity
      if( un .gt. 0.0d0)then
         S   = p/r**gamma
         utf = ut
      else
         S   = pinf/rinf**gamma
         utf = ut_inf
      endif

      unf = R1 + R2
      af  = 0.5d0*(R2 - R1)*GAMMA1








      rf  = (af**2/S/gamma)**(1.0d0/gamma1)
      pf  = S*rf**gamma
      uf  = unf*nx - utf*ny
      vf  = unf*ny + utf*nx
      ef  = pf/gamma1 + 0.5d0*rf*(uf**2 + vf**2)

      res(1) = res(1) + dr*(rf*unf)
      res(2) = res(2) + dr*(pf*nx + rf*uf*unf)
      res(3) = res(3) + dr*(pf*ny + rf*vf*unf)
      res(4) = res(4) + dr*((ef + pf)*unf)

      return
      end
