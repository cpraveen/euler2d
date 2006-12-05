C-----------------------------------------------------------------------------
C Flux for a farfield edge: Steger-Warming
C-----------------------------------------------------------------------------
      subroutine farfield2_flux(x1, x2, qc, cl, cd, res)
      implicit none
      include 'common.h'
      include 'inf.h'
      double precision x1(2), x2(2), qc(nvar), cl, cd, res(nvar)

      double precision dx, dy, dr, nx, ny, r, u, v, p, q2, a,
     +                 un, l2, l3
      double precision dref, theta, circ, fact1, fact2, fact3, fact4,
     +                 fact, uinf, vinf, pinf, rinf, ainf,
     +                 q2inf, l1p, l2p, l3p, l1n, l2n, l3n, fp(4), 
     +                 fn(4), f1, f2, a1, a2

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
         q2inf = uinf**2 + vinf**2
         fact4 = 1.0d0 + 0.5d0*GAMMA1*(q_inf**2 - q2inf)/a_inf**2
         rinf  = r_inf*fact4**(1.0d0/GAMMA1)
         pinf  = p_inf*(rinf/r_inf)**GAMMA
         ainf  = dsqrt(GAMMA*pinf/rinf)
      else
         uinf  = u_inf
         vinf  = v_inf
         q2inf = uinf**2 + vinf**2
         pinf  = p_inf
         rinf  = r_inf
         ainf  = a_inf
      endif

      dx =  x2(1) - x1(1)
      dy =  x2(2) - x1(2)
      dr =  dsqrt(dx**2 + dy**2)
      nx =  dy/dr
      ny = -dx/dr

c Positive flux
      r = qc(1)
      u = qc(2)/r
      v = qc(3)/r
      q2= u**2 + v**2
      p = gamma1*( qc(4) - 0.5d0*r*q2 )
      a = dsqrt(gamma*p/r)

      un = u*nx + v*ny
      l2 = un + a
      l3 = un - a
      l1p= dmax1(un, 0.0d0)
      l2p= dmax1(l2, 0.0d0)
      l3p= dmax1(l3, 0.0d0)
      a1 = 2.0d0*(GAMMA-1.0d0)*l1p + l2p + l3p
      f1 = 0.5d0*r/GAMMA

      fp(1) = f1*a1
      fp(2) = f1*( a1*u + a*(l2p - l3p)*nx )
      fp(3) = f1*( a1*v + a*(l2p - l3p)*ny )
      fp(4) = f1*( 0.5d0*a1*q2 + a*un*(l2p - l3p) +
     +             a**2*(l2p+l3p)/GAMMA1 )

c Negative flux
      un = uinf*nx + vinf*ny
      l2 = un + ainf
      l3 = un - ainf
      l1n= dmin1(un, 0.0d0)
      l2n= dmin1(l2, 0.0d0)
      l3n= dmin1(l3, 0.0d0)
      a2 = 2.0d0*(GAMMA-1.0d0)*l1n + l2n + l3n
      f2 = 0.5d0*rinf/GAMMA

      fn(1) = f2*a2
      fn(2) = f2*( a2*uinf + ainf*(l2n - l3n)*nx )
      fn(3) = f2*( a2*vinf + ainf*(l2n - l3n)*ny )
      fn(4) = f2*( 0.5d0*a2*q2inf + ainf*un*(l2n - l3n) +
     +             ainf**2*(l2n+l3n)/GAMMA1 )

      res(1) = res(1) + dr*( fp(1) + fn(1) )
      res(2) = res(2) + dr*( fp(2) + fn(2) )
      res(3) = res(3) + dr*( fp(3) + fn(3) )
      res(4) = res(4) + dr*( fp(4) + fn(4) )

      return
      end
