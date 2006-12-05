c Roe flux function
      subroutine hcusp_flux(x1, x2, qcl, qcr, qvl, qvr, resl, resr) 
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar),
     +                 qvr(nvar), resl(nvar), resr(nvar)

      integer          i
      double precision rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr,
     &                 ua, va, qa2, aa2, aa, ha,
     &                 ql2, qr2, rl12, rr12, rd,
     &                 unl, unr, una,
     &                 ct, st, nl, Fl(4), Fr(4), Fp(4), flux,
     &                 gg, l1, l2, lneg, lpos, mach,
     &                 dl, dr, li, LIMIT,
     &                 ql(4), qr(4), dqh(4), aqh(4), m0, m1,
     &                 alpha, beta, a1, a2, a3, b1, b2, b3
      intrinsic        dmax1

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      nl = dsqrt(ct**2 + st**2)

      do i=1,4
         dl    = qcl(i) - qvl(i)
         dr    = qvr(i) - qcr(i)
         li    = LIMIT(dl, dr)
         ql(i) = qcl(i) + 0.5d0*li
         qr(i) = qcr(i) - 0.5d0*li
      enddo

C     Left state
      rl = ql(1)
      ul = ql(2)/rl
      vl = ql(3)/rl
      ql2= ul**2 + vl**2
      pl = gamma1*( ql(4) - 0.5d0*rl*ql2 )
      al2= GAMMA*pl/rl
      hl = al2/GAMMA1 + 0.5d0*ql2

C     Right state
      rr = qr(1)
      ur = qr(2)/rr
      vr = qr(3)/rr
      qr2= ur**2 + vr**2
      pr = gamma1*( qr(4) - 0.5d0*rr*qr2 )
      ar2= GAMMA*pr/rr
      hr = ar2/GAMMA1 + 0.5d0*qr2

C     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

C     Centered flux
      Fl(1) = rl*unl
      Fl(2) = pl*ct + rl*ul*unl
      Fl(3) = pl*st + rl*vl*unl
      Fl(4) = rl*hl*unl

      Fr(1) = rr*unr
      Fr(2) = pr*ct + rr*ur*unr
      Fr(3) = pr*st + rr*vr*unr
      Fr(4) = rr*hr*unr

      Fp(1) = 0.0d0
      Fp(2) = (pr - pl)*ct
      Fp(3) = (pr - pl)*st
      Fp(4) = 0.0d0

      dqh(1) = qr(1) - ql(1)
      dqh(2) = qr(2) - ql(2)
      dqh(3) = qr(3) - ql(3)
      dqh(4) = rr*hr - rl*hl

      aqh(1) = 0.5d0*( qr(1) + ql(1) )
      aqh(2) = 0.5d0*( qr(2) + ql(2) )
      aqh(3) = 0.5d0*( qr(3) + ql(3) )
      aqh(4) = 0.5d0*( rr*hr + rl*hl )

C     Roe average
      rl12 = dsqrt(rl)
      rr12 = dsqrt(rr)
      rd   = 1.0d0/(rl12 + rr12)

      ua   = (ul*rl12 + ur*rr12)*rd
      va   = (vl*rl12 + vr*rr12)*rd
      ha   = (hl*rl12 + hr*rr12)*rd
      qa2  = ua**2 + va**2
      aa2  = GAMMA1*(ha - 0.5d0*qa2)

#ifdef DEBUG
      if(aa2 .le. 0.0d0)then
            print*,'Sonic speed is negative'
            print*,qcl(1),qcl(2),qcl(3),qcl(4)
            print*,qcr(1),qcr(2),qcr(3),qcr(4)
            print*
            print*,qvl(1),qvl(2),qvl(3),qvl(4)
            print*,qvr(1),qvr(2),qvr(3),qvr(4)
            print*
            print*,rl,ul,vl,pl
            print*,rr,ur,vr,pr
            print*,li
            stop
      endif
#endif

      aa  = dsqrt(aa2)
      una = ua*ct + va*st
      mach= una/aa/nl

C     Eigenvalues with entropy fix
      gg   = 0.5d0*(gamma + 1.0d0)/gamma
      l1   = (gg*una)**2 + ( (aa*nl)**2 - una**2 )/gamma
      l2   = dsqrt(l1)
      lpos = gg*una + l2
      lneg = gg*una - l2

      m0      = 0.01d0
      m1      = dabs(mach)
      if(m1 .gt. m0)then
         alpha= m1
      else
         alpha= 0.5d0*( m0 + mach**2/m0 )
      endif

      if(mach .le. -1.0d0)then
         beta = -1.0d0
      elseif(mach .gt. -1.0d0 .and. mach .le. 0.0d0)then
         a1   =  una + lpos
         a2   =  una - lpos
         a3   =  a1/a2
         beta = -dmax1(0.0d0, a3)
      elseif(mach .gt.  0.0d0 .and. mach .le. 1.0d0)then
         a1   =  una + lneg
         a2   =  una - lneg
         a3   =  a1/a2
         beta =  dmax1(0.0d0, a3)
      else
         beta =  1.0d0
      endif
         
c     Total flux
      b1 = 0.5d0*alpha*aa*nl
      b2 = 0.5d0*beta*(unr - unl)
      b3 = 0.5d0*beta
      do i=1,4
         flux    = 0.5d0*(Fl(i) + Fr(i)) - b1*dqh(i) - b2*aqh(i) -
     +             b3*Fp(i)
         resl(i) = resl(i) + flux
         resr(i) = resr(i) - flux
      enddo

      return
      end
