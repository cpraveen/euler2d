c Roe flux function
      subroutine roe_flux(x1, x2, qcl, qcr, qvl, qvr, resl, resr) 
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar),
     +                 qvr(nvar), resl(nvar), resr(nvar)

      integer          i
      double precision rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr,
     &                 ua, va, qa2, aa2, aa, ha,
     &                 ql2, qr2, rl12, rr12, rd,
     &                 unl, unr, una, vna, ct, st, Fc(4), Fd(4),
     &                 m1, m2, a1, a2, a3, a4, l1, l2, l3, l4,
     &                 a1l1, a2l2, a3l3, a4l4, aact, aast,
     &                 du1, du2, du3, du4, flux, dl, dr, li, limit,
     &                 ql(4), qr(4), lent, uabs, mabs, ETOL
      parameter(ETOL=0.01d0)
      intrinsic        dmin1

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)
      ct = ct/lent
      st = st/lent

      do i=1,4
         dl    = qcl(i) - qvl(i)
         dr    = qvr(i) - qcr(i)
c        li    = LIMIT(dl, dr)
         ql(i) = qcl(i) + 0.5d0*dl
         qr(i) = qcr(i) - 0.5d0*dr
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
      Fc(1) = rl*unl            + rr*unr
      Fc(2) = pl*ct + rl*ul*unl + pr*ct + rr*ur*unr
      Fc(3) = pl*st + rl*vl*unl + pr*st + rr*vr*unr
      Fc(4) = rl*hl*unl         + rr*hr*unr

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
      vna =-ua*st + va*ct

C     Eigenvalues with entropy fix
      l1 = dabs(una - aa)
      uabs = dabs(una)
      mabs = uabs/aa
      if(mabs .gt. ETOL)then
         l2 = uabs
         l3 = uabs
      else
         l2 = 0.5d0*( ETOL + mabs**2/ETOL)*aa
         l3 = l2
      endif

      l4 = dabs(una + aa)

c     Difference of conserved variables
      du1 = qr(1) - ql(1)
      du2 = qr(2) - ql(2)
      du3 = qr(3) - ql(3)
      du4 = qr(4) - ql(4)

c     Amplitudes
      m1 = (ct*du2 + st*du3 - una*du1)/aa
      m2 = GAMMA1*(du4 - ua*du2 - va*du3 + qa2*du1)/aa**2

      a4 = 0.5d0*(m1 + m2)
      a1 = 0.5d0*(m2 - m1)
      a3 = du1 - a1 - a4
      a2 = ( st*du2 - ct*du3 + vna*du1 )/aa

c     Diffusive flux
      a1l1  = a1*l1
      a2l2  = a2*l2
      a3l3  = a3*l3
      a4l4  = a4*l4
      aact  = aa*ct
      aast  = aa*st

      Fd(1) = a1l1               +               a3l3           + a4l4
      Fd(2) = a1l1*(ua - aact)   + a2l2*aa*st  + a3l3*ua        +
     &        a4l4*(ua + aact)
      Fd(3) = a1l1*(va - aast)   - a2l2*aa*ct  + a3l3*va        +
     &        a4l4*(va + aast)
      Fd(4) = a1l1*(ha - una*aa) + a2l2*aa*vna + a3l3*0.5d0*qa2 +
     &        a4l4*(ha + una*aa)

c     Total flux
      do i=1,4
         flux    = 0.5d0*lent*( Fc(i) - Fd(i) )
         resl(i) = resl(i) + flux
         resr(i) = resr(i) - flux
      enddo

      return
      end
