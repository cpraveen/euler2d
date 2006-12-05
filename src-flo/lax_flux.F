c Roe flux function
      subroutine lax_flux(x1, x2, qcl, qcr, qvl, qvr, resl, resr) 
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar),
     +                 qvr(nvar), resl(nvar), resr(nvar)

      integer          i
      double precision rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr,
     &                 ql2, qr2, al, ar, unl, unr, ct, st, Fc(4),
     &                 l1, l2, l3, dq(4), flux, dl, dr, li, limit,
     &                 ql(4), qr(4), lent

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)

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
      al = dsqrt(al2)
      hl = al2/GAMMA1 + 0.5d0*ql2

C     Right state
      rr = qr(1)
      ur = qr(2)/rr
      vr = qr(3)/rr
      qr2= ur**2 + vr**2
      pr = gamma1*( qr(4) - 0.5d0*rr*qr2 )
      ar2= GAMMA*pr/rr
      ar = dsqrt(ar2)
      hr = ar2/GAMMA1 + 0.5d0*qr2

C     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

C     Centered flux
      Fc(1) = 0.5d0*(rl*unl            + rr*unr)
      Fc(2) = 0.5d0*(pl*ct + rl*ul*unl + pr*ct + rr*ur*unr)
      Fc(3) = 0.5d0*(pl*st + rl*vl*unl + pr*st + rr*vr*unr)
      Fc(4) = 0.5d0*(rl*hl*unl         + rl*hr*unr)

C     Eigenvalues with entropy fix
      l1 = dabs(unl) + al*lent
      l2 = dabs(unr) + ar*lent
      l3 = dmax1(l1, l2)

c     Difference of conserved variables
      dq(1) = qr(1) - ql(1)
      dq(2) = qr(2) - ql(2)
      dq(3) = qr(3) - ql(3)
      dq(4) = qr(4) - ql(4)

c     Total flux
      do i=1,4
         flux    = Fc(i) - 0.5d0*l3*dq(i)
         resl(i) = resl(i) + flux
         resr(i) = resr(i) - flux
      enddo

      return
      end
