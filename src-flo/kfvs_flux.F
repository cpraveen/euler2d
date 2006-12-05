C.....Kinetic split fluxes
      subroutine kfvs_flux(x1, x2, qcl, qcr, qvl, qvr, resl, resr) 
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar),
     +                 qvr(nvar), resl(nvar), resr(nvar)

      integer          i
      double precision rl, ul, vl, pl, el, rr, ur, vr, pr, er,
     &                 ql2, qr2, unl, unr, vnl, vnr,
     &                 lent, ct, st, Al, Bl, Ar, Br,
     &                 sl, betal, sr, betar,
     &                 Fp(4), Fm(4), Ff(4), dl, dr, li, ql(nvar),
     &                 qr(nvar), limit, flux(nvar)

      ct   =  (x2(2) - x1(2))
      st   = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)
      ct   = ct/lent
      st   = st/lent

      do i=1,4
         dl = qcl(i) - qvl(i)
         dr = qvr(i) - qcr(i)
c        li = LIMIT(dl, dr)
         ql(i) = qcl(i) + 0.5d0*dl
         qr(i) = qcr(i) - 0.5d0*dr
      enddo

C     Left state
      rl = ql(1)
      ul = ql(2)/rl
      vl = ql(3)/rl
      ql2= ul**2 + vl**2
      pl = gamma1*( ql(4) - 0.5d0*rl*ql2 )
      el = pl/GAMMA1 + 0.5d0*rl*ql2

C     Right state
      rr = qr(1)
      ur = qr(2)/rr
      vr = qr(3)/rr
      qr2= ur**2 + vr**2
      pr = gamma1*( qr(4) - 0.5d0*rr*qr2 )
      er = pr/GAMMA1 + 0.5d0*rr*qr2

C     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

      vnl =-ul*st + vl*ct
      vnr =-ur*st + vr*ct

c     Positive flux
      betal = 0.5d0*rl/pl
      sl    = unl*dsqrt(betal)
      Al    = 0.5d0*(1.0d0 + DERF(sl))
      Bl    = 0.5d0*dexp(-sl**2)/dsqrt(M_PI*betal)

      Fp(1) = rl*(unl*Al + Bl)
      Fp(2) = (pl + rl*unl**2)*Al + rl*unl*Bl
      Fp(3) = rl*(unl*Al + Bl)*vnl
      Fp(4) = (el + pl)*unl*Al + (el + 0.5d0*pl)*Bl

c     Negative flux
      betar = 0.5d0*rr/pr
      sr    = unr*dsqrt(betar)
      Ar    = 0.5d0*(1.0d0 - DERF(sr))
      Br    = 0.5d0*dexp(-sr**2)/dsqrt(M_PI*betar)

      Fm(1) = rr*(unr*Ar - Br)
      Fm(2) = (pr + rr*unr**2)*Ar - rr*unr*Br
      Fm(3) = rr*(unr*Ar - Br)*vnr
      Fm(4) = (er + pr)*unr*Ar - (er + 0.5d0*pr)*Br

c     Total flux
      do i=1,4
         Ff(i) = Fp(i) + Fm(i)
      enddo

      flux(1) = Ff(1)*lent
      flux(2) = (ct*Ff(2) - st*Ff(3))*lent
      flux(3) = (st*Ff(2) + ct*Ff(3))*lent
      flux(4) = Ff(4)*lent

c     Total flux
      do i=1,4
         resl(i) = resl(i) + flux(i)
         resr(i) = resr(i) - flux(i)
      enddo

      return
      end
