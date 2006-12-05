      subroutine lusgs(elem, esue, coord, qcold, qc, res, dt, 
     +                 tarea)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), esue(3,ntmax)
      double precision coord(2,npmax), qcold(nvar,ntmax), 
     +                 qc(nvar,ntmax), res(nvar,ntmax),
     +                 dt(ntmax), tarea(ntmax)

      integer          it, iv, p1, p2, p3, t1, t2, t3
      double precision u, v, q2, a, lam, sx1, sx2, sx3, sy1, sy2, sy3,
     +                 ds1, ds2, ds3, D(nt), dqc(nvar,nt),
     +                 cres(nvar), flux1(nvar), flux2(nvar),
     +                 lam1, lam2, lam3, omega

      omega = 1.5d0

c Forward loop
      do it=1,nt
         p1  = elem(1,it)
         p2  = elem(2,it)
         p3  = elem(3,it)

         u   = qc(2,it)/qc(1,it)
         v   = qc(3,it)/qc(1,it)
         q2  = u*u + v*v
         a   = GAMMA*GAMMA1*( qc(4,it) - 0.5d0*qc(1,it)*q2)/qc(1,it)
         a   = dsqrt(a)

         lam =  0.0d0

         sx1 =  coord(2,p2) - coord(2,p1)
         sy1 = -(coord(1,p2) - coord(1,p1))
         ds1 =  dsqrt(sx1*sx1 + sy1*sy1)
         lam =  lam + dabs(u*sx1 + v*sy1) + a*ds1

         sx2 =  coord(2,p3) - coord(2,p2)
         sy2 = -( coord(1,p3) - coord(1,p2) )
         ds2 =  dsqrt(sx2*sx2 + sy2*sy2)
         lam =  lam + dabs(u*sx2 + v*sy2) + a*ds2

         sx3 =  coord(2,p1) - coord(2,p3)
         sy3 = -( coord(1,p1) - coord(1,p3) )
         ds3 =  dsqrt(sx3*sx3 + sy3*sy3)
         lam =  lam + dabs(u*sx3 + v*sy3) + a*ds3

         D(it) = tarea(it)/dt(it) + 0.5d0*lam*omega

         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         do iv=1,nvar
            cres(iv) = 0.0d0
         enddo

         if(t1 .lt. it .and. t1 .ne. 0)then
            call normalflux(sx1, sy1, qcold(1,t1),  flux1)
            call normalflux(sx1, sy1, qc(1,t1), flux2)
            call maxeig(sx1, sy1, ds1, qcold(1,t1), lam1)
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - 
     +                    omega*lam1*dqc(iv,t1)
            enddo
         endif

         if(t2 .lt. it .and. t2 .ne. 0)then
            call normalflux(sx2, sy2, qcold(1,t2),  flux1)
            call normalflux(sx2, sy2, qc(1,t2), flux2)
            call maxeig(sx2, sy2, ds2, qcold(1,t2), lam2)
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - 
     +                    omega*lam2*dqc(iv,t2)
            enddo
         endif

         if(t3 .lt. it .and. t3 .ne. 0)then
            call normalflux(sx3, sy3, qcold(1,t3),  flux1)
            call normalflux(sx3, sy3, qc(1,t3), flux2)
            call maxeig(sx3, sy3, ds3, qcold(1,t3), lam3)
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - 
     +                    omega*lam3*dqc(iv,t3)
            enddo
         endif

         do iv=1,nvar
            dqc(iv,it) = ( -res(iv,it) - 0.5d0*cres(iv) )/D(it)
            qc(iv,it)  = qcold(iv,it) + dqc(iv,it)
         enddo

      enddo

c Reverse loop
      do it=nt,1,-1
         p1  = elem(1,it)
         p2  = elem(2,it)
         p3  = elem(3,it)

         sx1 =   coord(2,p2) - coord(2,p1)
         sy1 = -(coord(1,p2) - coord(1,p1))
         ds1 =  dsqrt(sx1*sx1 + sy1*sy1)

         sx2 =   coord(2,p3) - coord(2,p2)
         sy2 = -(coord(1,p3) - coord(1,p2))
         ds2 =  dsqrt(sx2*sx2 + sy2*sy2)

         sx3 =   coord(2,p1) - coord(2,p3)
         sy3 = -(coord(1,p1) - coord(1,p3))
         ds3 =  dsqrt(sx3*sx3 + sy3*sy3)

         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         do iv=1,nvar
            cres(iv) = 0.0d0
         enddo

         if(t1 .gt. it .and. t1 .ne. 0)then
            call normalflux(sx1, sy1, qcold(1,t1),  flux1)
            call normalflux(sx1, sy1, qc(1,t1), flux2)
            call maxeig(sx1, sy1, ds1, qcold(1,t1), lam1)
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - 
     +                    omega*lam1*dqc(iv,t1)
            enddo
         endif

         if(t2 .gt. it .and. t2 .ne. 0)then
            call normalflux(sx2, sy2, qcold(1,t2),  flux1)
            call normalflux(sx2, sy2, qc(1,t2), flux2)
            call maxeig(sx2, sy2, ds2, qcold(1,t2), lam2)
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - 
     +                    omega*lam2*dqc(iv,t2)
            enddo
         endif

         if(t3 .gt. it .and. t3 .ne. 0)then
            call normalflux(sx3, sy3, qcold(1,t3),  flux1)
            call normalflux(sx3, sy3, qc(1,t3), flux2)
            call maxeig(sx3, sy3, ds3, qcold(1,t3), lam3)
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - 
     +                    omega*lam3*dqc(iv,t3)
            enddo
         endif

         do iv=1,nvar
            dqc(iv,it) = ( D(it)*dqc(iv,it) - 0.5d0*cres(iv) )/D(it)
            qc(iv,it)  = qcold(iv,it) + dqc(iv,it)
         enddo

      enddo

      return
      end

C Computes flux along (sx,sy)
      subroutine normalflux(sx, sy, qc, flux)
      implicit none
      include 'param.h'
      double precision sx, sy, qc(nvar), flux(nvar)

      double precision d, u, v, p, un, e

      call contoprim(qc, d, u, v, p)
      e       = p/GAMMA1 + 0.5d0*d*(u*u + v*v)
      un      = u*sx + v*sy
      flux(1) = d*un
      flux(2) = p*sx + d*u*un
      flux(3) = p*sy + d*v*un
      flux(4) = (e + p)*un

      return
      end

c Computes maximum eigenvalue normal to a face with normal (sx, sy)
c and ds = sqrt(sx*sx + sy*sy) is face length
      subroutine maxeig(sx, sy, ds, qc, lam)
      implicit none
      include 'param.h'
      double precision sx, sy, ds, qc(nvar), lam

      double precision d, u, v, p, a

      call contoprim(qc, d, u, v, p)
      a   = dsqrt(GAMMA*p/d)
      lam = dabs(u*sx + v*sy) + a*ds

      return
      end
