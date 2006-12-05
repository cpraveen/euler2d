      subroutine lusgs(elem, esue, coord, qp, qcold, qc, res, dt, 
     +                 tarea)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), esue(3,ntmax)
      double precision coord(2,npmax), qp(nvar,ntmax),
     +                 qcold(nvar,ntmax), qc(nvar,ntmax), 
     +                 res(nvar,ntmax), dt(ntmax), tarea(ntmax)

      integer          it, iv, p1, p2, p3, t1, t2, t3
      double precision u, v, q2, a, lam, sx1, sx2, sx3, sy1, sy2, sy3,
     +                 ds1, ds2, ds3, D(nt), dqc(nvar,nt),
     +                 cres(nvar), flux(nvar),
     +                 lam1, lam2, lam3, omega

c Over-relaxation factor; value > one increases stability but reduces
c convergence rate
      omega = 1.5d0

c Forward loop
      do it=1,nt
         p1  = elem(1,it)
         p2  = elem(2,it)
         p3  = elem(3,it)

         u   = qp(2,it)/qp(1,it)
         v   = qp(3,it)/qp(1,it)
         q2  = u*u + v*v
         a   = GAMMA*GAMMA1*( qp(4,it) - 0.5d0*qp(1,it)*q2)/qp(1,it)
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
            call jacobprod(sx1, sy1, qp(1,it), dqc(1,t1), flux)
            call maxeig(sx1, sy1, ds1, qp(1,it), lam1)
            do iv=1,nvar
               cres(iv) = cres(iv) - flux(iv) - omega*lam1*dqc(iv,t1)
            enddo
         endif

         if(t2 .lt. it .and. t2 .ne. 0)then
            call jacobprod(sx2, sy2, qp(1,it), dqc(1,t2), flux)
            call maxeig(sx2, sy2, ds2, qp(1,it), lam2)
            do iv=1,nvar
               cres(iv) = cres(iv) - flux(iv) - omega*lam2*dqc(iv,t2)
            enddo
         endif

         if(t3 .lt. it .and. t3 .ne. 0)then
            call jacobprod(sx3, sy3, qp(1,it), dqc(1,t3), flux)
            call maxeig(sx3, sy3, ds3, qp(1,it), lam3)
            do iv=1,nvar
               cres(iv) = cres(iv) - flux(iv) - omega*lam3*dqc(iv,t3)
            enddo
         endif

         do iv=1,nvar
            dqc(iv,it) = ( -res(iv,it) - 0.5d0*cres(iv) )/D(it)
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
            call jacobprod(sx1, sy1, qp(1,it), dqc(1,t1), flux)
            call maxeig(sx1, sy1, ds1, qp(1,it), lam1)
            do iv=1,nvar
               cres(iv) = cres(iv) - flux(iv) - omega*lam1*dqc(iv,t1)
            enddo
         endif

         if(t2 .gt. it .and. t2 .ne. 0)then
            call jacobprod(sx2, sy2, qp(1,it), dqc(1,t2), flux)
            call maxeig(sx2, sy2, ds2, qp(1,it), lam2)
            do iv=1,nvar
               cres(iv) = cres(iv) - flux(iv) - omega*lam2*dqc(iv,t2)
            enddo
         endif

         if(t3 .gt. it .and. t3 .ne. 0)then
            call jacobprod(sx3, sy3, qp(1,it), dqc(1,t3), flux)
            call maxeig(sx3, sy3, ds3, qp(1,it), lam3)
            do iv=1,nvar
               cres(iv) = cres(iv) - flux(iv) - omega*lam3*dqc(iv,t3)
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
      subroutine jacobprod(sx, sy, qp, dq, flux)
      implicit none
      include 'param.h'
      double precision sx, sy, qp(nvar), dq(nvar), flux(nvar)

      integer          i, j
      double precision d, u, v, q2, p, un, h, An(nvar,nvar)

      call contoprim(qp, d, u, v, p)
      q2      = u*u + v*v
      h       = GAMMA*p/d/GAMMA1 + 0.5d0*q2
      un      = u*sx + v*sy

C Normal jacobian A*sx + B*sy
      An(1,1) = 0.0d0
      An(2,1) = -u*un + 0.5d0*GAMMA1*q2*sx
      An(3,1) = -v*un + 0.5d0*GAMMA1*q2*sy
      An(4,1) = -(h - 0.5d0*GAMMA1*q2)*un

      An(1,2) = sx
      An(2,2) = (3.0d0-GAMMA)*u*sx + v*sy
      An(3,2) = v*sx - GAMMA1*u*sy
      An(4,2) = h*sx - GAMMA1*u*un

      An(1,3) = sy
      An(2,3) = -GAMMA1*v*sx + u*sy
      An(3,3) = u*sx + (3.0d0-GAMMA)*v*sy
      An(4,3) = h*sy - GAMMA1*v*un

      An(1,4) = 0.0d0
      An(2,4) = GAMMA1*sx
      An(3,4) = GAMMA1*sy
      An(4,4) = GAMMA*un

      do i=1,nvar
         flux(i) = 0.0d0
         do j=1,nvar
            flux(i) = flux(i) + An(j,i)*dq(j)
         enddo
      enddo

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
