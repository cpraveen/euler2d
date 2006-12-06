C Main program of vertex-centroid scheme of Jameson
      program main
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 esue(3,ntmax), vedge(2,nemax), spts(nspmax),
     +                 bdedge(2,nbpmax), esubp(mesubp,nbpmax),
     +                 ptype(npmax)
      double precision coord(2,npmax), qc(nvar,ntmax), 
     +                 qcold(nvar,ntmax), dt(ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax),
     +                 drmin(ntmax), res(nvar,ntmax)
      double precision cl, cd

      integer          i, j, ie, v1, v2, v3, e1, e2, c1, c2, irk

      call math
      call read_input
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, esubp, coord, drmin, carea, af)

C Set initial condition
      call initialize(qc, cl, cd)

      iter = 0
      fres = 1.0d20
      call system('rm -f FLO.RES')
      open(unit=99, file='FLO.RES')
      do while(iter .lt. MAXITER .and. fres .gt. MINRES)
         call time_step(drmin, qc, dt)
c        call time_step2(edge, tedge, carea, coord, qc, dt)
         call save_old(qc, qcold)

         do irk=1,nirk

            do i=1,nt
               do j=1,nvar
                  res(j,i) = 0.0d0
               enddo
            enddo

C Compute area averaged value at vertices
            do i=1,np
               do j=1,nvar
                  qv(j,i) = 0.0d0
               enddo
            enddo
            do i=1,nt
               v1 = elem(1,i)
               v2 = elem(2,i)
               v3 = elem(3,i)
               call vaverage(coord(1,v1), coord(1,v2), coord(1,v3),
     +                       af(1,v1), af(1,v2), af(1,v3), carea(i), 
     +                       qc(1,i), qv(1,v1), qv(1,v2), qv(1,v3))
            enddo
            do i=1,np
               do j=1,nvar
                  qv(j,i) = qv(j,i)/af(3,i)
               enddo
            enddo
            do i=1,nsp
               j = spts(i)
               e1= bdedge(1,i)
               e2= bdedge(2,i)
               v1= edge(1,e1)
               v2= j
               v3= edge(2,e2)
               call killnormalvel(coord(1,v1), coord(1,v2), coord(1,v3),
     +                            qv(1,j))
            enddo

C Edge count index
            ie = 0

C Compute flux for interior edges
            do i=1,nin
               ie = ie + 1
               e1 = edge(1,ie)
               e2 = edge(2,ie)
               c1 = tedge(1,ie)
               c2 = tedge(2,ie)
               v1 = vedge(1,ie)
               v2 = vedge(2,ie)
               call kfvs_flux(coord(1,e1), coord(1,e2),
     +                       qc(1,c1), qc(1,c2), qv(1,v1), qv(1,v2), 
     +                       res(1,c1), res(1,c2))
            enddo

C Compute flux for solid wall edges
            do i=1,nsw
               ie = ie + 1
               e1 = edge(1,ie)
               e2 = edge(2,ie)
               c1 = tedge(1,ie)
               call solid_flux(coord(1,e1), coord(1,e2),
     +                         qc(1,c1), res(1,c1))       
            enddo

C Flux for far-field points
            do i=1,nff
               ie = ie + 1
               e1 = edge(1,ie)
               e2 = edge(2,ie)
               c1 = tedge(1,ie)
               call farfield_flux(coord(1,e1), coord(1,e2), qc(1,c1),
     +                            cl, cd, res(1,c1))
            enddo

C Update the solution
            if(explicit .eq. yes)then
               do i=1,nt
                  do j=1,nvar
                     qc(j,i) = airk(irk)*qcold(j,i) + 
     +                       birk(irk)*(qc(j,i)-dt(i)*res(j,i)/carea(i))
                  enddo
               enddo
            else
               call lusgs(elem, esue, coord, qcold, qc, res, dt, carea)
            endif

         enddo

         iter = iter + 1
         call residue(res)
         call clcd(edge, tedge, coord, qc, cl, cd)
         write(99,'(i6,4e16.6)') iter, fres, fresi, cl, cd
         if(mod(iter,saveinterval) .eq. 0)then
            call write_result(coord, elem, edge, qc, qv, cl, cd)
         endif

      enddo
      close(99)

      call write_result(coord, elem, edge, qc, qv, cl, cd)

      stop
      end
