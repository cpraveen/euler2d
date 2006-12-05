C Main program of vertex-centroid scheme of Jameson
      program main
      implicit none
      include 'param.h'
      include 'param2.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 esue(3,ntmax), vedge(2,nemax), spts(nspmax),
     +                 bdedge(2,nbpmax), esubp(mesubp,nbpmax),
     +                 ptype(npmax)
      double precision coord(2,npmax), qc(nvar,ntmax), 
     +                 dt(ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax),
     +                 drmin(ntmax), res1(nvar), res2(nvar)
      double precision cl, cd
      double precision resb(nvar,ntmax), resbold(nvar,ntmax),
     +                 qcb(nvar,ntmax), qcb1(nvar,ntmax), 
     +                 qvb(nvar,npmax), coordb(2,npmax)
      double precision cost

      integer          i, j, ie, v1, v2, v3, e1, e2, c1, c2, irk

      call math
      call read_input
      call read_input_adj
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, esubp, coord, drmin, carea, af)

C Read the flow solution from file
      istart = restart
      call initialize(qc, cl, cd)
      call time_step(drmin, qc, dt)

C Compute averaged value at vertices
      print*,'Computing vertex average values ...'
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
     +                 af(1,v1), af(1,v2), af(1,v3), carea(i), 
     +                 qc(1,i), qv(1,v1), qv(1,v2), qv(1,v3))
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
     +                      qv(1,j))
      enddo

C Compute lift and drag coefficients
      call clcd(edge, tedge, coord, qc, cl, cd)

C Cost func derivative wrt flow solution. This is stored in qcb1 and is
C used in adjoint solution. It is the rhs of adjoint solution, dJ/dq. For
C boundary cost function most of the entries of qcb1 are zero. Can we do
C something to save memory ?
      print*,'Computing derivative of cost function wrt flow ...'
      call cost_q(spts, elem, edge, tedge, bdedge, coord, carea,
     +                  af, qc, qcb1, qv, qvb, cost)

C FV residual, no need to calculate, just set it to zero
      do i=1,nvar
         res1(i) = 0.0d0
         res2(i) = 0.0d0
      enddo

C Initialize adjoint variable
      call init_adj(resb)

C Start of adjoint iterations
      print*,'Beginning of adjoint iterations ...'
      iter = 0
      fres = 1.0d20
      call system('rm -f ADJ.RES')
      open(unit=99, file='ADJ.RES')
      do while(iter .lt. MAXITER .and. fres .gt. MINRES)
         call save_old(resb, resbold)

         do irk=1,nirk

            do i=1,nt
               do j=1,nvar
                  qcb(j,i) = 0.0d0
               enddo
            enddo
            do i=1,np
               do j=1,nvar
                  qvb(j,i) = 0.0d0
               enddo
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
               call kfvs_flux_bq(coord(1,e1), coord(1,e2),
     +                       qc(1,c1), qcb(1,c1), qc(1,c2), qcb(1,c2),
     +                       qv(1,v1), qvb(1,v1), qv(1,v2), qvb(1,v2),
     +                       res1, resb(1,c1), res2, resb(1,c2))
            enddo

C Compute flux for solid wall edges
            do i=1,nsw
               ie = ie + 1
               e1 = edge(1,ie)
               e2 = edge(2,ie)
               c1 = tedge(1,ie)
               call solid_flux_bq(coord(1,e1), coord(1,e2),
     +                            qc(1,c1), qcb(1,c1),
     +                            res1, resb(1,c1))       
            enddo

C Flux for far-field points
            do i=1,nff
               ie = ie + 1
               e1 = edge(1,ie)
               e2 = edge(2,ie)
               c1 = tedge(1,ie)
               call farfield_flux_bq(coord(1,e1), coord(1,e2), 
     +                               qc(1,c1), qcb(1,c1), cl, cd, 
     +                               res1, resb(1,c1))
            enddo

C Contribution from vertex averaging
            do i=1,nsp
               j = spts(i)
               e1= bdedge(1,i)
               e2= bdedge(2,i)
               v1= edge(1,e1)
               v2= j
               v3= edge(2,e2)
               call killnormalvel_bq(coord(1,v1), coord(1,v2), 
     +                               coord(1,v3), qv(1,j), qvb(1,j))
            enddo
            do i=1,np
               do j=1,nvar
                  qvb(j,i) = qvb(j,i)/af(3,i)
               enddo
            enddo
            do i=1,nt
               v1 = elem(1,i)
               v2 = elem(2,i)
               v3 = elem(3,i)
               call vaverage_bq(coord(1,v1), coord(1,v2), coord(1,v3),
     +                         af(1,v1), af(1,v2), af(1,v3), carea(i), 
     +                         qc(1,i), qcb(1,i), qv(1,v1), qvb(1,v1),
     +                         qv(1,v2), qvb(1,v2), qv(1,v3), qvb(1,v3))
            enddo

C Adjoint residual
            do i=1,nt
               do j=1,nvar
                  qcb(j,i)  = qcb(j,i) + qcb1(j,i)
               enddo
            enddo

C Update the solution
            if(explicit2 .eq. yes)then 
               do i=1,nt
                  do j=1,nvar
                     resb(j,i) = airk(irk)*resbold(j,i) + 
     +                           birk(irk)*(resb(j,i) - 
     +                           dt(i)*qcb(j,i)/carea(i))
                  enddo
               enddo
            else
               call lusgs(elem, esue, coord, qc, resbold, resb, qcb, dt,
     +                    carea)
            endif

         enddo

         iter = iter + 1
         call residue(qcb)
         write(99,'(i6,2e16.6)') iter, fres, fresi
         if(mod(iter,saveinterval) .eq. 0)then
            call cost_x(elem, edge, tedge, vedge, spts, bdedge, esubp,
     +                  ptype, coord, qc, qv, qvb, carea, af, resb,
     +                  coordb, cl, cd)
            call result(coord, spts, elem, edge, bdedge, carea, af, 
     +                  resb, qvb, coordb, cl, cd)
         endif

      enddo
      close(99)

      call cost_x(elem, edge, tedge, vedge, spts, bdedge, esubp,
     +            ptype, coord, qc, qv, qvb, carea, af, resb,
     +            coordb, cl, cd)

      call result(coord, spts, elem, edge, bdedge, carea, af, resb, qvb,
     +            coordb, cl, cd)

      stop
      end
