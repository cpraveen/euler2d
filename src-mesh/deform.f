C Given the boundary perturbation, deform the grid using spring model
      program main
      implicit none
      include 'param.h'
      include 'param2.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 esue(3,ntmax), vedge(2,nemax), spts(nspmax),
     +                 bdedge(2,nbpmax), esubp(mesubp,nbpmax),
     +                 ptype(npmax)
      double precision coord(2,npmax), af(3,npmax), carea(ntmax), 
     +                 drmin(ntmax)

      integer          n, ind(nbpmax)
      double precision dx(nbpmax), dy(nbpmax)

      call read_input
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, esubp, coord, drmin, carea, af)

      call bd_dx(n, ind, dx, dy)
      call smooth_grid(n, ind, ptype, edge, dx, dy, coord)
      call write_new_grid(coord, elem, ptype)

      stop
      end

C Read boundary point movement from file
      subroutine bd_dx(n, ind, dx, dy)
      implicit none
      include 'param.h'
      integer          n, ind(nbpmax)
      double precision dx(nbpmax), dy(nbpmax)

      integer          i

      open(20, file='def.dat')
      read(20,*) n
      do i=1,n
         read(20,*) ind(i), dx(i), dy(i)
      enddo
      close(20)

      return
      end

C Propagate shape deformation to interior grid points using spring model
      subroutine smooth_grid(n, ind, ptype, edge, dx, dy, coord)
      implicit none
      include 'param.h'
      include 'param2.h'
      integer          n, ind(nbpmax), ptype(npmax), edge(2,nemax)
      double precision dx(nbpmax), dy(nbpmax), coord(2,npmax)

      integer          i, giter, n1, n2
      double precision residue, residue1, dcoord(2,npmax), 
     +                 gres(2,npmax), eps(npmax), ddx, ddy

      print*,'Deforming the grid ...'

      do i=1,np
         dcoord(1,i) = 0.0d0
         dcoord(2,i) = 0.0d0
      enddo

      residue = 1.0d0
      residue1= 1.0d0
      giter   = 0
      call grid_eps(edge, ptype, coord, eps)
      do while(residue .gt. RTOL)
         giter = giter + 1

         do i=1,np
            gres(1,i) = 0.0d0
            gres(2,i) = 0.0d0
         enddo

         do i=1,ne
            n1 = edge(1,i)
            n2 = edge(2,i)
            call grid_residue(coord(1,n1), coord(1,n2), dcoord(1,n1),
     +                        dcoord(1,n2), gres(1,n1), gres(1,n2))
         enddo

         do i=1,n
            gres(1,ind(i)) = dx(i) - dcoord(1,ind(i))
            gres(2,ind(i)) = dy(i) - dcoord(2,ind(i))
         enddo

         do i=1,np
            if(ptype(i).ne.interior .and. ptype(i).ne.solid)then
               gres(1,i) = -dcoord(1,i)
               gres(2,i) = -dcoord(2,i)
            endif
         enddo

         residue = 0.0d0
         do i=1,np
            ddx         = efact*eps(i)*gres(1,i)
            ddy         = efact*eps(i)*gres(2,i)
            residue     = residue + ddx**2 + ddy**2
            dcoord(1,i) = dcoord(1,i) + ddx
            dcoord(2,i) = dcoord(2,i) + ddy
         enddo
         residue = dsqrt(residue)/np

         if(giter .eq. 1) residue1 = residue
         if(residue1 .ne. 0.0) residue = residue/residue1
         if(mod(giter,100) .eq. 0) print*,giter,residue
      enddo

      if(residue .gt. RTOL)then
         print*,'============== WARNING ====================='
         print*,'Grid smoothing iterations have not converged'
         print*,'Do you still want to continue ?'
         print*,'============== WARNING ====================='
      endif

C Update the grid coordinates
      do i=1,np
         coord(1,i) = coord(1,i) + dcoord(1,i)
         coord(2,i) = coord(2,i) + dcoord(2,i)
      enddo

      return
      end
