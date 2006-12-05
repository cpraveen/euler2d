C-----------------------------------------------------------------------------
C Differentiates grid deformation to get gradient wrt boundary
C coordinates. Requires dJ/dx and dJ/dy which are given by adjoint
C solver.
C-----------------------------------------------------------------------------
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

      integer          i, n, ind(nbpmax)
      double precision dx(nbpmax), dy(nbpmax)
      double precision dxb(nbpmax), dyb(nbpmax), coordb(2,npmax)

      call read_input
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, esubp, coord, drmin, carea, af)

      n = 0
      do i=1,np
         if(ptype(i) .eq. solid)then
            n      = n + 1
            ind(n) = i
         endif
      enddo

      do i=1,n
         dx(i) = 0.0d0
         dy(i) = 0.0d0
      enddo

      print*,'Reading gradients from ADJ.DXY ...'
      open(20, file='ADJ.DXY', status='OLD')
      do i=1,np
         read(20,*) coordb(1,i), coordb(2,i)
      enddo
      close(20)

      call smooth_grid_bx(n, ind, ptype, edge, dx, dxb, dy, dyb,
     +                    coord, coordb)

      return
      stop
      end

C-----------------------------------------------------------------------------
C Propagate shape deformation to interior grid points using spring model
C coordb is the gradient of cost function wrt grid point coordinates;
C this is given by the adjoint solver.
C Input:
C Output:
C-----------------------------------------------------------------------------
      subroutine smooth_grid_bx(n, ind, ptype, edge, dx, dxb, 
     +                         dy, dyb, coord, coordb)
      implicit none
      include 'param.h'
      include 'param2.h'
      integer n, ind(nbpmax), ptype(npmax), edge(2,nemax)
      double precision dx(nbpmax), dy(nbpmax), coord(2,npmax)
      double precision dxb(nbpmax), dyb(nbpmax), coordb(2,npmax)

      integer          i, giter, n1, n2
      double precision residue, residue1, dcoordb(2,npmax), 
     +                 gres1(2), gres2(2), gresb(2,npmax), eps(npmax), 
     +                 ddxb, ddyb, dcoord1(2), dcoord2(2), gresb1(2),
     +                 gresb2(2)

      print*,'Adjoint of grid deformation ...'

C Forward sweep
C dcoord is not actually used anywhere; at this stage we dont know dx
C and dy; they are not required
C     do i=1,np
C        dcoord(1,i) = 0.0d0
C        dcoord(2,i) = 0.0d0
C     enddo

C     do i=1,n
C        dcoord(1,ind(i)) = dx(i)
C        dcoord(2,ind(i)) = dy(i)
C     enddo

C Dummy variables; in place of dcoord since it is not required
      dcoord1(1) = 0.0d0
      dcoord1(2) = 0.0d0
      dcoord2(1) = 0.0d0
      dcoord2(2) = 0.0d0

C Iterative solution

C Dummy variables; in place of gres, which is anyway zero
      gres1(1) = 0.0d0
      gres1(2) = 0.0d0
      gres2(1) = 0.0d0
      gres2(2) = 0.0d0

C Initialize adjoint variable to zero
      do i=1,np
         gresb(1,i)   = 0.0d0
         gresb(2,i)   = 0.0d0
      enddo

      residue = 1.0d0
      residue1= 1.0d0
      giter   = 0
      call grid_eps(edge, ptype, coord, eps)
      do while(residue .gt. RTOL)
         giter = giter + 1

         do i=1,np
            dcoordb(1,i) = coordb(1,i)
            dcoordb(2,i) = coordb(2,i)
         enddo

         do i=1,np
            if(ptype(i).ne.interior)then
               dcoordb(1,i) = dcoordb(1,i) - gresb(1,i)
               dcoordb(2,i) = dcoordb(2,i) - gresb(2,i)
            endif
         enddo

         do i=1,ne
            n1 = edge(1,i)
            n2 = edge(2,i)

            if(ptype(n1).ne.interior)then
               gresb1(1) = 0.0d0
               gresb1(2) = 0.0d0
            else
               gresb1(1) = gresb(1,n1)
               gresb1(2) = gresb(2,n1)
            endif

            if(ptype(n2).ne.interior)then
               gresb2(1) = 0.0d0
               gresb2(2) = 0.0d0
            else
               gresb2(1) = gresb(1,n2)
               gresb2(2) = gresb(2,n2)
            endif

            call grid_residue_bx(coord(1,n1), coord(1,n2), 
     +                           dcoord1,     dcoordb(1,n1),
     +                           dcoord2,     dcoordb(1,n2),
     +                           gres1,       gresb1,
     +                           gres2,       gresb2)
         enddo

         residue = 0.0d0
         do i=1,np
            ddxb       = dcoordb(1,i)
            ddyb       = dcoordb(2,i)
            gresb(1,i) = gresb(1,i) + efact*eps(i)*ddxb
            gresb(2,i) = gresb(2,i) + efact*eps(i)*ddyb
            residue    = residue + ddxb**2 + ddyb**2
         enddo
         residue = dsqrt(residue)/np

         if(giter .eq. 1)then
            residue1 = residue
            print*,'Residue in first iteration = ', residue1
         endif
         if(residue1 .ne. 0.0) residue = residue/residue1
         if(mod(giter,100) .eq. 0) print*,giter,residue
      enddo

      if(residue .gt. RTOL)then
         print*,'Adjoint grid smoothing iterations have not converged'
         print*,'Do you still want to continue ?'
         pause
      endif

C Finally, gradient wrt boundary coordinates
      do i=1,n
         dxb(i) = gresb(1,ind(i))
         dyb(i) = gresb(2,ind(i))
      enddo

      print*,'Boundary gradients written into DEF.OUT'
      open(20, file='DEF.OUT')
      do i=1,n
         write(20,*) coord(1,ind(i)), coord(2,ind(i)), dxb(i), dyb(i)
      enddo
      close(20)


      return
      end
