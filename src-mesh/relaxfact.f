C-----------------------------------------------------------------------------
C.....Time-step or relaxation factor
C-----------------------------------------------------------------------------
      subroutine grid_eps(edge, ptype, coord, eps)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), ptype(npmax)
      double precision coord(2,npmax), eps(npmax)

      integer          i, n1, n2, count(npmax)
      double precision dx, dy, ds

      do i=1,np
         eps(i)   = 1.0e10
         count(i) = 0
      enddo

      do i=1,ne
         n1        = edge(1,i)
         n2        = edge(2,i)
         dx        = coord(1,n2) - coord(1,n1)
         dy        = coord(2,n2) - coord(2,n1)
         ds        = dsqrt(dx**2 + dy**2)
         eps(n1)   = dmin1(eps(n1), ds)
         eps(n2)   = dmin1(eps(n2), ds)
         count(n1) = count(n1) + 1
         count(n2) = count(n2) + 1
      enddo

      do i=1,np
         eps(i) = eps(i)/count(i)
         if(ptype(i).ne.interior) eps(i) = 1.0d0
      enddo

      return
      end
