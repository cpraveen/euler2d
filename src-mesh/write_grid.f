C-----------------------------------------------------------------------------
C.....Read grid data from a file
C.....Currently supports only triangular elements
C-----------------------------------------------------------------------------
      subroutine write_new_grid(coord, elem, ptype)
      implicit none
      include 'param.h'
      double precision coord(2,npmax)
      integer          elem(3,ntmax), ptype(npmax)

      integer          ngrid, ip, i, j

      print*,'Write new grid to file MESH'

      ngrid = 10
      open(ngrid, file="MESH")
      write(ngrid,*) np, nt

      do ip=1,np
         write(ngrid,*) ip, coord(1,ip), coord(2,ip), ptype(ip)
      enddo

      j = 0
      do ip=1,nt
         write(ngrid,*) ip, elem(1, ip), elem(2, ip), elem(3, ip), j
      enddo

      close(ngrid)

      return
      end
