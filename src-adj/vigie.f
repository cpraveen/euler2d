      subroutine vigie(coord, elem, qv, coordb)
      implicit none
      include 'param.h'
      integer           elem(3,ntmax)
      double precision  coord(2,npmax), qv(nvar,npmax), coordb(2,npmax)

      integer           i, vig

      vig = 55
      open(unit=vig, file='ADJ.VIG')

      write(vig,*)'points',np
      do i=1,np
         write(vig,*) coord(1,i), coord(2,i)
      enddo

      write(vig,*)'triangles',nt
      do i=1,nt
         write(vig,*) elem(1,i), elem(2,i), elem(3,i)
      enddo

      write(vig,*)'scalars  Adjoint_Density'
      do i=1,np
         write(vig,*) qv(1,i)
      enddo

      write(vig,*)'scalars  Adjoint_x_Momentum'
      do i=1,np
         write(vig,*) qv(2,i)
      enddo

      write(vig,*)'scalars  Adjoint_y_Momentum'
      do i=1,np
         write(vig,*) qv(3,i)
      enddo

      write(vig,*)'scalars  Adjoint_Energy'
      do i=1,np
         write(vig,*) qv(4,i)
      enddo

      write(vig,*)'scalars  J_x'
      do i=1,np
         write(vig,*) coordb(1,i)
      enddo

      write(vig,*)'scalars  J_y'
      do i=1,np
         write(vig,*) coordb(2,i)
      enddo

      write(vig,*)'vectors  vel  u  v  1e-02'
      do i=1,np
         write(vig,*) qv(2,i), qv(3,i)
      enddo
      write(vig,*)'end_block'

      close(vig)

      return
      end
