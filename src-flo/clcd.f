C.....Calculate lift and drag coefficients
      subroutine clcd(edge, tedge, coord, qc, cl, cd)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision coord(2,npmax), qc(nvar,ntmax), cl, cd

      integer          i, e1, e2, c1
      double precision xf, yf, d, u, v, p, dx, dy, nx, ny

      xf = 0.0d0
      yf = 0.0d0
      do i=nsw1,nsw2
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         call contoprim(qc(1,c1), d, u, v, p)
         dx = coord(1,e2) - coord(1,e1)
         dy = coord(2,e2) - coord(2,e1)
         nx = dy
         ny =-dx
         xf = xf + p*nx
         yf = yf + p*ny
      enddo

      cl =(-dsin(aoa)*xf + dcos(aoa)*yf)/(0.5d0*r_inf*q_inf**2)
      cd =( dcos(aoa)*xf + dsin(aoa)*yf)/(0.5d0*r_inf*q_inf**2)

      return
      end
