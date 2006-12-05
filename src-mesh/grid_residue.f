C-----------------------------------------------------------------------------
C Find rhs for spring model
C-----------------------------------------------------------------------------
      subroutine grid_residue(x1, x2, dx1, dx2, gres1, gres2)
      implicit none
      double precision x1(2), x2(2), dx1(2), dx2(2), gres1(2), gres2(2)

      double precision dx, dy, ds, k, ddx, ddy

      dx = x2(1) - x1(1)
      dy = x2(2) - x1(2)
      ds = dsqrt(dx**2 + dy**2)
      k  = 1.0d0/ds

      ddx= dx2(1) - dx1(1)
      ddy= dx2(2) - dx1(2)

      gres1(1) = gres1(1) + k*ddx
      gres1(2) = gres1(2) + k*ddy

      gres2(1) = gres2(1) - k*ddx
      gres2(2) = gres2(2) - k*ddy

      return
      end
