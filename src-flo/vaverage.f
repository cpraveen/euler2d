      subroutine vaverage(x1, x2, x3, af1, af2, af3, tarea, qc, 
     +                    qv1, qv2, qv3)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), x3(2), af1(3), af2(3), af3(3), 
     +                 tarea, qc(nvar), qv1(nvar), qv2(nvar), qv3(nvar)

      integer          i
      double precision xt, yt, dx1, dx2, dx3, dy1, dy2, dy3, dr1, dr2,
     +                 dr3, w1, w2, w3

      xt = ( x1(1) + x2(1) + x3(1) )/3.0d0
      yt = ( x1(2) + x2(2) + x3(2) )/3.0d0

      dx1= xt - x1(1)
      dy1= yt - x1(2)
      dr1= dsqrt(dx1**2 + dy1**2)
      w1 = 1.0d0 - (af1(1)*dx1 + af1(2)*dy1)/dr1

      dx2= xt - x2(1)
      dy2= yt - x2(2)
      dr2= dsqrt(dx2**2 + dy2**2)
      w2 = 1.0d0 - (af2(1)*dx2 + af2(2)*dy2)/dr2

      dx3= xt - x3(1)
      dy3= yt - x3(2)
      dr3= dsqrt(dx3**2 + dy3**2)
      w3 = 1.0d0 - (af3(1)*dx3 + af3(2)*dy3)/dr3

      do i=1,nvar
         qv1(i) = qv1(i) + w1*qc(i)/dr1
         qv2(i) = qv2(i) + w2*qc(i)/dr2
         qv3(i) = qv3(i) + w3*qc(i)/dr3
      enddo

      return
      end
