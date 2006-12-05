C----------------------------------------------------------------------------
C Add elements of matrix for computing averaging weights
C----------------------------------------------------------------------------
      subroutine afact1(x1, x2, x3, sax1, sax2, sax3, say1, say2,
     +                  say3, sax21, sax22, sax23, say21, say22, 
     +                  say23, saxy1, saxy2, saxy3)
      implicit none
      double precision x1(2), x2(2), x3(2), sax1, sax2, sax3, say1,
     +                 say2, say3, sax21, sax22, sax23, say21, say22, 
     +                 say23, saxy1, saxy2, saxy3

      double precision xt, yt, dx1, dy1, dr12, dr1, dx2, dy2, dr22, dr2,
     +                 dx3, dy3, dr32, dr3

      xt    = (x1(1) + x2(1) + x3(1))/3.0d0
      yt    = (x1(2) + x2(2) + x3(2))/3.0d0

      dx1   = xt     - x1(1)
      dy1   = yt     - x1(2)
      dr12  = dx1**2 + dy1**2
      dr1   = dsqrt(dr12)
      sax1  = sax1  + dx1/dr1
      say1  = say1  + dy1/dr1
      sax21 = sax21 + dx1**2/dr12
      say21 = say21 + dy1**2/dr12
      saxy1 = saxy1 + dx1*dy1/dr12

      dx2   = xt     - x2(1)
      dy2   = yt     - x2(2)
      dr22  = dx2**2 + dy2**2
      dr2   = dsqrt(dr22)
      sax2  = sax2  + dx2/dr2
      say2  = say2  + dy2/dr2
      sax22 = sax22 + dx2**2/dr22
      say22 = say22 + dy2**2/dr22
      saxy2 = saxy2 + dx2*dy2/dr22

      dx3   = xt     - x3(1)
      dy3   = yt     - x3(2)
      dr32  = dx3**2 + dy3**2
      dr3   = dsqrt(dr32)
      sax3  = sax3  + dx3/dr3
      say3  = say3  + dy3/dr3
      sax23 = sax23 + dx3**2/dr32
      say23 = say23 + dy3**2/dr32
      saxy3 = saxy3 + dx3*dy3/dr32

      return
      end

C----------------------------------------------------------------------------
C Finds normal for a boundary point
C NOTE: nx and ny must be initialized to zero before passing to this
C function. This is required for AD.
C----------------------------------------------------------------------------
      subroutine bd_normal(x1, xp, x2, nx, ny)
      implicit none
      double precision x1(2), xp(2), x2(2), nx, ny
      double precision dx1, dy1, ds1, nx1, ny1,
     +                 dx2, dy2, ds2, nx2, ny2,
     +                 nxa, nya, nsa
      dx1= xp(1) - x1(1)
      dy1= xp(2) - x1(2)
      ds1= dsqrt(dx1**2 + dy1**2)
      nx1=-dy1/ds1
      ny1= dx1/ds1
      dx2= x2(1) - xp(1)
      dy2= x2(2) - xp(2)
      ds2= dsqrt(dx2**2 + dy2**2)
      nx2=-dy2/ds2
      ny2= dx2/ds2
      nxa= 0.5d0*(nx1 + nx2)
      nya= 0.5d0*(ny1 + ny2)
      nsa= dsqrt(nxa**2 + nya**2)
      nx = nx + nxa/nsa
      ny = ny + nya/nsa
      return
      end

C----------------------------------------------------------------------------
C Add elements of matrix for computing averaging weights with boundary
C correction
C----------------------------------------------------------------------------
      subroutine afact2(x1, x2, x3, xp, nx, ny, sax, say, sax2,
     +                  say2, saxy)
      implicit none
      double precision x1(2), x2(2), x3(2), xp(2), nx, ny, sax, say,
     +                 sax2, say2, saxy

      double precision xt, yt, xn, xg, yg, dx1, dy1, dr12, dr1

      xt   = (x1(1) + x2(1) + x3(1))/3.0d0
      yt   = (x1(2) + x2(2) + x3(2))/3.0d0

      xn   = (xt - xp(1))*nx + (yt - xp(2))*ny
      xg   = xt - 2.0d0*xn*nx
      yg   = yt - 2.0d0*xn*ny

      dx1  = xg - xp(1)
      dy1  = yg - xp(2)
      dr12 = dx1**2 + dy1**2
      dr1  = dsqrt(dr12)

      sax  = sax  + dx1/dr1
      say  = say  + dy1/dr1
      sax2 = sax2 + dx1**2/dr12
      say2 = say2 + dy1**2/dr12
      saxy = saxy + dx1*dy1/dr12

      return
      end

C----------------------------------------------------------------------------
C Computes the averaging weights
C----------------------------------------------------------------------------
      subroutine afact3(sax, say, sax2, say2, saxy, afact)
      implicit none
      double precision sax, say, sax2, say2, saxy, afact(3)
      double precision det
      det      = sax2*say2 - saxy**2
      afact(1) = afact(1) + (say2*sax - saxy*say)/det
      afact(2) = afact(2) + (sax2*say - saxy*sax)/det
      return
      end

C----------------------------------------------------------------------------
C Computes denominator in vertex averaging formula
C----------------------------------------------------------------------------
      subroutine afact4(x1, x2, x3, tarea, af1, af2, af3)
      implicit none
      double precision x1(2), x2(2), x3(2), tarea, af1(3), af2(3),
     +                 af3(3)

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

      af1(3) = af1(3) + w1/dr1
      af2(3) = af2(3) + w2/dr2
      af3(3) = af3(3) + w3/dr3

      return
      end

C----------------------------------------------------------------------------
C Check min and max range of weights
C----------------------------------------------------------------------------
      subroutine checkweights(elem, coord, afact)
      implicit none
      include 'size.h'
      integer          elem(3,ntmax)
      double precision coord(2,npmax), afact(3,npmax)

      integer          i, v1, v2, v3
      double precision dx1, dx2, dx3, dy1, dy2, dy3, 
     +                 xt, yt, w1, w2, w3, wmin, wmax,
     +                 dr12, dr22, dr32, dr1, dr2, dr3

      wmin = 1.0d20
      wmax =-1.0d20
      open(20, file='wt.dat')
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)

         xt = (coord(1,v1) + coord(1,v2) + coord(1,v3))/3.0d0
         yt = (coord(2,v1) + coord(2,v2) + coord(2,v3))/3.0d0

         dx1      = xt       - coord(1,v1)
         dy1      = yt       - coord(2,v1)
         dr12     = dx1**2 + dy1**2
         dr1      = dsqrt(dr12)
         w1       = 1.0d0 - (afact(1,v1)*dx1 + afact(2,v1)*dy1)/dr1

         dx2      = xt       - coord(1,v2)
         dy2      = yt       - coord(2,v2)
         dr22     = dx2**2 + dy2**2
         dr2      = dsqrt(dr22)
         w2       = 1.0d0 - (afact(1,v2)*dx2 + afact(2,v2)*dy2)/dr2

         dx3      = xt       - coord(1,v3)
         dy3      = yt       - coord(2,v3)
         dr32     = dx3**2 + dy3**2
         dr3      = dsqrt(dr32)
         w3       = 1.0d0 - (afact(1,v3)*dx3 + afact(2,v3)*dy3)/dr3

         wmin     = dmin1(wmin, w1)
         wmin     = dmin1(wmin, w2)
         wmin     = dmin1(wmin, w3)

         wmax     = dmax1(wmax, w1)
         wmax     = dmax1(wmax, w2)
         wmax     = dmax1(wmax, w3)

         if(w1 .le. 0.0d0) then
            print*,'Averageing factor is non-positive'
            print*,'Triangle =',i
            print*,'Vertex   =',v1
            write(20,*) coord(1,v1), coord(2,v1)
            write(20,*) xt, yt
            write(20,*)
         endif

         if(w2 .le. 0.0d0) then
            print*,'Averageing factor is non-positive'
            print*,'Triangle =',i
            print*,'Vertex   =',v2
            write(20,*) coord(1,v2), coord(2,v2)
            write(20,*) xt, yt
            write(20,*)
         endif

         if(w3 .le. 0.0d0) then
            print*,'Averageing factor is non-positive'
            print*,'Triangle =',i
            print*,'Vertex   =',v3
            write(20,*) coord(1,v3), coord(2,v3)
            write(20,*) xt, yt
            write(20,*)
         endif

      enddo
      close(20)

      write(*,'(2x, "Minimum weight           =", f8.4)') wmin
      write(*,'(2x, "Maximum weight           =", f8.4)') wmax

      return
      end
