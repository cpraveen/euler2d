      subroutine smooth(ptype, psup1, psup2, elem, esup1, esup2, coord)
      implicit none
      include 'param.h'
      integer          ptype(npmax), psup1(mpsup*npmax), psup2(npmax+1),
     &                 elem(3,ntmax), esup1(mesup*npmax), esup2(npmax+1)
      double precision coord(2,npmax)

      integer          i, ipoin, ip, non(npmax), ismooth, nfil
      double precision res

      print*,'Smoothing the grid...'

      do ipoin = 1,np
         non(ipoin) = 0
         do ip=psup2(ipoin)+1, psup2(ipoin+1)
            non(ipoin) = non(ipoin) + 1
         enddo
      enddo

      ismooth = 20
      do i=1,ismooth
         call laplace(ptype, psup1, psup2, elem, esup1, esup2, non, 
     &                coord, res)
         print*,'Iteration = ',i,' Residue = ',res
      enddo

      nfil = 20
      open(unit=nfil, file='node')
      write(nfil,*) np, 2, 0, 1
      do i=1,np
         write(nfil,10) i, coord(1,i), coord(2,i), ptype(i)
      enddo
      close(nfil)
10    format(I8, 2X, 2F18.8, 2X, I8)
      return
      end

c Laplacian smoothing
      subroutine laplace(ptype, psup1, psup2, elem, esup1, esup2, non, 
     &                   coord, res)
      implicit none
      include 'param.h'
      integer          ptype(npmax), psup1(mpsup*npmax), psup2(npmax+1),
     &                 non(npmax), elem(3,ntmax), esup1(mesup*npmax),
     &                 esup2(npmax+1)
      double precision coord(2,npmax), res

      integer ipoin, ip, neigh, stat, convex
      double precision dx, dy, sfact, xx, yy

      sfact = 0.5d0
      res   = 0.0d0

      do ipoin = 1,np
         if(ptype(ipoin) .eq. interior)then
            dx = 0.0d0
            dy = 0.0d0
            do ip=psup2(ipoin)+1, psup2(ipoin+1)
                  neigh = psup1(ip)
                  dx = dx + ( coord(1,neigh) - coord(1,ipoin) )
                  dy = dy + ( coord(2,neigh) - coord(2,ipoin) )
            enddo
            dx = dx/non(ipoin)
            dy = dy/non(ipoin)
            xx = coord(1,ipoin) + sfact*dx
            yy = coord(2,ipoin) + sfact*dy
            res = res + dabs(coord(1,ipoin)-xx) +
     &                  dabs(coord(2,ipoin)-yy)
            stat = convex(elem, esup1, esup2, coord, xx, yy)
            if(stat .eq. 1)then
                  coord(1,ipoin) = xx
                  coord(2,ipoin) = yy
            else
                  print*,'Convexity violated for node =', ipoin
                  print*,'x =',coord(1,ipoin)
                  print*,'y =',coord(2,ipoin)
            endif
         endif
      enddo

      res = res/np
      
      return
      end

c Check if new point lies inside the convex hull - NOT COMPLETE
      integer function convex(elem, esup1, esup2, coord, xx, yy)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), esup1(mesup*npmax), esup2(npmax+1)
      double precision coord(2,npmax), xx, yy

      integer          ipoin, iesup, ie, n1, n2, n3

      convex = 1
      do ipoin=1,np
         do iesup=esup2(ipoin)+1, esup2(ipoin+1)
            ie = esup1(iesup)
            n1 = elem(1,ie)
            n2 = elem(2,ie)
            n3 = elem(3,ie)
         enddo
      enddo

      return
      end
