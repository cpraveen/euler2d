C----------------------------------------------------------------------------
C Computes weights for vertex averaging using my new method
C----------------------------------------------------------------------------
      subroutine avgfact(ptype, elem, edge, bdedge, esubp, spts, coord,
     +                   tarea, afact)
      implicit none
      include 'param.h'
      integer          ptype(npmax), elem(3,ntmax), edge(2,nemax),
     +                 bdedge(2,nbpmax), esubp(mesubp,nbpmax), 
     +                 spts(nspmax)
      double precision coord(2,npmax), tarea(ntmax), afact(3,npmax)

      integer          i, j, ip, it, e1, e2, p1, p2, v1, v2, v3
      double precision sax(npmax), say(npmax), sax2(npmax), say2(npmax),
     +                 saxy(npmax), nx(nspmax), ny(nspmax)

      print*,'Finding weights for area averaging ...'

      do i=1,np
         sax(i)  = 0.0d0
         say(i)  = 0.0d0
         sax2(i) = 0.0d0
         say2(i) = 0.0d0
         saxy(i) = 0.0d0
         afact(1,i)= 0.0d0
         afact(2,i)= 0.0d0
         afact(3,i)= 0.0d0
      enddo

      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call afact1(coord(1,v1), coord(1,v2), coord(1,v3), 
     +               sax(v1),  sax(v2),  sax(v3), 
     +               say(v1),  say(v2),  say(v3), 
     +               sax2(v1), sax2(v2), sax2(v3), 
     +               say2(v1), say2(v2), say2(v3), 
     +               saxy(v1), saxy(v2), saxy(v3))
      enddo

C Add ghost cell contributions for solid wall points
c We dont need this as it does not still give linear consistency for wall points
c     do i=1,nsp
c        ip = spts(i)
c        e1 = bdedge(1,i)
c        e2 = bdedge(2,i)
c        p1 = edge(1,e1)
c        p2 = edge(2,e2)
c        nx(i) = 0.0d0
c        ny(i) = 0.0d0
c        call bd_normal(coord(1,p1), coord(1,ip), coord(1,p2), 
c    +                  nx(i), ny(i))
c        do j=1,esubp(1,i)
c           it = esubp(j+1,i)
c           v1 = elem(1,it)
c           v2 = elem(2,it)
c           v3 = elem(3,it)
c           call afact2(coord(1,v1), coord(1,v2), coord(1,v3),
c    +                  coord(1,ip), nx(i), ny(i), sax(ip), say(ip), 
c    +                  sax2(ip), say2(ip), saxy(ip))
c        enddo
c     enddo

C Compute weights by inverting least squares matrix
      do i=1,np
         call afact3(sax(i), say(i), sax2(i), say2(i), saxy(i),
     +               afact(1,i))
c        if(ptype(i) .ne. interior .and. ptype(i) .ne. solid)then
         if(ptype(i) .ne. interior)then
            afact(1,i)= 0.0d0
            afact(2,i)= 0.0d0
         endif
      enddo

C Compute denominator in vertex averaging formula
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call afact4(coord(1,v1), coord(1,v2), coord(1,v3),
     +               tarea(i), afact(1,v1), afact(1,v2), 
     +               afact(1,v3))
      enddo

C Check min and max range of weights
      call checkweights(elem, coord, afact)

      return
      end
