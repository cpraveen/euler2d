C Computes weights for vertex averaging using my new method
      subroutine avgfact_x(ptype, elem, edge, bdedge, esubp, spts, 
     +                     coord, coordb, tarea, af, afb)
      implicit none
      include 'param.h'
      integer          ptype(npmax), elem(3,ntmax), edge(2,nemax),
     +                 bdedge(2,nbpmax), esubp(mesubp,nbpmax), 
     +                 spts(nspmax)
      double precision coord(2,npmax), coordb(2,npmax), tarea(ntmax), 
     +                 af(3,npmax), afb(3,npmax)

      integer          i, j, ip, it, e1, e2, p1, p2, v1, v2, v3
      double precision sax(npmax), say(npmax), sax2(npmax), say2(npmax),
     +                 saxy(npmax), nx(nspmax), ny(nspmax)
      double precision saxb(npmax), sayb(npmax), sax2b(npmax), 
     +                 say2b(npmax), saxyb(npmax), nxb(nspmax), 
     +                 nyb(nspmax), xpb(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Forward sweep
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,np
         sax(i)  = 0.0d0
         say(i)  = 0.0d0
         sax2(i) = 0.0d0
         say2(i) = 0.0d0
         saxy(i) = 0.0d0
      enddo

C Add elements of least squares matrix
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Reverse sweep
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call afact4_bx(coord(1,v1), coordb(1,v1),
     +                  coord(1,v2), coordb(1,v2),
     +                  coord(1,v3), coordb(1,v3),
     +                  tarea(i), af(1,v1), afb(1,v1),
     +                  af(1,v2), afb(1,v2), af(1,v3), afb(1,v3))
      enddo

      do i=1,np
         saxb(i)  = 0.0d0
         sayb(i)  = 0.0d0
         sax2b(i) = 0.0d0
         say2b(i) = 0.0d0
         saxyb(i) = 0.0d0
      enddo

C Compute weights by inverting least squares matrix
      do i=1,np
c        if(ptype(i) .ne. interior .and. ptype(i) .ne. solid)then
         if(ptype(i) .ne. interior)then
            afb(1,i)= 0.0d0
            afb(2,i)= 0.0d0
         endif
         call afact3_bx(sax(i), saxb(i), say(i), sayb(i), 
     +                  sax2(i), sax2b(i), say2(i), say2b(i),
     +                  saxy(i), saxyb(i), af(1,i), afb(1,i))
      enddo

C Add ghost cell contributions for solid wall points
c     do i=1,nsp
c        ip = spts(i)
c        e1 = bdedge(1,i)
c        e2 = bdedge(2,i)
c        p1 = edge(1,e1)
c        p2 = edge(2,e2)
c        nxb(i) = 0.0d0
c        nyb(i) = 0.0d0
c        xpb(1) = 0.0d0
c        xpb(2) = 0.0d0
c        do j=1,esubp(1,i)
c           it = esubp(j+1,i)
c           v1 = elem(1,it)
c           v2 = elem(2,it)
c           v3 = elem(3,it)
c           call afact2_bx(coord(1,v1), coordb(1,v1),
c    +                     coord(1,v2), coordb(1,v2),
c    +                     coord(1,v3), coordb(1,v3),
c    +                     coord(1,ip), xpb,
c    +                     nx(i), nxb(i), ny(i), nyb(i),
c    +                     sax(ip), saxb(ip), say(ip), sayb(ip),
c    +                     sax2(ip), sax2b(ip), say2(ip), say2b(ip),
c    +                     saxy(ip), saxyb(ip))
c        enddo
c        coordb(1,ip) = coordb(1,ip) + xpb(1)
c        coordb(2,ip) = coordb(2,ip) + xpb(2)
c        call bd_normal_bx(coord(1,p1), coordb(1,p1),
c    +                     coord(1,ip), coordb(1,ip),
c    +                     coord(1,p2), coordb(1,p2),
c    +                     nx(i), nxb(i), ny(i), nyb(i))
c     enddo

C Add elements of least squares matrix
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call afact1_bx(coord(1,v1), coordb(1,v1),
     +                  coord(1,v2), coordb(1,v2),
     +                  coord(1,v3), coordb(1,v3),
     +                  sax(v1), saxb(v1), 
     +                  sax(v2), saxb(v2), 
     +                  sax(v3), saxb(v3),
     +                  say(v1), sayb(v1), 
     +                  say(v2), sayb(v2), 
     +                  say(v3), sayb(v3),
     +                  sax2(v1), sax2b(v1), 
     +                  sax2(v2), sax2b(v2),
     +                  sax2(v3), sax2b(v3),
     +                  say2(v1), say2b(v1), 
     +                  say2(v2), say2b(v2),
     +                  say2(v3), say2b(v3),
     +                  saxy(v1), saxyb(v1), 
     +                  saxy(v2), saxyb(v2),
     +                  saxy(v3), saxyb(v3) )
      enddo

      return
      end
