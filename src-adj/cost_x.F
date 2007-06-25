C Gradient of cost function with respect to grid coordinates
      subroutine cost_x(elem, edge, tedge, vedge, spts, bdedge, esubp,
     +                  ptype, coord, qc, qv, qvb, carea, af, resb,
     +                  coordb, cl, cd)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 vedge(2,nemax), spts(nspmax), bdedge(2,nbpmax),
     +                 esubp(mesubp,nbpmax), ptype(npmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), qvb(nvar,npmax), carea(ntmax), 
     +                 afb(3,npmax), resb(nvar,ntmax), coordb(2,npmax)
      double precision cl, cd
      double precision res1(nvar), res2(nvar)
      double precision cost, costb, xd, yd

      integer          i, j, ie, v1, v2, v3, e1, e2, c1, c2

      do i=1,nvar
         res1(i) = 0.0d0
         res2(i) = 0.0d0
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Cost func derivative wrt grid.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      print*,'Computing derivative of cost function wrt grid ...'
#if defined COST1
      call cost1_x(edge, tedge, coord, coordb, qc, cost)
#elif defined COST2
      call cost2_x(ptype, elem, esubp, spts, edge, tedge, bdedge,
     +             coord, coordb, af, afb, carea, qc, qv, qvb, cost)
#elif defined COST3
      do i=1,np
         coordb(1,i) = 0.0d0
         coordb(2,i) = 0.0d0
      enddo
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Residual derivative wrt grid
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,np
         do j=1,nvar
            qvb(j,i)  = 0.0d0
         enddo
         afb(1,i)     = 0.0d0
         afb(2,i)     = 0.0d0
         afb(3,i)     = 0.0d0
      enddo

C Edge count index
      ie = 0

C Compute flux for interior edges
      do i=1,nin
         ie = ie + 1
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         c2 = tedge(2,ie)
         v1 = vedge(1,ie)
         v2 = vedge(2,ie)
         call kfvs_flux_bx(coord(1,e1), coordb(1,e1),
     +                    coord(1,e2), coordb(1,e2),
     +                    qc(1,c1), qc(1,c2),
     +                    qv(1,v1), qvb(1,v1), 
     +                    qv(1,v2), qvb(1,v2),
     +                    res1, resb(1,c1), res2, resb(1,c2))
      enddo

C Compute flux for solid wall edges
      do ie=nsw1,nsw2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call solid_flux_bx(coord(1,e1), coordb(1,e1),
     +                      coord(1,e2), coordb(1,e2),
     +                      qc(1,c1), res1, resb(1,c1))       
      enddo

C Flux for far-field points
      do i=1,nff
         ie = ie + 1
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call farfield_flux_bx(coord(1,e1), coordb(1,e1),
     +                         coord(1,e2), coordb(1,e2),
     +                         qc(1,c1), cl, cd, 
     +                         res1, resb(1,c1))
      enddo

C Contribution from vertex averaging
      do i=1,nsp
         j = spts(i)
         e1= bdedge(1,i)
         e2= bdedge(2,i)
         v1= edge(1,e1)
         v2= j
         v3= edge(2,e2)
         call killnormalvel_bx(coord(1,v1), coordb(1,v1),
     +                         coord(1,v2), coordb(1,v2),
     +                         coord(1,v3), coordb(1,v3),
     +                         qv(1,j), qvb(1,j))
      enddo
      do i=1,np
         do j=1,nvar
            afb(3,i) = afb(3,i) - qv(j,i)*qvb(j,i)/af(3,i)**2
            qvb(j,i) = qvb(j,i)/af(3,i)
         enddo
      enddo
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call vaverage_bx(coord(1,v1), coordb(1,v1),
     +                    coord(1,v2), coordb(1,v2), 
     +                    coord(1,v3), coordb(1,v3),
     +                    af(1,v1), afb(1,v1), 
     +                    af(1,v2), afb(1,v2), 
     +                    af(1,v3), afb(1,v3), 
     +                    carea(i), qc(1,i), 
     +                    qv(1,v1), qvb(1,v1),
     +                    qv(1,v2), qvb(1,v2), 
     +                    qv(1,v3), qvb(1,v3))
      enddo

      call avgfact_x(ptype, elem, edge, bdedge, esubp, spts, 
     +                     coord, coordb, carea, af, afb)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Final derivative
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      costb = 0.0d0
      do i=1,np
         xd    = coord(2,i)
         yd    =-coord(1,i)
         costb = costb + coordb(1,i)*xd + coordb(2,i)*yd
      enddo
      print*,'Gradient =',costb,
     +       costb*dsqrt(1.0d0-mach_inf**2)/2.0d0/M_PI

C Save boundary gradients
      open(15, file='ADJ.DER')
      do i=1,nsp
         j = spts(i)
         write(15, '(4e20.10)') coord(1,j), coord(2,j), coordb(1,j),
     +                          coordb(2,j)
      enddo
      close(15)

C Save all gradients
      open(15, file='ADJ.DXY')
      do i=1,np
         write(15, '(2e20.10)') coordb(1,i), coordb(2,i)
      enddo
      close(15)

      return
      end

C---------------------------------------------------------------------
C Cost func partial derivative wrt grid.
C Cost func depends on cell-center value
C---------------------------------------------------------------------
      subroutine cost1_x(edge, tedge, coord, coordb, qc, cost)
      implicit none
      include 'common.h'
      include 'size.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision coord(2,npmax), coordb(2,npmax), qc(nvar,ntmax),
     +                 cost

      integer          i, e1, e2, c1
      double precision costb

      do i=1,np
         coordb(1,i) = 0.0d0
         coordb(2,i) = 0.0d0
      enddo

      costb = 1.0d0
      do i=nsw1,nsw2
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         call costfunc_bx(coord(1,e1), coordb(1,e1), 
     +                    coord(1,e2), coordb(1,e2),
     +                    qc(1,c1), cost, costb)
      enddo

      return
      end

C---------------------------------------------------------------------
C Cost func partial derivative wrt grid.
C Cost func depends on vertex values
C---------------------------------------------------------------------
      subroutine cost2_x(ptype, elem, esubp, spts, edge, tedge, bdedge,
     +                   coord, coordb, af, afb, carea, qc, qv, qvb, 
     +                   cost)
      implicit none
      include 'common.h'
      include 'size.h'
      integer          ptype(npmax), elem(3,ntmax), 
     +                 esubp(mesubp,nbpmax), spts(nspmax), 
     +                 edge(2,nemax), tedge(2,nemax), bdedge(2,nbpmax)
      double precision coord(2,npmax), coordb(2,npmax), qc(nvar,ntmax),
     +                 qv(nvar,npmax), qvb(nvar,npmax), af(3,npmax), 
     +                 afb(3,npmax), carea(ntmax), cost

      integer          i, j, e1, e2, v1, v2, v3, c1
      double precision costb

      do i=1,np
         do j=1,nvar
            qvb(j,i)  = 0.0d0
         enddo
         coordb(1,i) = 0.0d0
         coordb(2,i) = 0.0d0
         afb(1,i)     = 0.0d0
         afb(2,i)     = 0.0d0
         afb(3,i)     = 0.0d0
      enddo

      costb = 1.0d0
      do i=nsw1,nsw2
         e1 = edge(1,i)
         e2 = edge(2,i)
         call costfunc_bx(coord(1,e1), coordb(1,e1), 
     +                    coord(1,e2), coordb(1,e2),
     +                    qv(1,e1), qvb(1,e1), qv(1,e2), qvb(1,e2),
     +                    cost, costb)
      enddo

C Contribution from vertex averaging
      do i=1,nsp
         j = spts(i)
         e1= bdedge(1,i)
         e2= bdedge(2,i)
         v1= edge(1,e1)
         v2= j
         v3= edge(2,e2)
         call killnormalvel_bx(coord(1,v1), coordb(1,v1),
     +                         coord(1,v2), coordb(1,v2),
     +                         coord(1,v3), coordb(1,v3),
     +                         qv(1,j), qvb(1,j))
      enddo
      do i=1,np
         do j=1,nvar
            afb(3,i) = afb(3,i) - qv(j,i)*qvb(j,i)/af(3,i)**2
            qvb(j,i) = qvb(j,i)/af(3,i)
         enddo
      enddo
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call vaverage_bx(coord(1,v1), coordb(1,v1),
     +                    coord(1,v2), coordb(1,v2), 
     +                    coord(1,v3), coordb(1,v3),
     +                    af(1,v1), afb(1,v1), 
     +                    af(1,v2), afb(1,v2), 
     +                    af(1,v3), afb(1,v3), 
     +                    carea(i), qc(1,i), 
     +                    qv(1,v1), qvb(1,v1),
     +                    qv(1,v2), qvb(1,v2), 
     +                    qv(1,v3), qvb(1,v3))
      enddo

      call avgfact_x(ptype, elem, edge, bdedge, esubp, spts, 
     +                     coord, coordb, carea, af, afb)

      return
      end
