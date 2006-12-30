C Cost func derivative wrt flow solution. This is stored in qcb1 and is
C used in adjoint solution. It is the rhs of adjoint solution, dJ/dq. For
C boundary cost function most of the entries of qcb1 are zero. Can we do
C something to save memory ?

      subroutine cost_q(spts, elem, edge, tedge, bdedge, coord, carea,
     +                  af, qc, qcb1, qv, qvb, cost)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), spts(nspmax),
     +                 tedge(2,nemax), bdedge(2,nbpmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax)
      double precision qcb1(nvar,ntmax), qvb(nvar,npmax), cost

#if defined COST1
      call cost1_q(edge, tedge, coord, qc, qcb1, cost)
#elif defined COST2
      call cost2_q(spts, elem, edge, bdedge, coord, carea, af, qc, qcb1,
     +                   qv, qvb, cost)
#endif

      return
      end


#ifdef COST1
C------------------------------------------------------------------------------
C Cost function depends on cell-center values
C------------------------------------------------------------------------------
      subroutine cost1_q(edge, tedge, coord, qc, qcb1, cost)
      implicit none
      include 'param.h'
      include 'param2.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision coord(2,npmax), qc(nvar,ntmax)
      double precision qcb1(nvar,ntmax)
      double precision cost, costb

      integer          i, j, e1, e2, c1

      do i=1,nt
         do j=1,nvar
            qcb1(j,i) = 0.0d0
         enddo
      enddo

      costb = 1.0d0
      do i=nsw1,nsw2
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         call costfunc_bq(coord(1,e1), coord(1,e2), 
     +                    qc(1,c1), qcb1(1,c1), cost, costb)
      enddo

      return
      end
#endif

#ifdef COST2
C------------------------------------------------------------------------------
C Cost function depends on vertex values
C------------------------------------------------------------------------------
      subroutine cost2_q(spts, elem, edge, bdedge, coord, carea, af, qc,
     +                   qcb1, qv, qvb, cost)
      implicit none
      include 'param.h'
      include 'param2.h'
      integer          elem(3,ntmax), edge(2,nemax), spts(nspmax),
     +                 bdedge(2,nbpmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax)
      double precision qcb1(nvar,ntmax), qvb(nvar,npmax)
      double precision cost, costb

      integer          i, j, v1, v2, v3, e1, e2

      do i=1,nt
         do j=1,nvar
            qcb1(j,i) = 0.0d0
         enddo
      enddo
      do i=1,np
         do j=1,nvar
            qvb(j,i) = 0.0d0
         enddo
      enddo

      costb = 1.0d0
      do i=nsw1,nsw2
         e1 = edge(1,i)
         e2 = edge(2,i)
         call costfunc_bq(coord(1,e1), coord(1,e2), 
     +                    qv(1,e1), qvb(1,e1), qv(1,e2), qvb(1,e2),
     +                    cost, costb)
      enddo

      do i=1,nsp
         j = spts(i)
         e1= bdedge(1,i)
         e2= bdedge(2,i)
         v1= edge(1,e1)
         v2= j
         v3= edge(2,e2)
         call killnormalvel_bq(coord(1,v1), coord(1,v2), 
     +                         coord(1,v3), qv(1,j), qvb(1,j))
      enddo
      do i=1,np
         do j=1,nvar
            qvb(j,i) = qvb(j,i)/af(3,i)
         enddo
      enddo
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call vaverage_bq(coord(1,v1), coord(1,v2), coord(1,v3),
     +                   af(1,v1), af(1,v2), af(1,v3), carea(i), 
     +                   qc(1,i), qcb1(1,i), qv(1,v1), qvb(1,v1),
     +                   qv(1,v2), qvb(1,v2), qv(1,v3), qvb(1,v3))
      enddo

      return
      end
#endif
