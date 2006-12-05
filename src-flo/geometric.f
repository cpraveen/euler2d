C-----------------------------------------------------------------------------
      subroutine geometric(elem, edge, edneigh, esue, vedge, spts, 
     +                     ptype, bdedge, esubp, coord, drmin, elarea, 
     +                     afact)
      implicit none
      include 'param.h'

      integer          ptype(npmax), elem(3,ntmax),
     &                 esup1(mesup*npmax), esup2(npmax+1),
     &                 psup1(mpsup*npmax), psup2(npmax+1),
     &                 edge(2,nemax), edneigh(2,nemax), spts(nspmax),
     &                 opts(nopmax), bpts(nbpmax), fpts(nfpmax), 
     &                 vedge(2,nemax), esue(3,ntmax), bdedge(2,nbpmax),
     &                 esubp(mesubp,nbpmax)
      double precision coord(2, npmax), elarea(ntmax), drmin(ntmax),
     &                 afact(3,npmax)

c     Read grid from file
      call read_grid(coord, elem, ptype, spts, fpts, opts, bpts)

c     Generate .gnu files for visualization
      call prep_gnuplot

c     Make sure elements are oriented ccw
      call tri_orient(elem, coord)

c     Find elements surrounding a point
      call el_surr_point(elem, esup1, esup2)

c     Find elements surrounding solid boundary point
      call el_surr_bp(spts, esup1, esup2, esubp)

c     Find points surrounding a point
      call pt_surr_pt(esup1, esup2, elem, psup1, psup2)

c     Create edges
      call create_edge(psup1, psup2, edge)

c     Put edge list in a particular order
      call order_edges(ptype, edge)

c     Find element adjoining each edge
      call el_surr_edge(esup1, esup2, elem, edge, edneigh, vedge)

c For each solid wall point find two edges adjoining it
      call bdedges(spts, ptype, edge, bdedge)

c     Element surrounding element
      call el_surr_el(elem, edge, edneigh, esue)

c     Smooth the grid using Laplacian smoothing
c     call smooth(ptype, psup1, psup2, elem, esup1, esup2, coord)
      
      if(explicit .eq. no)then
         call renumber(elem, esue, edneigh, esubp)
      endif

c     Calculate triangle areas
      call tri_area(coord, elem, elarea)

c     Length scale for time-step calculation
      call dtlength(coord, elarea, elem, drmin)

c     Write grid in gnuplot format for visualization
      call write_grid(coord, edge, edneigh)

c     Area averaging factors
      call avgfact(ptype, elem, edge, bdedge, esubp, spts, coord, 
     &             elarea, afact)

      return
      end

C-----------------------------------------------------------------------------
C.....Read grid data from a file
C.....Currently supports only triangular elements
C-----------------------------------------------------------------------------
      subroutine read_grid(coord, elem, ptype, spts, fpts, opts, bpts)
      implicit none
      include 'param.h'
      double precision coord(2,npmax)
      integer          elem(3,ntmax), ptype(npmax), spts(nspmax),
     &                 fpts(nfpmax), opts(nopmax), bpts(nbpmax)

      integer          ngrid, ip, i, j

      print*,'Reading grid from file ',gridfile

      ngrid = 10
      open(ngrid, file=gridfile, status="old")
      rewind(ngrid)
      read(ngrid,*) np, nt
      write(*, '( " Number of points    :", i8)') np
      write(*, '( " Number of triangles :", i8)') nt

      if(np.gt.npmax) then
         print*, 'Increase the size of npmax'
         stop
      endif

      if(nt.gt.ntmax) then
         print*, 'Increase the size of ntmax'
         stop
      endif

      do ip=1,np
         read(ngrid,*) i, coord(1,ip), coord(2,ip), ptype(ip)
      enddo

      do ip=1,nt
         read(ngrid,*) i, elem(1, ip), elem(2, ip), elem(3, ip), j
      enddo

      close(ngrid)

      call rotate_grid(coord)

c     Find bounding box
      xmin = 1000000.0d0
      ymin = 1000000.0d0
      xmax =-1000000.0d0
      ymax =-1000000.0d0
      do ip=1,np
         xmin = dmin1(xmin, coord(1,ip))
         ymin = dmin1(ymin, coord(2,ip))

         xmax = dmax1(xmax, coord(1,ip))
         ymax = dmax1(ymax, coord(2,ip))
      enddo

      write(*,'(" Bounding box:")')
      write(*,'(10x, "xmin =", f8.3)') xmin
      write(*,'(10x, "xmax =", f8.3)') xmax
      write(*,'(10x, "ymin =", f8.3)') ymin
      write(*,'(10x, "ymax =", f8.3)') ymax

      nsp = 0
      nfp = 0
      nop = 0
      nbp = 0

      do ip=1,np
         if(ptype(ip) .eq. solid)then
            nsp       = nsp + 1
            spts(nsp) = ip
         endif
         if(ptype(ip) .eq. farfield)then
            nfp       = nfp + 1
            fpts(nfp) = ip
         endif
         if(ptype(ip) .eq. outflow)then
            nop       = nop + 1
            opts(nop) = ip
         endif
         if(ptype(ip) .ne. interior)then
            nbp       = nbp + 1
            bpts(nbp) = ip
         endif
      enddo

      write(*,'(" Number of solid    points =", i8)') nsp
      write(*,'(" Number of farfield points =", i8)') nfp
      write(*,'(" Number of outflow  points =", i8)') nop
      write(*,'(" Number of boundary points =", i8)') nbp

      if( nsp+nfp+nop .ne. nbp )then
         print*,'There seem to be some unrecognized point types'
         stop
      endif

      return
      end

C-----------------------------------------------------------------------------
C Rotate grid by aoa
C-----------------------------------------------------------------------------
      subroutine rotate_grid(x)
      implicit none
      include 'param.h'
      double precision x(2,npmax)

      integer          i
      double precision x1, y1

      print*,'========Translating and rotating the grid================'

C Translate grid to ref point of vortex model
      do i=1,np
         x(1,i) = x(1,i) - xref
         x(2,i) = x(2,i) - yref
      enddo

C Rotate grid
      aoa = aoa_deg*M_PI/180.0d0
      do i=1,np
         x1 = x(1,i)
         y1 = x(2,i)
         x(1,i) = dcos(aoa)*x1 + dsin(aoa)*y1
         x(2,i) = dcos(aoa)*y1 - dsin(aoa)*x1
      enddo

      aoa     = 0.0d0
      aoa_deg = 0.0d0
      xref    = 0.0d0
      yref    = 0.0d0

      return
      end

C-----------------------------------------------------------------------------
C Generate gnuplot script file for plotting grid
C-----------------------------------------------------------------------------
      subroutine prep_gnuplot
      implicit none
      include 'param.h'
      integer gnu

      gnu = 10
      open(unit=gnu, file='grid.gnu')
      write(gnu,*)"set xrange[",xmin,":",xmax,"]"
      write(gnu,*)"set yrange[",ymin,":",ymax,"]"
      write(gnu,*)"set size ratio -1"
      write(gnu,*)"set nokey"
      write(gnu,*)"p \'BD.DAT\' w l"
      write(gnu,*)"pause 5"
      write(gnu,*)"p \'GRID.DAT\' w l"
      write(gnu,*)"pause 5"
      write(gnu,*)"p \'DUAL.DAT\' w l,\'BD.DAT\' w l"
      write(gnu,*)"pause 5"
      close(gnu)

      return
      end

C-----------------------------------------------------------------------------
C.....Check whether ordering of triangle is counter-clockwise
C.....Otherwise correct it
C-----------------------------------------------------------------------------
      subroutine tri_orient(elem, coord)
      implicit none
      include 'param.h'
      integer           elem(3,ntmax)
      double precision  coord(2,npmax)

      integer           cw, ccw, tmp, i, p1, p2, p3
      double precision  dx1, dy1, dx2, dy2, cross

      cw = 0
      ccw= 0

      do i=1,nt
         p1    = elem(1,i)
         p2    = elem(2,i)
         p3    = elem(3,i)

         dx1   = coord(1,p2) - coord(1,p1)
         dy1   = coord(2,p2) - coord(2,p1)

         dx2   = coord(1,p3) - coord(1,p2)
         dy2   = coord(2,p3) - coord(2,p2)

         cross = dx1*dy2 - dx2*dy1

         if(cross .eq. 0.0d0)then
            print*,'Fatal: triangle',i,' is degenerate'
            stop
         endif

         if(cross .lt. 0.0d0)then
            cw        = cw + 1
            tmp       = elem(2,i)
            elem(2,i) = elem(3,i)
            elem(3,i) = tmp
         else
            ccw       = ccw + 1
         endif
      enddo

      write(*,'(" No of cw  triangles       =", i8)') cw
      write(*,'(" No of cww triangles       =", i8)') ccw

      return
      end

C-----------------------------------------------------------------------------
C.....Calculate element and control volume areas for median cell
C-----------------------------------------------------------------------------
      subroutine tri_area(coord, elem, elarea)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax)
      double precision coord(2,npmax), elarea(ntmax)

      double precision dx1, dy1, dx2, dy2
      integer          i, n1, n2, n3

      write(*,'(" Finding triangle areas")')

      maxelarea = 0.0d0
      minelarea = 1.0d8
      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)

c Triangle area
         dx1= coord(1,n2) - coord(1,n1)
         dy1= coord(2,n2) - coord(2,n1)

         dx2= coord(1,n3) - coord(1,n1)
         dy2= coord(2,n3) - coord(2,n1)

         elarea(i) = 0.5d0*( dx1*dy2 - dx2*dy1 )
         maxelarea = dmax1(maxelarea, elarea(i))
         minelarea = dmin1(minelarea, elarea(i))

      enddo

      write(*,'(2x,"Minimum triangle area    =",e12.4)')minelarea
      write(*,'(2x,"Maximum triangle area    =",e12.4)')maxelarea

      return
      end

C-----------------------------------------------------------------------------
C.....Finds elements surrounding a point
C.....Taken from Lohner
C.....esup1 stores the elements
C.....ordering is such that the elements surrounding point ipoin are stored in
C.....locations esup2(ipoin)+1 to esup2(ipoin+1)
C-----------------------------------------------------------------------------
      subroutine el_surr_point(elem, esup1, esup2)
      implicit none
      include 'param.h'
      integer esup1(mesup*npmax), esup2(npmax+1), elem(3,ntmax)
      integer i, ie, inode, ipoi1, ipoin, istor

      do i=1,np+1
            esup2(i) = 0
      enddo

      do ie=1,nt
         do inode=1,3
            ipoi1        = elem(inode, ie) + 1
            esup2(ipoi1) = esup2(ipoi1) + 1
         enddo
      enddo

      do ipoin=2, np+1
         esup2(ipoin) = esup2(ipoin) + esup2(ipoin-1)
      enddo

      do ie=1, nt
         do inode=1,3
            ipoin        = elem(inode, ie)
            istor        = esup2(ipoin) + 1
            esup2(ipoin) = istor
            esup1(istor) = ie
         enddo
      enddo

      do ipoin=np+1, 2, -1
         esup2(ipoin) = esup2(ipoin-1)
      enddo

      esup2(1) = 0

      return
      end

C-----------------------------------------------------------------------------
C.....Finds points surrounding a point
C.....Taken from Lohner
C.....psup1 contains the points
C.....Neighbours of ipoin are between psup2(ipoin)+1 to psup2(ipoin+1)
C-----------------------------------------------------------------------------
      subroutine el_surr_bp(spts, esup1, esup2, esubp)
      implicit none
      include 'param.h'
      integer spts(nspmax), esup1(mesup*npmax), esup2(npmax+1), 
     +        esubp(mesubp,nbpmax)

      integer ip, ipoin, count, ie, iesup

      do ip=1,nsp
         ipoin = spts(ip)
         count = 1
         do iesup=esup2(ipoin)+1, esup2(ipoin+1)
            ie = esup1(iesup)
            esubp(count+1, ip) = ie
            count = count + 1
         enddo
         esubp(1,ip) = count - 1
         if( count .gt. mesubp )then
            print*,'el_surr_bp: Memory error, increase mesubp to'
            print*,'            at least ', count
            stop
         endif
      enddo

      return
      end

C-----------------------------------------------------------------------------
C.....Finds points surrounding a point
C.....Taken from Lohner
C.....psup1 contains the points
C.....Neighbours of ipoin are between psup2(ipoin)+1 to psup2(ipoin+1)
C-----------------------------------------------------------------------------
      subroutine pt_surr_pt(esup1, esup2, elem, psup1, psup2)
      implicit none
      include 'param.h'
      integer esup1(mesup*npmax), esup2(npmax+1)
      integer psup1(mpsup*npmax), psup2(npmax+1), elem(3,ntmax)

      integer ipoin, jpoin, inode, istor, ie, iesup, lpoin(np)

      do ipoin=1,np
         lpoin(ipoin) = 0
      enddo

      psup2(1) = 0
      istor     = 0

      do ipoin=1,np
         do iesup=esup2(ipoin)+1, esup2(ipoin+1)
            ie = esup1(iesup)
            do inode=1,3
               jpoin = elem(inode, ie)
               if(jpoin.ne.ipoin .and. lpoin(jpoin).ne.ipoin) then
                  istor = istor + 1
                  psup1(istor) = jpoin
                  lpoin(jpoin) = ipoin
               endif
            enddo
         enddo
         psup2(ipoin+1) = istor
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Create edges of triangles
C-----------------------------------------------------------------------------
      subroutine create_edge(psup1, psup2, edge)
      implicit none
      include 'param.h'
      integer psup1(mpsup*npmax), psup2(npmax+1), edge(2,nemax)
      integer ipoin, ip, neigh

      ne = 0

      do ipoin = 1,np
         do ip=psup2(ipoin)+1, psup2(ipoin+1)
            neigh = psup1(ip)
            if(neigh.gt.ipoin) then
               ne          = ne + 1
               edge(1,ne) = ipoin
               edge(2,ne) = neigh
            endif
         enddo
      enddo

      write(*,'(" Number of edges           =", i8)') ne

      return
      end

C-----------------------------------------------------------------------------
C Sort edges based on their type
C-----------------------------------------------------------------------------
      subroutine order_edges(ptype, edge)
      implicit none
      include 'param.h'
      integer ptype(npmax), edge(2,nemax)

      integer swedge(2,nbemax), inedge(2,nemax), ffedge(2,nbemax)
      integer n1, n2, i

      write(*,'(" Sorting edges based on their type ...")')

      nsw = 0
      nin = 0
      nff = 0
      nif = 0
      nof = 0

c Put edges into different arrays based on type
      do i=1,ne
         n1 = edge(1,i)
         n2 = edge(2,i)
         if(ptype(n1).eq.solid .and. ptype(n2).eq.solid)then
            nsw = nsw + 1
            swedge(1,nsw) = edge(1,i)
            swedge(2,nsw) = edge(2,i)
         elseif(ptype(n1).eq.farfield .and. ptype(n2).eq.farfield)then
            nff = nff + 1
            ffedge(1,nff) = edge(1,i)
            ffedge(2,nff) = edge(2,i)
         else
            nin = nin + 1
            inedge(1,nin) = edge(1,i)
            inedge(2,nin) = edge(2,i)
         endif
      enddo

c Now put them back into the common array "edge"
c Order is important
      ne = 0

      nin1 = ne + 1
      do i=1,nin
         ne = ne + 1
         edge(1,ne) = inedge(1,i)
         edge(2,ne) = inedge(2,i)
      enddo
      nin2 = ne

      nsw1 = ne + 1
      do i=1,nsw
         ne = ne + 1
         edge(1,ne) = swedge(1,i)
         edge(2,ne) = swedge(2,i)
      enddo
      nsw2 = ne

      nff1 = ne + 1
      do i=1,nff
         ne = ne + 1
         edge(1,ne) = ffedge(1,i)
         edge(2,ne) = ffedge(2,i)
      enddo
      nff2 = ne

      if(nin+nsw+nff .ne. ne)then
         print*,'order_edges: Error -> number of edges dont add up'
         stop
      endif

      write(*,'(2x, "Number of solid    edges =", i8)')nsw
      write(*,'(2x, "Number of interior edges =", i8)')nin
      write(*,'(2x, "Number of farfield edges =", i8)')nff
      write(*,'(2x, "Total no of edges        =", i8)')ne

      return
      end

C-----------------------------------------------------------------------------
C.....For each solid wall point find two edges that share the point
C-----------------------------------------------------------------------------
      subroutine bdedges(spts, ptype, edge, bdedge)
      implicit none
      include 'param.h'
      integer spts(nspmax), ptype(npmax), edge(2,nemax), 
     +        bdedge(2,nbpmax)

      integer i, ip, p1, p2, ed1(npmax), ed2(npmax)

      do i=nsw1,nsw2
         p1 = edge(1,i)
         p2 = edge(2,i)
         if(ptype(p1) .eq. solid)then
            ed2(p1) = i
         endif
         if(ptype(p2) .eq. solid)then
            ed1(p2) = i
         endif
      enddo

      do i=1,nsp
         ip = spts(i)
         bdedge(1,i) = ed1(ip)
         bdedge(2,i) = ed2(ip)
      enddo

      return
      end

C-----------------------------------------------------------------------------
C.....For each edge find the elements to its right and left
C-----------------------------------------------------------------------------
      subroutine el_surr_edge(esup1, esup2, elem, edge, edneigh, vedge)
      implicit none
      include 'param.h'
      integer esup1(mesup*npmax), esup2(npmax+1), edneigh(2,nemax),
     +        vedge(2,nemax), elem(3,ntmax), edge(2,nemax)

      integer i, jt, n1, n2, el, tmp

      do i=1,ne
            edneigh(1,i) = 0
            edneigh(2,i) = 0
            n1 = edge(1,i)
            n2 = edge(2,i)
            do jt=esup2(n1)+1, esup2(n1+1)
                  el = esup1(jt)
                  if( (n1.eq.elem(1,el) .and. n2.eq.elem(2,el)) .or.
     &                (n1.eq.elem(2,el) .and. n2.eq.elem(3,el)) .or.
     &                (n1.eq.elem(3,el) .and. n2.eq.elem(1,el)) )
     &                  edneigh(1,i) = el
                  if( (n2.eq.elem(1,el) .and. n1.eq.elem(2,el)) .or.
     &                (n2.eq.elem(2,el) .and. n1.eq.elem(3,el)) .or.
     &                (n2.eq.elem(3,el) .and. n1.eq.elem(1,el)) )
     &                  edneigh(2,i) = el
            enddo

            if(edneigh(1,i) .eq. 0)then
                  edneigh(1,i) = edneigh(2,i)
                  edneigh(2,i) = 0
                  tmp          = edge(1,i)
                  edge(1,i)    = edge(2,i)
                  edge(2,i)    = tmp
                  if(edneigh(1,i) .eq. 0)then
                     print*,'Fatal: No edge neighbour'
                     stop
                  endif
            endif
      enddo

C For each edge, find vertex of adjacent triangle which is opposite to
C the edge. This is used for limited reconstruction of inviscid fluxes
      do i=1,ne
         vedge(1,i) = 0
         vedge(2,i) = 0
         n1 = edge(1,i)
         n2 = edge(2,i)
         if(edneigh(1,i) .ne. 0)then
            el = edneigh(1,i)
            if(elem(1,el) .ne. n1 .and. 
     +         elem(1,el) .ne. n2) vedge(1,i) = elem(1,el)
            if(elem(2,el) .ne. n1 .and. 
     +         elem(2,el) .ne. n2) vedge(1,i) = elem(2,el)
            if(elem(3,el) .ne. n1 .and. 
     +         elem(3,el) .ne. n2) vedge(1,i) = elem(3,el)
         endif
         
         if(edneigh(2,i) .ne. 0)then
            el = edneigh(2,i)
            if(elem(1,el) .ne. n1 .and. 
     +         elem(1,el) .ne. n2) vedge(2,i) = elem(1,el)
            if(elem(2,el) .ne. n1 .and. 
     +         elem(2,el) .ne. n2) vedge(2,i) = elem(2,el)
            if(elem(3,el) .ne. n1 .and. 
     +         elem(3,el) .ne. n2) vedge(2,i) = elem(3,el)
         endif
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Find triangles surrounding a triangle, ie, face neighbours
C-----------------------------------------------------------------------------
      subroutine el_surr_el(elem, edge, edneigh, esue)
      implicit none
      include 'param.h'
      integer elem(3,ntmax), edge(2,nemax), edneigh(2,nemax),
     +        esue(3,ntmax)

      integer it, ie, e1, e2, c1, c2, p1, p2, p3, q1, q2, q3

      do it=1,nt
         esue(1,it) = 0
         esue(2,it) = 0
         esue(3,it) = 0
      enddo

      do ie=1,ne
         e1 = edge(1,ie)
         e2 = edge(2,ie)

         c1 = edneigh(1,ie)
         c2 = edneigh(2,ie)

         p1 = elem(1,c1)
         p2 = elem(2,c1)
         p3 = elem(3,c1)

         if(c2 .ne. 0)then
            if    (p1 .eq. e1)then
               esue(1,c1) = c2
            elseif(p2 .eq. e1)then
               esue(2,c1) = c2
            elseif(p3 .eq. e1)then
               esue(3,c1) = c2
            else
               print*,'esue: Fatal error'
               stop
            endif

            q1 = elem(1,c2)
            q2 = elem(2,c2)
            q3 = elem(3,c2)
            if    (q1 .eq. e2)then
               esue(1,c2) = c1
            elseif(q2 .eq. e2)then
               esue(2,c2) = c1
            elseif(q3 .eq. e2)then
               esue(3,c2) = c1
            else
               print*,'esue: Fatal error'
               stop
            endif

         endif

      enddo

      return

      do ie=nsw1,nsw2
         e1 = edge(1,ie)
         e2 = edge(2,ie)

         c1 = edneigh(1,ie)
         c2 = edneigh(2,ie)
         if    (p1 .eq. e1)then
            esue(1,c1) = -1
         elseif(p2 .eq. e1)then
            esue(2,c1) = -1
         elseif(p3 .eq. e1)then
            esue(3,c1) = -1
         else
            print*,'esue: Fatal error'
            stop
         endif
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Renumber cells for lusgs
C-----------------------------------------------------------------------------
      subroutine renumber(elem, esue, edneigh, esubp)
      implicit none
      include 'size.h'
      integer elem(3,ntmax), esue(3,ntmax), edneigh(2,nemax), 
     +        esubp(mesubp,nbpmax)

      integer oldnum(ntmax), newnum(ntmax), telem(3,ntmax), 
     +        tesue(3,ntmax)
      integer i, j, it, t1, t2, t3, tcount

      print*,'Renumbering cells for LUSGS...'

      do i=1,nt
         oldnum(i) = 0
         newnum(i) = 0
      enddo

      tcount    = 1
      oldnum(1) = 1
      newnum(1) = 1

      do i=1,nt
         it = oldnum(i)
         if(it .eq. 0)then
            print*,'renumber: Fatal error. it is zero for i=',i
            stop
         endif

         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         if(newnum(t1) .eq. 0 .and. t1 .gt. 0)then
            tcount         = tcount + 1
            oldnum(tcount) = t1
            newnum(t1)     = tcount
         endif

         if(newnum(t2) .eq. 0 .and. t2 .gt. 0)then
            tcount         = tcount + 1
            oldnum(tcount) = t2
            newnum(t2)     = tcount
         endif

         if(newnum(t3) .eq. 0 .and. t3 .gt. 0)then
            tcount         = tcount + 1
            oldnum(tcount) = t3
            newnum(t3)     = tcount
         endif
      enddo

      if(nt .ne. tcount)then
         print*,'renumber: tcount does not match nt.'
         print*,'          Possible bug'
         stop
      endif

c Now rearrange some other element data structure
      do i=1,nt
         telem(1,i) = elem(1,i)
         telem(2,i) = elem(2,i)
         telem(3,i) = elem(3,i)
         tesue(1,i) = esue(1,i)
         tesue(2,i) = esue(2,i)
         tesue(3,i) = esue(3,i)
      enddo

      do i=1,nt
         it        = oldnum(i)
         elem(1,i) = telem(1,it)
         elem(2,i) = telem(2,it)
         elem(3,i) = telem(3,it)
         do j=1,3
            t1 = tesue(j,it)
            if(t1 .gt. 0)then
               esue(j,i) = newnum(t1)
            else
               esue(j,i) = t1
            endif
         enddo
      enddo

      do i=1,ne
         t1 = edneigh(1,i)
         t2 = edneigh(2,i)
         if(t1 .gt. 0)then
            edneigh(1,i) = newnum(t1)
         else
            edneigh(1,i) = t1
         endif
         if(t2 .gt. 0)then
            edneigh(2,i) = newnum(t2)
         else
            edneigh(2,i) = t2
         endif
      enddo

      do i=1,nsp
         do j=1,esubp(1,i)
            it           = esubp(j+1,i)
            esubp(j+1,i) = newnum(it)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------------
C.....Calculate length used for time step
C.....For each triangle find the minimum altitude, or
C.....area/perimeter
C-----------------------------------------------------------------------------
      subroutine dtlength(coord, elarea, elem, drmin)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax)
      double precision coord(2,npmax), elarea(ntmax), drmin(ntmax)

      integer          it, n1, n2, n3
      double precision dx1, dx2, dx3, dy1, dy2, dy3, dr1, dr2, dr3, 
     &                 h1, h2, h3, perim


      do it=1,nt
         n1  = elem(1,it)
         n2  = elem(2,it)
         n3  = elem(3,it)

         dx1 = coord(1,n2) - coord(1,n3)
         dy1 = coord(2,n2) - coord(2,n3)
         dr1 = dsqrt(dx1**2 + dy1**2)
         h1  = 2.0d0*elarea(it)/dr1

         dx2 = coord(1,n3) - coord(1,n1)
         dy2 = coord(2,n3) - coord(2,n1)
         dr2 = dsqrt(dx2**2 + dy2**2)
         h2  = 2.0d0*elarea(it)/dr2

         dx3 = coord(1,n1) - coord(1,n2)
         dy3 = coord(2,n1) - coord(2,n2)
         dr3 = dsqrt(dx3**2 + dy3**2)
         h3  = 2.0d0*elarea(it)/dr3

         perim     = dr1 + dr2 + dr3
         drmin(it) = elarea(it)/perim
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Write grid into a file for visualization with gnuplot
C-----------------------------------------------------------------------------
      subroutine write_grid(coord, edge, edneigh)
      implicit none
      include 'param.h'
      double precision coord(2,npmax)
      integer          edge(2,nemax), edneigh(2,nemax)

      integer          gfile, i, n1, n2

c     Write boundary edges to BD.DAT
      open(unit=10, file='BD.DAT')
      do i=1,ne
      n1 = edge(1,i)
      n2 = edge(2,i)
      if( edneigh(1,i)*edneigh(2,i) .eq. 0)then
         write(10,*)coord(1,n1), coord(2,n1)
         write(10,*)coord(1,n2), coord(2,n2)
         write(10,*)
      endif
      enddo
      close(10)

c     Write grid into file GRID.DAT
      gfile=15
      open(unit=gfile, file='GRID.DAT')

      do i=1,ne
         n1 = edge(1,i)
         n2 = edge(2,i)
         write(gfile,*) coord(1,n1), coord(2,n1)
         write(gfile,*) coord(1,n2), coord(2,n2)
         write(gfile,*)
      enddo

      close(gfile)

c     call system('gnuplot -noraise grid.gnu &')

      return
      end
