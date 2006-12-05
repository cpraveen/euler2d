      subroutine result(coord, spts, elem, edge, bdedge, carea, af, 
     &                      qc, qv, coordb, cl, cd)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), spts(nspmax),
     &                 bdedge(2,nbpmax)
      double precision coord(2,npmax), carea(ntmax), af(3,npmax),
     &                 qc(nvar,ntmax), qv(nvar,npmax),
     &                 coordb(2,npmax), cl, cd

      integer          ifile, ii, is, it, i, j, e1, e2, v1, v2, v3
      double precision ft1(3), ft2(3), ft3(3), ft4(3),
     &                 val1(100), val2(100), val3(100), val4(100),
     &                 delro1, delro2, delro3, delro4, coort(2,3), 
     &                 deltat, xx0, yy0, xx1, yy1,
     &                 u, v
 
C Compute vertex average values of adjoint variable
      do i=1,np
         do j=1,nvar
            qv(j,i) = 0.0d0
         enddo
      enddo
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call vaverage(coord(1,v1), coord(1,v2), coord(1,v3),
     +                 af(1,v1), af(1,v2), af(1,v3), carea(i),
     +                 qc(1,i), qv(1,v1), qv(1,v2), qv(1,v3))
      enddo
      do i=1,np
         do j=1,nvar
            qv(j,i) = qv(j,i)/af(3,i)
         enddo
      enddo
      do i=1,nsp
         j = spts(i)
         e1= bdedge(1,i)
         e2= bdedge(2,i)
         v1= edge(1,e1)
         v2= j
         v3= edge(2,e2)
         call killnormalvel(coord(1,v1), coord(1,v2), coord(1,v3),
     +                      qv(1,j))
      enddo

C Compute range of adjoint solution
      rmin =  1.0d10
      rmax = -1.0d10
      umin =  1.0d10
      umax = -1.0d10
      vmin =  1.0d10
      vmax = -1.0d10
      emin =  1.0d10
      emax = -1.0d10
      do i=1,np
          rmin = dmin1(rmin, qv(1,i))
          rmax = dmax1(rmax, qv(1,i))
          umin = dmin1(umin, qv(2,i))
          umax = dmax1(umax, qv(2,i))
          vmin = dmin1(vmin, qv(3,i))
          vmax = dmax1(vmax, qv(3,i))
          emin = dmin1(emin, qv(4,i))
          emax = dmax1(emax, qv(4,i))
      enddo

      write(*,9)('-',i=1,70)
9     format(70a)
      write(*,10)mach_inf,aoa_deg,Rey,CFL
10    format(' Mach =',f6.3,', AOA =',f6.2, ', Rey = ',e10.4,
     &       ', CFL =',f8.2)
      write(*,11)iflux,ilimit,gridfile
11    format(' Flux =',i2, ',     Lim = ',i2,',    Grid= ',a30)
      write(*,9)('-',i=1,70)
      write(*,'(" Iterations        =",i12)')iter
      write(*,'(" Global dt         =",e16.6)') dtglobal
      write(*,'(" L2 residue        =",e16.6)')fres
      write(*,'(" Linf residue      =",e16.6)')fresi
      write(*,'(" Linf triangle     =",i12)') iresi
      write(*,*)
      write(*,'(" Cl, Cd            =",2f12.6)')cl,cd
      write(*,*)
      write(*,'(27x,"Min",8x,"Max")')          
      write(*,'(" Density           =",2f12.6)')rmin, rmax
      write(*,'(" u momentum        =",2f12.6)')umin, umax
      write(*,'(" v momentum        =",2f12.6)')vmin, vmax
      write(*,'(" Energy            =",2f12.6)')emin, emax
      call flush(6)

c iso-pressure/mach/visc for gnuplot
      open(22,file='ADJ.D',status='unknown')
      open(23,file='ADJ.U',status='unknown')
      open(24,file='ADJ.V',status='unknown')
      open(25,file='ADJ.E',status='unknown')
      rewind(22)
      rewind(23)
      rewind(24)
      rewind(25)
         
      delro1 = (rmax-rmin)/niso
      delro2 = (umax-umin)/niso
      delro3 = (vmax-vmin)/niso
      delro4 = (emax-emin)/niso
      do ii=1,niso+1
         val1(ii) = rmin       + (ii-1)*delro1
         val2(ii) = umin       + (ii-1)*delro2
         val3(ii) = vmin       + (ii-1)*delro3
         val4(ii) = emin       + (ii-1)*delro4
      enddo

      do it=1,nt
         do i=1,3
            is         = elem(i,it)
            coort(1,i) = coord(1,is)
            coort(2,i) = coord(2,is)
            ft1(i)     = qv(1,is)
            ft2(i)     = qv(2,is)
            ft3(i)     = qv(3,is)
            ft4(i)     = qv(4,is)
         enddo
         call isocont(22,ft1,coort,niso,val1)
         call isocont(23,ft2,coort,niso,val2)
         call isocont(24,ft3,coort,niso,val3)
         call isocont(25,ft4,coort,niso,val4)
      enddo
      close(22)
      close(23)
      close(24)
      close(25)

cc velocity vector for gnuplot
      ifile  = 26
      open(ifile,file='ADJ.VECT',status='unknown')
      rewind(ifile)
      deltat = 1.0d-2
      do is=1,np
         xx0 = coord(1,is)
         yy0 = coord(2,is)
         u   = qv(2,is)
         v   = qv(3,is)
         xx1 = u*deltat
         yy1 = v*deltat
         write(ifile,*) xx0, yy0, xx1, yy1
      enddo
      close(ifile)

c Run gnuplot and Start xv if not already started
      if(xvstatus .eq. no .and. iterlast .eq. 0)then
         call system('gnuplot adj.gnu')
         call system('xv -poll adj.png &')
         xvstatus = yes
      else
         call system('gnuplot adj.gnu &')
      endif

      call vigie(coord, elem, qv, coordb)
      call mayavi(coord, elem, qv, coordb)
      call save_adj(qc)

      return
      end

C-----------------------------------------------------------------------------
C Save adjoint solution into a file.
C-----------------------------------------------------------------------------
      subroutine save_adj(qc)
      implicit none
      include 'common.h'
      include 'size.h'
      double precision qc(nvar,ntmax)

      integer          i, j

      open(unit=50, file='ADJ.DAT')
      do i=1,nt
         write(50,'(4e20.10)') (qc(j,i), j=1,nvar)
      enddo
      close(50)

      return
      end

C-----------------------------------------------------------------------------
C Read adjoint solution from file
C-----------------------------------------------------------------------------
      subroutine read_adj(qc)
      implicit none
      include 'common.h'
      include 'size.h'
      double precision qc(nvar,ntmax)

      integer          i, j

      open(unit=50, file='ADJ.DAT', status='OLD')
      do i=1,nt
         read(50,'(4e20.10)') (qc(j,i), j=1,nvar)
      enddo
      close(50)

      return
      end

