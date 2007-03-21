      subroutine write_result(coord, elem, edge, qc, qv, cl, cd)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax)
      double precision coord(2,npmax), qc(nvar,ntmax), qv(nvar,npmax),
     &                 cl, cd

      integer          ifile, ii, is, it, i
      double precision ft1(3), ft2(3), coort(2,3), val1(100), 
     &                 val2(100), delro1, delro2,
     &                 deltat, xx0, yy0, xx1, yy1, q2, mach, cp, ent, 
     &                 r, u, v, p, prim(nvar,npmax)
 
      do i=1,np
         call contoprim(qv(1,i), r, u, v, p)
         prim(1,i) = r
         prim(2,i) = u
         prim(3,i) = v
         prim(4,i) = p
      enddo

      call vigie(coord, elem, prim)
      call mayavi(coord, elem, prim)
      call save_flow(qc)

      pmin =  1.0d10
      pmax = -1.0d10
      mmin =  1.0d10
      mmax = -1.0d10
      rmin =  1.0d10
      rmax = -1.0d10
      umin =  1.0d10
      umax = -1.0d10
      vmin =  1.0d10
      vmax = -1.0d10
      emin =  1.0d10
      emax = -1.0d10
      do is=1,nt
          r    = qc(1,is)
          u    = qc(2,is)/r
          v    = qc(3,is)/r
          q2   = u**2 + v**2
          p    = gamma1*( qc(4,is) - 0.5d0*r*q2 )
          rmin = dmin1(rmin, r) 
          rmax = dmax1(rmax, r) 
          umin = dmin1(umin, u) 
          umax = dmax1(umax, u) 
          vmin = dmin1(vmin, v) 
          vmax = dmax1(vmax, v) 
          pmin = dmin1(pmin, p) 
          pmax = dmax1(pmax, p) 
          mach = dsqrt(q2*r/(GAMMA*p))
          mmin = dmin1(mmin, mach)
          mmax = dmax1(mmax, mach)
          ent  = dlog10(p/r**GAMMA/ent_inf)
          emin = dmin1(emin, ent)
          emax = dmax1(emax, ent)
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
      write(*,'(" Pressure          =",2f12.6)')pmin, pmax
      write(*,'(" Mach number       =",2f12.6)')mmin, mmax
      write(*,'(" x velocity        =",2f12.6)')umin, umax
      write(*,'(" y velocity        =",2f12.6)')vmin, vmax
      write(*,'(" Entropy           =",2f12.6)')emin, emax
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      call flush(6)

C Write pressure coefficient
      open(unit=10, file='WALL.DAT')
      do i=nsw1,nsw2
         is = edge(1,i)
         r  = prim(1,is)
         u  = prim(2,is)
         v  = prim(3,is)
         q2 = u**2 + v**2
         p  = prim(4,is)
         cp = -(p - p_inf)/(0.5d0*r_inf*q_inf**2)
         ent= p/r**gamma/ent_inf - 1.0d0
         write(10,'(3e16.6)') coord(1,is), cp, ent
      enddo
      close(10)

c iso-pressure/mach/visc for gnuplot
      open(22,file='FLO.P',status='unknown')
      open(23,file='FLO.M',status='unknown')
      rewind(22)
      rewind(23)
         
      delro1 = (pmax-pmin)/niso
      delro2 = (mmax-mmin)/niso
      do ii=1,niso+1
         val1(ii) = pmin       + (ii-1)*delro1
         val2(ii) = mmin       + (ii-1)*delro2
      enddo

      do it=1,nt
         do i=1,3
            is         = elem(i,it)
            coort(1,i) = coord(1,is)
            coort(2,i) = coord(2,is)
            r          = prim(1,is)
            u          = prim(2,is)
            v          = prim(3,is)
            q2         = u**2 + v**2
            p          = prim(4,is)
            ft1(i)     = p
            ft2(i)     = dsqrt(q2*r/(GAMMA*p))
         enddo
         call isocont(22,ft1,coort,niso,val1)
         call isocont(23,ft2,coort,niso,val2)
      enddo
      close(22)
      close(23)

cc velocity vector for gnuplot
      ifile  = 26
      open(ifile,file='FLO.VECT',status='unknown')
      rewind(ifile)
      deltat = 1.0d-2
      do is=1,np
         xx0 = coord(1,is)
         yy0 = coord(2,is)
         u   = prim(2,is)
         v   = prim(3,is)
         xx1 = u*deltat
         yy1 = v*deltat
         write(ifile,*) xx0, yy0, xx1, yy1
      enddo
      close(ifile)

c Run gnuplot and Start xv if not already started
      if(xvstatus .eq. no .and. iterlast .eq. 0)then
         call system('gnuplot flo.gnu')
         call system('xv -poll flo.png &')
         xvstatus = yes
      else
         call system('gnuplot flo.gnu &')
      endif

      return
      end

c------------------------------------------------------------------------------
c     Write some solution values: Used for optimization
c------------------------------------------------------------------------------
      subroutine write_sol(iter, residue, cl, cd)
      implicit none
      integer          iter
      double precision residue, cl, cd

      integer          fid

      fid = 20
      open(fid, file="FLO.OUT")
      write(fid,*) iter, residue, cl, cd
      close(fid)

      return
      end
