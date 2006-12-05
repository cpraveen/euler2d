C-----------------------------------------------------------------------------
C.....Definition of some constants
C-----------------------------------------------------------------------------
      subroutine math
      implicit none
      include 'param.h'

      GAMMA        = 1.4d0
      GAMMA1       = GAMMA-1.0d0
      GAS_CONST    = 1.0d0
      M_PI         = 4.0d0*datan2(1.0d0, 1.0d0)

      prandtl      = 0.72d0
      prandtl_turb = 0.9d0

      ALBADA11     = 2.0d0/3.0d0
      ALBADA12     = 1.0d0 - ALBADA11

      ALBADA21     = 4.0d0/3.0d0
      ALBADA22     = 1.0d0 - ALBADA21

      xvstatus     = no

c Constants for Spalart-Allmaras model
      Cb1          = 0.1355d0
      Cb2          = 0.622d0
      sigma_sa     = 2.0d0/3.0d0
      kolm         = 0.41d0
      Cw1          = Cb1/kolm**2 + (1.0d0 + Cb2)/sigma_sa
      Cw2          = 0.3d0
      Cw3          = 2.0d0
      Cv1          = 7.1d0
      Cv2          = 5.0d0

      Cv11         = Cv1**3
      Cw31         = 1.0d0 + Cw3**6
      Cw32         = Cw3**6
      kolm2        = kolm**2
      Cb2Sig1      = (1.0d0 + Cb2)/sigma_sa
      Cb2Sig2      = Cb2/sigma_sa

      return
      end

C-----------------------------------------------------------------------------
C.....Read parameters from an input file and set freestream values
C-----------------------------------------------------------------------------
      subroutine read_input
      implicit none
      include 'param.h'
      integer inp, iargc, n, inpstatus
      character sdummy*32

      n = iargc()
      if(n .eq. 0)then
            print*,'You must specify an input file.'
            stop
      endif

      call getarg(1,inpfile)

      inp = 11
      open(unit=inp, file=inpfile, status='old')
      print*,'Reading parameters from ',inpfile
      read(inp,*)sdummy, istart
      read(inp,*)sdummy, iflow
      read(inp,*)sdummy, mach_inf
      read(inp,*)sdummy, aoa_deg
      read(inp,*)sdummy, Rey
      read(inp,*)sdummy, cfl
      read(inp,*)sdummy, explicit
      read(inp,*)sdummy, iterlast
      read(inp,*)sdummy, maxiter
      read(inp,*)sdummy, minres
      read(inp,*)sdummy, saveinterval
      read(inp,*)sdummy, niso
      read(inp,*)sdummy, iflux
      read(inp,*)sdummy, ILIMIT
      read(inp,*)sdummy, vortex, xref, yref
      read(inp,*)sdummy, gridfile
      close(inp)

      inpstatus = yes

      if(istart .ne. scratch .and. istart .ne. restart)then
            print*,'Unknown start option',istart
            print*,'Possible values: 1=scratch or 2=restart'
            inpstatus = no
      endif

      if(iflow .ne. inviscid .and. iflow .ne. laminar .and.
     &   iflow .ne. turbulent)then
            print*,'Unknown flow type',iflow
            print*,'Possible values: 1=inviscid, 2=laminar, 3=turbulent'
            inpstatus = no
      endif

      if(iflux .ne. iroe .and. iflux .ne. ikfvs .and. 
     &   iflux .ne. ihllc)then
            print*,'Unknown flux',iflux
            print*,'Possible values: 1=roe, 2=kfvs'
            inpstatus = no
      endif

      if(ilimit .ne. no .and. ilimit .ne. yes)then
            print*,'Unknown limiter option',ilimit
            print*,'Possible values: 0=no, 1=yes'
            inpstatus = no
      endif

      if(vortex .ne. yes .and. vortex .ne. no)then
            print*,'Unknown vortex option',vortex
            print*,'Possible values: 0=no, 1=yes'
            inpstatus = no
      endif

      if(inpstatus .eq. no) stop

      return
      end

C-----------------------------------------------------------------------------
C.....Variables stored are primitive - density, u, v, pressure
C.....Initialize primitive variables to free stream values
C-----------------------------------------------------------------------------
      subroutine initialize(qc, cl, cd)
      implicit none
      include 'param.h'
      double precision qc(nvar, ntmax), cl, cd
      
      integer          j

      q_inf   = 1.0d0
      r_inf   = 1.0d0
      p_inf   = 1.0d0/(GAMMA*mach_inf**2)
      T_inf   = p_inf/(r_inf*GAS_CONST)
      aoa     = aoa_deg*M_PI/180.0d0
      u_inf   = q_inf*dcos(aoa)
      v_inf   = q_inf*dsin(aoa)
      ent_inf = p_inf/r_inf**GAMMA
      a_inf   = dsqrt(GAMMA*p_inf/r_inf)
      H_inf   = a_inf**2/GAMMA1 + 0.5d0*q_inf

c Required by Sutherland law
      T_infd  = 300.0d0
      SCONST  = 110.4d0*T_inf/T_infd

c Primitive variables in free-stream
      prim_inf(1) = r_inf
      prim_inf(2) = u_inf
      prim_inf(3) = v_inf
      prim_inf(4) = p_inf

      if(iflow.eq.inviscid) print*,'Euler computation'
      if(iflow.eq.laminar)  print*,'Laminar Navier-Stokes computation'
      if(iflow.eq.turbulent)print*,'Turbulent Navier-Stokes computation'
      print*,'Free-stream values:'
      write(*,'(5x, " Mach number =", f8.4)')mach_inf
      write(*,'(5x, " AOA         =", f8.4)')aoa_deg
      write(*,'(5x, " u velocity  =", f8.4)')u_inf
      write(*,'(5x, " v velocity  =", f8.4)')v_inf
      write(*,'(5x, " Pressure    =", f8.4)')p_inf

      if(vortex .eq. yes)then
            print*,'Using point vortex correction for far-field points'
            write(*,'(" Vortex center =", 2e12.4)') xref, yref
      endif

C Runge-Kutta time stepping
      NIRK    = 3
      airk(1) = 0.0d0
      airk(2) = 3.0d0/4.0d0
      airk(3) = 1.0d0/3.0d0
      birk(1) = 1.0d0
      birk(2) = 1.0d0/4.0d0
      birk(3) = 2.0d0/3.0d0

      if(explicit .eq. no) NIRK = 1

      cl = 0.0d0
      cd = 0.0d0

      if(istart .eq. scratch)then
         print*,'Initializing solution to free stream values'
         do j=1,nt
            qc(1,j) = prim_inf(1)
            qc(2,j) = prim_inf(1)*prim_inf(2)
            qc(3,j) = prim_inf(1)*prim_inf(3)
            qc(4,j) = prim_inf(4)/gamma1 + 0.5d0*prim_inf(1)*(
     +                prim_inf(2)**2 + prim_inf(3)**2 )
         enddo

      else
         print*,'Initializing solution to old values'
         call read_flow(qc)
      endif

      return
      end

C-----------------------------------------------------------------------------
C Save flow solution into a file. Conserved variables are saved.
C-----------------------------------------------------------------------------
      subroutine save_flow(qc)
      implicit none
      include 'common.h'
      include 'size.h'
      double precision qc(nvar,ntmax)

      integer          i, j

      open(unit=50, file='FLO.DAT')
      do i=1,nt
         write(50,'(4e20.10)') (qc(j,i), j=1,nvar)
      enddo
      close(50)

      return
      end

C-----------------------------------------------------------------------------
C Read flow solution from file
C-----------------------------------------------------------------------------
      subroutine read_flow(qc)
      implicit none
      include 'common.h'
      include 'size.h'
      double precision qc(nvar,ntmax)

      integer          i, j

      open(unit=50, file='FLO.DAT', status='OLD')
      do i=1,nt
         read(50,'(4e20.10)') (qc(j,i), j=1,nvar)
      enddo
      close(50)

      return
      end

C-----------------------------------------------------------------------------
C Time-step from cfl condition
C-----------------------------------------------------------------------------
      subroutine time_step(drmin, q, dt)
      implicit none
      include 'param.h'
      double precision drmin(ntmax), q(nvar,ntmax), dt(ntmax)

      integer          i
      double precision d, u, v, p, s, a

      dtglobal = 1.0d20
      do i=1,nt
         call contoprim(q(1,i), d, u, v, p)
         s = dsqrt(u*u + v*v)
         a = dsqrt(gamma*p/d)
         dt(i) = cfl*drmin(i)/(s + a)
         dtglobal = dmin1(dtglobal, dt(i))
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Time-step from cfl condition
C-----------------------------------------------------------------------------
      subroutine time_step2(edge, tedge, tarea, coord, q, dt)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision tarea(ntmax), coord(2,npmax), q(nvar,ntmax), 
     +                 dt(ntmax)

      integer          i, j, t1, t2, e1, e2
      double precision sumspeed(nt), qa(nvar), d, u, v, p, a, nx, ny, 
     +                 nl, un, ll, lr

      do i=1,nt
         sumspeed(i) = 0.0d0
      enddo

      do i=1,ne
         t1 = tedge(1,i)
         t2 = tedge(2,i)
         e1 =  edge(1,i)
         e2 =  edge(2,i)
         if(t2 .eq. 0)then
            do j=1,nvar
               qa(j) = q(j,t1)
            enddo
         else
            do j=1,nvar
               qa(j) = 0.5d0*( q(j,t1) + q(j,t2) )
            enddo
         endif
         call contoprim(qa, d, u, v, p)

         nx =  ( coord(2,e2) - coord(2,e1) )
         ny = -( coord(1,e2) - coord(1,e1) )
         nl =  dsqrt(nx*nx + ny*ny)

         un = u*nx + v*ny
         a  = dsqrt(gamma*p/d)*nl
         ll = dmax1(0.0d0,  un + a)
         lr = dmax1(0.0d0, -un + a)

         sumspeed(t1) = sumspeed(t1) + ll
         if(t2 .ne. 0) sumspeed(t2) = sumspeed(t2) + lr
      enddo

      dtglobal = 1.0d20
      do i=1,nt
         dt(i)    = cfl*tarea(i)/sumspeed(i)
         dtglobal = dmin1(dtglobal, dt(i))
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Save old solution
C-----------------------------------------------------------------------------
      subroutine save_old(q, qold)
      implicit none
      include 'param.h'
      double precision q(nvar,ntmax), qold(nvar,ntmax)

      integer          i, j

      do i=1,nt
         do j=1,nvar
            qold(j,i) = q(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------------
C convert conserved variable to primitive variable
C-----------------------------------------------------------------------------
      subroutine contoprim(q, d, u, v, p)
      implicit none
      include 'param.h'
      double precision d, u, v, p, q(nvar)

      d = q(1)
      u = q(2)/d
      v = q(3)/d
      p = gamma1*( q(4) - 0.5d0*d*(u*u + v*v) )

      return
      end

C-----------------------------------------------------------------------------
C L2 and Linf norm of the finite volume residual
C-----------------------------------------------------------------------------
      subroutine residue(res)
      implicit none
      include 'param.h'
      double precision res(nvar,ntmax)

      integer          i, j
      double precision fr

      fres  = 0.0d0
      fresi = 0.0d0
      iresi = 0

      do i=1,nt
         fr = 0.0d0
         do j=1,nvar
            fr = fr + res(j,i)**2
         enddo
         fres = fres + fr
         fr   = dsqrt(fr)
         if(fr .gt. fresi)then
            fresi = fr
            iresi = i
         endif
      enddo

      fres = dsqrt(fres/nt)

      if(iter .eq. 1)then
         fres1 = fres
         print*,'Residue in first iteration =',fres1
      endif

      if(fres1 .ne. 0.0d0) fres = fres/fres1

      return
      end
