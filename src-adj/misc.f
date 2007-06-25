C-----------------------------------------------------------------------------
C.....Read parameters from an input file and set freestream values
C-----------------------------------------------------------------------------
      subroutine read_input_adj
      implicit none
      include 'param.h'
      include 'param2.h'
      integer   i, inp, iargc, n, inpstatus
      character sdummy*32, adjfile*32

      n = iargc()
      if(n .ne. 2)then
            print*,'You must specify flow and adjoint input files'
            stop
      endif

      call getarg(2,adjfile)

      inp = 11
      open(unit=inp, file=adjfile, status='old')
      print*,'Reading parameters from ',adjfile
      read(inp,*)sdummy, istart2
      read(inp,*)sdummy, cfl
      read(inp,*)sdummy, explicit2
      read(inp,*)sdummy, maxiter
      read(inp,*)sdummy, minres
      read(inp,*)sdummy, saveinterval
      read(inp,*)sdummy, niso
      read(inp,*)sdummy, iflux
      read(inp,*)sdummy, costtype
      close(inp)

      inpstatus = yes

      if(istart2 .ne. scratch .and. istart2 .ne. restart)then
            print*,'Unknown start option',istart2
            print*,'Possible values: 1=scratch or 2=restart'
            inpstatus = no
      endif

      if(iflux .ne. iroe .and. iflux .ne. ikfvs .and. 
     &   iflux .ne. ihllc)then
            print*,'Unknown flux',iflux
            print*,'Possible values: 1=roe, 2=kfvs'
            inpstatus = no
      endif

      if(inpstatus .eq. no) stop

      if(explicit2 .eq. no)then
         nirk = 1
         print*,'Implicit time timestepping...'
      else
         nirk = 3
         print*,'Explicit time timestepping...'
      endif

c     for pressure matching problem, read target cp
      if(costtype.eq.1)then
         print*,'Reading target pressure from cp0.dat'
         inp =12
         open(unit=inp, file='cp0.dat')
         do i=1,nsp
            read(inp,*) xcp0(i), cp0(i)
         enddo
         close(inp)
      endif

      return
      end

C-----------------------------------------------------------------------------
C Initialize adjoint solution to zero or from previous solution
C-----------------------------------------------------------------------------
      subroutine init_adj(q)
      implicit none
      include 'param.h'
      include 'param2.h'
      double precision q(nvar,ntmax)

      integer          i, j

      if(istart2 .eq. scratch)then
         print*,'Setting adjoint solution to zero ...'
         do i=1,nt
            do j=1,nvar
               q(j,i) = 0.0d0
            enddo
         enddo
      else
         print*,'Setting adjoint solution from file ...'
         call read_adj(q)
      endif

      return
      end
