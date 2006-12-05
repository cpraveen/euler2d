      include 'common.h'
      include 'inf.h'
      include 'size.h'
      include 'visc.h'

      integer istart, scratch, restart
      parameter(scratch=1, restart=2)
      common/starttype/istart

C Range of bounding box
      double precision xmin, xmax, ymin, ymax
      common/range/xmin, xmax, ymin, ymax

c Type of grid
c any other value implies hybrid grid
      character gridfile*32, inpfile*32
      common/files/gridfile,inpfile

      double precision CFL, MINRES, dtglobal
      integer iter, ITERLAST, MAXITER, saveinterval, explicit
      common/itparam/CFL,MINRES,dtglobal,iter,ITERLAST,
     &               MAXITER,saveinterval, explicit

      double precision airk(3), birk(3)
      integer NIRK
      common/timeintegration/airk,birk,NIRK

C     Number of contours
      integer niso
      common/contour/niso

C Range of some variables
C rmin,rmax = density
C pmin,pmax = pressure
C mmin,mmax = mach number
      double precision rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, 
     &                 mmin, mmax, emin, emax, nmin, nmax
      common/minmaxprim/rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, 
     &                  mmin, mmax, emin, emax, nmin, nmax

      double precision fres, fres1, fresi
      integer          iresi
      common/resparam/fres,fres1,fresi,iresi

      double precision wd1(nspmax), wdmin, wdmax
      common/wdparam/wd1, wdmin, wdmax

C     Define tags for point types
      integer interior, solid, farfield, outflow
      parameter(interior=0)
      parameter(solid=3)
      parameter(outflow=4)
      parameter(farfield=5)

C Small number
      double precision EPSILON
      parameter(EPSILON=1.0d-16)

C Limiter factor for MUSCL
      integer          GRADTYPE, ILIMIT
      double precision LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22
      common/lim/LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22, 
     &           GRADTYPE, ILIMIT

C Size of connectivity list; required by mayavi
      integer lsize
      common/maya/lsize

      integer iflux, iroe, ikfvs, ihllc
      parameter(iroe=1, ikfvs=2, ihllc=3)
      common/flux/iflux

      integer inviscid, laminar, turbulent
      parameter(inviscid=1, laminar=2, turbulent=3)

      integer xvstatus
      common/xv/xvstatus

      double precision minelarea, maxelarea, mincvarea, maxcvarea
      double precision minflen, maxflen
      common/minmaxarea/minelarea, maxelarea, mincvarea, maxcvarea,
     &                  minflen, maxflen
