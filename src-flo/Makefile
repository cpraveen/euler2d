#
# fortran flags
#
 
FC      = ${F77}
GCC     = gcc
RM      = /bin/rm
FFLAGS  = ${CFLAGS}

TARGETS = flo

#
# program dependencies
#

flo   =   geometric.o \
          common.o \
          main.o \
          vaverage.o \
	       killnormalvel.o \
	       solid_flux.o \
	       farfield_flux.o \
	       farfield2_flux.o \
	       kfvs_flux.o \
	       lax_flux.o \
	       roe_flux.o \
	       hcusp_flux.o \
	       limit.o \
	       result.o \
	       mayavi.o \
	       vigie.o \
	       clcd.o \
	       avgfact.o \
	       avgfact2.o \
	       lusgs.o \
	       smooth.o \
			 math.o \
			 isocont.o

#
# programs
#

ALL: $(TARGETS)

flo: ${flo}
	  ${FC} -o flo ${flo}

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

#Flux subroutines: They must be named as somename_flux.F
#Example: roe_flux.F
%_flux.o:  %_flux.F
		${GCC} -DDEBUG -E -C -P $*_flux.F > $*_flux.f
		${FC} ${FFLAGS} -c $*_flux.f
		${RM} -f $*_flux.f

##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o $(TARGETS)
