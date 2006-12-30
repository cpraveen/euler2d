#Include common variables
include ../Makefile.in

#Tapenade command
TPN   = ${TAPENADE_HOME}/bin/tapenade -b -O ${SRCMSH}

#
# fortran flags
#
FFLAGS = ${CFLAGS} -I${SRCFLO}

TARGETS =  deform adjoint

#
# program dependencies
#

pri     = deform.o \
          relaxfact.o \
			 grid_residue.o \
          ${SRCFLO}/geometric.o \
			 ${SRCFLO}/avgfact.o \
			 ${SRCFLO}/avgfact2.o \
			 ${SRCFLO}/common.o

adj     = adjoint.o \
          relaxfact.o \
			 grid_residue_bx.o \
          ${SRCFLO}/geometric.o \
			 ${SRCFLO}/avgfact.o \
			 ${SRCFLO}/avgfact2.o \
			 ${SRCFLO}/common.o

#
# programs
#

ALL:		$(TARGETS)

deform:	${pri}
		   ${FC} -o deform ${pri}

adjoint:	${adj}
		   ${FC} -o adjoint ${adj}

# Compile only, dont run Tapenade
only:
		touch *_bq.f *_bx.f
		make adj

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

${SRCFLO}/%.o: ${SRCFLO}/%.f
		${CDFLO}; ${FC} ${FFLAGS} -c $*.f

#This is to prevent these files from being automatically deleted by make
.SECONDARY: grid_residue_bx.f

##############################################################################
# Adjoint versions of nonlinear routines
##############################################################################
grid_residue_bx.f:	grid_residue.f
		${TPN} \
		       -head         grid_residue   \
		       -vars         "dx1 dx2 gres1 gres2"  \
		       -outvars      "dx1 dx2 gres1 gres2"  \
		       -difffuncname "_bx"      \
		        grid_residue.f;
#-----------------------------------------------------------------------------
${TPDIR}/adStack.o: ${ADSTACK}
		${CC} ${CFLAGS} -c ${ADSTACK}
#-----------------------------------------------------------------------------
${TPDIR}/adBuffer.o: ${ADBUFFER}
		${FC} ${FFLAGS} -c ${ADBUFFER}
##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o *.msg *~ $(TARGETS)

cleanall:	
	${RM} -f *.o *.msg *~ *_bx.f $(TARGETS)
