#This requires a fairly recent version of make, preferably GNU's make.
# gcc is also required for some pre-processing actions. 
#The variables F77 and CC are assumed to be set to the fortran and C compilers. 
#Otherwise set them here
#F77=
#CC=

#Location of flow solver files
SRC=../src-flo

#Location of mesh adjoint files
SRCADJ=../src-mesh

# Location of tapenade executable
# TAPENADE_HOME must be defined, probably in your bashrc
# Else define it first
#TAPENADE_HOME=
GenLib= ${SRCADJ}/MyGeneralLib
ADLib = ${SRCADJ}/MyADLib
TPN   = ${TAPENADE_HOME}/bin/tapenade -b -O ${SRCADJ}
ADSTACK=../tapenade/adStack.c

#Used by tapenade to change to source file directory
CDSRC = cd ${SRC}

#
# fortran flags
#
 
FC     = ${F77}
GCC    = gcc
RM     = /bin/rm
FFLAGS = ${CFLAGS} -I${SRC}

TARGETS =  deform adjoint

#
# program dependencies
#

pri     = deform.o \
          relaxfact.o \
			 grid_residue.o \
          ${SRC}/geometric.o \
			 ${SRC}/avgfact.o \
			 ${SRC}/avgfact2.o \
			 ${SRC}/common.o

adj     = adjoint.o \
          relaxfact.o \
			 grid_residue_bx.o \
          ${SRC}/geometric.o \
			 ${SRC}/avgfact.o \
			 ${SRC}/avgfact2.o \
			 ${SRC}/common.o

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

${SRC}/%.o: ${SRC}/%.f
		${CDSRC}; ${FC} ${FFLAGS} -c $*.f

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
