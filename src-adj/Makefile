#This requires a fairly recent version of make, preferably GNU's make.
# gcc is also required for some pre-processing actions. 
#The variables F77 and CC are assumed to be set to the fortran and C compilers. 
#Otherwise set them here
#F77=
#CC=

#Location of flow solver files
SRC=../src-flo

#Location of adjoint files
SRCADJ=../src-adj

#Which cost function to use
#COST1 depends only on cell-center values
#COST2 depends only on vertex values
COST=COST1

# Location of tapenade executable
# TAPENADE_HOME must be defined, probably in your bashrc
# Else define it first
#TAPENADE_HOME=
GenLib= ${SRCADJ}/MyGeneralLib
ADLib = ${SRCADJ}/MyADLib
TPN   = ${TAPENADE_HOME}/bin/tapenade \
        -ext ${GenLib} -extAD ${ADLib} -b -O ${SRCADJ}
ADSTACK= ${TAPENADE_HOME}/lib/adStack.c
ADBUFFER= ${TAPENADE_HOME}/lib/adBuffer.f

#Used by tapenade to change to source file directory
CDSRC = cd ${SRC}

#
# fortran flags
#
 
FC     = ${F77}
GCC    = gcc
RM     = /bin/rm
FFLAGS = ${CFLAGS} -I${SRC}

TARGETS =  adj

#
# program dependencies
#

adj     = main.o \
			 ${SRC}/avgfact.o \
			 ${SRC}/avgfact2.o \
			 ${SRC}/clcd.o \
			 ${SRC}/common.o \
          ${SRC}/geometric.o \
			 ${SRC}/isocont.o \
			 ${SRC}/killnormalvel.o \
			 ${SRC}/limit.o \
          ${SRC}/math.o \
			 ${SRC}/vaverage.o \
			 avgfact_x.o \
			 afact1_bx.o \
			 afact2_bx.o \
			 afact3_bx.o \
			 afact4_bx.o \
			 bd_normal_bx.o \
			 cost_q.o \
			 cost_x.o \
			 costfunc_bq.o \
			 costfunc_bx.o \
			 farfield_flux_bq.o \
			 farfield_flux_bx.o \
			 kfvs_flux_bq.o \
			 kfvs_flux_bx.o \
			 killnormalvel_bq.o \
			 killnormalvel_bx.o \
			 lax_flux_bq.o \
			 lax_flux_bx.o \
			 limit_bq.o \
			 limit_bx.o \
			 lusgs.o \
			 mayavi.o \
			 misc.o \
			 result.o \
          roe_flux_bq.o \
          roe_flux_bx.o \
			 solid_flux_bq.o \
			 solid_flux_bx.o \
			 vaverage_bq.o \
			 vaverage_bx.o \
			 vigie.o \
			 adStack.o \
			 adBuffer.o

#
# programs
#

ALL:		$(TARGETS)

adj:	${adj}
		${FC} -o adj ${adj}

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
.SECONDARY: roe_flux_bq.f   roe_flux_bx.f \
            kfvs_flux_bq.f  kfvs_flux_bx.f \
            hcusp_flux_bq.f hcusp_flux_bx.f \
            lax_flux_bq.f lax_flux_bx.f

##############################################################################
# Adjoint versions of nonlinear routines
##############################################################################
# Differentiation of flux routines
# This takes care of all the flux routines
# They must be named "name_flux" where name is the name of the flux function
# for example ausm_flux.F, kfvs_flux.F, roe_flux.F, etc.
%_flux_bq.f:  ${SRC}/%_flux.F
		${CDSRC}; ${GCC} -E -C -P $*_flux.F > $*_flux.f
		${CDSRC}; ${TPN} \
		       -head       $*_flux                     \
		       -vars       "qcl qcr qvl qvr resl resr" \
		       -outvars    "qcl qcr qvl qvr resl resr" \
		       -difffuncname "_bq"                     \
		        $*_flux.f
		${CDSRC}; ${RM} -f $*_flux.f

%_flux_bx.f:  ${SRC}/%_flux.F
		${CDSRC}; ${GCC} -E -C -P $*_flux.F > $*_flux.f
		${CDSRC}; ${TPN} \
		       -head       $*_flux                   \
		       -vars       "x1 x2 qvl qvr resl resr" \
		       -outvars    "x1 x2 qvl qvr resl resr" \
		       -difffuncname "_bx"                   \
		        $*_flux.f
		${CDSRC}; ${RM} -f $*_flux.f
#-----------------------------------------------------------------------------
# Flux for solid edges
solid_flux_bq.f:	${SRC}/solid_flux.f
		${CDSRC}; ${TPN} \
		       -head       solid_flux \
		       -vars       "qc res"   \
		       -outvars    "qc res"   \
		       -difffuncname "_bq"    \
		        solid_flux.f;

solid_flux_bx.f:	${SRC}/solid_flux.f
		${CDSRC}; ${TPN} \
		       -head       solid_flux   \
		       -vars       "x1 x2 res"  \
		       -outvars    "x1 x2 res"  \
		       -difffuncname "_bx"      \
		        solid_flux.f;
#-----------------------------------------------------------------------------
# Flux for farfield face
farfield_flux_bq.f: ${SRC}/farfield_flux.F
		${CDSRC}; ${GCC} -E -C -P farfield_flux.F > farfield_flux.f
		${CDSRC}; ${TPN} \
		       -head       farfield_flux \
		       -vars       "qc res"      \
		       -outvars    "qc res"      \
		       -difffuncname "_bq"       \
		        farfield_flux.f;

farfield_flux_bx.f: ${SRC}/farfield_flux.F
		${CDSRC}; ${GCC} -E -C -P farfield_flux.F > farfield_flux.f
		${CDSRC}; ${TPN} \
		       -head       farfield_flux \
		       -vars       "x1 x2 res"   \
		       -outvars    "x1 x2 res"   \
		       -difffuncname "_bx"       \
		        farfield_flux.f;
#-----------------------------------------------------------------------------
# Vertex averaging
vaverage_bq.f: ${SRC}/vaverage.f
		${CDSRC}; ${TPN} \
		       -head       vaverage         \
		       -vars       "qc qv1 qv2 qv3" \
		       -outvars    "qc qv1 qv2 qv3" \
		       -difffuncname "_bq"          \
		        vaverage.f;

vaverage_bx.f: ${SRC}/vaverage.f
		${CDSRC}; ${TPN} \
		       -head       vaverage                           \
		       -vars       "x1 x2 x3 af1 af2 af3 qv1 qv2 qv3" \
		       -outvars    "x1 x2 x3 af1 af2 af3 qv1 qv2 qv3" \
		       -difffuncname "_bx"          \
		        vaverage.f;
#-----------------------------------------------------------------------------
# Differentiation of averaging weights
afact1_bx.f: ${SRC}/avgfact2.f
		${CDSRC}; ${TPN} \
		       -head       afact1           \
		       -vars       "x1 x2 x3 sax1 sax2 sax3 say1 say2 say3  sax21 sax22 sax23 say21 say22 say23 saxy1 saxy2 saxy3" \
		       -outvars    "x1 x2 x3 sax1 sax2 sax3 say1 say2 say3  sax21 sax22 sax23 say21 say22 say23 saxy1 saxy2 saxy3" \
		       -difffuncname "_bx"          \
		        ${SRC}/avgfact2.f;

afact2_bx.f: ${SRC}/avgfact2.f
		${CDSRC}; ${TPN} \
		       -head       afact2           \
		       -vars       "x1 x2 x3 xp nx ny sax say sax2 say2 saxy" \
		       -outvars    "x1 x2 x3 xp nx ny sax say sax2 say2 saxy" \
		       -difffuncname "_bx"          \
		        ${SRC}/avgfact2.f;

afact3_bx.f: ${SRC}/avgfact2.f
		${CDSRC}; ${TPN} \
		       -head       afact3 \
		       -vars       "sax say sax2 say2 saxy afact" \
		       -outvars    "sax say sax2 say2 saxy afact" \
		       -difffuncname "_bx"          \
		        ${SRC}/avgfact2.f;

afact4_bx.f: ${SRC}/avgfact2.f
		${CDSRC}; ${TPN} \
		       -head       afact4 \
		       -vars       "x1 x2 x3 af1 af2 af3" \
		       -outvars    "x1 x2 x3 af1 af2 af3" \
		       -difffuncname "_bx"          \
		        ${SRC}/avgfact2.f;

bd_normal_bx.f: ${SRC}/avgfact2.f
		${CDSRC}; ${TPN} \
		       -head       bd_normal \
		       -vars       "x1 x2 xp nx ny" \
		       -outvars    "x1 x2 xp nx ny" \
		       -difffuncname "_bx"          \
		        ${SRC}/avgfact2.f;
#-----------------------------------------------------------------------------
killnormalvel_bq.f: ${SRC}/killnormalvel.f
		${CDSRC}; ${TPN} \
		       -head       killnormalvel \
		       -vars       "qv"          \
		       -outvars    "qv"          \
		       -difffuncname "_bq"       \
		        killnormalvel.f;

killnormalvel_bx.f: ${SRC}/killnormalvel.f
		${CDSRC}; ${TPN} \
		       -head       killnormalvel \
		       -vars       "x1 x2 x3 qv" \
		       -outvars    "x1 x2 x3 qv" \
		       -difffuncname "_bx"       \
		        killnormalvel.f;
#-----------------------------------------------------------------------------
# Limiter function
limit_bq.f:	${SRC}/limit.f
		${CDSRC}; ${TPN} \
		       -head       limit             \
		       -vars       "dl dr"           \
		       -outvars    "dl dr limit"     \
		       -difffuncname "_bq"           \
		        ${SRC}/limit.f;

limit_bx.f:	${SRC}/limit.f
		${CDSRC}; ${TPN} \
		       -head       limit             \
		       -vars       "dl dr"           \
		       -outvars    "dl dr limit"     \
		       -difffuncname "_bx"           \
		        ${SRC}/limit.f;
#-----------------------------------------------------------------------------
cost_q.o:	cost_q.F Makefile
		${GCC} -E -C -P -I${SRC} -D${COST} cost_q.F > cost_q.f
		${FC} ${FFLAGS} -c cost_q.f
		${RM} -f cost_q.f
#-----------------------------------------------------------------------------
cost_x.o:	cost_x.F Makefile
		${GCC} -E -C -P -I${SRC} -D${COST} cost_x.F > cost_x.f
		${FC} ${FFLAGS} -c cost_x.f
		${RM} -f cost_x.f
#-----------------------------------------------------------------------------
costfunc.o:	costfunc.F
		${GCC} -E -C -P -I${SRC} costfunc.F > costfunc.f
		${FC} ${FFLAGS} -c costfunc.f
		${RM} -f costfunc.f
#-----------------------------------------------------------------------------
costfunc_bq.f:	costfunc.F
		${GCC} -E -C -P -I${SRC} costfunc.F > costfunc.f
		${TPN} \
		       -head       costfunc        \
		       -vars       "qc cost"  \
		       -outvars    "qc cost"  \
		       -difffuncname "_bq"         \
		        costfunc.f;
		${RM} -f costfunc.f

costfunc_bx.f:	costfunc.F
		${GCC} -E -C -P -I${SRC} costfunc.F > costfunc.f
		${TPN} \
		       -head       costfunc              \
		       -vars       "x1 x2 cost"  \
		       -outvars    "x1 x2 cost"  \
		       -difffuncname "_bx"               \
		        costfunc.f;
		${RM} -f costfunc.f
#-----------------------------------------------------------------------------
adStack.o: ${ADSTACK}
		${CC} ${CFLAGS} -c ${ADSTACK}
#-----------------------------------------------------------------------------
adBuffer.o: ${ADBUFFER}
		${FC} ${FFLAGS} -c ${ADBUFFER}
##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o *.msg *~ $(TARGETS)

cleanall:	
	${RM} -f *.o *.msg *~ *_bq.f *_bx.f $(TARGETS)
