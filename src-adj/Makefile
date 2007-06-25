#Include common variables
include ../Makefile.in

#Which cost function to use
#COST1 depends only on cell-center values
#COST2 depends only on vertex values
#COST3 for pressure matching, also set costtype=1 in adj.inp
COST=COST1

#Tapenade command
TPN     = ${TAPENADE_HOME}/bin/tapenade \
        -ext ${GenLib} -extAD ${ADLib} -b -O ${SRCADJ}

#
# fortran flags
#
FFLAGS = ${CFLAGS} -I${SRCFLO}

TARGETS =  adj

#
# program dependencies
#

adj     = main.o \
			 ${SRCFLO}/avgfact.o \
			 ${SRCFLO}/avgfact2.o \
			 ${SRCFLO}/clcd.o \
			 ${SRCFLO}/common.o \
          ${SRCFLO}/geometric.o \
			 ${SRCFLO}/isocont.o \
			 ${SRCFLO}/killnormalvel.o \
			 ${SRCFLO}/limit.o \
          ${SRCFLO}/math.o \
			 ${SRCFLO}/vaverage.o \
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
			 costfunc_pre_bq.o \
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
			 ${TPDIR}/adStack.o \
			 ${TPDIR}/adBuffer.o

#
# programs
#

ALL:		$(TARGETS)

adj:	${adj}
		${FC} -o adj ${adj}
		${RM} -f *.msg *.msg~

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
%_flux_bq.f:  ${SRCFLO}/%_flux.F
		${CDFLO}; ${GCC} -E -C -P $*_flux.F > $*_flux.f
		${CDFLO}; ${TPN} \
		       -head       $*_flux                     \
		       -vars       "qcl qcr qvl qvr resl resr" \
		       -outvars    "qcl qcr qvl qvr resl resr" \
		       -difffuncname "_bq"                     \
		        $*_flux.f
		${CDFLO}; ${RM} -f $*_flux.f

%_flux_bx.f:  ${SRCFLO}/%_flux.F
		${CDFLO}; ${GCC} -E -C -P $*_flux.F > $*_flux.f
		${CDFLO}; ${TPN} \
		       -head       $*_flux                   \
		       -vars       "x1 x2 qvl qvr resl resr" \
		       -outvars    "x1 x2 qvl qvr resl resr" \
		       -difffuncname "_bx"                   \
		        $*_flux.f
		${CDFLO}; ${RM} -f $*_flux.f
#-----------------------------------------------------------------------------
# Flux for solid edges
solid_flux_bq.f:	${SRCFLO}/solid_flux.f
		${CDFLO}; ${TPN} \
		       -head       solid_flux \
		       -vars       "qc res"   \
		       -outvars    "qc res"   \
		       -difffuncname "_bq"    \
		        solid_flux.f;

solid_flux_bx.f:	${SRCFLO}/solid_flux.f
		${CDFLO}; ${TPN} \
		       -head       solid_flux   \
		       -vars       "x1 x2 res"  \
		       -outvars    "x1 x2 res"  \
		       -difffuncname "_bx"      \
		        solid_flux.f;
#-----------------------------------------------------------------------------
# Flux for farfield face
farfield_flux_bq.f: ${SRCFLO}/farfield_flux.F
		${CDFLO}; ${GCC} -E -C -P farfield_flux.F > farfield_flux.f
		${CDFLO}; ${TPN} \
		       -head       farfield_flux \
		       -vars       "qc res"      \
		       -outvars    "qc res"      \
		       -difffuncname "_bq"       \
		        farfield_flux.f;
		${CDFLO}; ${RM} -f farfield_flux.f

farfield_flux_bx.f: ${SRCFLO}/farfield_flux.F
		${CDFLO}; ${GCC} -E -C -P farfield_flux.F > farfield_flux.f
		${CDFLO}; ${TPN} \
		       -head       farfield_flux \
		       -vars       "x1 x2 res"   \
		       -outvars    "x1 x2 res"   \
		       -difffuncname "_bx"       \
		        farfield_flux.f;
		${CDFLO}; ${RM} -f farfield_flux.f
#-----------------------------------------------------------------------------
# Vertex averaging
vaverage_bq.f: ${SRCFLO}/vaverage.f
		${CDFLO}; ${TPN} \
		       -head       vaverage         \
		       -vars       "qc qv1 qv2 qv3" \
		       -outvars    "qc qv1 qv2 qv3" \
		       -difffuncname "_bq"          \
		        vaverage.f;

vaverage_bx.f: ${SRCFLO}/vaverage.f
		${CDFLO}; ${TPN} \
		       -head       vaverage                           \
		       -vars       "x1 x2 x3 af1 af2 af3 qv1 qv2 qv3" \
		       -outvars    "x1 x2 x3 af1 af2 af3 qv1 qv2 qv3" \
		       -difffuncname "_bx"          \
		        vaverage.f;
#-----------------------------------------------------------------------------
# Differentiation of averaging weights
afact1_bx.f: ${SRCFLO}/avgfact2.f
		${CDFLO}; ${TPN} \
		       -head       afact1           \
		       -vars       "x1 x2 x3 sax1 sax2 sax3 say1 say2 say3  sax21 sax22 sax23 say21 say22 say23 saxy1 saxy2 saxy3" \
		       -outvars    "x1 x2 x3 sax1 sax2 sax3 say1 say2 say3  sax21 sax22 sax23 say21 say22 say23 saxy1 saxy2 saxy3" \
		       -difffuncname "_bx"          \
		        ${SRCFLO}/avgfact2.f;

afact2_bx.f: ${SRCFLO}/avgfact2.f
		${CDFLO}; ${TPN} \
		       -head       afact2           \
		       -vars       "x1 x2 x3 xp nx ny sax say sax2 say2 saxy" \
		       -outvars    "x1 x2 x3 xp nx ny sax say sax2 say2 saxy" \
		       -difffuncname "_bx"          \
		        ${SRCFLO}/avgfact2.f;

afact3_bx.f: ${SRCFLO}/avgfact2.f
		${CDFLO}; ${TPN} \
		       -head       afact3 \
		       -vars       "sax say sax2 say2 saxy afact" \
		       -outvars    "sax say sax2 say2 saxy afact" \
		       -difffuncname "_bx"          \
		        ${SRCFLO}/avgfact2.f;

afact4_bx.f: ${SRCFLO}/avgfact2.f
		${CDFLO}; ${TPN} \
		       -head       afact4 \
		       -vars       "x1 x2 x3 af1 af2 af3" \
		       -outvars    "x1 x2 x3 af1 af2 af3" \
		       -difffuncname "_bx"          \
		        ${SRCFLO}/avgfact2.f;

bd_normal_bx.f: ${SRCFLO}/avgfact2.f
		${CDFLO}; ${TPN} \
		       -head       bd_normal \
		       -vars       "x1 x2 xp nx ny" \
		       -outvars    "x1 x2 xp nx ny" \
		       -difffuncname "_bx"          \
		        ${SRCFLO}/avgfact2.f;
#-----------------------------------------------------------------------------
killnormalvel_bq.f: ${SRCFLO}/killnormalvel.f
		${CDFLO}; ${TPN} \
		       -head       killnormalvel \
		       -vars       "qv"          \
		       -outvars    "qv"          \
		       -difffuncname "_bq"       \
		        killnormalvel.f;

killnormalvel_bx.f: ${SRCFLO}/killnormalvel.f
		${CDFLO}; ${TPN} \
		       -head       killnormalvel \
		       -vars       "x1 x2 x3 qv" \
		       -outvars    "x1 x2 x3 qv" \
		       -difffuncname "_bx"       \
		        killnormalvel.f;
#-----------------------------------------------------------------------------
# Limiter function
limit_bq.f:	${SRCFLO}/limit.f
		${CDFLO}; ${TPN} \
		       -head       limit             \
		       -vars       "dl dr"           \
		       -outvars    "dl dr limit"     \
		       -difffuncname "_bq"           \
		        ${SRCFLO}/limit.f;

limit_bx.f:	${SRCFLO}/limit.f
		${CDFLO}; ${TPN} \
		       -head       limit             \
		       -vars       "dl dr"           \
		       -outvars    "dl dr limit"     \
		       -difffuncname "_bx"           \
		        ${SRCFLO}/limit.f;
#-----------------------------------------------------------------------------
cost_q.o:	cost_q.F Makefile
		${GCC} -E -C -P -I${SRCFLO} -D${COST} cost_q.F > cost_q.f
		${FC} ${FFLAGS} -c cost_q.f
		${RM} -f cost_q.f
#-----------------------------------------------------------------------------
cost_x.o:	cost_x.F Makefile
		${GCC} -E -C -P -I${SRCFLO} -D${COST} cost_x.F > cost_x.f
		${FC} ${FFLAGS} -c cost_x.f
		${RM} -f cost_x.f
#-----------------------------------------------------------------------------
costfunc.o:	costfunc.F
		${GCC} -E -C -P -I${SRCFLO} costfunc.F > costfunc.f
		${FC} ${FFLAGS} -c costfunc.f
		${RM} -f costfunc.f
#-----------------------------------------------------------------------------
costfunc_bq.f:	costfunc.F
		${GCC} -E -C -P -I${SRCFLO} costfunc.F > costfunc.f
		${TPN} \
		       -head       costfunc        \
		       -vars       "qc cost"  \
		       -outvars    "qc cost"  \
		       -difffuncname "_bq"         \
		        costfunc.f;
		${RM} -f costfunc.f

costfunc_bx.f:	costfunc.F
		${GCC} -E -C -P -I${SRCFLO} costfunc.F > costfunc.f
		${TPN} \
		       -head       costfunc              \
		       -vars       "x1 x2 cost"  \
		       -outvars    "x1 x2 cost"  \
		       -difffuncname "_bx"               \
		        costfunc.f;
		${RM} -f costfunc.f
#-----------------------------------------------------------------------------
costfunc_pre_bq.f:	costfunc_pre.F
		${GCC} -E -C -P -I${SRCFLO} costfunc_pre.F > costfunc_pre.f
		${TPN} \
		       -head       costfunc_pre    \
		       -vars       "qv cost"  \
		       -outvars    "qv cost"  \
		       -difffuncname "_bq"         \
		        costfunc_pre.f;
		${RM} -f costfunc_pre.f

#-----------------------------------------------------------------------------
${TPDIR}/adStack.o: ${ADSTACK}
		cd ${TPDIR}; ${CC} ${CFLAGS} -c ${ADSTACK}
#-----------------------------------------------------------------------------
${TPDIR}/adBuffer.o: ${ADBUFFER}
		cd ${TPDIR}; ${FC} ${FFLAGS} -c ${ADBUFFER}
##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o *.msg *~ $(TARGETS)

cleanall:	
	${RM} -f *.o *.msg *~ *_bq.f *_bx.f $(TARGETS)
