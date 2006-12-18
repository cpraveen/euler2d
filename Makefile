RM=rm -f

#Location of flow solver files
SRCFLO=./src-flo

#Location of adjoint files
SRCADJ=./src-adj

#Mesh deformation and adjoint files
SRCMESH=./src-mesh

#Directory with some tapenade-related stuff
TPDIR=./tapenade

TARGETS =  flo adj mesh

ALL:		$(TARGETS)

flo:
		cd ${SRCFLO}; make

adj:
		cd ${SRCADJ}; make

mesh:
		cd ${SRCMESH}; make

##############################################################################
# clean things up
##############################################################################

clean:	
	cd ${SRCFLO};  ${RM} *.o *.msg *~ flo 
	cd ${SRCADJ};  ${RM} *.o *.msg *~ adj
	cd ${SRCMESH}; ${RM} *.o *.msg *~ deform adjoint
	cd ${TPDIR};   ${RM} *.o
