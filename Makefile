# This makefile doesnt actually do anything. It just runs the Makefile in the 
# specific directories.
RM=rm -f

#Location of flow solver files
SRCFLO=./src-flo

#Location of adjoint files
SRCADJ=./src-adj

#Mesh deformation and adjoint files
SRCMSH=./src-mesh

#Directory with some tapenade-related stuff
TPDIR=./tapenade

TARGETS =  flo adj mesh

ALL:		$(TARGETS)

flo:
		cd ${SRCFLO}; make

adj:
		cd ${SRCADJ}; make

mesh:
		cd ${SRCMSH}; make

##############################################################################
# clean things up
##############################################################################

clean:	
	cd ${SRCFLO};  ${RM} *.o *.msg *~ flo 
	cd ${SRCADJ};  ${RM} *.o *.msg *~ adj
	cd ${SRCMSH};  ${RM} *.o *.msg *~ deform adjoint
	cd ${TPDIR};   ${RM} *.o
