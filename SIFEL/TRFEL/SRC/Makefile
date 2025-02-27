####################################################################
#
#  Makefile for transport problems TRFEL
#
####################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../../Makefile.in



################################################################### 
#  List of include directories
#
INCLUDES        = ../../GEFEL/:../../GEFEL/PARSER/:../../GEFEL/OPTIM:../../GEFEL/RANDY:../../GEFEL/CMLFILE



################################################################### 
#  List of libraries
#
LOC_LIBS        = -lgef



####################################################################
#  List of source files
#
SRCS  = aepointst.cpp backupsolt.cpp baroghelB.cpp baroghel_reten.cpp bazant_reten.cpp bazped.cpp bnodvalt.cpp boundfluxes.cpp carbmat1.cpp cemhydmat.cpp 
SRCS += cerny_concrete.cpp climatcond.cpp climatcond2.cpp concreteB.cpp consol_awf1.cpp consol_wf1.cpp consol_awf2.cpp consol_hawf3.cpp consol_hwf2.cpp consol_wf2.cpp constrel.cpp cpnnpsolvert.cpp 
SRCS += cpnpsolvert.cpp crsection1d.cpp crsection2d.cpp crsection3d.cpp C30baroghel.cpp 
SRCS += C60baroghel.cpp C30bazant.cpp C60bazant.cpp damisotrmat.cpp dampermeability.cpp devriesmat.cpp discisotrmat.cpp discmat.cpp 
SRCS += dnnpsolvert.cpp dnpsolvert.cpp edget.cpp elementt.cpp elemswitcht.cpp endnodet.cpp febex_granit_reten.cpp fourmedia.cpp gardner.cpp gen2delem.cpp
SRCS += glasgowmat.cpp globalt.cpp globmatt.cpp gmultiphase.cpp grunewaldmat.cpp homogt.cpp homogmat.cpp homogtrans.cpp hydrationheat.cpp
SRCS += interfacematt.cpp interfacequadrilat.cpp intpointst.cpp isotherm.cpp isotrmat.cpp 
SRCS += kunzel.cpp kunzel2.cpp lhsrhst.cpp lewis_schrefler.cpp linbart.cpp linbart3d.cpp linbartax.cpp lincoupmat.cpp linhext.cpp lintett.cpp linwedget.cpp 
SRCS += loadcaset.cpp loadelt.cpp masin_reten.cpp millymat.cpp moistheat.cpp multiphase.cpp nlisotrmat.cpp npglvec.cpp nnpsolvert.cpp nodet.cpp nonlinmant.cpp 
SRCS += npsolvert.cpp nspsolvert.cpp onemedium.cpp outdiagt.cpp outdrivert.cpp o30bazant.cpp pedersen.cpp 
SRCS += probdesct.cpp pvalt.cpp potts.cpp quadbart.cpp quadbartax.cpp quadhext.cpp quadlinaxisym.cpp quadlineart.cpp 
SRCS += quadquadrilatt.cpp quadquadrilattax.cpp radiationmat.cpp richards.cpp richardscontam.cpp saltmat1.cpp saltmat2.cpp saltmat3.cpp saltmat4.cpp seebeckheat.cpp
SRCS += sejtkrmat.cpp soil1mat.cpp solvert.cpp sorpisohansen.cpp sorpisorootgrun.cpp sourcet.cpp spsolvert.cpp stochdrivert.cpp tdisotrmat.cpp
SRCS += threemedia.cpp transbclc.cpp transcrsec.cpp transmat.cpp transprint.cpp transtop.cpp trfelinit.cpp 
SRCS += trlinaxisym.cpp trlineart.cpp twomedia.cpp van_genuchten.cpp adaptivityt.cpp

HERMTRFSRCS = hermtrf.cpp



####################################################################
#  List of object files
#
OBJS         := $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
HERMTRFOBJS  := $(HERMTRFSRCS:%.cpp=$(OUTPUTPATH)%.o)
MOBJS        := $(OUTPUTPATH)trfel.o
TRFEL        := $(OUTPUTPATH)trfel
LIBTRFEL     := $(OUTPUTPATH)libtrf.a
LIBHERMTRF   := $(OUTPUTPATH)libhermtrf.a
DEPS         := $(SRCS:%.cpp=$(OUTPUTPATH)%.d) $(HERMTRFSRCS:%.cpp=$(OUTPUTPATH)%.d) $(MOBJS:%.o=%.d)



####################################################################
#  List of lib files which TRFEL depends on
#

# list of file names of dependent libraries
LIBDEPS  = $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a

# list of paths to dependent libraries
DIR_LIBS = -L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : gefel



####################################################################
#  Targets
#
all:  $(TRFEL)

$(TRFEL): $(LIBDEPS) $(MOBJS) $(OBJS) $(HERMTRFOBJS)
	@($(RM) $@)
	@($(RM) ../$@)
	@(echo "##### Creating TRFEL library . . .")
	@(echo $(AR) $(ARFLAGS) $(LIBTRFEL) $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $(LIBTRFEL) $(OBJS))
	@(echo "##### TRFEL library was successfully created")
	@(echo "##### Creating HERMTRF library . . .")
	@(echo $(AR) $(ARFLAGS) $(LIBHERMTRF) $(HERMTRFSRCS:%.cpp=$(OUTPUTDIR)%.o) $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $(LIBHERMTRF) $(HERMTRFOBJS) $(OBJS))
	@(echo "##### HERMTRF library was successfully created")
	@(echo "##### Creating trfel executable . . .")
	$(CC) -o $@ $(MOBJS) $(EXECFLAGS) -L$(OUTPUTPATH) -ltrf $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)
	@(echo "##### trfel was successfully created")

$(MOBJS):       | $(OUTPUTPATH)
$(OBJS):        | $(OUTPUTPATH)
$(HERMTRFOBJS): | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS): gefel;

clean:
	@($(RM) $(TRFEL) $(LIBTRFEL) $(OBJS) $(MOBJS) $(HERMTRFOBJS) $(LIBHERMTRF) $(DEPS) *~ core)
	@(echo " trfel was successfully cleaned")

cleanall: $(LIBDEPS)
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " trfel was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " trfel was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
gefel:
	+@(cd ../../GEFEL; $(MAKE) $(MAKECMDGOALS))



####################################################################
#  Dependencies generated by compilers
#
-include $(DEPS)



####################################################################
#  Suffix rules
#
$(OUTPUTPATH)%.o : %.cpp
	@export CPATH=$(INCLUDES); \
	echo $(CC) $(CFLAGS) -c $(<D)/$(<F); \
	$(CC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) -c $(<D)/$(<F) -o $(@D)/$(@F)
