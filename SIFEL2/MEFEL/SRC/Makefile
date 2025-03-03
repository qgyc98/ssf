####################################################################
#
#  Makefile for mechanical part of code MEFEL
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
#  List of local libraries
#
LOC_LIBS  = -lgef



####################################################################
#  List of source files
#
SRCS  =  aci78.cpp adaptivity.cpp aepoints.cpp anisodam.cpp anisodamrot.cpp arclength.cpp argyrisplate.cpp argyristriangle.cpp
SRCS +=  axisrotrec.cpp axisymlq.cpp axisymlqiface.cpp axisymlt.cpp axisymqq.cpp axisymqt.cpp axisymcq.cpp backupsol.cpp barel2d.cpp barel3d.cpp barelq2d.cpp barelq3d.cpp 
SRCS +=  beamel2d.cpp beamel3d.cpp beamgen3d.cpp bbm.cpp boermat.cpp camclay.cpp camclaycoup.cpp cct.cpp 
SRCS +=  cebfip78.cpp cebfipcontactmat.cpp chen.cpp consol.cpp contactmat.cpp cpsolver.cpp creep.cpp creep_b3.cpp creep_dpl.cpp 
SRCS +=  creepb.cpp creepdam.cpp creep_effym.cpp creep_rspec.cpp crsec2dbar.cpp crsec2dbeam.cpp 
SRCS +=  crsec3d.cpp crsec3dbeam.cpp crseclayer.cpp crsecnod.cpp crsecplstr.cpp cusatismaterial.cpp damplast.cpp damplastifacemat.cpp dkq.cpp dkt.cpp dloadcase.cpp 
SRCS +=  dloadn.cpp dloadel.cpp dloadpd.cpp drprag.cpp doubdp.cpp drprag2.cpp dsm.cpp dst.cpp edgem.cpp edsolver.cpp effstress.cpp eigvalsol.cpp elastgmat2d.cpp elastgmat3d.cpp 
SRCS +=  elastisomat.cpp elastisopdmat.cpp elastortomat.cpp elasttime.cpp element.cpp elemparticle.cpp elemswitch.cpp endnodem.cpp epsolver.cpp fdsolver.cpp flsubdom.cpp fixortodam.cpp 
SRCS +=  generalmod.cpp
SRCS +=  geoelast.cpp glasgmech.cpp glasgowdam.cpp global.cpp globmat.cpp gnodvalvm.cpp
SRCS +=  graphmat.cpp hardsoft.cpp hissplas.cpp hexiface.cpp homog.cpp homogmatm.cpp homogmech.cpp hvisolver.cpp inicd.cpp intpoints.cpp isoviscous.cpp hypoplunsatexptherm2.cpp

SRCS +=  j2flow.cpp layplate.cpp lbsolver.cpp lemaitre.cpp lenjonesmat.cpp lfssolver.cpp lhsrhs.cpp 
SRCS +=  linhex.cpp linhexrot.cpp lintet.cpp lintetrot.cpp linwedge.cpp llssolver.cpp loadcase.cpp 
SRCS +=  loadel.cpp loadn.cpp lssolver.cpp mechbclc.cpp mechcrsec.cpp mechmat.cpp mechprint.cpp 
SRCS +=  mechtop.cpp mefelinit.cpp meshtransfer.cpp microfiber.cpp microM4.cpp microSIM.cpp 
SRCS +=  mohrc.cpp mohrcparab.cpp mtglvec.cpp mtsolver.cpp ndsolver.cpp newtonraph.cpp nfssolver.cpp node.cpp 
SRCS +=  nonlocdamg2.cpp nonlinman.cpp nonlocmicroM4.cpp nonlocplast.cpp nssolver.cpp ortodam.cpp ortodam2.cpp ortodamrot.cpp 
SRCS +=  outdiagm.cpp outdriverm.cpp outquantm.cpp outresfilem.cpp plastifacemat.cpp plelemlq.cpp plelemlt.cpp plelemqq.cpp plelemqt.cpp plelemrotlq.cpp
SRCS +=  plelemrotlt.cpp plelemsubqt.cpp plquadiface.cpp probdesc.cpp problem.cpp q4plate.cpp quadhex.cpp quadrilatkirch.cpp
SRCS +=  quadtet.cpp quadwedge.cpp relaxeuroc.cpp scaldam.cpp scaldamcc.cpp seismtool.cpp sequent.cpp 
SRCS +=  shefplast.cpp shellq.cpp shelltr.cpp shelltrm.cpp shrinkmat.cpp simviscous.cpp slipsurf.cpp slopesol.cpp slsolver.cpp 
SRCS +=  soilbeam.cpp soilplateq.cpp soilplatetr.cpp solverm.cpp splas1d.cpp springel.cpp stochdriver.cpp strretalg.cpp tensorcname.cpp
SRCS +=  tetralattice.cpp therisomat.cpp thervolisomat.cpp timeswmat.cpp varelastisomat.cpp vectorcnamem.cpp vecttens.cpp viselast.cpp visplast.cpp
SRCS +=  winpast.cpp z2_smoothing.cpp umatunsat2.cpp



####################################################################
#  List of object files
#
OBJS         = $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
ifeq ($(FC_DEF), enabled)
#  OBJS         += $(OBJS) $(OUTPUTPATH)umatunsat.o
endif
MOBJS        = $(OUTPUTPATH)mefel.o
MEFEL        = $(OUTPUTPATH)mefel

#GENMODTEST   = $(OUTPUTPATH)genmod-test
#GENMODOBJS   = $(OUTPUTPATH)genmod-test.o

LIBMEFEL     = $(OUTPUTPATH)libmef.a
DEPS         = $(SRCS:%.cpp=$(OUTPUTPATH)%.d) $(MOBJS:%.o=%.d)



####################################################################
#  List of lib files
#

# libraries which MEFEL depends on
LIBDEPS  = $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a

# path to directories with needed libraries for the MEFEL target
DIR_LIBS = -L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : gefel



####################################################################
#  Targets
#
all:  $(MEFEL)

#all:  $(MEFEL) $(GENMODTEST)

$(MEFEL): $(LIBDEPS) $(OBJS) $(MOBJS)
	@($(RM) $@)
	@($(RM) ../$@)
	@(echo "##### Creating MEFEL library . . .")
	@(echo $(AR) $(ARFLAGS) $(OUTPUTPATH)libmef.a $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $(OUTPUTPATH)libmef.a $(OBJS))
	@(echo "##### MEFEL library was successfully created")
	@(echo "##### Creating mefel executable . . .")
	@(echo "Runtime support for tracing of source files was $(TRACE_SRC_F)")
	$(CC) -o $@ $(MOBJS) $(EXECFLAGS) -L$(OUTPUTPATH) -lmef $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)
	@(echo "##### mefel was successfully created")

$(MOBJS): | $(OUTPUTPATH)

$(OBJS):  | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS): gefel;

clean:    
	@($(RM)  $(MEFEL) $(LIBMEFEL) $(MOBJS) $(OBJS) $(DEPS) *~ core)
	@(echo " mefel was successfully cleaned")

cleanall: $(LIBDEPS)
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " mefel was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " mefel was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
gefel:
	+@(cd ../../GEFEL;  $(MAKE) $(MAKECMDGOALS))



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

$(OUTPUTPATH)%.o : %.for
	@echo $(FC) -c $(CFLAGS) $(FCFLAGS) $(<D)/$(<F); \
	$(FC) -c $(CFLAGS) $(FCFLAGS) $(<D)/$(<F) -o $(@D)/$(@F)
