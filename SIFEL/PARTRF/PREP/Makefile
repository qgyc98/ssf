################################################################
#
#  Makefile for preprocessor of parallel mechcode
#  
################################################################



# the following variable indicates the makefile as the parallel part of code
PARSIFEL = yes



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../../Makefile.in



################################################################### 
#  List of include directories
#
INCLUDES         = ../../PARTRF/SRC:../../TRFEL/PREP/:../../TRFEL/SRC/:../../PARGEF/:../../GEFEL/:../../GEFEL/PARSER/:../../GEFEL/OPTIM:../../GEFEL/RANDY:../../GEFEL/CMLFILE:



################################################################### 
#  List of libraries
#

# list of local libraries
LOC_LIBS         = -ltrf -ltprep -lptrf -lpargef -lgef



####################################################################
# List of source files
#
SRCS = ptransprep.cpp



####################################################################
# List of object files
OBJS          := $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
MOBJS         := $(OUTPUTPATH)ptransprep.o
PTRANSPREP    := $(OUTPUTPATH)ptransprep
LIBPTRANSTPEP := $(OUTPUTPATH)libptransprep.a
DEPS          := $(SRCS:%.cpp=$(OUTPUTPATH)%.d) $(MOBJS:%.o=%.d)



####################################################################
#  List of lib files
#

# list of dependent libraries
LIBDEPS  = $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)TRFEL/SRC/$(OUTPUTDIR)libtrf.a
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)TRFEL/PREP/$(OUTPUTDIR)libtprep.a
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)PARGEF/$(OUTPUTDIR)libpargef.a
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)PARTRF/SRC/$(OUTPUTDIR)libptrf.a

# list of paths for dependent libraries
DIR_LIBS  = -L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)TRFEL/SRC/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)TRFEL/PREP/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)PARGEF/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)PARTRF/SRC/$(OUTPUTDIR)

#list of all necessary libraries
LIBS =  $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)
#${PETSC_C_SH_LIB_PATH} ${PETSC_LIB}



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : transprep partrf 



#########################################################################
#  Targets
#
all:  $(PTRANSPREP)

$(PTRANSPREP): $(LIBDEPS) $(MOBJS)
	@($(RM) $@)
	@($(RM) ../$@)
	@(echo "##### Creating ptransprep executable . . .")
	$(PCC) -o $@ $(MOBJS) $(EXECFLAGS) $(LIBS)
	@(echo "##### ptransprep was successfully created")

$(MOBJS): | $(OUTPUTPATH)
$(OBJS):  | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS) : transprep partrf;

clean:    
	@($(RM)  $(PTRANSPREP) $(OBJS) $(DEPS) *~ core)
	@(echo " ptransprep was successfully cleaned")

cleanall: $(LIBDEPS)
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " ptransprep was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " ptransprep was successfully cleaned")



################################################################
#  List of targets with dependent modules
#
transprep:
	+@(cd ../../TRFEL/PREP; $(MAKE) $(MAKECMDGOALS))

partrf: transprep
	+@(cd ../SRC; $(MAKE) $(MAKECMDGOALS))



####################################################################
#  Dependencies generated by compilers
#
-include $(DEPS)



####################################################################
#  Suffix rules
#
$(OUTPUTPATH)%.o : %.cpp
	@export CPATH=$(INCLUDES); \
	echo $(PCC) $(CFLAGS) -c $(<D)/$(<F); \
	$(PCC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) -c $(<D)/$(<F) -o $(@D)/$(@F)
