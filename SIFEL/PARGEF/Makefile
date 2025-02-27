###################################################################
#
#  Makefile for PARGEF - general parallel code library
# 
###################################################################



# the following variable indicates the makefile as the parallel part of code
PARSIFEL = yes



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../Makefile.in



################################################################### 
#  List of include directories
#
INCLUDES        = ../GEFEL/:../GEFEL/PARSER/:../GEFEL/OPTIM:../GEFEL/RANDY:../GEFEL/CMLFILE
#INCLUDES        = ../GEFEL/:../GEFEL/PARSER/:../GEFEL/OPTIM:../GEFEL/RANDY:../GEFEL/CMLFILE:$(MPI_HOME)/include:$(PETSC_DIR)/include:$(PETSC_DIR)/bmake/linux-gnu-cxx-debug



################################################################### 
#  List of libraries
#
LOC_LIBS	= -lgef

# list of all libraries - SYS_LIBS are defined in the default Makefile.in
LIBS 		= $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)



####################################################################
# List of source files
#
SRCS = fixnodesel.cpp parcongrad.cpp partop.cpp partopjk.cpp dpfeti.cpp feti1.cpp lplate.cpp psolver.cpp schurcompl.cpp selnodes.cpp

GEFEL_DIR=../GEFEL/




####################################################################
# List of object files
#
OBJS := $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
LIBPARGEF := $(OUTPUTPATH)libpargef.a	
DEPS      := $(SRCS:%.cpp=$(OUTPUTPATH)%.d)



####################################################################
#  List of libraries needed by PARGEF
#

# file names of dependent libraries
LIBDEPS = $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a

# path to dependent libraries
DIR_LIBS=-L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : gefel



#########################################################################
#  Targets
#
all: $(LIBPARGEF)

$(LIBPARGEF): $(LIBDEPS) $(OBJS)
	@($(RM) $@)
	@(echo "##### Creating PARGEF library . . .")
	@(echo $(AR) $(ARFLAGS) $@  $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $@ $(OBJS))
	@(echo "##### PARGEF library was successfully created")

$(OBJS):  | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS) : gefel;

clean:
	$(RM) $(LIBPARGEF) $(OBJS) $(DEPS) *~ core
	@(echo " pargef was successfully cleaned")

cleanall: $(LIBDEPS)
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " pargef was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " pargef BIN was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
gefel:
	+@(cd $(GEFEL_DIR); $(MAKE) $(MAKECMDGOALS))



####################################################################
# Dependencies generated by compilers
#
-include $(DEPS)



####################################################################
#  Suffix rules
#
$(OUTPUTPATH)%.o : %.cpp
	@export CPATH=$(INCLUDES); \
	echo $(PCC) $(CFLAGS) -c $(<D)/$(<F); \
	$(PCC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) -c $(<D)/$(<F) -o $(@D)/$(@F)
