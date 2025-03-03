####################################################################
#
#  Makefile for parallel mechanical part of code PMEFEL
#
####################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../Makefile.in



###################################################################
#  Process targets given on the command line 
#
ALLTGTS := all pmefel pmechprep
TGT     := $(filter $(ALLTGTS), $(MAKECMDGOALS))
TGTTYPE := $(filter-out $(ALLTGTS), $(MAKECMDGOALS))

ifeq ($(strip $(TGT)),)  
  # if no target has been specified then set default one 'all'
  TGT := all
endif

DEPTGT  := pmefel
PATTERN := dep
ifeq ($(findstring $(PATTERN),$(TGTTYPE)),dep)
  DEPTGT := 
endif
$(info Making '$(TGT)' target with option '$(TGTTYPE)')



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : pmefel pmechprep

all: pmechprep

# cleans MEFEL and MECHPREP
clean:             $(TGT)

# cleans PMEFEL PMECHPREP and dependent targets such as GEFEL
cleanall:          $(TGT)

# cleans all debug and optimized versions of PMEFEL and PMECHPREP and also dependent targets such as PARGEF, MEFEL, MECHPREP and GEFEL including module directories in BIN
cleandepall:       $(TGT)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " pmefel BIN was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
pmechprep: $(DEPTGT)
	+@(cd ./PREP/;  $(MAKE) $(TGTTYPE))

pmefel:
	+@(cd ./SRC/;   $(MAKE) $(TGTTYPE))
