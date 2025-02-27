#
# Generated Automatically From Makefile.In By Configure.
####################################################################
#
#  Makefile for mechanical part of code MEFEL
#
####################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../Makefile.in



ALLTGTS := all mefel mechprep
TGT     := $(filter $(ALLTGTS), $(MAKECMDGOALS))
TGTTYPE := $(filter-out $(ALLTGTS), $(MAKECMDGOALS))

ifeq ($(strip $(TGT)),)  
  # if no target has been specified then set default one 'all'
  TGT := all
endif

DEPTGT  := mefel
PATTERN := dep
ifeq ($(findstring $(PATTERN),$(TGTTYPE)),dep)
  DEPTGT := 
endif
$(info Making '$(TGT)' target with option '$(TGTTYPE)')


###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : mefel mechprep

all: mechprep

# cleans MEFEL and MECHPREP
clean:             $(TGT)

# cleans MEFEL MECHPREP and dependent targets such as GEFEL
cleanall:          $(TGT)

# cleans all debug and optimized versions of MEFEL and MECHPREP and also dependent targets such as GEFEL including module directories in BIN
cleandepall:       $(TGT)
	$(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH)
	@(echo " mefel BIN was successfully cleaned")




####################################################################
#  List of targets with dependent modules
#
mechprep: $(DEPTGT)
	+@(cd ./PREP/;  $(MAKE) $(TGTTYPE))

mefel:
	+@(cd ./SRC/;   $(MAKE) $(TGTTYPE))
