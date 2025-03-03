####################################################################
#
#  Makefile for parallel transport part of code PTRFEL
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
ALLTGTS := all ptrfel ptransprep
TGT     := $(filter $(ALLTGTS), $(MAKECMDGOALS))
TGTTYPE := $(filter-out $(ALLTGTS), $(MAKECMDGOALS))

ifeq ($(strip $(TGT)),)  
  # if no target has been specified then set default one 'all'
  TGT := all
endif

DEPTGT  := ptrfel
PATTERN := dep
ifeq ($(findstring $(PATTERN),$(TGTTYPE)),dep)
  DEPTGT := 
endif
$(info Making '$(TGT)' target with option '$(TGTTYPE)')



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : ptrfel ptransprep



####################################################################
#  Targets
#
all:            ptransprep

# cleans TRFEL and TRANSPREP
clean:             $(TGT)

# cleans PTRFEL, PTRANSPREP and dependent targets such as GEFEL
cleanall:          $(TGT)

# cleans all debug and optimized versions of TRFEL and TRANSPREP and also dependent targets such as GEFEL including modules directories in BIN
cleandepall:       $(TGT)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " ptrfel BIN was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
ptransprep: $(DEPTGT)
	+@(cd ./PREP/;  $(MAKE) $(TGTTYPE))

ptrfel:
	+@(cd ./SRC/;   $(MAKE) $(TGTTYPE))
