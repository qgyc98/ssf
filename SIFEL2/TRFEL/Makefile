####################################################################
#
#  Makefile for transport part of code TRFEL
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
ALLTGTS := all trfel transprep
TGT     := $(filter $(ALLTGTS), $(MAKECMDGOALS))
TGTTYPE := $(filter-out $(ALLTGTS), $(MAKECMDGOALS))

ifeq ($(strip $(TGT)),)  
  # if no target has been specified then set default one 'all'
  TGT := all
endif

DEPTGT  := trfel
PATTERN := dep
ifeq ($(findstring $(PATTERN),$(TGTTYPE)),dep)
  DEPTGT := 
endif
$(info Making '$(TGT)' target with option '$(TGTTYPE)')



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : trfel transprep



####################################################################
#  Targets
#
all:               transprep

# cleans TRFEL and TRANSPREP
clean:             $(TGT)

# cleans TRFEL, TRANSPREP and dependent targets such as GEFEL
cleanall:          $(TGT)

# cleans all debug and optimized versions of TRFEL and TRANSPREP and also dependent targets such as GEFEL including modules directories in BIN
cleandepall:       $(TGT)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " trfel BIN was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
transprep: $(DEPTGT)
	+@(cd ./PREP/;  $(MAKE) $(TGTTYPE))

trfel:
	+@(cd ./SRC/;   $(MAKE) $(TGTTYPE))
