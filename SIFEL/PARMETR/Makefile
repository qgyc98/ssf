###############################################################################
#
#  Makefile for parallel coupled mechanical and transport parts of code PMETR
#
###############################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../Makefile.in



################################################################### 
#  Process targets given on the command line
#
ALLTGTS := all pmetr noprep
TGT     := $(filter $(ALLTGTS), $(MAKECMDGOALS))
TGTTYPE := $(filter-out $(ALLTGTS), $(MAKECMDGOALS))

ifeq ($(strip $(TGT)),)  
  # if no target has been specified then set default one 'all'
  TGT := all
endif

$(info Making '$(TGT)' target with option '$(TGTTYPE)')



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : pmefel ptrfel pmetr pmetrprep noprep



####################################################################
#  Targets
#
all:            pmetrprep

# cleans PMETR
clean:             $(TGT)

# cleans PMETR and dependent targets such as GEFEL, TRFEL and MEFEL and their dependent targets
cleanall:          $(TGT)

# cleans all debug and optimized versions of PMETR and also dependent targets such as PMEFEL, PTRFEL, PARGEF, MEFEL, MECHPREP, TRFEL, TRANSPREP and GEFEL including module directories in BIN
cleandepall:       $(TGT)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " pmetr BIN was successfully cleaned")



####################################################################
#  List of targets with dependent modules
#
noprep:
	+@(cd ./SRC/;  $(MAKE) $(TGTTYPE))

pmetrprep:
	+@(cd ./PREP/;  $(MAKE) $(TGTTYPE))

pmetr:  pmefel ptrfel
	+@(cd ./SRC/;  $(MAKE) $(TGTTYPE))

pmefel:
	+@(cd ../PARMEF/;   $(MAKE) $(TGTTYPE))

ptrfel: |pmefel;
	+@(cd ../PARTRF/;   $(MAKE) $(TGTTYPE))
