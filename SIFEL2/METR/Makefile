####################################################################
#
#  Makefile for coupled mechanical and transport part of code METR
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
ALLTGTS := all metr noprep
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
.PHONY : mefel trfel metr noprep



####################################################################
#  Targets
#
all:            metr

# cleans METR
clean:             $(TGT)

# cleans METR and dependent targets such as GEFEL, TRFEL and MEFEL and their dependent targets
cleanall:          $(TGT)

# cleans all debug and optimized versions of METR and also dependent targets such as MEFEL, MECHPREP, TRFEL, TRANSPREP and GEFEL including module directories in BIN
cleandepall:       $(TGT)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " metr BIN was successfully cleaned")



####################################################################
#  Module targets
#

# no preprocessor compilation in MEFEL and TRFEL 
noprep:
	+@(cd ./SRC/;  $(MAKE) $(TGTTYPE))

# full compilation including MEFEL and TRFEL preprocessors
metr:  mefel trfel
	+@(cd ./SRC/;  $(MAKE) $(TGTTYPE))

# full compilation of MEFEL including preprocessor
mefel:
	+@(cd ../MEFEL/;   $(MAKE) $(TGTTYPE))

# full compilation of TRFEL including preprocessor
trfel: |mefel;
	+@(cd ../TRFEL/;   $(MAKE) $(TGTTYPE))
