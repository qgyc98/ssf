####################################################################
#
#  Makefile for complete package SIFEL
#
####################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ./Makefile.in



################################################################### 
#  Process targets given on the command line
#
ALLTGTS := all metr mefel trfel pmefel ptrfel pmetr
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
.PHONY : all mefel trfel metr pmefel ptrfel pmetr



####################################################################
#  Targets
#
all:            pmetr

# makes optimized version of code without debugging symbols of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR
opt:            $(TGT)

# makes version of code with debug symbols and strict code checking of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR
deb:            $(TGT)

# makes optimized version of code with debugging symbols of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR
opt:            $(TGT)

# makes OpenMP version of code with debugging symbols of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR and their parallel versions
omp:               $(TGT)

# makes optimized OpenMP version of code without debugging symbols of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR and their parallel versions
ompopt:            $(TGT)

# makes OpenMP version of code with debug symbols and strict code checking of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR and their parallel versions
ompdeb:            $(TGT)

# makes OpenMP optimized version of code with debugging symbols of MEFEL, MECHPREP, TRFEL, TRANSPREP and METR and their parallel versions
ompoptdeb:         $(TGT)

# cleans debug versions of PMETR
clean:          $(TGT)

# cleans OpenMP debug versions of PMETR
cleanomp:          $(TGT)

# cleans debug versions of PMETR
cleandeb:       $(TGT)

# cleans optimized versions of PMETR
cleanopt:       $(TGT)

# cleans optimized versions with debugging symbols of PMETR
cleanoptdeb:    $(TGT)

# cleans OpenMP debug versions of PMETR
cleanompdeb:       $(TGT)

# cleans OpenMP optimized versions of PMETR
cleanompopt:       $(TGT)

# cleans OpenMP optimized versions with debugging symbols of PMETR
cleanompoptdeb:    $(TGT)

# cleans all debug and optimized versions of PMETR and also dependent targets such as MEFEL, MECHPREP, TRFEL, TRANSPREP, GEFEL and their parallel versions
cleandep:       $(TGT)

# cleans debug versions of PMETR and also dependent targets such as MEFEL, MECHPREP, TRFEL, TRANSPREP, GEFEL and their parallel versions
cleandepdeb:    $(TGT)

# cleans all debug and optimized versions of PMETR and also dependent targets such as MEFEL, MECHPREP, TRFEL, TRANSPREP, GEFEL and their parallel versions including module directories in BIN
cleandepall:       $(TGT)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT))
	@(echo " SIFEL BIN was successfully cleaned")

# cleans optimized versions of PMETR and also dependent targets such as MEFEL, MECHPREP, TRFEL, TRANSPREP, GEFEL and their parallel versions
cleandepopt:    $(TGT)

# cleans optimized versions with debugging symbols of PMETR and also dependent targets such as MEFEL, MECHPREP, TRFEL, TRANSPREP and GEFEL and their parallel versions
cleandepoptdeb: $(TGT)

# cleans OpenMP debug versions of PMETR and also dependent targets such as PMEFEL, PTRFEL, PARGEF, MEFEL, MECHPREP, TRFEL, TRANSPREP and GEFEL
cleandepompdeb:    $(TGT)

# cleans OpenMP optimized versions of PMETR and also dependent targets such as PMEFEL, PTRFEL, PARGEF, MEFEL, MECHPREP, TRFEL, TRANSPREP and GEFEL
cleandepompopt:    $(TGT)

# cleans OpenMP optimized versions with debugging symbols of PMETR and also dependent targets such as PMEFEL, PTRFEL, PARGEF, MEFEL, MECHPREP, TRFEL, TRANSPREP and GEFEL
cleandepompoptdeb: $(TGT)



####################################################################
#  List of targets with dependent modules
#
pmetr:
	+@(cd ./PARMETR/; $(MAKE) $(TGTTYPE))

pmefel:
	+@(cd ./PARMEF/;  $(MAKE) $(TGTTYPE))

ptrfel:
	+@(cd ./PARTRF/;  $(MAKE) $(TGTTYPE))

metr:
	+@(cd ./METR/;    $(MAKE) $(TGTTYPE))

mefel:
	+@(cd ./MEFEL/;   $(MAKE) $(TGTTYPE))

trfel:
	+@(cd ./TRFEL/;   $(MAKE) $(TGTTYPE))
