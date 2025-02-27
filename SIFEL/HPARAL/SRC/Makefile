################################################################################
#
#  Makefile for hybrid parallel mechanical, transport, and coupled  code HPARAL
#
################################################################################


# Default values for make command line targets
PARAM    = $(MAKECMDGOALS)
TGTHPMET = -f Makefilec
TGTHPMEF = -f Makefilem
TGTHPTRF = -f Makefilet
PARDEP = all

ifneq ($(and $(findstring hpmetr,$(MAKECMDGOALS)),  $(findstring hpmefel,$(MAKECMDGOALS)), $(findstring hptrfel,$(MAKECMDGOALS))),)
#
# set parameters in the case that hpmetr, hpmefel, and hptrfel were specified on the command line
#

PARAM = $(filter-out hpmefel,$(MAKECMDGOALS)) # remove 'hpmefel' from command line string
AUX = $(filter-out hptrfel,$(PARAM))          # remove 'hptrfel' from command line string
PARAM := $(AUX)                               # store remaining target options
PARDEP = ; @(echo)                            # avoiding make message 'Nothing to be done for...'
TRFDEP = |hpmefel;                            # if hpmefel and hptrfel are specified, hptrfel waits on finish of hpmefel and hpmetr waits on finish of hptrfel
METDEP = |hptrfel;                            # if hpmefel and hptrfel are specified, hptrfel waits on finish of hpmefel and hpmetr waits on finish of hptrfel

else
#
# othwerwise set parameters for individual targets
#

# hpmefel was specified on the command line
ifeq ($(findstring hpmefel,$(MAKECMDGOALS)), hpmefel)
TGTHPTRF = 
PARAM = $(filter-out hpmefel,$(MAKECMDGOALS)) # remove hptrfel from make parameters and left the remaining target options
ifneq ($(PARAM),)
PARDEP = ; @(echo)                            # avoiding the make message 'Nothing to be done for...'
endif
endif

# hptrfel was specified on the command line
ifeq ($(findstring hptrfel,$(MAKECMDGOALS)), hptrfel)
TGTHPMEF = 
PARAM = $(filter-out hptrfel,$(MAKECMDGOALS)) # remove hptrfel from make parameters and left the rest target options
ifneq ($(PARAM),)
PARDEP = ; @(echo)                            # avoiding the make message 'Nothing to be done for...'
endif
endif

# hpmetr was specified on the command line
ifeq ($(findstring hpmetr,$(MAKECMDGOALS)), hpmetr)
TGTHPMET = 
PARAM = $(filter-out hpmetr,$(MAKECMDGOALS)) # remove hpmetr from make parameters and left the remaining target options
ifneq ($(PARAM),)
PARDEP = ; @(echo)                            # avoiding the make message 'Nothing to be done for...'
endif
endif

endif


#########################################################################
#  Targets
#

.PHONY : all opt deb optdeb depend
.PHONY : omp ompopt ompdeb ompoptdeb 
.PHONY : clean cleandeb cleanopt cleanoptdeb cleandeb dependclean
.PHONY : cleanomp cleanompdeb cleanompopt cleanompoptdeb 
.PHONY : cleandep cleandepdeb clenadepopt cleandepoptdeb 
.PHONY : cleandepompdeb cleandepompopt cleandepompoptdeb 
.PHONY : hpmefel hptrfel hpmetr


all: hpmefel hptrfel hpmetr


hpmefel:
ifneq ($(TGTHPMEF),) 
	+@($(MAKE) $(TGTHPMEF) $(PARAM))
endif



hptrfel: $(TRFDEP)
ifneq ($(TGTHPTRF),) 
	+@($(MAKE) $(TGTHPTRF) $(PARAM))
endif


hpmetr: $(METDEP)
ifneq ($(TGTHPMET),) 
	+@($(MAKE) $(TGTHPMET) $(PARAM))
endif



opt:    $(PARDEP)
deb:   	$(PARDEP)
optdeb: $(PARDEP)

omp:       $(PARDEP)
ompdeb:    $(PARDEP)
ompopt:    $(PARDEP)
ompoptdeb: $(PARDEP)

clean:       $(PARDEP)
cleandeb:    $(PARDEP)
cleanopt:    $(PARDEP)
cleanoptdeb: $(PARDEP)

cleanomp:       $(PARDEP)
cleanompdeb:    $(PARDEP)
cleanompopt:    $(PARDEP)
cleanompoptdeb: $(PARDEP)

cleandep:       $(PARDEP)
cleandepdeb:    $(PARDEP)
cleandepopt:    $(PARDEP)
cleandepoptdeb: $(PARDEP)

cleandepompdeb:    $(PARDEP)
cleandepompopt:    $(PARDEP)
cleandepompoptdeb: $(PARDEP)
