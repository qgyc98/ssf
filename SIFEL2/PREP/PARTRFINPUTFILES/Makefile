# Generated Automatically From Makefile.In By Configure.
####################################################################
#
#  Makefile for mechcode
#
####################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../../Makefile.in



################################################################### 
#  List of include directories
#
INCLUDES        =../../GEFEL/:../../GEFEL/PARSER/:../../GEFEL/OPTIM:../../GEFEL/RANDY:../../TRFEL/PREP/:../../TRFEL/SRC

LOC_LIBS        = -ltprep -ltrf -lgef 
SYS_LIBS        = -lm -rdynamic -ldl -lbfd



####################################################################
#  List of source files
#
SRCS  = partrfinputfiles.cpp



####################################################################
#  List of object files
#
OBJS             := $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
PARTRFINPUTFILES := $(OUTPUTPATH)partrfinputfiles



####################################################################
#  List of lib files
#
# libraries which preprocessors depends on
LIBDEPS  := $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a
LIBDEPS  += $(SIFEL_ROOT)$(BIN_ROOT)TRFEL/SRC/$(OUTPUTDIR)libtrf.a
LIBDEPS  += $(SIFEL_ROOT)$(BIN_ROOT)TRFEL/PREP/$(OUTPUTDIR)libtprep.a

# path to directories with needed libraries for the preprocessors target
DIR_LIBS := -L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)TRFEL/SRC/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)TRFEL/PREP/$(OUTPUTDIR)



####################################################################
#  Targets
#
.PHONY : transprep

all:  transprep $(PARTRFINPUTFILES)

$(PARTRFINPUTFILES): $(OBJS) $(LIBDEPS) | $(OUTPUTPATH) 
	@(echo "##### Creating partrfinputfiles executable . . .")
	$(CC) -o $@ $(OBJS) $(EXECFLAGS) $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)
	@(echo "##### partrfinputfiles was successfully created")

$(OBJS):  | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS): transprep;

clean:
	@($(RM) $(PARTRFINPUTFILES) $(OBJS) *~ core)
	@(echo " partrfinputfiles was successfully cleaned")

cleanall: $(LIBDEPS)
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " partrfinputfiles was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " partrfinputfiles was successfully cleaned")



####################################################################
#  Module targets
#

transprep:
	@(cd ../../TRFEL/PREP; make $(MAKECMDGOALS))



####################################################################
#  Suffix rules
#
$(OUTPUTPATH)%.o : %.cpp
	@export CPATH=$(INCLUDES); \
	echo $(CC) $(CFLAGS) $(DEFS) -c $(<D)/$(<F); \
	$(CC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) -c $(<D)/$(<F) -o $(@D)/$(@F)



####################################################################
#  Dependencies
#

depend:
	@(cd ../../GEFEL; make depend)
	@(cd ../../MEFEL/SRC; make depend)
	@(cd ../../TRFEL/SRC; make depend)
	$(DEPEND)  $(DEFS) $(INCLUDES) $(SRCS)


dependclean:
	(t="/tmp/t_$$" ; sed 1,/"DO\ NOT\ DELETE"/p Makefile -n > t ; mv t Makefile)





####################################################################
# DO NOT DELETE
