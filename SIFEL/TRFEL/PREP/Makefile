####################################################################
#
#  Makefile for TRFEL preprocessor
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
INCLUDES        = ../SRC/:../../GEFEL/:../../GEFEL/PARSER/:../../GEFEL/OPTIM:../../GEFEL/RANDY:../../GEFEL/CMLFILE



################################################################### 
#  List of libraries
#

# list of local libraries
LOC_LIBS        = -ltrf -lgef



####################################################################
#  List of source files
#
SRCS  = advectvel.cpp bocont.cpp dbcrst.cpp dbmatt.cpp descript.cpp entitybocon.cpp hangnodet.cpp inputt.cpp intdomt.cpp
SRCS += outputt.cpp



####################################################################
#  List of object files
#
OBJS         := $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
MOBJS        := $(OUTPUTPATH)transprep.o
TRANSPREP    := $(OUTPUTPATH)transprep
LIBTPREP     := $(OUTPUTPATH)libtprep.a
DEPS         := $(SRCS:%.cpp=$(OUTPUTPATH)%.d) $(MOBJS:%.o=%.d)



####################################################################
#  List of lib files which transprep depends on
#

# list of file names to dependent libraries
LIBDEPS  = $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a $(SIFEL_ROOT)$(BIN_ROOT)TRFEL/SRC/$(OUTPUTDIR)libtrf.a

# list of paths to dependent libraries
DIR_LIBS = -L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR) -L$(SIFEL_ROOT)$(BIN_ROOT)TRFEL/SRC/$(OUTPUTDIR)



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : trfel



####################################################################
#  Targets
#
all:  $(TRANSPREP)

$(TRANSPREP): $(LIBDEPS) $(MOBJS) $(OBJS)
	@($(RM) $@)
	@($(RM) ../$@)
	@(echo "\n##### Creating TRANSPREP library . . .")
	@(echo $(AR) $(ARFLAGS) $(OUTPUTPATH)libtprep.a $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $(OUTPUTPATH)libtprep.a $(OBJS))
	@(echo "##### TRANSPREP library was successfully created")
	@(echo "##### Creating transprep executable . . .")
	$(CC) -o $@ $(MOBJS) $(EXECFLAGS) -L$(OUTPUTPATH) -ltprep $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)
	@(echo "##### transprep was successfully created")

$(MOBJS):    | $(OUTPUTPATH)
$(OBJS):     | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS): trfel;

clean:    
	$(RM) $(TRANSPREP) $(LIBTPREP) $(OBJS) $(MOBJS) $(DEPS) *~ core
	@(echo " transprep was successfully cleaned")

cleanall: $(LIBDEPS)
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " transprep was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " transprep was successfully cleaned")


####################################################################
#  List of targets with dependent modules
#
trfel:
	+@(cd ../SRC; $(MAKE) $(MAKECMDGOALS))



####################################################################
#  Dependencies generated by compilers
#
-include $(DEPS)



####################################################################
#  Suffix rules
#
$(OUTPUTPATH)%.o : %.cpp
	@export CPATH=$(INCLUDES); \
	echo $(CC) $(CFLAGS) -c $(<D)/$(<F); \
	$(CC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) -c $(<D)/$(<F) -o $(@D)/$(@F)
