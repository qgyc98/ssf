################################################################
#
#  Makefile for GEFEL - general code library
# 
################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../Makefile.in



################################################################### 
#  List of include directories
#
INCLUDES  := ./OPTIM/:./RANDY/:./CMLFILE/:./PARSER/:./SPARSE/:
ifeq ($(INC_PERMON), yes)
 INCLUDES := $(INCLUDES)$(PERMON_DIR)/include/:$(PETSC_DIR)/include/:$(PETSC_DIR)/$(PETSC_ARCH)/include/:
endif



################################################################### 
#  List of libraries
#

# list of local libraries
LOC_LIBS	= -lgef

# path to directories with needed libraries for the GEFEL target
DIR_LIBS=-L$(OUTPUTPATH)

# list of all libraries - SYS_LIBS are defined in the default Makefile.in
LIBS 		= $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)



####################################################################
#  List of source files and directories
#
SRCS  = aggregator.cpp basefun.cpp bmatrix.cpp boundbox.cpp cc.cpp cr.cpp csv.cpp densemat.cpp diagmat.cpp difcalc.cpp 
SRCS += dskyline.cpp eigsol.cpp elemmat.cpp elemnode.cpp elemtools.cpp endnode.cpp
SRCS += floatsub.cpp fllopsif.cpp fuzzygen.cpp fuzzynum.cpp gadaptivity.cpp gedge.cpp gelement.cpp gfunct.cpp gmatrix.cpp gfmatrix.cpp
SRCS += gnode.cpp gphase.cpp gsurf.cpp gtopology.cpp hdbcontr.cpp import_T3d.cpp intools.cpp intp.cpp 
SRCS += iotools.cpp ipmap.cpp itablefunct.cpp kwdset.cpp lcoordsys.cpp lgnode.cpp libtrace.cpp list.cpp mathem.cpp matrix.cpp nodmap.cpp obj_funct.cpp 
SRCS += obj_funct_gefel.cpp ordering.cpp precond.cpp saddpoint.cpp scc.cpp scr.cpp selection.cpp seqfeti.cpp seqschur.cpp
SRCS += seqselnodes.cpp seqtop.cpp siftop.cpp skyline.cpp slesolv.cpp fllopsif.cpp stacktrace.cpp tablefunct.cpp 
SRCS += tensor.cpp timecontr.cpp vector.cpp xfile.cpp

CMLFILE_DIR  = CMLFILE/
RANDY_DIR    = RANDY/
OPTIM_DIR    = OPTIM/
PARSER_DIR   = PARSER/
SPARSE_DIR   = SPARSE/



####################################################################
#  List of created object files
#
OBJS         = $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
GOBJS        = $(OUTPUTPATH)gefel.o
GEFEL        = $(OUTPUTPATH)gefel
LIBGEFEL     = $(OUTPUTPATH)libgef.a
DEPS         = $(SRCS:%.cpp=$(OUTPUTPATH)%.d) $(GOBJS:%.o=%.d)



####################################################################
#  List of dependent libraries with directories
#
LIBDEPS  = $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(CMLFILE_DIR)$(OUTPUTDIR)libcmlfile.a
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(RANDY_DIR)$(OUTPUTDIR)librandy.a 
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OPTIM_DIR)$(OUTPUTDIR)liboptim.a 
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(PARSER_DIR)$(OUTPUTDIR)libparser.a
LIBDEPS += $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(SPARSE_DIR)$(OUTPUTDIR)libsparse.a
ifeq ($(INC_PERMON), yes) # PERMON library support should be included
 LIBDEPS += $(PERMON_DIR)/$(PETSC_ARCH)/lib/libfllop.a
endif



###########################################################################################
#  Following targets will be made regardless of whether there are files with the same names
#
.PHONY : randy optim cmlfile parser sparse



####################################################################
#  Targets
#
all:  $(GEFEL)

$(GEFEL): $(LIBGEFEL) $(GOBJS)
	@($(RM) $@)
	@(echo "##### Creating gefel executable . . .")
	@(echo "Runtime support for tracing of source files was $(TRACE_SRC_F)")
	$(CC) -o $@ $(GOBJS) $(EXECFLAGS) $(LIBS)
	@(echo "##### gefel executable was successfully created")

$(LIBGEFEL): $(LIBDEPS) $(OBJS)
	@(echo "##### Creating GEFEL base library . . .")
	@(echo $(AR) $(ARFLAGS) $@ $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $@ $(OBJS))
	@(echo "##### GEFEL base library was successfully created")
	@(echo "##### Merging libraries . . .")
	$(MERGELIB) $@ $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(CMLFILE_DIR)$(OUTPUTDIR)libcmlfile.a   ||  $(RM) $@
	@(test -e $@)
	$(MERGELIB) $@ $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(RANDY_DIR)$(OUTPUTDIR)librandy.a       ||  $(RM) $@
	@(test -e $@)
	$(MERGELIB) $@ $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OPTIM_DIR)$(OUTPUTDIR)liboptim.a       ||  $(RM) $@
	@(test -e $@)
	$(MERGELIB) $@ $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(PARSER_DIR)$(OUTPUTDIR)libparser.a     ||  $(RM) $@
	@(test -e $@)
	$(MERGELIB) $@ $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(SPARSE_DIR)$(OUTPUTDIR)libsparse.a     ||  $(RM) $@
	@(test -e $@)
	@(echo "##### GEFEL library was successfully created")

$(OBJS):  | $(OUTPUTPATH)

$(GOBJS): | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS): cmlfile randy optim parser sparse;

clean:  $(LIBDEPS)
	@($(RM) $(LIBGEFEL) $(GEFEL) $(GOBJS) $(OBJS) $(DEPS) *~ core)
	@(echo " gefel was successfully cleaned")

cleanall:  $(LIBDEPS)      
	@($(RM) -r $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " gefel was successfully cleaned")

cleandepall: $(LIBDEPS)
	@($(RM) -r $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " gefel BIN was successfully cleaned")


####################################################################
#  List of targets with dependent modules
#
cmlfile:
	+@(cd $(CMLFILE_DIR);  $(MAKE) $(MAKECMDGOALS))

randy:  
	+@(cd $(RANDY_DIR);    $(MAKE) $(MAKECMDGOALS))

optim:  
	+@(cd $(OPTIM_DIR);    $(MAKE) $(MAKECMDGOALS))

parser:
	+@(cd $(PARSER_DIR);   $(MAKE) $(MAKECMDGOALS))

sparse:  
	+@(cd $(SPARSE_DIR);   $(MAKE) $(MAKECMDGOALS))




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
	$(CC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) -I/usr/include/suitesparse/ -c $(<D)/$(<F) -o $(@D)/$(@F)
