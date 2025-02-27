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
INCLUDES        = ../../GEFEL/:../../MEFEL/PREP/:../../MEFEL/SRC/:../../GEFEL/PARSER/:../../GEFEL/OPTIM/:../../GEFEL/RANDY/:../../GEFEL/CMLFILE/:./metis-5.0pre2/include/:./metis-5.0pre2/libmetis/:metis-5.0pre2/programs/:metis-5.0pre2/GKlib/trunk:./
#Chaco-2.2/code/main/:PARTY/:scotch_5.1/include/



################################################################### 
#  List of local libraries
#
LOC_LIBS        =  -lmprep -lmef -lgef  -lmetis 



####################################################################
#  List of source files
#
SRCS = meshdecomp.cpp 



####################################################################
#  List of object files
#
OBJS         := $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
MESHDECOMP   := $(OUTPUTPATH)meshdecomp
DEPS         := $(SRCS:%.cpp=$(OUTPUTPATH)%.d)



####################################################################
#  List of lib files
#

# libraries which preprocessors depends on
LIBDEPS  := $(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)libgef.a
LIBDEPS  += $(SIFEL_ROOT)$(BIN_ROOT)MEFEL/SRC/$(OUTPUTDIR)libmef.a
LIBDEPS  += $(SIFEL_ROOT)$(BIN_ROOT)MEFEL/PREP/$(OUTPUTDIR)libmprep.a
LIBDEPS  += ./metis-5.0pre2/build/Linux-x86_64/libmetis.a
# path to directories with needed libraries for the preprocessors target
DIR_LIBS := -L$(SIFEL_ROOT)$(BIN_ROOT)GEFEL/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)MEFEL/SRC/$(OUTPUTDIR)
DIR_LIBS += -L$(SIFEL_ROOT)$(BIN_ROOT)MEFEL/PREP/$(OUTPUTDIR)
DIR_LIBS += -L./metis-5.0pre2/build/Linux-x86_64/



####################################################################
#  Targets
#

all:  mechprep metispart $(MESHDECOMP)

$(MESHDECOMP): $(OBJS) $(LIBDEPS)
	@(echo "##### Creating meshdecomp executable . . .")
	$(CC) -o $@ $(OBJS) $(EXECFLAGS) $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)
	@(echo "##### meshdecomp was successfully created")

$(OBJS):  | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

$(LIBDEPS): mechprep metispart

clean:    
	$(RM) $(MESHDECOMP) $(OBJS) $(DEPS) *~ core
	@(echo " meshdecomp was successfully cleaned")

cleanall: mechprep
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(cd metis-5.0pre2; make realclean)	
	@(echo " meshdecomp was successfully cleaned")

cleandepall: mechprep
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(cd metis-5.0pre2; make realclean)	
	@(echo " Preprocessors were successfully cleaned")



####################################################################
#  List of dependent module targets
#
mechprep:
	@(cd ../../MEFEL/PREP; make $(MAKECMDGOALS))

metispart:
	@(cd metis-5.0pre2/; make all)


####################################################################
#  Dependencies generated by compilers
#
-include $(DEPS)



####################################################################
#  Suffix rules
#
$(OUTPUTPATH)%.o : %.cpp
	@export CPATH=$(INCLUDES); \
	echo $(CC) $(CFLAGS) $(DEFS) $(INCLMETIS) -c $(<D)/$(<F); \
	$(CC) $(CFLAGS) $(COMPATFLAGS) $(DEFS) $(INCLMETIS) -c $(<D)/$(<F) -o $(@D)/$(@F)
