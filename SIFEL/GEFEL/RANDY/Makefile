########################
#                      #
#  Makefile for RANDY  #
#                      #
########################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../../Makefile.in



### DEFs
#DEFS = -D__make_fast__



################################################################### 
#  List of include directories
#
INCLUDES = 



################################################################### 
#  List of libraries
#

# list of local libraries
LOC_LIBS	=

# list of paths to dependent libraries
DIR_LIBS 	=

# list of all libraries - SYS_LIBS are defined in the default Makefile.in
LIBS = $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)



####################################################################
#  List of files
#
SRCS =	randy.cpp



####################################################################
#  List of object files
#
OBJS      = $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
LIBRANDY  = $(OUTPUTPATH)librandy.a
DEPS      = $(SRCS:%.cpp=$(OUTPUTPATH)%.d)



####################################################################
#  Targets
#
all: $(LIBRANDY)

$(LIBRANDY) : $(OBJS)
	$(RM) $@
	@(echo "##### Creating RANDY library . . .")
	@(echo $(AR) $(ARFLAGS) $@ $(SRCS:%.cpp=$(OUTPUTDIR)%.o))
	@($(AR) $(ARFLAGS) $@ $(OBJS))
	@(echo "##### RANDY library was successfully created")

$(OBJS) : | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

clean:
	@($(RM) $(OBJS) $(DEPS) *~ core)
	@($(RM) $(LIBRANDY))
	@(echo " randy was successfully cleaned")

cleanall:       
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " randy was successfully cleaned")

cleandepall:       
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " randy BIN was successfully cleaned")



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
