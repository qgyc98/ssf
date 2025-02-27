################################################################
#
#  Makefile for SPARSE direct solver
# 
################################################################



################################################################### 
#  Include default part of Makefile - 
#  compiler definitions, flags, libraries ...
#
include ../../Makefile.in



################################################################### 
#  List of include directories
#
INCLUDES        = ../GEFEL/:../PARSER:.



################################################################### 
#  List of libraries
#

# list of local libraries
LOC_LIBS	= 

# list of paths to dependnet libraries
DIR_LIBS 	= 

# list of all libraries - SYS_LIBS are defined in the default Makefile.in
LIBS 		= $(DIR_LIBS) $(LOC_LIBS) $(SYS_LIBS)



####################################################################
#  List of source files and directories
#
SRCS = BigMatrix.cpp BiSection.cpp DSSolver.cpp  DenseMatrixArithmeticsNN.cpp IntArrayList.cpp 
SRCS += MathTracer.cpp SparseConectivityMtx.cpp SparseGridMtx.cpp SparseGridMtxLDL.cpp SparseGridMtxLL.cpp 
SRCS += SparseGridMtxLU.cpp SparseGridMtxPD.cpp SparseMatrixF.cpp



####################################################################
#  List of object files
#
OBJS = $(SRCS:%.cpp=$(OUTPUTPATH)%.o)
LIBSPARSE = $(OUTPUTPATH)libsparse.a
DEPS      = $(SRCS:%.cpp=$(OUTPUTPATH)%.d)



####################################################################
#  Targets
#
all:  $(LIBSPARSE)

$(LIBSPARSE): $(OBJS) 
	$(RM) $@
	@(echo "##### Creating SPARSE library . . .")
	@(echo $(AR) $(ARFLAGS) $@ $(SRCS:%.cpp=$(OUTPUTPATH)%.o))
	@($(AR) $(ARFLAGS) $@ $(OBJS))
	@(echo "##### SPARSE library was successfully created")

$(OBJS): | $(OUTPUTPATH)

$(OUTPUTPATH):
	@(mkdir -p $(OUTPUTPATH))

clean: 
	@($(RM) $(LIBSPARSE) $(OBJS) $(DEPS) *~ core)
	@(echo " sparse was successfully cleaned")

cleanall:       
	@($(RM) -r -f $(addprefix $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH),$(REMOVEDIR)))
	@(echo " sparse was successfully cleaned")

cleandepall:
	@($(RM) -r -f $(SIFEL_ROOT)$(BIN_ROOT)$(MODULE_PATH))
	@(echo " sparse BIN was successfully cleaned")



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
