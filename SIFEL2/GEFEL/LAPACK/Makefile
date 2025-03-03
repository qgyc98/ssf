# Generated Automatically From Makefile.In By Configure.
################################################################
#
#  Makefile for mechcode
# 
################################################################

.NOEXPORT:
SHELL = /bin/sh

srcdir = .

prefix = .
exec_prefix = $(prefix)

################################################################
#  If you want to install, rename these variables appropriately
################################################################
incdir = $$HOME/include
libdir = $$HOME/lib

INSTALL         = @INSTALL@
INSTALL_PROGRAM = install.sh
INSTALL_DATA    = install.sh

#CC              = xlc
#CC              = mpiCC
CC              = g++
#CC              = gcc
#CC              = cc

DEFS            = 
CFLAGS          = -g -Wall
#CFLAGS          = -O3
#CFLAGS          = -g  -I/usr/local/qt/include
#CFLAGS          = -g  -I/usr/lib/qt-2.0.1/include
#CFLAGS          = -g  -I/usr/lib/qt-2.2.3/include
#CFLAGS          = -g  -I/usr/lib/qt-2.3.1/include
#CFLAGS          = -O3 -L/usr/lib/qt-2.3.1/lib -lqt
#CFLAGS          = -O3 -I/usr/lib/qt-2.3.1/include
#CFLAGS          = -g
LDFLAGS         = 
LDOPTIONS       = $(CDEBUGFLAGS) $(LDFLAGS)

LIBS            =  -lm
X_LIBS          =

LOC_LIBS        =
SYS_LIBS        = -lm

INCLUDES        = 

AR              = ar
AR_FLAGS        = rc
RANLIB          = ranlib

RM              = rm -f
DEPEND          = makedepend


#
# You probably don't need to change anything below here.
#




####################################################################
# List of object files 

OBJS = spallocate.o spbuild.o spfactor.o spoutput.o spsolve.o sputils.o

####################################################################
# List of source files
# 

SRCS = spallocate.cpp spbuild.cpp spfactor.cpp spoutput.cpp spsolve.cpp sputils.cpp



#########################################################################
#  Targets
#

all: lapack

lapack: $(OBJS) 
	$(RM) $@

 
clean:
#	$(RM) pars core
	$(RM) lapack *.o core

distclean: clean
	$(RM) makefile config.status config.cache config.log

################################################################

makefile: makefile.in configure
	cd $(srcdir); 
	./configure

configure: configure.in
	cd $(srcdir); 
	autoconf

################################################################
# Dependencies

.cpp.o:
	$(CC) -c $(CFLAGS) $(COMMONFLAGS) $(DEFS) $(INCLUDES) $<

depend: $(SRCS)
	$(DEPEND) -DMAKEDEPEND $(DEFS) $(INCLUDES) $(SRCS)
