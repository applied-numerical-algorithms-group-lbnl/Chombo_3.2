# -*- Mode: Makefile; -*-

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../../../Chombo/lib


ebase = extract3DncTo2dAMR

MINDIM = 2
MAXDIM = 3

# this is the local GNUmakefile which contains this example's multidim
# build info -- libraries and source directory information
MULTIDIM_MAKEFILE = GNUmakefile.multidim

# this is the default target, which just invokes the all-multidim target
all: all-multidim

noLibs: all-multidim-nolibs


#include ../../../mk/Make.defs
include ./Make.defs



#if NETCDF_INC was not defined, try to work it out from nc-config or NETCDF_HOME
ifeq ($(NETCDF_INC),)

#no NETCDF_HOME? try nc-config
ifeq ($(NETCDF_HOME), )
NETCDF_HOME := $(shell nc-config --prefix)

ifneq ($(NETCDF_HOME), )
NETCDF_INC := -I$(shell nc-config --includedir)

ifeq ($(shell nc-config --has-f90), yes)
NETCDF_LIBS := $(shell nc-config --flibs) #does no harm to link the fortran libs if present
else ifeq ($(shell nc-config --has-fortran), yes)
NETCDF_LIBS := $(shell nc-config --flibs) #does no harm to link the fortran libs if present
else
NETCDF_LIBS := $(shell nc-config --libs)
endif

endif

else
NETCDF_INC :=-I$(NETCDF_HOME)/include
NETCDF_LIBS := -L$(NETCDF_HOME)/lib -lnetcdf
endif
endif


$(info NetCDF: $(NETCDF_INC) $(NETCDF_LIBS) )
ifneq ($(NETCDF_INC), )
CPPFLAGS += -DHAVE_NETCDF $(NETCDF_INC)
#compileOptions="-DHAVE_NETCDF $${NETCDF_INC}  $${compileOptions}"; 
endif

#compileOptions ="-DTEST"

# this file contains the basic rules used to build multidim codes (using the 
# GNUmakefile.multidim in this directory), including the shell script
# which orchestrates the make process 
include $(CHOMBO_HOME)/mk/Make.multidim.basic
