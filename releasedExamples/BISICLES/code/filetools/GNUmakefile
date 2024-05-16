
# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase :=  nctoamr amrtotxt amrtoplot flatten extract merge addbox amrtocf stats glfaces glfacesNew faces rescale sum pythonf diagnostics subregion gia-stats ppismip6
src_dirs = src
all:  clean libbike 

.NOTPARALLEL: all

include ../mk/Make.defs
 

#if NETCDF_INC was not defined, try to work it out from nc-config or NETCDF_HOME
ifeq ($(NETCDF_INC),)

#no NETCDF_HOME? try nc-config
ifeq ($(NETCDF_HOME), )
NETCDF_HOME := $(shell nc-config --prefix)

ifneq ($(NETCDF_HOME), )
NETCDF_INC := -I$(shell nc-config --includedir)

ifeq ($(shell nc-config --has-f90), yes)
NETCDF_LIBS := $(shell nc-config --flibs) $(shell nc-config --libs) #does no harm to link the fortran libs if present
else ifeq ($(shell nc-config --has-fortran), yes)
NETCDF_LIBS := $(shell nc-config --flibs) $(shell nc-config --libs) #does no harm to link the fortran libs if present
else
NETCDF_LIBS := $(shell nc-config --libs)
endif

endif

else
NETCDF_INC := -I$(NETCDF_HOME)/include
NETCDF_LIBS := -L$(NETCDF_HOME)/lib -lnetcdf
endif
endif


$(info NetCDF: $(NETCDF_INC) $(NETCDF_LIBS) )
ifneq ($(NETCDF_INC), )
CPPFLAGS += -DHAVE_NETCDF $(NETCDF_INC)
endif


##if GDAL_INC was not defined, try to work it out from gdal-config or GDAL_HOME
#ifeq ($(GDAL_INC),)

##no GDAL_HOME? try gdal-config
#ifeq ($(GDAL_HOME), )
#GDAL_HOME := $(shell gdal-config --prefix)
#ifneq ($(GDAL_HOME), )
#GDAL_INC := $(shell gdal-config --cflags)
#GDAL_LIBS := $(shell gdal-config --libs)
#endif

##GDAL_HOME was define, so assume the usual layout
#else
#GDAL_INC :=-I$(GDAL_HOME)/include
#GDAL_LIBS := -L$(GDAL_HOME)/lib -lgdal
#endif
#endif


$(info gdal: $(GDAL_INC) $(GDAL_LIBS) )
ifneq ($(GDAL_INC),)
CPPFLAGS+=$(GDAL_INC) -DHAVE_GDAL
endif


libbike:
	cd ../lib; $(MAKE) lib


LIBFLAGS := -L../lib/ -lBisicles$(config) $(NETCDF_LIBS) $(GDAL_LIBS)  $(LIBFLAGS)
