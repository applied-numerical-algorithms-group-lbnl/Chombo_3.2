FC = gfortran

include ../../code/mk/Make.defs

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
NETCDF_LIBS := $(shell nc-config --flibs) $(shell nc-config --libs)  #does no harm to link the fortran libs if present
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

#ifneq ($(NETCDF_INC), )
#CPPFLAGS += -DHAVE_NETCDF $(NETCDF_INC)
#endif


FFLAGS = -fPIC -g -I. $(NETCDF_INC) $(HDFINCFLAGS)   
LIBFLAGS = $(NETCDF_LIBS) -lhdf5_hl $(HDFLIBFLAGS)
LDFLAGS = -g

OBJS = mgrelax.o friction.o

friction: $(OBJS)
	$(FC) $(FFLAGS) -o friction mgrelax.o friction.o $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FFLAGS) -c  $< -o $@

#clean: 
#	$(RM) $(OBJS) a.out

#realclean:
#	$(RM) $(OBJS) a.out *.mod
