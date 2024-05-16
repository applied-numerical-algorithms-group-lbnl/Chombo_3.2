FC = gfortran

include ../../code/mk/Make.defs

FFLAGS = -g -I. -I$(NETCDF_HOME)/include $(HDFINCFLAGS)   
LIBFLAGS = -L$(NETCDF_HOME)/lib -L$(NETCDF_HOME)/lib $(NETCDF_LIBS) -lhdf5_hl $(HDFLIBFLAGS)
LDFLAGS = -g

OBJS = mgrelax.o friction.o

all: $(OBJS)
	$(FC) $(FFLAGS) -o friction mgrelax.o friction.o $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FFLAGS) -c  $< -o $@

clean: 
	$(RM) $(OBJS) a.out

realclean:
	$(RM) $(OBJS) a.out *.mod
