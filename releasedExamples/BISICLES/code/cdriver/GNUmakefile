# -*- Mode: Makefile -*- 

include ../mk/Make.defs

lib:
	cd ../lib; $(MAKE) lib



fcc := $(subst FALSE,$(FC),$(subst TRUE,mpif90,$(MPI)))

fcflag := $(subst FALSE,,$(subst TRUE,-DCH_MPI,$(MPI)))
#not sure if this only works with open mpi : we need the C++ mpi lib as well as the fortran lib
lmpicxx := $(subst FALSE,,$(subst TRUE,-lmpi_cxx,$(MPI)))
ftestwrapper.$(config).ex : lib
	$(fcc) $(fcflag) -I../src/ -g -DCH_LANG_CC -o  ftestwrapper.$(config).ex testwrapper.F90 -L ../lib -lBisicles$(config)  $(LIBFLAGS) -lstdc++ $(lmpicxx)



ftestwrapper:  ftestwrapper.$(config).ex

all: ftestwrapper

.PHONY : ctestwrapper ftestwrapper lib
