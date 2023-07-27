MPI = TRUE
HDF5 = TRUE
USE_PETSC=FALSE

ebase = comm_split_test

CHOMBO_HOME = ../../../lib

LibNames = AMRElliptic AMRTools BoxTools BaseTools

base_dir = .
src_dirs = 

INPUT = split.inputs

include $(CHOMBO_HOME)/mk/Make.example

pristine:
	\rm -rf d f o p *.ex *.hdf5 pout* time.* ../../../lib/*.a ../../../lib/src/*/d ../../../lib/src/*/f ../../../lib/src/*/o ../../../lib/src/*/p
