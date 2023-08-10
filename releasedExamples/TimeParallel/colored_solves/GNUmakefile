MPI = TRUE
HDF5 = TRUE
USE_PETSC=FALSE

ebase = colored_solves

CHOMBO_HOME = ../../../lib

LibNames = AMRElliptic AMRTools BoxTools BaseTools

base_dir = .
src_dirs = 

INPUT = example.inputs

include $(CHOMBO_HOME)/mk/Make.example

pristine:
	\rm -rf LB.txt d f o p *.ex *.hdf5 pout* time.* ../../../lib/*.a ../../../lib/src/*/d ../../../lib/src/*/f ../../../lib/src/*/o ../../../lib/src/*/p
