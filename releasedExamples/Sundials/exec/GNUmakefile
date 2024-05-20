# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_Sundials_exec

## Define the variables needed by Make.example
USE_SUNDIALS=TRUE
MPI=TRUE

# Point to extra Sundials libraries needed for this example
XTRALIBFLAGS += -lsundials_arkode -lsundials_sunnonlinsolnewton -lsundials_ida -lsundials_sunlinsolpcg

# the base name(s) of the application(s) in this directory
ebase = testSundialsPDE

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = Sundials AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src

# input file for 'run' target
# INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
