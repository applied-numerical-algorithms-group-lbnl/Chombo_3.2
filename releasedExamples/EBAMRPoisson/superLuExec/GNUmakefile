# -*- Mode: Makefile -*-

# This is not an embedded boundary example, just a simple test
USE_EB=FALSE

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = testCooToCsc

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTools \
		BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src ../superLuSrc

# input file for 'run' target
INPUT = 

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets (not needed)
XTRACXXFLAGS = -DCH_USE_SUPERLU

# XTRALIBFLAGS is not needed because I am using the /usr/lib installation 
# from ubuntu
# XTRALIBFLAGS = /home/mando/Development/superlu/lib/libsuperlu_5.1.a
