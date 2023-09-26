# the base name(s) of the application(s) in this directory
USE_PROTO=TRUE
NAMESPACE=TRUE
PROTO_HOME= $(HOME)/_proto/proto

ebase = data_transfer

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../common

# input file for 'run' target
INPUT = regression.inputs
# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example
XTRACPPFLAGS += -I$(PROTO_HOME)/include/transition
XTRACXXFLAGS += -I$(PROTO_HOME)/include/transition

# application-specific variables

# application-specific targets
