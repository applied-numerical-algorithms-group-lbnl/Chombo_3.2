USE_PROTO=TRUE
NAMESPACE=TRUE
PROTO_HOME= $(HOME)/_proto/proto
ebase = old_viscous_tensor_solve
CHOMBO_HOME = ../../../lib
LibNames = AMRElliptic AMRTools BoxTools
base_dir = .
src_dirs = ../common
INPUT = _inputs/standard.inputs


include $(CHOMBO_HOME)/mk/Make.example
XTRACPPFLAGS += -I$(PROTO_HOME)/include/transition
XTRACXXFLAGS += -I$(PROTO_HOME)/include/transition

pristine:
	make realclean
	cd $(CHOMBO_HOME); make realclean

main:
	make realclean
	make all
	mv $(ebase)*.ex main.exe
