USE_PROTO=TRUE
NAMESPACE=TRUE
PROTO_HOME = $(SCRATCH)/proto

ebase = amr_helmholtz
CHOMBO_HOME = ../../../lib
LibNames = AMRElliptic AMRTools BoxTools
base_dir = .
src_dirs = ../common
INPUT = _inputs/standard.inputs


include $(CHOMBO_HOME)/mk/Make.example

pristine:
	make realclean
	cd $(CHOMBO_HOME); make realclean

main:
	make realclean
	make -j32 all
	mv $(ebase)*.ex main.exe
