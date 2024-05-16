# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := solverBenchmark
all:  clean libbike  solverBenchmark
solverBenchmark: clean libbike

.NOTPARALLEL: all

include ../mk/Make.defs

reallyclean: realclean libclean

libbike:
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)



