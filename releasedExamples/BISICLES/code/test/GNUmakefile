# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += BISICLES_test

# the base name(s) of the application(s) in this directory
ebase = testPetsc testLevelSigmaCS testConstitutive testL1L2 testFortranInterfaceIBC testAVCPoisson testCGOptimize

# general stuff....
include ../mk/Make.defs

lib:
	cd ../lib; $(MAKE) lib

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)

