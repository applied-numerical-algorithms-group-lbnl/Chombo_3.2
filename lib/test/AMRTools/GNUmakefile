# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_AMRTools

ebase :=  testPWL quadCFInterpTest testPeriodic fluxRegisterTest \
          testRegionGather testCoarseAverage testPeriodicFillPatch \
	testComputeSum  FineInterpEdgeTest refluxEdgeTest testPeriodicFR \
	testFourthOrderFineInterp testFineInterp  fourthOrderCFInterpTest\
	nwoQuadCFInterpTest CoDimCopierTest testFourthOrderTensorFineInterp\


LibNames := AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test

LIBFLAGS += $(LAPACKLIBS)

