# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_BoxTools

ebase = transformTest ldIVSFABCopyTest ldIVSFABCopyTestInt interiorExchangeTest copy2Test \
  broadcastTest copyTest   threadTest domainSplitTest fabTest gatherTest         \
  HDF5attributes  HDF5boxIO HDF5data newIVSTest testBox               \
  testIntVectSet testBaseFabMacros testLoadBalance testMeshRefine     \
  testPeriodic ivsfabTest testRealVect codimensionBoundaryTest        \
  testTreeIntVectSet scopingTest reductionTest testRealTensor         \
  testCHArray mortonTest testIndicesTransformation matrixTest stdIVSTest \
  boxCountThreadTest edgeAndCellTest FaceSumOpTest testMDArrayMacros \
  testWriteCenterings


LibNames = BoxTools

include $(CHOMBO_HOME)/mk/Make.test

LIBFLAGS += $(LAPACKLIBS)
