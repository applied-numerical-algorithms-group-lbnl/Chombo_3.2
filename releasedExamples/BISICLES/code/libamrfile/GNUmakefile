# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
DIM=2



ebase := testlibamrfile
src_dirs := src
include ../mk/Make.defs


libamrfile.so: $(objects) testlibamrfile
	$(CXX) -o libamrfile.so -shared  $(objects)  $(CHOMBO_HOME)/libamrtools$(config).a $(CHOMBO_HOME)/libboxtools$(config).a $(CHOMBO_HOME)/libbasetools$(config).a $(HDFLIBFLAGS)


