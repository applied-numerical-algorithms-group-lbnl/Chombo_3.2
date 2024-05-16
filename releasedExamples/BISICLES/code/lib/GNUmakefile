# -*- Mode: Makefile -*- 
src_dirs =../src ../util
include ../mk/Make.defs

libChomboLibs$(config).a: $(CHOMBO_HOME)/libamrelliptic$(config).a $(CHOMBO_HOME)/libamrtimedependent$(config).a $(CHOMBO_HOME)/libamrtools$(config).a $(CHOMBO_HOME)/libboxtools$(config).a $(CHOMBO_HOME)/libbasetools$(config).a
	mkdir tmp_lib_objects; \
	cd tmp_lib_objects; \
		pwd; \
		ar -x $(CHOMBO_HOME)/libamrelliptic$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtimedependent$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libboxtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libbasetools$(config).a; 
	ar -cvq libChomboLibs$(config).a  tmp_lib_objects/*.o
	rm tmp_lib_objects/*.o
	rmdir tmp_lib_objects

libBisicles$(config).a:  $(objects)  
	if test -f libBisicles$(config).a; then rm libBisicles$(config).a; fi 
	ar -cvq libBisicles$(config).a o/$(config)/*.o

lib: libBisicles$(config).a libChomboLibs$(config).a

all: lib

.PHONY: all lib
