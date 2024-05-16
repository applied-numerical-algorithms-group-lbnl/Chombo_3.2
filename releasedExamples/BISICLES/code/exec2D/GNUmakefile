# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := driver

#all:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
#driver:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
all: clean libbike driver
driver: clean libbike

.NOTPARALLEL: all

include ../mk/Make.defs


tttt: tttt.txt
	echo "a" > tttt.txt

reallyclean: realclean libclean

libbike: 
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)

driver: clean libbike ../lib/libBisicles$(config).a 

alltest: libbike ../lib/libBisicles$(config).a 



.PHONY: all lib libbike libclean alltest


