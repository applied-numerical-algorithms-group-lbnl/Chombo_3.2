
# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
CHOMBO_HOME=../../../../../../Chombo/lib/


BISICLES_HOME=../../../code
BISICLES_BRANCH=../../../


DIM=2

#include $(BISICLES_HOME)/mk/Make.defs

SVN_REV=$(shell cd $(BISICLES_BRANCH) ; svnversion -n)
ifeq ($(SVN_REV),)
SVN_REV := unknown
SVN_REP := unknown
SVN_URL := unknown
else

SVN_REP := $(shell cd $(BISICLES_BRANCH) ; svn info | awk '/Repository Root:/ {print $$3}')
ifeq ($(SVN_REP),)
SVN_REP := unknown
endif

SVN_URL := $(shell cd $(BISICLES_BRANCH) ; svn info | awk '/^URL:/ {print $$2}')
ifeq ($(SVN_URL),)
SVN_URL := unknown
endif

endif

ebase :=  createSectorsFromVel
#ebase :=  createMaskFromVel

#ifneq ($(PYTHON_INC),)
#ebase += pythonf
#endif

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRElliptic AMRTimeDependent AMRTools BoxTools 

base_dir = .

SRC_DIR = $(BISICLES_HOME)/src

UTIL_DIR = $(BISICLES_HOME)/util

BASE_DIR = ../

EXEC_DIR = .

src_dirs = $(SRC_DIR) $(UTIL_DIR)


include $(CHOMBO_HOME)/mk/Make.example

#ifneq ($(PYTHON_INC),)
#CPPFLAGS+=-I$(PYTHON_INC)	
#LIBFLAGS+=$(PYTHON_LIBS)
#endif

CPPFLAGS += -DSVN_REV='"$(SVN_REV)"'
CPPFLAGS += -DSVN_REP='"$(SVN_REP)"'
CPPFLAGS += -DSVN_URL='"$(SVN_URL)"'

