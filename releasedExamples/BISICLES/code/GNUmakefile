# -*- Mode: Makefile -*- 

### This makefile requires GNU make.
### It builds and run the IAMR examples.
### This file is essentially the same as Chombo/example/GNUmakefile.

### User-definable variables are in UPPER case.
### Internal variables are in lower case.
###
### Usage:
###   target:    purpose:
###   -------    ----------
###   all        builds all the executables, including the libraries they depend on
###   test       builds all the executables, but not the libraries
###   run        runs all the executables
###   clean      deletes files for this configuration
###   realclean  deletes files for all configurations that weren't in the original distribution
###
### Method:
###   The way this makefile works is the target specified by the user
###   on the command line is treated as the `action' to perform on all
###   the subdirectories in $dirs.  In the rules, the subdirectories are
###   treated like targets and the action is specified by a variable.
###   The makefile builds each directory target by invoking make
###   recursively on the directory with the action as the recursive target.
###   This assumes that the makefiles in the subdirectories implement all
###   the actions in the $(_all_actions) variable.
###
### Notes:
###
##################################################################################

export makefiles+=bisicles

## This causes all commands to be executed with `csh' instead of the
## default `sh'.  The `-f' option makes csh load faster by skipping
## the user's init files (we don't want them anyway).
##[NOTE: do _not_ change this unless you really know what you're doing!]
SHELL = /bin/csh -f

## This variable specifies where to find the Chombo libraries and source dirs
##export CHOMBO_HOME := ../../../../Chombo/lib

##
## Default values of configuration variables.
##
##include $(CHOMBO_HOME)/mk/Make.defs.defaults


## This variable defines the names of the subdirectories to be built by default.
## One or more can be used as a target on the command line.
dirs = lib exec2D benchmark filetools cdriver test 


##
## These are the actions (targets for the recursive make) supported by this
## makefile.  It assumes the makefiles in the subdirectories support them too.
##
_all_actions = all example run clean

.PHONY: $(_all_actions) $(dirs) doc doxygen usage vars

##
## Default target
##
usage:
	@echo "This makefile supports targets:"
	@echo "     $(_all_actions) $(dirs) doc doxygen realclean usage vars"

##
## For each action, set the action variable to tell the recursive
## makefile what to do.  There should be a separate dependency for
## every value in $(_all_actions).
##[NOTE: the `NODEPENDS=TRUE' variable disables automatic
##       dependency analysis in the recursive makefile, which is only
##       needed for actions that compile or link code.]
##
_action = all
all       : _action = all
example   : _action = all
run       : _action = run       NODEPENDS=TRUE
clean     : _action = clean     NODEPENDS=TRUE

# this dependency ensures the action is done for every directory.
$(_all_actions) : $(dirs)

# this rule performs the action
$(dirs):
	$(MAKE) --no-print-directory --directory $@ -f GNUmakefile $(_action)

realclean:
	-find . -type d \( -name d -o -name p -o -name f -o -name o \) -depth -exec $(RM) -r {} \;
	-find . \( -name \*.ex -o -name \*.a -o -name \*_F.H -o -name core \) -exec $(RM) {} \;
	-find . \( -name \*.hdf5 -o -name \*.h5 -o -name \*~ -o -name work.pc\* -o -name gmon.out -o -name pout.\[0-9-\]\* -o -name leak.out\* -o -name .\#\* \) -exec $(RM) {} \;
	-$(RM) -r doc/doxygen
	$(MAKE) --no-print-directory --directory doc $@

# this rule prints out various 'make' variables to the screen for debugging
vars:
	@$(MAKE) --no-print-directory -f $(CHOMBO_HOME)/mk/Make.printVariables $@

doc:
	cd doc ; ./mkdoc

ifeq ($(DOXYGEN),)
DOXYGEN=doxygen
endif

doxygen:
	mkdir -p doc/doxygen
	$(DOXYGEN) doxygen.config ; echo "" ; echo "point browser at doc/doxygen/html/index.html" ;

.NOTPARALLEL: all
