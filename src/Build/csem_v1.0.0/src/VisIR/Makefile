############################################################################################
#
# GNU Makefile 
# This is a top-level GNU makefile of the JCSDA Community Surface Emissivity Modeling (CSEM)
# system. It is used as a general driver of the lower-tire makefiles to compile all 
# the F90/F03 module files of the Infrared surface emissivity modling subsystem. 
#
# Before running this makefile,  set the proper compiler
# e.g., source ifort.setup
#
# CREATION HISTORY:
#       Written by:     Ming Chen, 28-Feb-2016
#                       ming.chen@noaa.gov
#  
# 
############################################################################################

CSEM_HOME = ../..
include $(CSEM_HOME)/Build/make.macros

#.SILENT :

PWD=$(shell pwd)
PREFIX=$(PWD)


DIRS	= Land Water Snow Ice
	
all: 
	@echo; echo Builting up Infrared library from:
	@echo; echo $(PWD); echo;                      \
	for d in $(DIRS); do (cd $$d; $(MAKE)); done

install :
	@echo; echo installing Infrared library files from:
	@echo; echo $(PWD); echo;  \
	for d in $(DIRS); do (cd $$d; \
	 $(MAKE) install PREFIX=$(PREFIX) ); done

copy :
	@echo; echo installing Infrared library files from:
	@echo; echo $(PWD); echo;  \
	for d in $(DIRS); do (cd $$d; \
	 $(MAKE) copy PREFIX=$(PREFIX) ); done

clean :
	@echo; echo cleaning up Infrared library files ....;\
	for d in $(DIRS); do (cd $$d; $(MAKE) clean PREFIX=$(PREFIX) ); done
