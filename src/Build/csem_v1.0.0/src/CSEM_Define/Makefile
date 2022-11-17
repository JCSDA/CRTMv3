############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03/F08 module CSEM_LifeCycle. 
#
# CREATION HISTORY:
#       Written by:     Ming Chen, 28-Feb-2015
#                       Ming.chen@noaa.gov
#  
# 
############################################################################################
CSEM_HOME = ../..
include $(CSEM_HOME)/Build/make.macros
PREFIX=$(CSEM_HOME)
# PATH to the shared source files
INC_DIR = $(PREFIX)/include
LIB_DIR = $(PREFIX)/lib
OBJ_FILE = CSEM_LifeCycle.o	   

all:  $(OBJ_FILE)


install: $(OBJ_FILE)  $(LIB_DIR) $(INC_DIR)
	@echo "=============================================="
	@echo "Install files ..."
	@echo "=============================================="
	@echo; if [ -f ${OBJ_FILE} ]; then \
	     echo "copy files to ${LIB_DIR}"; \
	     (cp -f $(OBJ_FILE)  $(LIB_DIR) ; rm -f $(OBJ_FILE));  \
	     ( mv *.mod  $(INC_DIR)); \
             fi

# A clean target
clean:
	@echo "=============================================="
	@echo "Cleaning up "
	@echo "=============================================="
	@echo; rm -f $(OBJ_FILE)  *.o *.mod

# ...Create directories if necessary
$(LIB_DIR) \
$(INC_DIR) :
	- test -d $@ || mkdir $@

# Suffix rule definitions
# -----------------------
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INC_DIR) $<
