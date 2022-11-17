############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03/F08 module files of Microwave Water surface 
# modeling system. When the system is expanded/improved, users/developers may add their
# files to the OBJ_FILE list and recompile the system.
#
# CREATION HISTORY:
#       Written by:     Ming Chen, 28-Feb-2016
#                       ming.chen@noaa.gov
#  
# 
############################################################################################
CSEM_HOME = ../../..

include $(CSEM_HOME)/Build/make.macros
FC_FLAGS=$(FCFLAGS) 
ifeq ($(strip $(FC)),gfortran)
  FC_FLAGS=$(FCFLAGS) -Wno-unused -Wno-uninitialized -Wno-unused-dummy-argument
endif

PWD = $(shell pwd)
PREFIX  = $(shell (cd $(CSEM_HOME); pwd))
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include
LIB_FILE = libLandMW.a


# PATH to the shared source files
CSEM_DEFINE = $(CSEM_HOME)/src/CSEM_Define
UTILITY   = $(CSEM_HOME)/src/CSEM_Utility
FASTEM    = ./CRTM_FASTEM ./RTTOV_FASTEM ./ANN_FASTEM ./RSS_RTM
VPATH     = $(CSEM_DEFINE) $(UTILITY)  $(FASTEM)


# Macros of individual models
COMMON =  CSEM_Type_Kinds.o \
          CSEM_Struct_Define.o  \
          CSEM_Model_Manager.o  \
          CSEM_Exception_Handler.o 

CRTM_FASTEM = \
          CRTM_FASTEM_Constants.o \
	  CRTM_FASTEM_Parameters.o \
	  CSEM_Search_Utility.o \
	  CSEM_Interpolation.o \
	  FASTEM_Fresnel.o  Ellison.o  Guillou.o  Liu.o\
	  CSEM_FitCoeff_Define.o CSEM_Hyperbolic_Step.o Slope_Variance.o \
	  CRTM_MWwaterLUT_Define.o \
	  CRTM_MWwaterCoeff_Define.o \
          Azimuth_Emissivity_Module.o Azimuth_Emissivity_F6_Module.o \
          Foam_Utility_Module.o  \
          Reflection_Correction_Module.o \
          Small_Scale_Correction_Module.o  Large_Scale_Correction_Module.o \
          CRTM_LowFrequency_MWSSEM.o \
	  CRTM_Fastem1.o CRTM_Fastem3.o  CRTM_FastemXX.o  FASTEM_Coeff_Reader.o\
          CRTM_FASTEM_MODULE.o  

RTTOV_FASTEM  = \
	  rttov_fastem6_coef.o rttov_fastem6.o \
	  rttov_fastem6_tl.o rttov_fastem6_ad.o \
	  rttov_fastem5r1_coef.o rttov_fastem5r1.o \
	  rttov_fastem5r1_tl.o rttov_fastem5r1_ad.o \
	  RTTOV_FASTEM_MODULE.o
	
# Module files involved	   
OBJ_FILE = ${COMMON} ${CRTM_FASTEM} ${RTTOV_FASTEM}  \
	   CSEM_WaterMW_SfcOptics.o 
	
LIB_FILE = libWaterMW.a
all:  $(LIB_FILE) 
	

# Target definitions
# ------------------
# The main target
$(LIB_FILE):  $(OBJ_FILE) 
	@echo; echo; echo; \
	echo "=============================================="; \
	echo "Building  $(LIB_FILE) using:"; \
	echo "  Compiler: $(FC)"; \
	echo "  Switches: $(FC_FLAGS)"; \
	echo "  MW Water library at: $(PWD)"; \
	echo "=============================================="; \
	$(AR) $(AR_FLAG) $(LIB_FILE) $(OBJ_FILE)


copy: $(LIB_FILE) $(LIB_DIR) $(INC_DIR)
	@echo "=============================================="
	@echo "Install files ..."
	@echo "=============================================="
	@echo; if [ -f ${LIB_FILE} ]; then \
	  if [ "$(PREFIX)" == "$(PWD)" ];  then \
	     echo "move ${LIB_FILE} to ${LIB_DIR}"; \
	     (mv $(LIB_FILE)  $(LIB_DIR) ; rm -f $(OBJ_FILE));  \
	     mv *.mod  $(INC_DIR) ;  else \
	     echo "copy files to ${LIB_DIR}"; \
	     (cp -f $(OBJ_FILE)  $(LIB_DIR) ; rm -f $(OBJ_FILE));  \
	     (mv $(LIB_FILE)  $(LIB_DIR) ; mv *.mod  $(INC_DIR)); \
	  fi; fi

install: $(LIB_FILE) $(LIB_DIR) $(INC_DIR)
	@echo "=============================================="
	@echo "Install files ..."
	@echo "=============================================="
	@echo; if [ -f ${LIB_FILE} ]; then \
	  if [ "$(PREFIX)" == "$(PWD)" ];  then \
	     echo "move ${LIB_FILE} to ${LIB_DIR}"; \
	     (cp $(LIB_FILE)  $(LIB_DIR) );  \
	     mv *.mod  $(INC_DIR) ;  else \
	     echo "copy files to ${LIB_DIR}"; \
	     (cp -f $(OBJ_FILE)  $(LIB_DIR) );  \
	     (cp $(LIB_FILE)  $(LIB_DIR) ; cp *.mod  $(INC_DIR)); \
	  fi; fi

# A clean target
clean:
	@echo "=============================================="
	@echo "Cleaning up MW Water library files"
	@echo "=============================================="
	@echo; rm -f  $(OBJ_FILE)  *.o *.mod $(LIB_FILE) 
	@echo; rm -f  $(LIB_DIR)/$(LIB_FILE) 

# ...Create directories if necessary
$(LIB_DIR) \
$(INC_DIR) :
	- test -d $@ || mkdir $@

# Suffix rule definitions
# -----------------------
.SUFFIXES:
.SUFFIXES: .f90 .F90 .o
.f90.o:
	$(FC) -c $(FCFLAGS)  $(NETCDF_INC) $<
.F90.o:
	$(FC) -c $(FC_FLAGS) $(NETCDF_INC) $<

