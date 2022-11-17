############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03/F08 module files of Vis-Infrared Snow surface 
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

PWD = $(shell pwd)
PREFIX  = $(shell (cd $(CSEM_HOME); pwd))
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include


LDFLAGS_HDF5 =  $(HDF5_LIB) 
FFLAGS_HDF5  =  $(HDF5_INC)


# PATH to the shared source files
CSEM_DEFINE = $(CSEM_HOME)/src/CSEM_Define
UTILITY   = $(CSEM_HOME)/src/CSEM_Utility
VPATH = $(CSEM_DEFINE) $(UTILITY) ../LUT_Atlas/CRTM_LUT 

# Macros of individual models
COMMON =  CSEM_Type_Kinds.o \
          CSEM_Struct_Define.o  \
          CSEM_Model_Manager.o  \
          CSEM_Exception_Handler.o \
          CSEM_Interpolation.o 
	  
NPOESS_LUT = \
          NPOESS_LUT_Reader.o \
	  NPOESS_LUT_Module.o 
	  
PHY_MODEL  = \
          NESDIS_SnowIR_PhyModel.o\
          NESDIS_SnowVIS_PhyModel.o\

# Module files involved
OBJ_FILE = $(COMMON) ${NPOESS_LUT} ${PHY_MODEL} \
	  CSEM_SnowIR_SfcOptics.o  \
	  CSEM_SnowVIS_SfcOptics.o 

RTTOV_LIB = $(RTTOV_DIR)/libRTTOV_VisIR.a
LIB_FILE = libSnowVisIR.a

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
	echo "  VisIR Snow library at: $(PWD)"; \
	echo "=============================================="; \
	$(AR) $(AR_FLAG) $(LIB_FILE) $(OBJ_FILE)


copy: $(LIB_FILE) $(LIB_DIR) $(INC_DIR)
	@echo "=============================================="
	@echo "Install files ..."
	@echo "=============================================="
	@echo; if [ -f ${LIB_FILE} ]; then \
	  if [ "$(PREFIX)" == "$(PWD)" ];  then \
	     rm -f  $(LIB_DIR)/*.o;\
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
	@echo "Cleaning up VisIR Snow library files"
	@echo "=============================================="
	@echo; rm -f  $(OBJ_FILE)  *.o *.mod $(LIB_FILE)
	@echo; rm -f  $(LIB_DIR)/$(LIB_FILE) 

# ...Create directories if necessary
$(INC_DIR) \
$(LIB_DIR) :
	- test -d $@ ||  mkdir -p $@


$(RTTOV_LIB):
	@echo; echo  library files ....;\
	cd $(RTTOV_DIR); $(MAKE) ; \

# Suffix rule definitions
# -----------------------
.SUFFIXES:
.SUFFIXES: .f90 .f .F90 .o
.f90.o:
	$(FC) -c $(FCFLAGS) $(FFLAGS_HDF5)  $(NETCDF_INC) $<
.f.o:
	$(FC) -c $(NETCDF_INC) $<
.F90.o:
	$(FC) -c $(FCFLAGS) $(FFLAGS_HDF5)  $(NETCDF_INC) $<

