############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03/F08 module files of Microwave Snow surface 
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

# PATH to the shared source files
CSEM_DEFINE = $(CSEM_HOME)/src/CSEM_Define
UTILITY   = $(CSEM_HOME)/src/CSEM_Utility
VPATH = $(CSEM_DEFINE) $(UTILITY) ../Soil 

# Macros of individual models
COMMON  = CSEM_Type_Kinds.f90 \
          CSEM_Struct_Define.f90  \
          CSEM_Model_Manager.f90  \
          CSEM_Exception_Handler.f90 \
          CSEM_Interpolation.f90 \
	  CSEM_String_Utility.f90

PHY_MODEL = \
          CSEM_Fresnel.f90 \
          MW_SoilWater_Permittivity.f90  \
          MW_Soil_Permittivity.f90 \
          SoilMW_Roughness_Correction.f90 \
          NESDIS_SnowEM_Parameters.f90 \
	  NESDIS_SnowMW_Optical_Model.f90 \
          NESDIS_SnowMW_PhyModel.f90  
 
SRC_SENSOR = \
	  NESDIS_SnowMW_LUT.f90 \
          NESDIS_SnowMW_Util.f90 \
	  NESDIS_SnowMW_AMSU_Module.f90 \
          NESDIS_SnowMW_MHS_Module.f90 \
          NESDIS_SnowMW_AMSRE_Module.f90 \
	  NESDIS_SnowMW_SSMI_Module.f90 \
	  NESDIS_SnowMW_SSMIS_Module.f90 \
          NESDIS_SnowEM_ATMS_Parameters.f90 \
	  NESDIS_SnowMW_ATMS_Module.f90 \
	  NESDIS_SnowMW_Sensors.f90 

# Module files involved	   
SRC_FILE = ${COMMON} ${PHY_MODEL} ${SRC_SENSOR} \
	   CSEM_SnowMW_SfcOptics.f90
	   
OBJ_FILE =  $(SRC_FILE:.f90=.o)

LIB_FILE = libSnowMW.a
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
	echo "  MW Snow library at: $(PWD)"; \
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
	@echo "Cleaning up MW Snow library files"
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
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FCFLAGS) $(NETCDF_INC)  $< 
	
