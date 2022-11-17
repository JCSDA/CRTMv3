############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03/F08 module files of Microwave ice surface 
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
VPATH = $(CSEM_DEFINE) $(UTILITY) 

# Macros of individual models
COMMON  = CSEM_Type_Kinds.o \
          CSEM_Struct_Define.o  \
          CSEM_Model_Manager.o  \
          CSEM_Exception_Handler.o \
          CSEM_Interpolation.o \
	  CSEM_String_Utility.o
	  
PHY_MODEL = CSEM_Fresnel.o \
          NESDIS_IceMW_PhyModel.o   

OBJ_SENSOR = \
 	  NESDIS_IceMW_LUT.o \
          NESDIS_IceMW_Util.o \
	  NESDIS_IceMW_AMSU_Module.o\
          NESDIS_IceMW_MHS_Module.o \
          NESDIS_IceMW_AMSRE_Module.o\
	  NESDIS_IceMW_SSMI_Module.o\
	  NESDIS_IceMW_SSMIS_Module.o\
          NESDIS_ATMS_SeaICE_LIB.o \
	  NESDIS_IceMW_ATMS_Module.o \
          NESDIS_IceMW_Sensors.o 

# Module files involved
OBJ_FILE = ${COMMON} ${PHY_MODEL} ${OBJ_SENSOR} \
	   CSEM_IceMW_SfcOptics.o 

LIB_FILE = libIceMW.a
all: $(LIB_FILE) 


# Target definitions
# ------------------
# The main target
$(LIB_FILE):  $(OBJ_FILE) 
	@echo; echo; echo; \
	echo "=============================================="; \
	echo "Building  $(LIB_FILE) using:"; \
	echo "  Compiler: $(FC)"; \
	echo "  Switches: $(FC_FLAGS)"; \
	echo "  MW Ice library at: $(PWD)"; \
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
	@echo "Cleaning up MW Ice library files"
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

