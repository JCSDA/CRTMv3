############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03 module files of Infrared Water surface 
# modeling system. When the system is expanded/improved, users/developers may add their
# files to the OBJ_FILE list and recompile the system.
#
# CREATION HISTORY:
#       Written by:     Ming Chen, 28-Feb-2015
#                       ming.chen@noaa.gov
#  
# 
############################################################################################

include ../Build/make.macros

# PATH to the shared source files
CSEM_DIR = ..

INC_FILE = -I$(CSEM_DIR)/include
LIB_FILE = $(CSEM_DIR)/lib/libcsem.a 

OBJ_FILE = 

LDFLAGS_HDF5 =  $(HDF5_LIB) 
FFLAGS_HDF5  =  $(HDF5_INC) 

# Main files
EXE_LandIR=TEST_CSEM_LandIR_SfcOptics
EXE_LandVis=TEST_CSEM_LandVis_SfcOptics
EXE_LandMW=TEST_CSEM_LandMW_SfcOptics
EXE_WaterIR=TEST_CSEM_WaterIR_SfcOptics
EXE_WaterIR_BRDF=TEST_CSEM_WaterIR_BRDF
EXE_WaterMW=TEST_CSEM_WaterMW_SfcOptics
EXE_SnowIR=TEST_CSEM_SnowIR_SfcOptics
EXE_SnowMW=TEST_CSEM_SnowMW_SfcOptics
EXE_IceMW=TEST_CSEM_IceMW_SfcOptics
EXE_IceIR=TEST_CSEM_IceIR_SfcOptics

EXE_FILE = $(EXE_SnowIR)

all:  $(EXE_FILE)

$(EXE_FILE): $ ${OBJ_FILE} ${EXE_FILE}.o
	$(FL) $(NETCDF_LIB)  $(FLFLAGS) -o $(EXE_FILE) \
	${OBJ_FILE} ${EXE_FILE}.o $(LIB_FILE) 



# A clean target
clean:
	@echo "=============================================="
	@echo "Cleaning up IR Water library files"
	@echo "=============================================="
	@echo; rm -f $(OBJ_FILE)  *.o *.mod ${EXE_FILE} 

# Suffix rule definitions
# -----------------------
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FCFLAGS) $(INC_FILE) $(NETCDF_INC)  $<
