#
# Makefile for TauCoeff_Inspect
#

# Define the source file
SRC_FILE = TauCoeff_Inspect.f90

# Define the object file
OBJ_FILE = $(SRC_FILE:.f90=.o)

# Define the executable file
EXE_FILE = TauCoeff_Inspect

# Define the include and library paths
# Assuming the standard CRTM build structure
INCLUDES = -I../../../../build/module/crtm/GNU/13.3.0
LIBS = -L../../../../build/lib -lcrtm

# Compiler
FC = gfortran
FFLAGS = -O2

# Build rule
all: $(EXE_FILE)

$(EXE_FILE): $(OBJ_FILE)
	$(FC) $(FFLAGS) -o $(EXE_FILE) $(OBJ_FILE) $(LIBS)

$(OBJ_FILE): $(SRC_FILE)
	$(FC) $(FFLAGS) $(INCLUDES) -c $(SRC_FILE)

clean:
	rm -f $(OBJ_FILE) $(EXE_FILE)
