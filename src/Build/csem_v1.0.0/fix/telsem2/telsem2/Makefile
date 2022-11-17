FC=gfortran
FCFLAGS=-g -fPIC -fbounds-check -finit-integer=-1 -finit-real=nan -finit-character=127 -std=f95
OBJS=mod_mwatlas_m2.o test_telsem2.o
MODS=mod_mwatlas_m2.mod

.SUFFIXES: .F90
%.o : %.mod

.F90.o :
	$(FC) -c $< -o $@ $(FCFLAGS)

all:$(OBJS)
	$(FC) $^ -o test_telsem2

clean:
	rm -f $(OBJS) $(MODS) test_telsem2
