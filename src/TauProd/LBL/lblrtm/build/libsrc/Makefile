# libsrc/Makefile.  Generated from Makefile.in by configure.

# libsrc/ level makefile template

# Package-specific substitution variables
package = LBLRTM I/O
version = v1.0.0
tarname = lblrtmio
distdir = $(tarname)_$(version)

# Prefix-specific substitution variables
prefix = /scratch/users/pstegmann/workspace/crtm/src/TauProd/LBL/lblrtm/build

# Mixed package and prefix substitution variables
installdir = $(prefix)/$(distdir)

# Tool-specific substitution variables
FC       = ifort
FCFLAGS  = -r8 -i4 -g -FR -I/opt/netcdf4/4.6.2-intel-18.0.3/include -c
FPPFLAGS = -DREAL_SIZE=8 -DINT_SIZE=4
LDFLAGS  =  -L/opt/netcdf4/4.6.2-intel-18.0.3/lib
LIBS     =  -lnetcdff
AR       = ar
ARFLAGS  = crvs
RANLIB   = ranlib
INSTALL      = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644

# The library name
LIBRARY = lib$(tarname).a

# The file definitions. This include must occur before targets.
include make.filelist

# The targets
all: library

$(LIBRARY): $(OBJ_FILES)
	$(AR) $(ARFLAGS) $@ $(OBJ_FILES)
	$(RANLIB) $@

library: $(LIBRARY)

clean:
	cd test && $(MAKE) $@
	-rm *.o *.mod *.a >/dev/null 2>&1

distclean:
	cd test && $(MAKE) $@
	-rm Makefile >/dev/null 2>&1

check: library
	cd test && $(MAKE) $@

install:
	$(INSTALL) -d $(DESTDIR)$(installdir)/include
	$(INSTALL_DATA) *.mod $(DESTDIR)$(installdir)/include
	$(INSTALL) -d $(DESTDIR)$(installdir)/lib
	$(INSTALL_DATA) $(LIBRARY) $(DESTDIR)$(installdir)/lib
	$(INSTALL_DATA) ../config.log $(DESTDIR)$(installdir)

uninstall:
	-rm -fr $(DESTDIR)$(installdir) >/dev/null 2>&1

Makefile: Makefile.in ../config.status
	cd .. && ./config.status libsrc/$@

../config.status: ../configure
	cd .. && ./config.status --recheck

# Specify targets that do not generate filesystem objects
.PHONY: all clean distclean check install uninstall

# Dependency include files
include make.dependencies

# Suffix rules
.SUFFIXES:
.SUFFIXES: .fpp .f90 .F90 .o
.fpp.o:
	$(FC) $(FCFLAGS) $(EXTRA_FCFLAGS) $(FPPFLAGS) $<

.f90.o:
	$(FC) $(FCFLAGS) $(EXTRA_FCFLAGS) $<

.F90.o:
	$(FC) $(FCFLAGS) $(EXTRA_FCFLAGS) $(FPPFLAGS) $<
