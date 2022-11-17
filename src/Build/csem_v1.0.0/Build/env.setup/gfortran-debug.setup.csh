#!/bin/csh
#-------------------------------------------------------------------------------#
# DEBUG build settings for Linux gfortran compiler
#-------------------------------------------------------------------------------#

setenv FC "gfortran"
setenv FCFLAGS "-fbounds-check -fimplicit-none -ffpe-trap=overflow,zero,invalid -ffree-form -fno-second-underscore -frecord-marker=4 -ggdb -Wall -Wno-conversion -std=f2008"
setenv LDFLAGS ""
setenv LIBS ""

setenv NETCDF_HOME "/data/jcsda/mchen/netcdf-fortran-4.4.4/gfortran"
setenv HDF5_HOME   "/data/jcsda/mchen/hdf5-1.8.20-gfortran"

echo "========================================="
echo " CSEM compilation environment variables:"
echo "   FC:       ${FC}"
echo "   FCFLAGS: ${FCFLAGS}"
echo "   FL:       ${FL}"
echo "   FLFLAGS: ${FLFLAGS}"
echo "========================================="
echo
