#!/bin/csh
#-------------------------------------------------------------------------------#
# DEBUG build settings for Linux ifort compiler
#-------------------------------------------------------------------------------#

setenv FC "ifort"
setenv FCFLAGS "-g -check bounds -e08 -traceback -free -assume byterecl,realloc_lhs -fp-stack-check -mieee-fp"
setenv LDFLAGS ""
setenv LIBS ""

setenv NETCDF_HOME "/data/jcsda/mchen/netcdf-fortran-4.4.4/ifort"
setenv HDF5_HOME   "/data/jcsda/mchen/hdf5-1.8.20-ifort"
echo "========================================="
echo " CRTM compilation environment variables:"
echo "   FC:       ${FC}"
echo "   FCFLAGS: ${FCFLAGS}"
echo "   FL:       ${FL}"
echo "   FLFLAGS: ${FLFLAGS}"
echo "========================================="
echo
