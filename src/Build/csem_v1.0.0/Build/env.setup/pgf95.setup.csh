#!/bin/csh
#-------------------------------------------------------------------------------#
# PRODUCTION build settings for Linux pgf95 compiler
#-------------------------------------------------------------------------------#

setenv FC "pgf95"
setenv FCFLAGS "-g -fast"
setenv LDFLAGS ""
setenv LIBS ""

setenv NETCDF_HOME ""
setenv HDF5_HOME   ""

echo "========================================="
echo " CSEM compilation environment variables:"
echo "   FC:       ${FC}"
echo "   FCFLAGS: ${FCFLAGS}"
echo "   FL:       ${FL}"
echo "   FLFLAGS: ${FLFLAGS}"
echo "========================================="
echo
