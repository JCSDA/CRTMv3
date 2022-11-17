#!/bin/csh
#-------------------------------------------------------------------------------#
# DEBUG build settings for Linux g95 compiler
#-------------------------------------------------------------------------------#

setenv FC "g95"
setenv FCFLAGS "-fbounds-check -ffree-form -fno-second-underscore -ftrace=frame -malign-double -Wall"
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
