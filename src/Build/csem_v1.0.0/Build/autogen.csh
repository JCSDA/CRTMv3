#!/bin/csh
#This script is used to generate the configure 
#and make files which are used to compile and
#install the CSEM package 
#
# CREATION HISTORY:
#       Written by:     Ming Chen, 28-Feb-2018
#                       ming.chen@noaa.gov
#  
# 
############################################################################################


setenv CSEM_HOME "/scratch4/NCEPDEV/jcsda/save/Ming.Chen/CRTM/CSEM"
setenv NETCDF_HOME "/apps/netcdf/4.3.0-intel"
setenv HDF5_HOME   "/apps/hdf5/1.8.14-intel"

if ( ! -d ${CSEM_HOME}/src ) then
  echo "No CSEM package is found under CSEM_HOME=${CSEM_HOME}..."
  exit 1
endif

source ${CSEM_HOME}/Build/env.setup/ifort.setup.csh.theia


if ( $#argv > 1 ) then
   echo "$0 [prefix=absolute path]"
   exit 1
endif

if ( $#argv == 0 ) then
  set prefix=`pwd`
else 
  if ( $1 =~ *"prefix="* ) then
    set prefix=`echo $1| cut -d'=' -f 2`
    set A=`echo $prefix |awk '{print substr($0, 0, 1)}'`
    if ( "$A" != "/" ) then
      echo "$0 prefix=[absolute path]"
      exit 1
    endif
  else
      echo "$0 prefix=[absolute path]"
      exit 1
  endif
endif

if ( -f 'aclocal.m4' ) then
  rm -f aclocal.m4
endif
if ( -d 'autom4te.cache' ) then
  rm -fR autom4te.cache
endif

aclocal
if ( -f configure ) then
  rm -f configure
endif

autoconf

./configure --prefix=$prefix
make clean
make
make install



