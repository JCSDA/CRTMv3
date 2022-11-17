#!/bin/bash
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


make_lib(){
  #source ./Build/ifort.setup
  #source ./Build/gfortran.setup
  
  source ../CSEM_OLD/configuration/gfortran.setup
  export CSEM_HOME=/path to csem src
  export NETCDF_DIR=/path to netcdf lib & include
  export HDF5_DIR=/path to hdf5 lib & include
  ./configure --prefix=$1
  make clean
  make install
}

aclocal
if [ ! -f install-sh ]; then
  automake --add-missing --copy >/dev/null 2>&1   
  rm -f missing
fi
autoconf

if [ "$#" == 1 ]; then
  if [[ "$1" =~ 'prefix=' ]]; then
    prefix=$(echo $1| cut -d'=' -f 2)
    if [ "${prefix:0:1}" = '/' ]; then
      echo 'make library ...'
      make_lib $prefix
    else
      echo "$0 prefix=[absolute path]"
    fi
  fi
fi



