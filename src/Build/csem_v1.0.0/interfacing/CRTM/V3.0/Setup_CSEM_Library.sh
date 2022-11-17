#!/bin/bash
############################################################################################
#
# Setup_CSEM_Library.sh 
# This bash shell script file is used to automate the modification/replacement
# of the related CRTM files for CRTM to use the CSEM libraray. 
#
# CRTM_ROOT and CSEM_HOME need to be specified accordingly.
#
#       Written by:     Ming Chen, 22-Feb-2022
#                       ming.chen@noaa.gov
#  
# 
############################################################################################

export CRTM_ROOT=/home/Ming.Chen/save/CRTM/REL-3.0.0/qliu_v3.0_merge.new
export CSEM_HOME=/home/Ming.Chen/save/CRTM/CSEM1.0.0

cd ${CRTM_ROOT} && source ${CRTM_ROOT}/Set_CRTM_Environment.sh
export export PATH=~/bin:$PATH
if [ ! -d $CRTM_SOURCE_ROOT ]; then
   echo "Error: $CRTM_SOURCE_ROOT not exists .."
   exit 1
fi

# backup the original crtm files
CRTM_BKUP=$CRTM_ROOT/src.ori
if [[ ! -d ${CRTM_BKUP} ]]; then
  mkdir -p ${CRTM_BKUP}
  mv $CRTM_SOURCE_ROOT/CRTM_LifeCycle.f90 ${CRTM_BKUP}/
  mv $CRTM_SOURCE_ROOT/SfcOptics/CRTM_SfcOptics.f90 ${CRTM_BKUP}/
  mv $CRTM_SOURCE_ROOT/Build/configure.ac ${CRTM_BKUP}/
  mv $CRTM_SOURCE_ROOT/Build/Makefile.in ${CRTM_BKUP}/
  mv $CRTM_SOURCE_ROOT/Build/libsrc/Makefile.in ${CRTM_BKUP}/Makefile.in.libsrc
  mv $CRTM_SOURCE_ROOT/Build/libsrc/make.filelist ${CRTM_BKUP}/
  mv $CRTM_SOURCE_ROOT/Build/libsrc/make.dependencies ${CRTM_BKUP}/
  mv $CRTM_SOURCE_ROOT/Build/libsrc/make.rules ${CRTM_BKUP}/
  # CMakefiles
  mv $CRTM_ROOT/CMakeLists.txt ${CRTM_BKUP}/top.CMakeLists.txt
  mv $CRTM_SOURCE_ROOT/CMakeLists.txt ${CRTM_BKUP}/src.CMakeLists.txt

fi

#CSEM interfacing files
ifile_dir=${CSEM_HOME}/interfacing/CRTM
interface_files=(\
 CRTM_LifeCycle.f90.csem CRTM_SfcOptics.f90.csem \
 make.dependencies.csem make.filelist.csem \
 configure.ac.csem Makefile.in.csem Makefile.in.libsrc \
 top.CMakeLists.csem src.CMakeLists.csem)
for file in "${interface_files[@]}"; do
  #echo $file
  ifile=${ifile_dir}/$file
  if [ ! -f $ifile ]; then
    echo "Error: $ifile not exists .."
    exit 2
  fi
done
# replace the crtm files
cp -f ${ifile_dir}/CRTM_LifeCycle.f90.csem $CRTM_SOURCE_ROOT/CRTM_LifeCycle.f90
cp -f ${ifile_dir}/CRTM_SfcOptics.f90.csem $CRTM_SOURCE_ROOT/SfcOptics/CRTM_SfcOptics.f90

# automake files
cp -f ${ifile_dir}/configure.ac.csem $CRTM_SOURCE_ROOT/Build/configure.ac
cp -f ${ifile_dir}/Makefile.in.csem $CRTM_SOURCE_ROOT/Build/Makefile.in
cp -f ${ifile_dir}/Makefile.in.libsrc  $CRTM_SOURCE_ROOT/Build/libsrc/Makefile.in
cp -f ${ifile_dir}/make.filelist.csem  $CRTM_SOURCE_ROOT/Build/libsrc/make.filelist
cp -f ${ifile_dir}/make.dependencies.csem $CRTM_SOURCE_ROOT/Build/libsrc/make.dependencies
# CMakefile
cp -f ${ifile_dir}/top.CMakeLists.csem $CRTM_ROOT/CMakeLists.txt
cp -f ${ifile_dir}/src.CMakeLists.csem $CRTM_SOURCE_ROOT/CMakeLists.txt

# remove CRTM_IRlandCoeff module 
sed -i '/CRTM_IRlandCoeff/d' $CRTM_SOURCE_ROOT/CRTM_Module.F90


# default csem_v1.0.0 is put/linked to the subdirectory
#  $CRTM_SOURCE_ROOT/Build 
csem_link=$CRTM_SOURCE_ROOT/Build/csem_v1.0.0
if [[ ! -d $csem_link || ! -L $csem_link ]]
then 
  ln -s $CSEM_HOME $csem_link
fi 
if [[ ! -f $csem_link/lib/libcsem.a  ]]
then 
  echo "$csem_link/lib/libcsem.a not exist .."
  echo "CSEM library needs to be pre-built .."
  exit 4
fi 

make_crtm_lib(){
 if [ "$#"  -lt 1 ]; then
   echo "usgae: make_crtm_lib CRTM_ROOT ..."
   exit 1 
 fi
 local CRTM_ROOT="$1"
  # specify compiler and netcdf
  compiler=$CRTM_ROOT/configuration/ifort.setup
  test -f $compiler && source $compiler || ("echo Error: $compiler" ; exit 4) 
  export NETCDF_ROOT=/apps/netcdf/4.7.0/intel/18.0.5.274
  export FCFLAGS="-I${NETCDF_ROOT}/include $FCFLAGS"
  export LDFLAGS="-L${NETCDF_ROOT}/lib -lnetcdff $LDFLAGS" 

  # make crtm library
  cd $CRTM_SOURCE_ROOT
  make
  cd $CRTM_SOURCE_ROOT/Build
  autogen.sh
  ./configure --prefix=$CRTM_SOURCE_ROOT/Build
  make install
}

#make_crtm_lib $CRTM_ROOT 
