####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl" )

if( HAVE_AUTOPROFILE )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finstrument-functions" )
endif()

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE
     "-O3 -unroll -no-heap-arrays -assume byterecl" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG
     "-O0 -g -check bounds -traceback -warn all -no-heap-arrays -fpe0 -ftz -check all -assume byterecl" )

####################################################################
# RELWITHDEBINFO FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO
     "-g -O0 -traceback" )
#     "-O2 -g -DNDEBUG -check bounds -traceback -no-heap-arrays -assume byterecl" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT
     "-O2 -no-heap-arrays -fp-model strict -assume byterecl" )

####################################################################
# OpenMP (gated on the top-level OPENMP option)
####################################################################

if( OPENMP )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -qopenmp" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -qopenmp" )
  set( CMAKE_Fortran_FLAGS_BIT     "${CMAKE_Fortran_FLAGS_BIT} -qopenmp" )
endif()

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS )

####################################################################
# Notes:
# - `-qopenmp` enables OpenMP for `ifx`, including linking the necessary runtime
# - For GPU offload: Add `-fopenmp-targets=spir64` (Intel GPU) or `-fopenmp-targets=nvptx64-nvidia-cuda` (NVIDIA)
# - Reproducibility: `-fp-model strict` ensures bitwise results but may hurt perf
# - Consider `-diag-disable=10448` if `ifx` emits unnecessary diagnostic noise
