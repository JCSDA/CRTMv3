# (C) Copyright 2026 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


# Set compiler flags for basic build types,
# for compilers where this is not provided by ecbuild.
include(build_type_compiler_flags)

# Set JEDI's common compiler flags
include(jedi_common_compiler_flags)

# Set CRTM-specific compiler flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  ecbuild_add_fortran_flags("-ffree-line-length-none")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  ecbuild_add_fortran_flags("-assume byterecl")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
  # to match historical behavior, disable RelWithDebugInfo optimization for IntelLLVM
  # TODO: determine source of problem and resolve in a more elegant way, then return
  #       to default compiler flags so that RELWITHDEBINFO does what it says
  set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O0 -g -DNDEBUG -traceback -heap-arrays 32" )
  ecbuild_add_fortran_flags("-assume byterecl")
endif()
