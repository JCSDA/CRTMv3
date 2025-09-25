
if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif()

#######################################################################################
# Fortran
#######################################################################################

message(STATUS "Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Compiler: ${CMAKE_Fortran_COMPILER}")

set(CMAKE_FORTRAN_STANDARD 08)
set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(CMAKE_FORTRAN_EXTENSIONS OFF)

if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  include( compiler_flags_GNU_Fortran )

elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" )
  include( compiler_flags_IntelLLVM_Fortran )

elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  include( compiler_flags_Intel_Fortran )

elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
  include( compiler_flags_NVHPC_Fortran )

elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "XL" )
  include( compiler_flags_XL_Fortran )

elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Cray" )
  include( compiler_flags_Cray_Fortran )

else()
  message( STATUS "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()
