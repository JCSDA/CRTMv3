####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")

if( HAVE_AUTOPROFILE )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finstrument-functions")
endif( )

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -unroll -inline -no-heap-arrays" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -check bounds -traceback -warn -heap-arrays -fpe-all=0 -fpe:0 -ftz -check all -no-inline" )

####################################################################
# RELEASE WITH DEBUG INFO FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -ip -unroll -inline -g -DNDEBUG -check bounds -traceback -no-heap-arrays -assume byterecl" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -ip -ipo -unroll -inline -no-heap-arrays" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################
# OpenMP (gated on the top-level OPENMP option)
####################################################################

if( OPENMP )
  set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE} -qopenmp" )
  set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG} -qopenmp" )
  set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -qopenmp" )
  set( CMAKE_Fortran_FLAGS_BIT            "${CMAKE_Fortran_FLAGS_BIT} -qopenmp" )
endif()

####################################################################

  
