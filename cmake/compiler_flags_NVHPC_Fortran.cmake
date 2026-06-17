####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_REAL8_ -Mfreeform -Mbackslash" )

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-fast -O3" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Mbounds -Ktrap=fp -traceback" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -Kieee" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################
# OpenMP (gated on the top-level OPENMP option)
####################################################################

if( OPENMP )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -mp" )
  set( CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -mp" )
endif()
