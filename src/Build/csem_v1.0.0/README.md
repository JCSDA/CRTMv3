The Community Surface Emissivity Model (CSEM) is a highly modularized Earth’s surface RT modeling system based on Object-Oriented Programming (OOP) design.  It evolved from the surface modules of the Community Radiative Transfer Model (CRTM), but with completely redesigned model structure to facilitate the implementation of various surface RT models. CSEM provides the surface emissivity and reflectivity simulations of diverse surface types in the spectral range from the ultraviolet, visible to microwave bands.  Enclosed in CSEM are not only the physical models based on sound radiative transfer equations, but also a variety of empirical and semi-empirical models, type-based emissivity lookup tables (LUT), and global emissivity atlases from satellite retrievals. The object-oriented design provides very flexible software interfaces for implementing and testing new model components.  Multiple model options of the same kind (e.g., microwave land models) may be easily implemented and accommodated in the CSEM framework. In practical applications, CSEM may be used as a standalone surface RT research tool, or used as a sub-system to provide the surface radiative conditions for CRTM, significantly leveraging the development and improvement efforts of CRTM. 


How to build CSEM:

1) Go to Build/env.setup and source the specific compiler config file

   e.g., source gfortran.setup

   note: make sure the following three environment variables are
   already defined or included in the setup file.
 
   B shell (sh, bash)


         export NETCDF_HOME=path to the netcdf

         export HDF5_HOME=path to the hdf5

    C shell (csh)


         setenv NETCDF_HOME  path to the netcdf

         setenv HDF5_HOME    path to the hdf5

    run step 1) for the first time fresh installation or as long as
    these three ENV variables are changed.

2) Generate the file "configure"
    ./autogen.sh
3) Generate the file "Makefile", you may specify where the CSEM library will be
   installed. The default is the current directory
    ./configure --prefix=path for the CSEM library to be installed
4) make
5) make install

