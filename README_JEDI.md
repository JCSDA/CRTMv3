README_JEDI.md

CRTM REL-3.1.0  Released October 31, 2023
CRTM REL-3.0.0  Released April 1, 2021



The README.md file contains a lot of general information about this repository and the legacy build system based on autotools.

CRTM REL-3.1.0 JEDI environment build instructions
=========================================================

Preamble
--------

CRTM v3.1.0 release (`REL-3.1.0`)  

This is a fully functional release of CRTM v3.1.0.

Basic requirements:  
(1) A Fortran 2003 compatible compiler.  
(2) A netCDF4 / HDF5 library.   
(3) A linux, macOS, or unix-style environment.  This has not been tested under any Windows Fortran environments.
(4) Bash shell is preferred.
(5) A suitable JCSDA JEDI environment: either HPC enabled, a JEDI container, spack-stack, or at a bare minimum ecbuild from ectools OR cmake.
(6) git and git-lfs > version 2.0 (tested on 2.10).

=========================================================

**Important Note**: If reading this, you're cloning the CRTM development repository.  The development repository is structured in a way that makes it less user friendly, but more amenable to development and testing.  You're also reading about the JEDI Environment instructions.  Please see README.md for all other uses.

In most cases, you'll be running CRTM inside of a JEDI bundle (e.g. fv3-bundle, ufo-bundle, etc.)  so you'll have no need to follow any directions here.  However, if you're interested in build CRTM stand-alone using a JEDI environment (i.e., for testing purposes, running the ctests, etc. ) continue reading.

=========================================================

Contents
========

1. Configuration  
2. Building the library  
3. Testing the library  
4. Installing the library  
5. Cleaning up  
6. Feedback and contact info  



Configuration, building, and testing the library
================================================  
JCSDA CRTM v3.1.0 Build Instructions
  
The CRTM **development** repository directory structure looks like:

<pre>
 .
  ├── LICENSE  (CC0 license)
  ├── COPYING  (CC0 legal document)
  ├── NOTES
  ├── README.md 
  ├── Get_CRTM_Binary_Files.sh (downloads the "fix" directory binary data from ftp, this is useful if you're doing out-of-jedi tests or if you want to override the default binary datasets. ) 
  ├── CMakeLists.txt           (top-level configuration file for ecbuild/cmake)
  ├── <b>configuration/</b>
  ├── <b>documentation/</b>
  ├── <b>fix/</b>
  │   ├── AerosolCoeff/
  │   ├── CloudCoeff/
  │   ├── EmisCoeff/
  │   ├── SpcCoeff/
  │   └── TauCoeff/
  ├── scripts/
  │   └── shell/
  ├── <b>src/</b>
  │   ├── Ancillary/
  │   ├── AntennaCorrection/
  │   ├── AtmAbsorption/
  │   ├── AtmOptics/
  │   ├── AtmScatter/
  │   ├── Atmosphere/
  │   ├── CRTM_Utility/
  │   ├── ChannelInfo/
  │   ├── Coefficients/
  │   ├── GeometryInfo/
  │   ├── InstrumentInfo/
  │   ├── Interpolation/
  │   ├── NLTE/
  │   ├── Options/
  │   ├── RTSolution/
  │   ├── SensorInfo/
  │   ├── SfcOptics/
  │   ├── Source_Functions/
  │   ├── Statistics/
  │   ├── Surface/
  │   ├── TauProd/
  │   ├── TauRegress/
  │   ├── Test_Utility/
  │   ├── User_Code/
  │   ├── Utility/
  │   ├── Validation/
  │   ├── Zeeman/
  └── <b>test/</b>
      └── Main/
</pre>

In the above list, the directories highlighted in bold (bold in markdown), are the key directories of interest to the casual developer.

JEDI Configuration
------------------
As of v3.1.0, binary data is obtained during the ecbuild/cmake step, it downloads a tarball from UCAR's GDEX service and unpacks it.  see `test/CMakeLists.txt`.   


**Configuration**
    git clone https://github.com/JCSDA/CRTMv3      (you've probably done this already)  
    cd CRTMv3
    git fetch 
    git pull

**Build Instructions**
<pre>
    mkdir build
    cd build
    cmake pathtocrtm 
</pre>
where `pathrocrtm` is where the `crtm/` diretory is located.  In this example if you're in the `crtm/build` directory, typing `cmake ..` will work.

<pre>
    make -j8     (-j8 means 8 parallel make processes, adjust the number to your machine)
    ctest -j8
</pre>
This should compile all of the source codes, create a libcrtm.so file, compile the tests, and finally run the various ctests.  If you're making changes to code, simply running the make command will detect your code changes and rebuild everything for you.  

Uninstalling the library
------------------------

To "uninstall" the library (assuming you haven't moved the installation directory contents somewhere else) you can type:
    cd build/
    rm -rf *  (make sure you do this in the build/ directory where you ran `cmake`)

Cleaning Up
-----------
<pre>
cd build/
make clean (removes compiled files, but not binary assets)
</pre>


**Additional options**
You can modify the various compiler flags, etc in the `CRTMv3/cmake/` directory.  There you will find several configuration files based on differen compilers.


**Feedback and Contact Information**

CRTM SUPPORT EMAIL: crtm-support@googlegroups.com OR visit https://forums.jcsda.org/

If you have problems building the library please include the generated "config.log" file in your email correspondence.

Known Issues
------------

(1) Any "Transmitance Coefficient" generation codes included in src/ are not functional.  Contact CRTM support above for details.  
(2) Testing was only done on modern gfortran compilers, with limited testing on intel fortran compilers.

Troubleshooting
---------------

TBD
