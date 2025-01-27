CRTM REL-3.1.1
====================

[![Build Status](https://app.travis-ci.com/JCSDA-internal/crtm.svg?token=r6aaq9P13fHcTi8yBgdM&branch=develop)](https://app.travis-ci.com/JCSDA-internal/crtm)

Preamble
--------

CRTM v3.1.1 release (`REL-3.1.1`)

v3.1.1 released August 12, 2024
v3.1.0 (alpha) Released October 31, 2023
v3.0.0 Released March, 2023  
v2.4.1-alpha Released on April 1, 2021 (internal realease only)
v2.4.0 Released on October 23, 2020

This is an experimental release of CRTM v3.0, some features may not be fully functional. Contact crtm-support@googlegroups.com.  
v3.x features will be rolled out in incremental updates. 

Basic requirements:  
(1) A Fortran 2008 compatible compiler
(2) A netCDF4 / HDF5 library that has been built with the compiler you're going to use (module environments are helpful here)
(3) A linux, macOS, or unix-style environment.  This has not been tested under any Windows Fortran environments.
(4) Bash shell is preferred. 
(5) git and git-lfs (minimum version TBD, but has been tested on git-lfs v2.10 and higher )
(6) cmake / make build system

=========================================================

**JEDI NOTE** This release branch is also designed to work directly in a JEDI container or JEDI environment. If you're doing JEDI things, you're probably in the right spot. However, you should stop reading right now and have a look at the README_JEDI.md file.   

If you're looking for an older version of CRTM (v2.3.0 or older) you should obtain the appropriate tarball from
https://bin.ssec.wisc.edu/pub/s4/CRTM/   OR https://github.com/JCSDA/crtm (old versions).   

If you're looking for version 2.4.0 or newer in a structure similar to older CRTM tarball releases, you should check out the appropriate release/ branch.
`git branch --remote | grep "release/"` to see a list of release branches OR you may checkout the appropriate tag on the master branch and build it yourself. 

Finally, you may follow the instructions here to build a "latest" release based on the most recent developments.

=========================================================

Contents
========

1. Configuration  
2. Building the library  
3. Testing the library  
4. Installing the library  
  a. GNU Install  
      - Linking to the library  
  b. Uninstalling the library  
5. Cleaning up  
6. Feedback and contact info  



Configuration, building, and testing the library
================================================  
JCSDA CRTM v3.1.1 Build Instructions

The CRTM repository directory structure looks (something) like:

<pre>
 .
  ├── LICENSE  (Public Domain)
  ├── COPYING  (Public Domain)
  ├── NOTES
  ├── README.md 
  ├── Get_CRTM_Binary_Files.sh  
  ├── <b>configuration/</b>
  ├── <b>documentation/</b>
  ├── <b>fix/</b>
  │   ├── AerosolCoeff/
  │   ├── CloudCoeff/
  │   ├── EmisCoeff/
  │   ├── SpcCoeff/
  │   └── TauCoeff/
  |── <b>src/</b>
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
A user is only likely to be interested in creating a "build" or use a previously created build (see releases/* on the github.com repository).

A typical "build release" of CRTM (what you would normally find in a tarball and see in libraries) is what will be contained under the `src/Build` directory after successful compilation.
But after a clean clone of the development repository, none of the links to source code have been created yet under `src/Build`.   To get there, follow the next steps.

Configuration
-------------
At present, the `fix/` directory is provided through ftp using the Get_CRTM_Binary_Files.sh script to obtain and unpack the dataset. 
If this directory doesn't exist during the `cmake` step, then cmake will download and install into `build/test_data/fix_REL-3.1.1.x/fix/`...

The fix/ directory (as of v3.1.1) contains most of the netCDF SpcCoeff and TauCoeff files, as part of our ongoing effort to transition toward netCDF-only CRTM.  We expect to deprecate the binary formats in v3.2.x 

As of CRTM v3.0.0, we no longer support legacy build system using autotools. (i.e., configure/make).  Only cmake / ecbuild (a cmake wrapper, but not required) is supported.   Many standalone Makefiles, make.dependencies, etc. have been removed, but not entirely.  Cleanup occurs as we work our way through the repository updating other things.  

**Build Step 1**
From the top level of the CRTM directory, e.g., `CRTMv3/` 
<pre>
mkdir build/
cd build/
cmake -D<cmake variables here, see below> ..
make clean
make -j8
make install (optional, see -DCMAKE_INSTALL_PREFIX below, default install location is `<build>/.`)
ctest -j8
</pre>

`-j8` runs 8 processes in parallel, adjust to your system. 
Now we have compiled the linked source codes that reside in the `src/` directory, and the ctests are built as well.

The CMake variables of interest are:
`-DCMAKE_BUILD_TYPE = RELEASE / DEBUG / RELWITHDEBINFO`  (default is `RELEASE` if not specified)
`-DBUILD_SHARED_LIBS = ON / OFF`   (build shared lib (`<build>/lib/libcrtm.so`) or static lib (`<build>/lib/libcrtm.a`) --  default is `ON` if not specified)
`-DCMAKE_INSTALL_PREFIX=<path-to-install>` (You have to run `make install` to install the libcrtm* into your desired directory `<build>/path-to-install`).


example:
```
cmake -DCMAKE_BUILD_TYPE=DEBUG -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=./install ..
```
this would make a debug build of CRTM, static library (`libcrtm.a`) and set the optional install location to `<build>/install/.` (or something similar, search for `libcrtm.*` and `*.mod`).  Custom Install only happens if you issue the `make install` command. 

The first time you run `cmake`, it will check for a `fix/` directory one level above, and if it does't find it, it will download the binary files (according to `test/CMakeLists.txt` file information), and store them in `<build>/test_data/**`.  

Linking to the library
----------------------

Let's assume the above install was moved into "/home/username/CRTMv3/", to use the library in this structure in your own application, the usual environment variables would need to be be modified something like:

<pre>
libroot="/home/username/CRTMv3/"
FCFLAGS="-I${libroot}/build/module/crtm/GNU/13.1.0 ${FCFLAGS}"  (as appropriate for your build environment)
LDFLAGS="-L${libroot}/src ${LDFLAGS}"
LIBS="-lcrtm ${LIBS}"
</pre>


**Feedback and Contact Information**

CRTM SUPPORT: visit https://forums.jcsda.org/ or visit https://github.com/JCSDA/CRTMv3 and post an issue (be sure to assign someone from the team).


Known Issues
------------

(1) Any "Transmitance Coefficient" generation codes included in src/ are not functional.  Contact CRTM support above for details.  
(2) No testing was done on PGI, XLF, or other less common compilers.  Feedback from users suggest that there's no major concerns though.  Please contact us with specifics.  Tested on GCC v5 and higher, and ifort v18 and higher.  Some specific compiler versions have issues, contact support if you run into problems.

  






