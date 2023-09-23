CRTM REL-3.0.0 
====================

[![Build Status](https://app.travis-ci.com/JCSDA-internal/crtm.svg?token=r6aaq9P13fHcTi8yBgdM&branch=develop)](https://app.travis-ci.com/JCSDA-internal/crtm)

Preamble
--------

CRTM v3.0.1 release (`REL-3.0.1`)  

v3.0.1 Released September 23, 2023
v3.0.0 Released March, 2023  
v2.4.1-alpha Released on April 1, 2021 (internal realease only)
v2.4.0 Released on October 23, 2020

This is an experimental/early release of CRTM v3.0.0, some features may not be fully functional. Contact crtm-support@googlegroups.com.  
v3.x features will be rolled out in incremental updates. 

Basic requirements:  
(1) A Fortran 2003 compatible compiler
(2) A netCDF4 / HDF5 library that has been built with the compiler you're going to use (module environments are helpful here)
(3) A linux, macOS, or unix-style environment.  This has not been tested under any Windows Fortran environments.
(4) Bash shell is preferred. 
(5) git and git-lfs (minimum version TBD, but has been tested on git-lfs v2.10 and higher )

=========================================================

**Important Note**: If reading this, you're cloning the CRTM development repository.  The development repository is structured in a way that makes it less user friendly, but more amenable to development and testing.

**JEDI NOTE** This develop (or master) branch is also designed to work directly in a JEDI container or JEDI environment. If you're doing JEDI things, you're probably in the right spot. However, you should stop reading right now and have a look at the README_JEDI.md file.   

If you're looking for an older version of CRTM (v2.3.0 or older) you should obtain the appropriate tarball from
https://ftp.emc.ncep.noaa.gov/jcsda/CRTM/

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
JCSDA CRTM v3.0.x Build Instructions

- Development Repository Build
- Note: the development repository build differs from a release build. 
  
The CRTM **development** repository directory structure looks like:

<pre>
 .
  ├── LICENSE  (CC0 license)
  ├── COPYING  (CC0 legal document)
  ├── NOTES
  ├── README.md 
  ├── Set_CRTM_Environment.sh
  ├── Get_CRTM_Binary_Files.sh  
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
  │   ├── <b>Build/</b>
  │   │   └── <b>libsrc/</b>
  │   │       └── <b>test/</b>
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
At present, the 'fix/' directory is provided through ftp, rather than git-LFS.  We're working on a way to make the process of getting binary data equitable and easy for all users.

The files therein are gzipped (*.gz extension).  Use the Get_CRTM_Binary_Files.sh script to obtain and unpack the dataset.    

At the top level (`crtm/`), the `configuration` directory contains the various compiler-specific configuration files.
<pre>
  ls  configuration/
  ftn.setup                 ftn.setup.csh
  g95-debug.setup           gfortran.setup.csh     pgf95.setup
  g95-debug.setup.csh       ifort-debug.setup      pgf95.setup.csh
  g95.setup                 ifort-debug.setup.csh  xlf2003-debug.setup
  g95.setup.csh             ifort.setup            xlf2003-debug.setup.csh
  gfortran-debug.setup      ifort.setup.csh        xlf2003.setup
  gfortran-debug.setup.csh  pgf95-debug.setup      xlf2003.setup.csh
  gfortran.setup            pgf95-debug.setup.csh

</pre>
[Note: as of the time of writing, October 2020, only `ifort.setup`, `ifort-debug.setup`, `gfortran.setup`, `gfortran-debug.setup` have been actively developed and tested.  It is strongly recommended that the user use one of these compilers until the remaining setup files are updated, or update them yourself.  Contact the support email address for specific compiler support requests.  The c-shell (.csh) extension files have not been updated.]

All of the above files define values for the environment variables `FC`, `FCFLAGS`, `LDFLAGS`, and `LIBS`.

It is also inside of these files that you should specify the path to your netCDF4 and HDF5 installations.  
Here's a trick to find these on your system:
find / -name "libnetcdff.a" -type f
find / -name "netcdf.mod" -type f
find / -name "libhdf5.a" -type f  

For example, if you're using the intel fortran compiler version 19.0.5:
```
find / -name "libhdf5.a" -type f  
<...>
/opt/hdf5/1.8.21-intel-19.0.5/lib/libhdf5.a
<...>
```
and you would modify the base path in the appropriate `ifort.setup` file, for example:
`export HDF_DIR="/opt/hdf5/1.8.21-intel-19.0.5"`
similarly for the `NC4_DIR` line.   

To use these files to define the CRTM build environment, you should source them. For example, if you use the sh/bash/ksh shells and you want to setup
for a build using the gfortran compiler using debug options you would type:

**Configuration Step 1**

    . configuration/gfortran-debug.setup

(note the `. ` -- for a detailed discussion of `.` vs. `source` see: https://unix.stackexchange.com/questions/58514/what-is-the-difference-between-and-source-in-shells)

**Configuration Step 2**

    . ./Set_CRTM_Environment.sh

Again noting the leading `. `.  This sets the required environment variables to identify various paths needed to build.  `CRTM_ROOT`, `CRTM_SOURCE_ROOT`, etc.

**Configuration Step 3**
<pre>
sh Uncompress_Binary_Files.sh
cd src/
cd Build/
make clean  
cd ..  
make realclean  
make  
</pre>
The commands `make clean` and `make realclean` ensures that the underlying links, compiled files, generated Makefiles. are removed to avoid conflicts with a clean build.
The command `make` at the `src/` level performs the linking process between the upper level `src/**` directories and the `src/Build/libsrc` directory.  

Note: After runnin `make`, you may see certain "nc4" files listed as missing, these are files that will be converted to netCDF4 format, but have not yet been added.  

Assuming no fatal error messages, continue to the Build steps below.

**Build Step 1**
<pre>
cd Build/
./configure --prefix=${PWD}
make clean
make -j4
</pre>

(See additional options for `configure` below.  `-j4` sets the number of parallel `make` processes to 4.)

Now we have finally compiled the linked source codes that reside in the `libsrc/` directory.    Please note that once the source codes are linked in the libsrc directory, all development and testing can occur at the `Build/` level.  In the `libsrc/` directory, the source codes link back to the version-controlled counterparts, so you'll want to answer "yes" to any queries about opening the version controlled codes when trying to edit them (this occurs in `emacs`, for example).

**Build Step 2**
<pre>
cd libsrc/
ls -l *.mod
ls -l *.a
make check
make install
</pre>

The `ls` commands are to verify that indeed the .mod files have been created and the library file (which external codes link against) has also been created.
The `make check` command builds and runs the default CRTM test `check_crtm`, located in the `src/Build/libsrc/test` directory.

`make install` installs the library in the directory defined by the command `configure --prefix=directory_name`.


Summary of Build Script
-----------------------
<pre>
. configuration/gfortran-debug.setup (or whatever compiler you're using)  
. ./Set_CRTM_Environment.sh  
sh Uncompress_Binary_Files.sh
cd src/  
make realclean  
make  
cd Build/  
./configure --prefix=${PWD}  
make clean  
make -j4  
cd libsrc/  
ls -l *.mod  
ls -l *.a  
rm -rf make_check.out  
make check > make_check.out  
head -n10 make_check.out  
tail -n10 make_check.out  
make install  
</pre>


Linking to the library
----------------------

Let's assume the above install was moved into "/home/username/CRTM/crtm_v3.0.0/", to use the library in this structure in your own application, the usual environment variables would need to be be modified something like:

<pre>
libroot="/home/username/CRTM/crtm_v3.0.0"
FCFLAGS="-I${libroot}/include ${FCFLAGS}"
LDFLAGS="-L${libroot}/lib ${LDFLAGS}"
LIBS="-lcrtm ${LIBS}"
</pre>

Uninstalling the library
------------------------

To uninstall the library (assuming you haven't moved the installation directory contents somewhere else) you can type:

    make uninstall

This will DELETE the created installation directory. So, for a library version, say, v3.0.0, if your configure script invocation was something like

    ./configure --prefix=${PWD} ...other command line arguments...

then the "uninstall" target will delete the "${PWD}/crtm_v3.0.0" directory.


Cleaning Up
-----------

Two targets are provided for cleaning up after the build. To remove all the build products type

<pre>
cd src/Build
make clean
</pre>

To also remove all the configuration products (i.e. the makefiles) type

<pre>
cd src/Build
make distclean
</pre>


(optional) "Build Release" Setup and Configuration:
--------------------------------------------------

Within the 'src/Build' directory, The legacy build system for the CRTM uses an autoconf-generated `configure` script, which depends on the existence of a few key files.  
(1) the `configure.ac` file, which contains instructions for how the `configure` file will be built when the `autoconf` command is executed, this is finicky, contact support for help.  
(2) The `Makefile.in` file, which contains instructions on how executing the `configure` script will generate `Makefile` in libsrc and test subdirectories.  

The Build `Makefile`s assume that environment variables (envars) will be defined that describe the compilation environment. The envars
that *must* be defined are:  
  FC:      the Fortran95/2003 compiler executable,  
  FCFLAGS: the flags/switches provided to the Fortran compiler.
These can be set (in the Build directory) by `. ./config-setup/<compiler>.setup` as described previously).

In src/Build:
<pre>
. ./config-setup/gfortran-debug.setup
./configure --prefix=${PWD}
make clean
make -j4
make check
make install
</pre>


**Additional options for `configure`**

`configure` sets an install path environment variable, among other things.  This, by default, will set the `lib/` and `include/` directory paths in the `/usr/local/crtm_v3.0.0/` (or whatever string in in `src/CRTM_Version.inc`).  

The `--prefix` switch sets the installation directory, make sure you have write access to that directory.  

You can override this by setting a different install directory as follows:  
  `   ./configure --prefix=<install directory>`  
For example, `./configure --prefix=${PWD}` will create the library in the directory in which you're currently in (e.g., crtm/src/Build/crtm_v3.0.0y/).

By default, the CRTM is built for big-endian I/O. The --disable-big-endian switch builds the library and test programs for little-endian I/O:

  `   ./configure --disable-big-endian --prefix=<install directory>`

If you need more flexibility in the library build you can specify the necessary information directly to the configure script that generates the Makefiles. For
example, for the Intel ifort compiler:
<pre>
  ./configure --prefix=${PWD} \
                --disable-big-endian \
                FC="ifort" \
                FCFLAGS="-O3 -qopenmp -g -traceback"  
</pre>
This overrides the FC and FCFLAGS variables that were set by "sourcing" the `configuration/` file earlier, it is strongly recommended that you use the provided configuration files since they contain flags that have been added after substantial debugging and testing.



**Feedback and Contact Information**

CRTM SUPPORT EMAIL: crtm-support@googlegroups.com OR visit https://forums.jcsda.org/

If you have problems building the library please include the generated "config.log" file in your email correspondence.

Known Issues
------------

(1) Any "Transmitance Coefficient" generation codes included in src/ are not functional.  Contact CRTM support above for details.  
(2) No testing was done on PGI, XLF, or other less popular compilers.  
(3) Compiler setup files do not contain "generic" ways to point to netCDF libraries - you need to edit those files and ensure that the paths point to the correct place.  This is the netCDF life.  Note: Building inside of a JEDI environment (e.g., singularity container) using ecbuild makes this part much easier. 

Troubleshooting
---------------
<pre>
Installing <path>/crtm based scripts...
  crtm_install_scripts.sh(INFORMATION): CRTM root directory is crtm
  crtm_install_scripts.sh(INFORMATION): /bin exists...
  crtm_install_scripts.sh(INFORMATION): Your $PATH does NOT contain /bin...
  crtm_install_scripts.sh(INFORMATION): Creating a crtmrc file with $PATH modification. For a permanent change modify your .bash_profile (or similar) file.
</pre>

This uncommon error message relates to the fact that you do not have a `$HOME` environment variable set.  You'll also need a `$HOME/bin` directory. Typically something like: `export HOME="/home/users/username/"` or `export HOME="~"` may work as well.    However, usually `$HOME` is set automatically by your system. If you're having this problem, you're likely to have even more problems later -- contact your Sysadmin first.  

<pre>
checking whether the Fortran compiler works... no
configure: error: in `<path>/src/Build':
configure: error: Fortran compiler cannot create executables
</pre>

Bash users, type `export | grep "FC"`, it should be set to the name of a compiler, e.g., `declare -x FC="ifort"`.  next simply type `ifort` (or whatever FC is trying to use)  at the command prompt to see if it's accessible.   If it says `command not found`, then you're missing the path to your compiler.  This could be a `module` command that needs to be run, or a valid compiler needs to be installed.  This varies based on operating system. 


<pre>
When issuing the `make` command in src/ :

File Type_Kinds.f90 not found in CRTM_Module.F90 hierarchy.
File File_Utility.f90 not found in CRTM_Module.F90 hierarchy.
<...> dozens of similar lines <...>
File FitCoeff_WriteFile.inc not found in CRTM_Module.F90 hierarchy.
File FitCoeff_Equal.inc not found in CRTM_Module.F90 hierarchy.

Returning to directory <directory>/crtm/src
</pre>
You forgot to `. ./Set_CRTM_Environment.sh` in the crtm/ directory, paying close attention to that leading `. `.
Then `cd src/` and `make`.



