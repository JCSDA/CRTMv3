# Tech Stack - CRTM v3.2.0

## Core Technologies
- **Fortran:** Primary programming language (Modern standards like 2008+ are encouraged for new development).
- **C:** Used for build system interactions and system-level operations.

## Build & Dependencies
- **Build System:** CMake (Minimum VERSION 3.20) and Make.
- **Data Formats:**
    - **NetCDF:** Primary library for scientific data handling and coefficient files.
    - **HDF5:** Core dependency for NetCDF.
- **Parallelization:** OpenMP for multi-threaded performance.

## Operating Environment
- **Platform:** Linux, macOS (Unix-style).
- **Compilers Supported:** Intel (ifort/ifx), GNU (gfortran), NVIDIA (nvfortran).

## Integration Context
- **JEDI:** Designed for seamless integration and performance within the Joint Effort for Data assimilation Integration environment.