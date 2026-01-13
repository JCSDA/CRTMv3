# Technology Stack - CRTMv3 Missing Jacobians

## Core Technologies
- **Language:** Fortran (2008 compatible) - The primary language for the CRTMv3 project.
- **Build System:** CMake (Minimum version 3.20) - Used for configuration, building, and testing.

## Libraries & Dependencies
- **NetCDF4 / HDF5:** Required for reading and writing satellite and coefficient data.
- **OpenMP:** Utilized for shared-memory parallelization of radiative transfer calculations.

## Development Tools
- **Testing:** CTest is used as the primary test runner for the project's verification suite.
- **Environment:** Linux/Unix-style environment is the target platform.
