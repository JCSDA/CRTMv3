# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CRTM (Community Radiative Transfer Model) v3.2.0 is a Fortran radiative transfer model for satellite data assimilation. It computes top-of-atmosphere radiances and brightness temperatures from atmospheric/surface states, and provides Jacobians (sensitivities) for data assimilation systems.

## Build Commands

```bash
mkdir build && cd build
cmake ..
make -j8           # Adjust -j to your CPU count
ctest -j8          # Run all tests in parallel
```

**CMake options:**
- `-DCMAKE_BUILD_TYPE=DEBUG|RELEASE|RELWITHDEBINFO` (default: RELEASE)
- `-DBUILD_SHARED_LIBS=ON|OFF` (default: ON, creates libcrtm.so)
- `-DFIX_FILE_PATH=<path>` - Path to coefficient files (auto-downloads if not found)
- `-DBUILD_TESTING=ON|OFF` (default: ON)
- `-DOPENMP=ON|OFF` (default: ON)

**Build outputs:**
- Library: `build/lib/libcrtm.so` or `libcrtm.a`
- Modules: `build/module/crtm/<compiler>/<version>/`
- Test data: `build/test_data/` (downloaded automatically)

## Running Tests

```bash
cd build
ctest -j4                      # Run all tests in parallel
ctest -VV -R <testname>        # Run specific test with verbose output
ctest --output-on-failure      # Show output only for failed tests
```

Test exit codes: `STOP 0` = success, `STOP 1` = failure.

## Code Architecture

### Core Modules (src/)

The four main CRTM entry points in order of complexity:
1. **CRTM_Forward_Module** - Forward radiative transfer (atmosphere/surface → radiances)
2. **CRTM_Tangent_Linear_Module** - Linearized forward model
3. **CRTM_K_Matrix_Module** - Jacobian/sensitivity matrix calculations
4. **CRTM_Adjoint_Module** - Adjoint (reverse mode) computations

**CRTM_LifeCycle** handles initialization and finalization.

### Component Libraries

| Directory | Purpose |
|-----------|---------|
| Atmosphere/ | Atmospheric profiles (pressure, temperature, gases, clouds, aerosols) |
| AtmAbsorption/ | Gas absorption (ODAS, ODPS, ODZeeman algorithms) |
| AtmOptics/ | Atmospheric optical properties |
| AtmScatter/ | Cloud, aerosol, and molecular scattering |
| SfcOptics/ | Surface emissivity/reflectivity (MW_Land, MW_Water, IR_*, VIS_*) |
| Surface/ | Surface type characterization |
| Coefficients/ | Coefficient I/O (SpcCoeff, TauCoeff, CloudCoeff, AerosolCoeff, EmisCoeff) |
| RTSolution/ | Radiative transfer solution structures |

### File Naming Conventions

- `*_Define.f90` - Derived type definitions
- `*_TL.f90` - Tangent-linear routines
- `*_AD.f90` - Adjoint routines
- `*_Binary_IO.f90` / `*_netCDF_IO.f90` - Coefficient I/O

### Test Structure (test/mains/)

- **application/** - Large end-to-end tests (check_crtm.F90)
- **regression/** - Functionality tests organized by mode:
  - forward/, k_matrix/, adjoint/, tangent_linear/
- **unit/** - Small-scope component tests

## Requirements

- Fortran 2008 compiler (GCC 5+, Intel 18+, Cray, NVHPC)
- netCDF4/HDF5 (built with same Fortran compiler)
- CMake 3.20+
- git-lfs 2.10+

## Key Technical Notes

- Coefficient files support both legacy binary (.bin) and netCDF (.nc4) formats
- OpenMP parallelization is enabled by default
- Most modules have parallel TL/AD implementations for sensitivity calculations
- Transmittance coefficient generation codes (TauProd/, TauRegress/) are not functional
