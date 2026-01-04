# Product Guide - CRTM v3.2.0

## Initial Concept
The Community Radiative Transfer Model (CRTM) is a fast radiative transfer model designed for coordinated development and releases of version 3.x. It serves as a critical component for satellite data assimilation and remote sensing applications.

## Target Users
- NWP (Numerical Weather Prediction) scientists and researchers.
- Satellite data assimilation engineers (e.g., JEDI/JCSDA community).
- Climate modelers and atmospheric physicists.

## Key Goals & Features
- **Modernization:** Complete the transition to netCDF-only coefficient files and deprecate binary formats.
- **Integration:** Enhance stability and performance within the JEDI (Joint Effort for Data assimilation Integration) environment.
- **Optimization:** Improve multi-threading efficiency and OpenMP support.
- **Advancement:** Focus on modern, netCDF-only architecture and support for new technologies like active sensors and updated physical models (aerosol/cloud properties).

## Development Philosophy
- **Stability First:** Rigorous testing and validation (ctest, regression tests) are mandatory before any feature is merged into the `develop` or `master` branches.

## Non-Functional Requirements
- **Portability:** Ensure reliable compilation across various platforms (Linux, macOS) and compilers (Intel, GNU, NVIDIA).
- **Performance:** Maintain high-speed radiative transfer calculations for real-time operational weather prediction.