# Specification - Aerosol Scattering Jacobians

## Overview
This track implements missing Jacobian calculations for the Aerosol modules in CRTMv3, as identified in issue #281. This involves updating Tangent Linear (TL) and Adjoint (AD) models and developing K-Matrix interfaces.

## Scope
- **Target Modules:** `src/AtmScatter/AerosolScatter` and related `AtmOptics` components for aerosols.
- **Mathematical Integrity:** New TL/AD code must be exact linearizations/adjoints of the Forward model.
- **Interfaces:** Ensure compatibility with the `CRTM_K_Matrix_Module`.

## Requirements
- Consistency checks between TL/AD and Forward models.
- Verification of K-Matrix outputs against reference sensitivities.
- Adherence to CRTMv3 Fortran coding standards.
