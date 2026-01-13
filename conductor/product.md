# Product Guide - CRTMv3 Missing Jacobians

## Initial Concept
Replacing missing jacobians in CRTM (Ref: JCSDA/CRTMv3 issue #281).

## Target Users
- Remote sensing researchers developing new sensor models who require complete Jacobian (K-matrix) functionality for their work.

## Project Goals
- Systematically replace or implement missing Jacobian calculations throughout the CRTMv3 codebase.
- Ensure all new Tangent Linear (TL), Adjoint (AD), and K-Matrix models are mathematically consistent and correctly integrated.

## Core Features
- **Missing Jacobian Implementation:** Development of missing Jacobian calculations for key modules, including Aerosol and Cloud scattering components.
- **TL/AD Development:** Implementation of the corresponding Tangent Linear and Adjoint models required for the Jacobians.
- **K-Matrix Models:** Development of K-Matrix models derived from the Adjoint implementations to provide sensitivity profiles.

## Validation & Testing
- **Unit Testing:** Comprehensive tests for individual Jacobian components to ensure mathematical correctness.
- **Consistency Verification:** Integration tests to verify that the Tangent Linear and Adjoint models are consistent with the Forward model (e.g., using finite difference checks).

## Coding Standards
- **CRTMv3 Standards:** Strict adherence to the established Fortran coding standards and architectural patterns used in the CRTMv3 project.
