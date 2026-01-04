# Track Specification: TL/AD/K-Matrix OMP Consistency

## Goal
Investigate and resolve consistency issues in Tangent Linear (TL), Adjoint (AD), and K-Matrix calculations when running with OpenMP enabled. The primary objective is to ensure that the multi-threaded execution produced identical (within numerical precision) results to single-threaded execution and passes all validation tests.

## Context
The branch `feature/btj_TL_AD_K_Matrix_Consistency_OMP` was created to improve consistency and performance under OpenMP. However, current tests (e.g., `test_k_matrix_Simple`, `test_adjoint_Simple`, `test_tangent_linear_Simple`) are failing, indicating either race conditions, incorrect data sharing, or unsynchronized updates in the parallel regions.

## Requirements
- **Functional:**
    - `test_k_matrix_Simple_*` must pass.
    - `test_adjoint_Simple_*` must pass.
    - `test_tangent_linear_Simple_*` must pass.
    - `test_tangent_linear_ClearSky_*` must pass.
    - Consistency between TL and AD must be verified (dot-product test).
    - Consistency between TL and K-Matrix must be verified.
- **Non-Functional:**
    - Adhere to JCSDA/JEDI coding standards.
    - Use `IMPLICIT NONE`.
    - Maintain portability across supported compilers (ifx, gfortran, etc.).
- **Validation:**
    - Tests must pass with `OMP_NUM_THREADS=1` and `OMP_NUM_THREADS > 1`.
