# RT Solver 4-Stokes Support for MW

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Update RT solvers to fully support 4-Stokes (I, Q, U, V) simulation for microwave channels. While infrastructure exists, several solvers have hardcoded limitations and the MW pathway defaults to scalar or 2-component calculations.

## Current Implementation

### Solver Capabilities (Nominal)

| Solver | Location | n_Stokes Support |
|--------|----------|------------------|
| ADA | `src/RTSolution/ADA/` | Up to 4 (VIS/UV) |
| AMOM | `src/RTSolution/AMOM/` | Up to 4 (VIS/UV) |
| SOI | `src/RTSolution/SOI/` | **Hardcoded to 2** |
| Emission | `src/RTSolution/Emission/` | Scalar |

### SOI Limitation
```fortran
! src/RTSolution/UWisc_SOI/SOI_CRTM_Forward_Module.f90 (line 22)
INTEGER, PARAMETER :: SFCOPTICS_N_STOKES = 2
```
This hardcodes SOI to 2-Stokes regardless of user setting.

### Default Behavior
```fortran
! src/Options/CRTM_Options_Define.f90 (line 149)
INTEGER :: n_Stokes = 1  ! Default is scalar
```

### Surface Boundary Condition
MW surface models only provide 2-component (V/H) data:
- Solvers expect diagonal 2×2 reflection, not full 4×4
- Need to handle expanded surface matrix (Issue #01)

## Files to Modify

### SOI Solver
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/UWisc_SOI/SOI_CRTM_Forward_Module.f90` | Remove 2-Stokes hardcode |
| `src/RTSolution/UWisc_SOI/SOI_Module.f90` | Verify 4-Stokes arrays |
| `src/RTSolution/UWisc_SOI/SOI_doubling.f90` | Check matrix dimensions |

### Main RT Interface
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/CRTM_RTSolution.f90` | Pass n_Stokes to all solvers |
| `src/RTSolution/RTV_Define.f90` | Verify array dimensions for n_Stokes=4 |
| `src/RTSolution/Common_RTSolution.f90` | Ensure MW path supports 4-Stokes |

### Emission Solver
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/Emission/Emission_Module.f90` | Add multi-Stokes capability |

### Options
| File | Changes Needed |
|------|----------------|
| `src/Options/CRTM_Options_Define.f90` | Document MW 4-Stokes usage |

## Technical Requirements

### 1. Remove SOI Hardcoding
```fortran
! BEFORE (hardcoded)
INTEGER, PARAMETER :: SFCOPTICS_N_STOKES = 2

! AFTER (from RTV)
n_Stokes = RTV%n_Stokes  ! Use user-specified value
```

### 2. Array Dimension Verification
Ensure all working arrays sized for 4-Stokes:
```fortran
! Phase matrices: (n_Angles*n_Stokes, n_Angles*n_Stokes, n_Layers)
REAL(fp) :: Pff(MAX_N_ANGLES*MAX_N_STOKES, MAX_N_ANGLES*MAX_N_STOKES+1, n_Layers)
REAL(fp) :: Pbb(MAX_N_ANGLES*MAX_N_STOKES, MAX_N_ANGLES*MAX_N_STOKES+1, n_Layers)

! Source terms: (n_Angles*n_Stokes)
REAL(fp) :: Source(MAX_N_ANGLES*MAX_N_STOKES)
```

### 3. Surface Boundary Interface
Update surface coupling for 4×4 matrix:
```fortran
! Current: 2×2 diagonal
DO i = 1, n_Angles
  Reflectivity_used(i,1,i,1) = SfcOptics%Reflectivity(i,1,i,1)
  Reflectivity_used(i,2,i,2) = SfcOptics%Reflectivity(i,2,i,2)
END DO

! Required: Full 4×4
DO i = 1, n_Angles
  DO is = 1, n_Stokes
    DO js = 1, n_Stokes
      Reflectivity_used(i,is,i,js) = SfcOptics%Reflectivity(i,is,i,js)
    END DO
  END DO
END DO
```

### 4. Emission Solver Extension
For non-scattering atmospheres, extend to 4-Stokes:
```fortran
! Stokes emission vector
Emission(1) = (1 - Emissivity_V - Emissivity_H) / 2  ! I
Emission(2) = (Emissivity_V - Emissivity_H) / 2      ! Q
Emission(3) = ...  ! U (from azimuth dependence)
Emission(4) = ...  ! V (if any)
```

### 5. MW-Specific Considerations
```fortran
! MW typically uses mth_Azi = 0 (azimuth-averaged)
! For full-pol, may need mth_Azi > 0 for U-Stokes
IF (RTV%n_Stokes > 2 .AND. SpcCoeff%Is_Microwave) THEN
  ! Consider azimuth treatment for U component
END IF
```

## Implementation Approach

### Phase 1: SOI Cleanup
1. Remove SFCOPTICS_N_STOKES hardcoding
2. Pass n_Stokes from RTV to SOI internals
3. Verify SOI works with n_Stokes = 1, 2, 4
4. Test against ADA for consistency

### Phase 2: Array Verification
1. Audit all RT arrays for dimension sufficiency
2. Identify any implicit 2-Stokes assumptions
3. Fix dimension issues found

### Phase 3: Surface Coupling
1. Update surface BC application for 4×4 matrix
2. Handle case when surface provides 2×2 (backward compat)
3. Test with mock 4×4 surface data

### Phase 4: Emission Solver
1. Extend emission module for Stokes vector
2. Add polarized emission capability
3. Test clear-sky MW polarization

### Phase 5: Integration Testing
1. End-to-end 4-Stokes MW forward model
2. Compare solvers (ADA vs SOI vs AMOM)
3. Verify against reference calculations

## Testing

### Unit Tests
- [ ] SOI with n_Stokes=1 (regression)
- [ ] SOI with n_Stokes=2 (regression)
- [ ] SOI with n_Stokes=4 (new capability)
- [ ] ADA/SOI/AMOM consistency at n_Stokes=4

### Integration Tests
- [ ] MW channel with n_Stokes=4, scattering atmosphere
- [ ] MW channel with n_Stokes=4, clear atmosphere
- [ ] Surface-only case (no atmosphere)

### Regression Tests
- [ ] All existing MW tests pass with n_Stokes=1
- [ ] V/H results unchanged when using n_Stokes=2

## Backward Compatibility

**Critical:** Must not break existing functionality
- n_Stokes=1 (default) produces identical results
- n_Stokes=2 for MW produces identical results
- Only n_Stokes=4 activates new code paths

## Relationship to Other Issues

This issue **enables** the other MW full-pol issues:
- Issue #01 (Surface): Provides 4×4 input to RT solver
- Issue #02 (Oriented): Provides phase matrix for 4-Stokes
- Issue #03 (Zeeman): Uses matrix transmittance in RT
- Issue #04 (Faraday): Applied after RT solution

**Should be worked in parallel with Issue #01** (Surface) for testing.

## Acceptance Criteria

1. SOI hardcoded 2-Stokes limitation removed
2. All solvers accept n_Stokes = 1, 2, 4
3. Surface boundary handles 4×4 reflection matrix
4. Emission solver supports Stokes vector
5. Regression tests pass (no change at n_Stokes=1,2)
6. ADA/SOI/AMOM produce consistent 4-Stokes results

## Estimated Complexity

- SOI cleanup: Low-Moderate
- Array verification: Low
- Surface coupling: Moderate
- Emission extension: Moderate
- Integration testing: Moderate

## Labels

`enhancement`, `MW`, `polarization`, `rt-solver`, `architecture`
