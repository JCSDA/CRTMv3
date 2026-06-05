# MW Full Polarization Validation Suite

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Create a comprehensive testing and validation suite for MW full polarization (4-Stokes) simulation capability. This includes unit tests, regression tests, and validation against observations and reference models.

## Current Testing State

### Existing MW Tests
Located in `test/mains/`:
- Forward model tests (scalar/2-component)
- K-matrix tests
- TL/AD consistency tests
- **No 4-Stokes MW tests**

### Test Infrastructure
- Uses STOP 0/1 convention
- CTest integration
- Coefficient data auto-downloads

## Test Categories Needed

### 1. Unit Tests (Component Level)

| Test | Purpose | Location |
|------|---------|----------|
| `test_MW_Surface_4x4` | Verify 4×4 emissivity/reflection | `test/mains/unit/` |
| `test_Phase_Matrix_4Stokes` | Verify 6-element → 4×4 mapping | `test/mains/unit/` |
| `test_Zeeman_Polarization` | Verify Zeeman V/H coupling | `test/mains/unit/` |
| `test_Faraday_Rotation` | Verify Q/U rotation | `test/mains/unit/` |
| `test_Stokes_Vector_RT` | Verify Stokes propagation | `test/mains/unit/` |

### 2. Regression Tests (Module Level)

| Test | Purpose | Location |
|------|---------|----------|
| `test_MW_Forward_4Stokes` | Full forward model 4-Stokes | `test/mains/regression/` |
| `test_MW_KMatrix_4Stokes` | Jacobians for all Stokes | `test/mains/regression/` |
| `test_MW_TLAD_4Stokes` | TL/AD consistency | `test/mains/regression/` |
| `test_MW_Backward_Compat` | n_Stokes=1,2 unchanged | `test/mains/regression/` |

### 3. Validation Tests (Science Level)

| Test | Purpose | Reference |
|------|---------|-----------|
| `test_WindSat_Validation` | Compare to WindSat obs | WindSat SDR |
| `test_SSMIS_Zeeman` | Zeeman at 50-60 GHz | SSMIS obs |
| `test_SMOS_Faraday` | L-band Faraday rotation | SMOS L1 |
| `test_Radar_Consistency` | ZDR/KDP forward model | Radar FO |

## Test Case Definitions

### Unit Test: MW Surface 4×4
```fortran
PROGRAM test_MW_Surface_4x4
  ! Test that MW water surface computes full 4x4 matrix

  ! Setup
  Surface%Water_Coverage = 1.0
  Surface%Wind_Speed = 10.0
  Options%n_Stokes = 4

  ! Call surface optics
  CALL CRTM_Compute_SfcOptics(...)

  ! Verify all 4 emissivity components
  IF (SfcOptics%Emissivity(1,1) <= 0) STOP 1  ! I
  IF (ABS(SfcOptics%Emissivity(1,2)) < TINY) STOP 1  ! Q (should be non-zero)
  ! U depends on azimuth
  ! V typically small

  ! Verify reflection matrix populated
  DO i = 1, 4
    DO j = 1, 4
      IF (i == j) THEN
        IF (SfcOptics%Reflectivity(1,i,1,j) <= 0) STOP 1
      END IF
    END DO
  END DO

  STOP 0
END PROGRAM
```

### Unit Test: Faraday Rotation
```fortran
PROGRAM test_Faraday_Rotation
  ! Test Faraday rotation angle and Stokes transformation

  ! Known values
  TEC = 10.0_fp        ! TECU
  B_parallel = 0.3_fp  ! Gauss
  Frequency = 1.4_fp   ! GHz

  ! Expected rotation (analytical)
  Omega_expected = 1.355e-5 * TEC * B_parallel / Frequency**2

  ! Compute
  CALL Compute_Faraday_Rotation(Ionosphere, GeometryInfo, Frequency, Omega)

  ! Verify
  IF (ABS(Omega - Omega_expected) > 1.0e-6) STOP 1

  ! Test Stokes rotation
  Stokes = [1.0, 0.5, 0.0, 0.0]  ! Partially polarized
  CALL Apply_Faraday_Rotation(Omega, Stokes)

  ! Q and U should rotate, I and V unchanged
  IF (ABS(Stokes(1) - 1.0) > TINY) STOP 1
  IF (ABS(Stokes(4) - 0.0) > TINY) STOP 1
  ! Q² + U² should be conserved
  IF (ABS(Stokes(2)**2 + Stokes(3)**2 - 0.25) > TINY) STOP 1

  STOP 0
END PROGRAM
```

### Regression Test: Backward Compatibility
```fortran
PROGRAM test_MW_Backward_Compat
  ! Verify n_Stokes=1 and n_Stokes=2 produce same results as before

  ! Reference values (from current CRTM)
  TB_V_reference = 250.5_fp
  TB_H_reference = 248.2_fp

  ! Test with n_Stokes = 1
  Options%n_Stokes = 1
  CALL CRTM_Forward(...)
  IF (ABS(RTSolution%Brightness_Temperature - TB_V_reference) > 0.01) STOP 1

  ! Test with n_Stokes = 2
  Options%n_Stokes = 2
  CALL CRTM_Forward(...)
  IF (ABS(RTSolution%Brightness_Temperature - TB_V_reference) > 0.01) STOP 1

  STOP 0
END PROGRAM
```

### Validation Test: WindSat
```fortran
PROGRAM test_WindSat_Validation
  ! Compare CRTM 4-Stokes to WindSat observations

  ! Load test cases with collocated WindSat + NWP profiles
  CALL Load_WindSat_Testcases(testcases)

  DO i = 1, n_testcases
    ! Run CRTM with n_Stokes = 4
    CALL CRTM_Forward(Atmosphere(i), Surface(i), ..., RTSolution)

    ! Compare all 4 Stokes
    bias_I = RTSolution%Stokes(1) - testcases(i)%obs_I
    bias_Q = RTSolution%Stokes(2) - testcases(i)%obs_Q
    bias_U = RTSolution%Stokes(3) - testcases(i)%obs_U
    bias_V = RTSolution%Stokes(4) - testcases(i)%obs_V

    ! Accumulate statistics
    ...
  END DO

  ! Check bias and RMSE within acceptable limits
  IF (ABS(mean_bias_I) > 2.0) STOP 1  ! 2K threshold
  IF (rmse_Q > 1.0) STOP 1  ! 1K threshold for Q

  STOP 0
END PROGRAM
```

## Test Data Requirements

### Reference Profiles
Create standardized test atmospheres:
- Clear ocean (various wind speeds/SST)
- Precipitating (stratiform/convective)
- Cold/warm atmospheres
- Various viewing geometries

### Observation Datasets
| Sensor | Data Type | Use |
|--------|-----------|-----|
| WindSat | SDR L1C | Full 4-Stokes validation |
| SSMIS | TDR | Zeeman validation |
| SMOS | L1C | Faraday validation |
| GPM/GMI | L1B | Rain depolarization |

### Reference Model Outputs
- Standalone Zeeman RT models
- T-matrix/Mie calculations for hydrometeors
- Analytical solutions for simple cases

## Test Infrastructure

### New Test Directory Structure
```
test/mains/
├── unit/
│   ├── test_MW_Surface_4x4.f90
│   ├── test_Phase_Matrix_4Stokes.f90
│   ├── test_Zeeman_Polarization.f90
│   ├── test_Faraday_Rotation.f90
│   └── test_Stokes_Vector_RT.f90
├── regression/
│   ├── test_MW_Forward_4Stokes.f90
│   ├── test_MW_KMatrix_4Stokes.f90
│   ├── test_MW_TLAD_4Stokes.f90
│   └── test_MW_Backward_Compat.f90
└── validation/
    ├── test_WindSat_Validation.f90
    ├── test_SSMIS_Zeeman.f90
    └── test_SMOS_Faraday.f90
```

### CMake Integration
```cmake
# Add new test category
add_test(NAME test_MW_Surface_4x4
         COMMAND test_MW_Surface_4x4)
set_tests_properties(test_MW_Surface_4x4
         PROPERTIES LABELS "unit;MW;polarization")

# Group tests
add_custom_target(test_MW_fullpol
         COMMAND ${CMAKE_CTEST_COMMAND} -L "polarization")
```

## Implementation Approach

### Phase 1: Unit Test Framework
1. Create test programs for each component
2. Define pass/fail criteria
3. Add to CMake/CTest

### Phase 2: Regression Tests
1. Generate reference results with current CRTM
2. Create comparison tests
3. Establish tolerance thresholds

### Phase 3: Validation Infrastructure
1. Acquire observation datasets
2. Create collocation tools
3. Implement statistical comparison

### Phase 4: Continuous Integration
1. Add tests to CI pipeline
2. Set up nightly validation runs
3. Track metrics over time

## Acceptance Criteria

### Unit Tests
- [ ] All unit tests pass
- [ ] Coverage of all new code paths
- [ ] Edge cases tested

### Regression Tests
- [ ] Backward compatibility verified
- [ ] TL/AD consistency < 10⁻⁶
- [ ] Solver consistency (ADA/SOI/AMOM)

### Validation Tests
- [ ] WindSat bias < 2K (I), < 1K (Q)
- [ ] Zeeman V/H ratio within 5%
- [ ] Faraday rotation within 1°

## Dependencies

- Requires completion of Issues #01-#07 for full testing
- Can create test stubs early, fill in as features implemented
- Validation requires external data access

## Estimated Effort

| Phase | Effort |
|-------|--------|
| Unit tests | Moderate |
| Regression tests | Moderate |
| Validation data | High (data acquisition) |
| Validation tests | Moderate |
| CI integration | Low |

## Labels

`testing`, `MW`, `polarization`, `validation`, `infrastructure`
