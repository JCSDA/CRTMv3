# Oriented Particle Scattering Capability

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Add support for scattering from preferentially oriented non-spherical particles (ice crystals, canted raindrops). The current implementation assumes randomly oriented particles with plane symmetry, which is insufficient for realistic MW polarimetric simulation of precipitation.

## Current Implementation

### Phase Matrix Structure
The RT solver uses 6 phase elements assuming random orientation symmetry:
```fortran
! src/CRTM_Parameters.f90 (line 265)
INTEGER, PUBLIC, PARAMETER :: MAX_N_PHASE_ELEMENTS = 6
```

### Symmetry-Based Parameterization
From `src/RTSolution/Common_RTSolution.f90`:
- Element 1 (α₁): P₁₁ - intensity scattering
- Element 2 (α₂): P₂₂+P₃₃ combination
- Element 3 (α₃): P₂₂-P₃₃ combination
- Element 4 (α₄): P₄₄ - circular polarization
- Element 5 (β₁): P₁₂, P₂₁ coupling
- Element 6 (β₂): P₃₄, P₄₃ coupling

**This is correct for randomly oriented particles but insufficient for oriented particles.**

## Physical Motivation

### Ice Crystals
- Horizontally aligned plates/dendrites in stratiform clouds
- Electric field alignment in convective clouds
- Orientation affects polarization signature significantly

### Raindrops
- Oblate spheroids with axis ratio dependent on size
- Canting angle distribution from turbulence/wind shear
- Differential extinction between H and V polarizations

### Polarimetric Radar Observations
- ZDR (differential reflectivity) requires oriented particles
- KDP (specific differential phase) requires oriented particles
- These are standard radar products that CRTM should support

## Files to Modify

### Scattering Modules
| File | Changes Needed |
|------|----------------|
| `src/AtmScatter/CRTM_CloudScatter.f90` | Add orientation distribution handling |
| `src/AtmScatter/CRTM_AerosolScatter.f90` | Add orientation for non-spherical aerosols |
| `src/AtmScatter/CRTM_AtmScatter_Define.f90` | Add orientation parameters to structure |

### RT Solver
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/Common_RTSolution.f90` | Handle asymmetric phase matrices |
| `src/RTSolution/RTV_Define.f90` | Add orientation state variables |

### Coefficient Files
| File | Changes Needed |
|------|----------------|
| `src/Coefficients/CloudCoeff/CloudCoeff_Define.f90` | Add oriented particle LUTs |
| `src/Coefficients/CloudCoeff/CloudCoeff_netCDF_IO.f90` | I/O for new coefficients |

### Atmosphere Definition
| File | Changes Needed |
|------|----------------|
| `src/Atmosphere/CRTM_Cloud_Define.f90` | Add canting angle, orientation parameters |

## Technical Requirements

### 1. Orientation Distribution Function
Add parameters to cloud/aerosol definitions:
```fortran
TYPE :: CRTM_Cloud_type
  ! ... existing fields ...
  REAL(fp) :: Mean_Canting_Angle      ! Mean canting angle (degrees)
  REAL(fp) :: Canting_Angle_Std       ! Standard deviation of canting
  INTEGER  :: Orientation_Type        ! 0=random, 1=horizontal, 2=vertical, 3=custom
END TYPE
```

### 2. T-Matrix or DDA Coefficients
For oriented particles, need pre-computed scattering properties:
- Extinction matrix (not just scalar)
- Full phase matrix as function of orientation
- Stored in CloudCoeff with orientation dimension

### 3. Modified Phase Matrix Calculation
When orientation is specified:
```fortran
IF (Cloud%Orientation_Type /= RANDOM_ORIENTATION) THEN
  ! Use full scattering matrix from oriented particle LUT
  ! Apply orientation distribution averaging
ELSE
  ! Use existing 6-element parameterization
END IF
```

### 4. Extinction Matrix
For oriented particles, extinction is polarization-dependent:
```
     | κ_I   κ_Q   0     0   |
K =  | κ_Q   κ_I   0     0   |
     | 0     0     κ_U   κ_V |
     | 0     0    -κ_V   κ_U |
```

## Implementation Approach

### Phase 1: Infrastructure
1. Add orientation parameters to Cloud_Define
2. Extend CloudCoeff structure for oriented particles
3. Create placeholder coefficient files

### Phase 2: Rain (Oblate Spheroids)
1. Implement Pruppacher-Pitter raindrop shape model
2. Add T-matrix computed coefficients for rain
3. Implement canting angle distribution averaging
4. Validate against polarimetric radar forward models

### Phase 3: Ice Crystals
1. Add ice crystal habit orientation models
2. Generate DDA/T-matrix coefficients for oriented ice
3. Implement horizontal alignment model
4. Validate against aircraft observations

### Phase 4: TL/AD
1. Linearize orientation-dependent calculations
2. Update adjoint code
3. Validate with finite difference tests

## Coefficient Data Requirements

New LUT dimensions needed:
- Particle size (existing)
- Frequency (existing)
- Temperature (existing)
- **Orientation angle** (new)
- **Canting angle distribution width** (new)

Estimated coefficient file size increase: 5-10x for oriented particles

## Testing

- [ ] Unit test: Random orientation reproduces current results
- [ ] Unit test: Horizontal alignment gives expected ZDR
- [ ] Regression test: No change when orientation=random
- [ ] Validation: Compare rain ZDR to disdrometer + T-matrix
- [ ] Validation: Compare ice depolarization to radar observations

## References

- Mishchenko, M. I. (2000). Calculation of the amplitude matrix for a nonspherical particle in a fixed orientation. Applied Optics.
- Ryzhkov, A. V. (2001). Interpretation of polarimetric radar covariance matrix for meteorological scatterers. JAMC.
- Hogan, R. J. (2012). Fast approximate calculation of multiply scattered lidar returns. Applied Optics.

## Acceptance Criteria

1. Cloud structure includes orientation parameters
2. CloudCoeff supports oriented particle data
3. RT solver handles non-symmetric phase matrices
4. Rain and ice orientation models implemented
5. Results validated against polarimetric observations
6. No regression in random-orientation cases

## Dependencies

- Requires coefficient generation tools (T-matrix, DDA)
- May require external collaboration for LUT generation

## Labels

`enhancement`, `MW`, `polarization`, `scattering`, `cloud-physics`
