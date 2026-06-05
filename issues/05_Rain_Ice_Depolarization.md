# Rain/Ice Depolarization Physics

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Implement polarization-dependent extinction and depolarization from precipitation (rain, snow, ice). Currently, CRTM treats hydrometeor extinction as scalar, missing the differential effects that are crucial for polarimetric observations.

## Current Implementation

### Scalar Extinction
All hydrometeor extinction treated identically for all polarizations:
```fortran
! Optical depth is scalar - same for V and H
AtmOptics%Optical_Depth(k) = ... ! Single value per layer
```

### Missing Physics
- No differential extinction (κ_V ≠ κ_H)
- No cross-polarization extinction
- No depolarization from multiple scattering

## Physical Background

### Differential Extinction
Non-spherical hydrometeors cause polarization-dependent extinction:
- **Rain**: Oblate drops → κ_H > κ_V (H attenuated more)
- **Ice**: Shape-dependent, often κ_V > κ_H for columns

### Depolarization Mechanisms
1. **Differential absorption**: ΔA = A_H - A_V
2. **Differential phase shift**: Creates rotation of polarization
3. **Multiple scattering**: Randomizes polarization state

### Specific Differential Attenuation (A_DP)
Key observable in polarimetric radar:
```
A_DP = A_H - A_V (dB/km)
```

### Relationship to Polarimetric Radar
| Parameter | Physical Meaning |
|-----------|------------------|
| ZDR | Differential reflectivity |
| KDP | Specific differential phase |
| A_DP | Specific differential attenuation |
| ρ_HV | Correlation coefficient |

CRTM should be consistent with these forward operators.

## Files to Modify

### Atmospheric Optics
| File | Changes Needed |
|------|----------------|
| `src/AtmOptics/CRTM_AtmOptics_Define.f90` | Add polarized extinction |
| `src/AtmOptics/CRTM_AtmOptics.f90` | Compute differential extinction |

### Cloud Scattering
| File | Changes Needed |
|------|----------------|
| `src/AtmScatter/CRTM_CloudScatter.f90` | Output polarized extinction |
| `src/Coefficients/CloudCoeff/CloudCoeff_Define.f90` | Store differential extinction coefficients |

### RT Solver
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/Common_RTSolution.f90` | Use polarized extinction |
| `src/RTSolution/ADA/ADA_Module.f90` | Matrix extinction |

## Technical Requirements

### 1. Polarized Extinction Structure
```fortran
TYPE :: Polarized_Extinction_type
  REAL(fp) :: Extinction_V        ! Vertical polarization (Np/km)
  REAL(fp) :: Extinction_H        ! Horizontal polarization
  REAL(fp) :: Extinction_Cross_Re ! Cross-polarization (real)
  REAL(fp) :: Extinction_Cross_Im ! Cross-polarization (imag)
  REAL(fp) :: Differential_Phase  ! Phase difference V-H (rad/km)
END TYPE
```

### 2. Extinction Matrix
For a layer with precipitation:
```
       | κ_I + κ_Q    κ_U - iκ_V |
K =    |                         |
       | κ_U + iκ_V   κ_I - κ_Q  |
```

Or in V/H basis:
```
       | κ_V    κ_VH |
K =    |             |
       | κ_HV   κ_H  |
```

### 3. Raindrop Model
Use Pruppacher-Pitter or similar axis ratio model:
```fortran
! Axis ratio as function of equivalent diameter
axis_ratio = 1.0 - 0.062 * D_eq  ! For D_eq in mm

! Then compute extinction from T-matrix or Mie
CALL Compute_Oblate_Extinction(D_eq, axis_ratio, frequency, &
                                extinction_V, extinction_H)
```

### 4. Integration with DSD
Integrate over drop size distribution:
```fortran
! Marshall-Palmer or gamma DSD
DO i = 1, n_size_bins
  kext_V = kext_V + N(D_i) * sigma_ext_V(D_i) * dD
  kext_H = kext_H + N(D_i) * sigma_ext_H(D_i) * dD
END DO
```

### 5. Layer Transmittance
Polarized transmittance for precipitation layer:
```fortran
! Simple diagonal approximation
T_V = EXP(-kext_V * path_length)
T_H = EXP(-kext_H * path_length)

! Full matrix if cross-pol significant
T_matrix = matrix_exponential(-K * path_length)
```

## Implementation Approach

### Phase 1: Rain Depolarization
1. Implement raindrop axis ratio model
2. Compute differential extinction from existing LUTs or T-matrix
3. Apply polarized extinction in RT solver
4. Validate with dual-pol radar forward models

### Phase 2: Ice Depolarization
1. Add ice crystal differential extinction
2. Handle multiple ice habits
3. Consider orientation effects (links to Issue #02)

### Phase 3: Multiple Scattering Depolarization
1. Track depolarization ratio through scattering
2. Implement ρ_HV-like diagnostic
3. Handle heavy precipitation regimes

### Phase 4: TL/AD
1. Linearize differential extinction
2. Implement adjoints
3. Enable retrievals of rain rate from polarimetric obs

## Coefficient Requirements

### New LUT Entries
For each hydrometeor type, add:
- `Extinction_V(size, frequency, temperature)`
- `Extinction_H(size, frequency, temperature)`
- `Differential_Phase(size, frequency, temperature)`

### Source Options
1. Pre-computed T-matrix results
2. DDA calculations for complex shapes
3. Parameterizations from literature

## Testing

- [ ] Unit test: Spherical particles → no differential extinction
- [ ] Unit test: Oblate rain → κ_H > κ_V
- [ ] Unit test: Extinction matrix symmetry
- [ ] Regression test: Depolarization off = current results
- [ ] Validation: Compare A_DP to radar forward model
- [ ] Validation: Brightness temperature depression in rain

## Sensors Affected

| Sensor | Frequencies | Rain Depol Impact |
|--------|-------------|-------------------|
| GMI | 10-183 GHz | Significant at 10-37 GHz |
| AMSR2 | 6-89 GHz | Significant at lower freq |
| SSMIS | 19-183 GHz | Moderate |
| WindSat | 6-37 GHz | Significant |

## References

- Battaglia, A., et al. (2010). Multiple scattering in observations of the GPM DPR. JGR.
- Bringi, V. N., & Chandrasekar, V. (2001). Polarimetric Doppler Weather Radar.
- Marzano, F. S., et al. (2010). Passive microwave modeling of rain at W band. IEEE TGRS.

## Acceptance Criteria

1. Differential extinction computed for rain
2. Differential extinction computed for ice
3. Polarized extinction used in RT solver
4. Results consistent with radar forward operators
5. TL/AD implemented
6. No regression when depolarization disabled

## Dependencies

- Related to Issue #02 (Oriented Particles) for ice crystals
- May share T-matrix/DDA infrastructure

## Labels

`enhancement`, `MW`, `polarization`, `precipitation`, `cloud-physics`
