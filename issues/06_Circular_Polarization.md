# Circular Polarization (Stokes V) Generation Mechanisms

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Identify and implement physical mechanisms that generate circular polarization (Stokes V) in the microwave. Currently, CRTM has no mechanism to produce non-zero Stokes V, limiting full 4-Stokes simulation capability.

## Current Implementation

**No Stokes V generation.** All sources produce V=0:
- Surface emission: V/H only, no circular component
- Scattering: Phase matrix element β₂ (P₃₄) exists but not exercised
- No coupling mechanisms implemented

## Physical Sources of Circular Polarization in MW

### 1. Zeeman Effect (Primary Source)
O₂ absorption lines in magnetic field create V/H asymmetry:
- Circular polarization from σ± transitions
- Significant at 50-60 GHz, 118 GHz
- **Covered by Issue #03 (Zeeman Integration)**

### 2. Oriented Hydrometeors with Canting
Tilted non-spherical particles can convert linear to circular:
- Canted raindrops
- Ice crystals with preferred tilt
- Requires non-zero mean canting angle
- **Related to Issue #02 (Oriented Particles)**

### 3. Surface Reflection (Small Effect)
Complex dielectric → slight circular component:
- Generally very small for natural surfaces
- Potentially measurable for specific conditions

### 4. Faraday Rotation + Differential Emission
Combination effect:
- Differential Q emission from surface
- Faraday rotation couples Q→U
- Does NOT create V directly
- **Related to Issue #04 (Faraday)**

### 5. Multiple Scattering
Successive scattering events can generate V:
- Phase matrix β₂ elements couple U↔V
- Requires significant optical depth
- Currently supported in RT solver but not exercised

## Priority Assessment

| Mechanism | MW Significance | Already Addressed |
|-----------|-----------------|-------------------|
| Zeeman | **High** (50-60 GHz) | Issue #03 |
| Oriented particles | **Moderate** | Issue #02 |
| Multiple scattering | Low-Moderate | RT solver ready |
| Surface reflection | Low | Not critical |

## Files to Modify

### For Multiple Scattering V Generation
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/Common_RTSolution.f90` | Verify β₂ elements used |
| `src/AtmScatter/CRTM_CloudScatter.f90` | Ensure P₃₄ coefficients computed |

### For Surface Circular Polarization (Low Priority)
| File | Changes Needed |
|------|----------------|
| `src/SfcOptics/CRTM_MW_Water_SfcOptics.f90` | Add V-component if significant |

## Technical Details

### Phase Matrix Element for V
The β₂ (element 6) couples Stokes U and V:
```fortran
! From Common_RTSolution.f90 (lines 2059-2075)
IF( RTV%n_Stokes == 4 ) THEN
  ! beta2 (3,4)
  RTV%Pff(i1+2,j1+3,k) = ... Phase_Coefficient(l,6,k) ...
  ! beta2 (4,3)
  RTV%Pff(i1+3,j1+2,k) = ... Phase_Coefficient(l,6,k) ...
  ! alpha4 (4,4)
  RTV%Pff(i1+3,j1+3,k) = ... Phase_Coefficient(l,4,k) ...
END IF
```

**The infrastructure exists** - need to ensure:
1. Phase coefficients (elements 4 and 6) are non-zero in CloudCoeff
2. n_Stokes = 4 is set
3. Source term includes V component

### Circular Polarization from Canting
Mean canting angle θ_c creates V from linear polarization:
```
V_scattered ∝ sin(2θ_c) × U_incident
```

For randomly canted particles (θ_c = 0), no net V is produced.
Non-zero mean canting → systematic V generation.

## Implementation Approach

### Phase 1: Verify Existing Infrastructure
1. Confirm β₂ phase elements computed in CloudCoeff
2. Test RT solver with n_Stokes=4 and artificial V source
3. Verify V propagates correctly through layers

### Phase 2: Enable Scattering-Generated V
1. Ensure CloudCoeff has non-zero P₃₄ elements
2. Connect to oriented particle model (Issue #02)
3. Test V generation from U scattering

### Phase 3: Surface V (If Needed)
1. Literature review of MW surface circular polarization
2. Implement if observationally significant
3. Low priority unless specific application

## Testing

- [ ] Unit test: V=0 for symmetric scattering
- [ ] Unit test: V≠0 for canted particles
- [ ] Unit test: V propagation through layers
- [ ] Integration test: Full 4-Stokes with scattering
- [ ] Validation: Compare to Zeeman V observations (Issue #03)

## Sensors Measuring Stokes V

| Sensor | V Capability | Primary V Source |
|--------|--------------|------------------|
| WindSat | Yes (all channels) | Surface + atmosphere |
| SSMIS | Partial | Zeeman (50-60 GHz) |
| Future pol radiometers | Planned | TBD |

## References

- Gasiewski, A. J. (1993). Microwave radiative transfer in hydrometeors. Chapter in Atmospheric Remote Sensing.
- Troitsky, A. V., et al. (2003). Circular polarization of thermal microwave radiation. Radio Science.
- Mätzler, C. (2006). Thermal Microwave Radiation. IET.

## Acceptance Criteria

1. Existing V infrastructure verified working
2. Scattering can generate V when conditions met
3. Integration with oriented particles (Issue #02)
4. Integration with Zeeman (Issue #03)
5. Documentation of V generation mechanisms

## Relationship to Other Issues

This issue is largely **dependent on other issues**:
- **Issue #02 (Oriented Particles)**: Canting → V generation
- **Issue #03 (Zeeman)**: Primary V source at 50-60 GHz
- **Issue #05 (Depolarization)**: Scattering effects

Consider this a **verification and integration** issue rather than new development.

## Estimated Complexity

- Infrastructure verification: Low
- Scattering V testing: Low
- Surface V (if needed): Low-Moderate
- Documentation: Low

## Labels

`enhancement`, `MW`, `polarization`, `stokes-V`, `integration`
