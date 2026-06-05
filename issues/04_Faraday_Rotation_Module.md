# Ionospheric Faraday Rotation Module

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Implement a new module to compute ionospheric Faraday rotation of polarized microwave radiation. This effect rotates the plane of linear polarization as radiation passes through the ionosphere and is significant at frequencies below ~10 GHz.

## Current Implementation

**Not implemented.** No ionosphere module exists in CRTM.

## Physical Background

### Faraday Rotation
When linearly polarized radiation passes through a magnetized plasma (ionosphere), the polarization plane rotates:

```
Ω = (e³ / 8π²ε₀m²c³) × (1/f²) × ∫ Nₑ B_∥ ds
```

Where:
- `Ω` = rotation angle (radians)
- `f` = frequency (Hz)
- `Nₑ` = electron density (m⁻³)
- `B_∥` = magnetic field component along propagation path
- Integration along slant path through ionosphere

### Simplified Form
```
Ω ≈ 1.355 × 10⁻⁵ × (TEC × B_∥) / f²
```

Where:
- `TEC` = Total Electron Content (TECU = 10¹⁶ electrons/m²)
- `B_∥` = in Gauss
- `f` = in GHz
- `Ω` = in radians

### Impact on Stokes Parameters
Faraday rotation transforms Stokes Q and U:
```
Q' = Q cos(2Ω) + U sin(2Ω)
U' = -Q sin(2Ω) + U cos(2Ω)
```
I and V are unchanged.

### Frequency Dependence
| Frequency | TEC=10 TECU | TEC=50 TECU |
|-----------|-------------|-------------|
| 1.4 GHz | ~5° | ~25° |
| 6.9 GHz | ~0.2° | ~1° |
| 10.7 GHz | ~0.1° | ~0.4° |
| 18.7 GHz | ~0.03° | ~0.15° |

Significant for L-band (1.4 GHz), decreasing rapidly with frequency.

## Files to Create/Modify

### New Module
| File | Purpose |
|------|---------|
| `src/Atmosphere/CRTM_Faraday_Module.f90` | Main Faraday rotation calculations |
| `src/Atmosphere/CRTM_Ionosphere_Define.f90` | Ionosphere input structure |

### Modify
| File | Changes Needed |
|------|----------------|
| `src/CRTM_Forward_Module.f90` | Apply Faraday rotation to output |
| `src/CRTM_Tangent_Linear_Module.f90` | TL of Faraday rotation |
| `src/CRTM_Adjoint_Module.f90` | AD of Faraday rotation |
| `src/CRTM_K_Matrix_Module.f90` | K-matrix for TEC Jacobian |
| `src/Options/CRTM_Options_Define.f90` | Add Faraday options |
| `src/RTSolution/CRTM_RTSolution_Define.f90` | Store rotation angle |

## Technical Requirements

### 1. Ionosphere Input Structure
```fortran
TYPE :: CRTM_Ionosphere_type
  REAL(fp) :: TEC                    ! Total Electron Content (TECU)
  REAL(fp) :: Be                     ! Earth magnetic field (Gauss)
  REAL(fp) :: Magnetic_Dip_Angle     ! Magnetic dip angle (degrees)
  REAL(fp) :: Magnetic_Declination   ! Magnetic declination (degrees)
  LOGICAL  :: Use_IGRF               ! Use built-in IGRF model
END TYPE
```

### 2. Faraday Rotation Calculation
```fortran
SUBROUTINE Compute_Faraday_Rotation( &
  Ionosphere,    &  ! Input: ionosphere parameters
  GeometryInfo,  &  ! Input: viewing geometry
  Frequency,     &  ! Input: channel frequency (GHz)
  Rotation_Angle )  ! Output: Faraday rotation (radians)

  ! Compute B parallel to propagation direction
  B_parallel = Ionosphere%Be * COS(dip_angle - zenith_angle)

  ! Compute rotation angle (f² dependence)
  Rotation_Angle = FARADAY_CONSTANT * Ionosphere%TEC * B_parallel / Frequency**2

END SUBROUTINE
```

### 3. Stokes Rotation Matrix
```fortran
SUBROUTINE Apply_Faraday_Rotation( &
  Rotation_Angle, &  ! Input
  Stokes_Vector   )  ! In/Out: [I, Q, U, V]

  cos2Omega = COS(TWO * Rotation_Angle)
  sin2Omega = SIN(TWO * Rotation_Angle)

  Q_new = Stokes_Vector(2) * cos2Omega + Stokes_Vector(3) * sin2Omega
  U_new = -Stokes_Vector(2) * sin2Omega + Stokes_Vector(3) * cos2Omega

  Stokes_Vector(2) = Q_new
  Stokes_Vector(3) = U_new
  ! I and V unchanged

END SUBROUTINE
```

### 4. Integration with RT Solution
Apply Faraday rotation **after** atmospheric RT, before returning to user:
```fortran
! In CRTM_Forward_Module.f90
IF (Options%Apply_Faraday_Rotation) THEN
  CALL Compute_Faraday_Rotation(Ionosphere, GeometryInfo, &
                                 SpcCoeff%Frequency(l), Omega)
  CALL Apply_Faraday_Rotation(Omega, RTSolution%Stokes)
  RTSolution%Faraday_Rotation_Angle = Omega
END IF
```

### 5. Magnetic Field Model (Optional)
Consider including simplified IGRF for magnetic field:
- Could use lookup table vs lat/lon/altitude
- Or require user to provide B field
- IGRF coefficients updated every 5 years

## Implementation Approach

### Phase 1: Basic Module
1. Create `CRTM_Ionosphere_Define.f90` with TEC input
2. Create `CRTM_Faraday_Module.f90` with rotation calculation
3. User provides TEC and B field

### Phase 2: Integration
1. Add Faraday option to `CRTM_Options`
2. Apply rotation in Forward module
3. Store rotation angle in RTSolution

### Phase 3: TL/AD/K
1. Implement tangent-linear (simple - just chain rule)
2. Implement adjoint
3. Add TEC to K-matrix for retrievals

### Phase 4: Magnetic Field Model (Optional)
1. Add IGRF lookup table
2. Compute B field from lat/lon
3. Allow user override

## Data Requirements

### User Input
- TEC: From GPS-derived maps, climatology, or NWP
- Magnetic field: User-provided or computed from IGRF

### Potential Data Sources
- IGS Global TEC maps
- IRI (International Reference Ionosphere) model
- IGRF magnetic field coefficients

## Testing

- [ ] Unit test: Rotation angle calculation vs analytic
- [ ] Unit test: Stokes rotation matrix correctness
- [ ] Unit test: f² frequency dependence
- [ ] Regression test: Faraday off = no change
- [ ] Validation: SMOS/Aquarius L-band polarization
- [ ] Validation: Compare to COSMIC TEC + forward model

## Sensors Affected

| Sensor | Frequency | Faraday Impact |
|--------|-----------|----------------|
| SMOS | 1.4 GHz | **Critical** (degrees) |
| Aquarius | 1.4 GHz | **Critical** |
| SMAP | 1.4 GHz | **Critical** |
| WindSat | 6.8-37 GHz | Moderate at 6.8 GHz |
| AMSR2 | 6.9-89 GHz | Small at 6.9 GHz |

## References

- Yueh, S. H. (2000). Estimates of Faraday rotation with passive microwave polarimetry. IEEE TGRS.
- Le Vine, D. M., & Abraham, S. (2002). The effect of the ionosphere on remote sensing of sea surface salinity. IEEE TGRS.
- Meissner, T., & Wentz, F. J. (2006). Polarization rotation and the third Stokes parameter. IEEE TGRS.

## Acceptance Criteria

1. Ionosphere input structure defined
2. Faraday rotation computed from TEC and B field
3. Rotation applied to Stokes Q and U
4. TL/AD implemented for TEC Jacobians
5. Option to enable/disable Faraday correction
6. Validated against L-band observations

## Estimated Complexity

- Module creation: Low-Moderate
- Forward integration: Low
- TL/AD: Low (analytic derivatives)
- IGRF model: Moderate (optional)
- Validation: Moderate (need L-band data)

## Labels

`enhancement`, `MW`, `polarization`, `ionosphere`, `new-module`
