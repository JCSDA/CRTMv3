# MW Surface Emissivity/Reflection: Extend to Full 4x4 Stokes Matrix

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Extend MW surface emissivity and reflection models to compute the full 4x4 Stokes matrix instead of the current 2-component (V/H) implementation. This is the primary blocker for MW full-polarization simulation.

## Current Implementation

### Surface Emissivity
Only vertical (V) and horizontal (H) components computed:
```fortran
! src/SfcOptics/CRTM_MW_Water_SfcOptics.f90
SfcOptics%Emissivity(i,1)  ! Vertical only
SfcOptics%Emissivity(i,2)  ! Horizontal only
```

### Surface Reflection Matrix
Only diagonal elements populated, all off-diagonal terms = 0:
```fortran
! src/SfcOptics/CRTM_MW_Water_SfcOptics.f90 (lines 235-237)
SfcOptics%Reflectivity(i,1,i,1)  ! V-V
SfcOptics%Reflectivity(i,2,i,2)  ! H-H
! Cross-polarization terms NOT computed
```

## Files to Modify

### Primary Files
| File | Changes Needed |
|------|----------------|
| `src/SfcOptics/CRTM_MW_Water_SfcOptics.f90` | Add Q, U, V emissivity; populate off-diagonal reflection |
| `src/SfcOptics/MW_Water/CRTM_FastemX.f90` | Extend FASTEM model for full Stokes |
| `src/SfcOptics/CRTM_MW_Land_SfcOptics.f90` | Land surface full-pol emissivity |
| `src/SfcOptics/CRTM_MW_Snow_SfcOptics.f90` | Snow surface full-pol emissivity |
| `src/SfcOptics/CRTM_MW_Ice_SfcOptics.f90` | Ice surface full-pol emissivity |

### Supporting Files
| File | Changes Needed |
|------|----------------|
| `src/SfcOptics/CRTM_SfcOptics_Define.f90` | Verify structure supports 4x4 (it does) |
| `src/SfcOptics/MW_Water/FASTEM_MWSSEM/` | Low-level FASTEM modules |

## Technical Requirements

### 1. Stokes Emissivity Vector (4 components)
```
e = [e_I, e_Q, e_U, e_V]^T
```
Where:
- `e_I` = (e_V + e_H) / 2 (intensity, currently computed)
- `e_Q` = (e_V - e_H) / 2 (linear polarization difference)
- `e_U` = cross-polarization component (wind azimuth dependent)
- `e_V` = circular polarization (typically small for ocean)

### 2. Reflection Matrix (4x4)
Full Mueller matrix for surface reflection:
```
     | R_II  R_IQ  R_IU  R_IV |
R =  | R_QI  R_QQ  R_QU  R_QV |
     | R_UI  R_UQ  R_UU  R_UV |
     | R_VI  R_VQ  R_VU  R_VV |
```

For ocean surface, key non-zero terms include:
- Diagonal: R_II, R_QQ, R_UU, R_VV
- Off-diagonal coupling from wind-roughened surface

### 3. Physical Effects to Include

**Ocean Surface:**
- [ ] Fresnel reflection (V/H basis → Stokes basis transformation)
- [ ] Wind-induced surface roughness cross-polarization
- [ ] Foam coverage effects on all Stokes components
- [ ] Azimuthal dependence of U component

**Land Surface:**
- [ ] Vegetation depolarization
- [ ] Soil moisture polarization signature

**Ice/Snow:**
- [ ] Volume scattering depolarization
- [ ] Surface roughness effects

## Implementation Approach

### Phase 1: Ocean Surface (FASTEM)
1. Modify `CRTM_FastemX.f90` to compute full Stokes emissivity
2. Add wind-azimuth dependent U-Stokes calculation
3. Implement 4x4 reflection matrix from Fresnel coefficients
4. Add cross-polarization from roughness model

### Phase 2: Other Surfaces
1. Extend land emissivity model
2. Extend ice/snow models
3. Ensure consistent interface across all surface types

### Phase 3: Tangent-Linear and Adjoint
1. Update `*_TL` routines for all modified forward code
2. Update `*_AD` routines for all modified forward code
3. Verify TL/AD consistency with finite difference tests

## Testing

- [ ] Unit test: Verify 4x4 matrix populated correctly
- [ ] Unit test: Energy conservation (emissivity + reflectivity = 1)
- [ ] Unit test: Reciprocity of reflection matrix
- [ ] Regression test: Compare V/H components to current implementation
- [ ] Validation: Compare against Meissner & Wentz model
- [ ] Validation: Compare against WindSat observations

## References

- Meissner, T., & Wentz, F. J. (2012). The emissivity of the ocean surface between 6 and 90 GHz. IEEE TGRS.
- Yueh, S. H. (1997). Modeling of wind direction signals in polarimetric sea surface brightness temperatures. JGR.
- Johnson, J. T. (2006). Third Stokes parameter emission from a periodic water surface. IEEE TGRS.

## Acceptance Criteria

1. `SfcOptics%Emissivity(i,1:4)` populated for all MW surface types
2. `SfcOptics%Reflectivity(i,j,k,l)` full 4x4 matrix for MW
3. TL/AD routines updated and validated
4. No regression in existing 2-component (V/H) results
5. Unit tests passing

## Labels

`enhancement`, `MW`, `polarization`, `surface-optics`
