# Zeeman Effect Integration with Operational RT Solver

**Parent Issue:** Microwave Full Polarization (4-Stokes) Simulation Capability

## Summary

Integrate the existing Zeeman line-by-line (LBL) code with the operational CRTM forward model to enable polarization coupling from O₂ Zeeman splitting in the microwave. The LBL code exists but is disconnected from the main RT solver.

## Current Implementation

### Zeeman LBL Code (EXISTS)
Located in `src/Zeeman/src_lbl/`:
```fortran
! O2_zeeman_lbl.f90 - Computes Zeeman-split O2 absorption
! Outputs:
!   tau12_Re : Real part of coherent (off-diagonal) transmittance
!   tau12_Im : Imaginary part of coherent transmittance
```

### ODZeeman Module (PARTIALLY IMPLEMENTED)
Located in `src/AtmAbsorption/ODZeeman/`:
- Parameterized Zeeman absorption coefficients
- Used for fast forward model
- **Does NOT output polarization coupling terms**

### Zeeman Input Structure (EXISTS)
```fortran
! src/Options/Zeeman_Input/Zeeman_Input_Define.f90
TYPE :: Zeeman_Input_type
  REAL(fp) :: Be           ! Earth magnetic field strength (Gauss)
  REAL(fp) :: Cos_ThetaB   ! Cosine of angle between B and propagation
  REAL(fp) :: Cos_PhiB     ! Azimuth of B field
END TYPE
```

### Gap: No Connection to RT Solver
The off-diagonal transmittance terms computed in LBL are **not propagated** through the radiative transfer equations. The RT solver treats Zeeman as scalar absorption only.

## Physical Background

### O₂ Zeeman Effect
- Earth's magnetic field splits O₂ rotational lines
- Creates polarization-dependent absorption
- Significant at 50-60 GHz and 118 GHz O₂ bands
- Causes cross-polarization coupling (V↔H mixing)

### Transmittance Matrix
For Zeeman-affected layers, transmittance is a 2x2 matrix (circular basis):
```
T = | τ₊₊   τ₊₋  |
    | τ₋₊   τ₋₋  |
```
Or in linear basis (V/H):
```
T = | τ_VV   τ_VH |
    | τ_HV   τ_HH |
```

Off-diagonal terms cause polarization rotation/mixing.

## Files to Modify

### Zeeman Modules
| File | Changes Needed |
|------|----------------|
| `src/Zeeman/src_lbl/O2_zeeman_lbl.f90` | Create interface for RT integration |
| `src/AtmAbsorption/ODZeeman/ODZeeman_AtmAbsorption.f90` | Output polarization terms |
| `src/AtmAbsorption/CRTM_AtmAbsorption_Define.f90` | Add off-diagonal absorption |

### RT Solver
| File | Changes Needed |
|------|----------------|
| `src/RTSolution/CRTM_RTSolution.f90` | Use matrix transmittance for Zeeman layers |
| `src/RTSolution/RTV_Define.f90` | Add Zeeman transmittance matrix storage |
| `src/RTSolution/Common_RTSolution.f90` | Matrix transmittance propagation |

### Forward Model
| File | Changes Needed |
|------|----------------|
| `src/CRTM_Forward_Module.f90` | Pass Zeeman options to RT |
| `src/Options/CRTM_Options_Define.f90` | Enable Zeeman polarization flag |

## Technical Requirements

### 1. Atmospheric Absorption with Polarization
Extend `CRTM_AtmAbsorption_type` or create new type:
```fortran
TYPE :: Zeeman_Absorption_type
  REAL(fp) :: Optical_Depth_VV    ! V-polarized absorption
  REAL(fp) :: Optical_Depth_HH    ! H-polarized absorption
  REAL(fp) :: Optical_Depth_VH_Re ! Cross-pol (real)
  REAL(fp) :: Optical_Depth_VH_Im ! Cross-pol (imaginary)
END TYPE
```

### 2. Layer Transmittance Matrix
For each Zeeman-affected layer:
```fortran
! Compute 2x2 transmittance matrix
T_matrix(1,1) = EXP(-tau_VV)
T_matrix(2,2) = EXP(-tau_HH)
T_matrix(1,2) = ... ! From tau_VH
T_matrix(2,1) = ... ! From tau_HV
```

### 3. RT Equation Modification
Standard scalar RT:
```
I_up(k) = I_up(k+1) * T(k) + B(k) * (1 - T(k))
```

With Zeeman polarization:
```
I_up(k) = T_matrix(k) @ I_up(k+1) + (I - T_matrix(k)) @ B(k)
```
Where `I` is identity matrix, `@` is matrix multiply.

### 4. Magnetic Field Input
Already exists in `Zeeman_Input_type`. Need to:
- Pass through to RT solver
- Compute at each layer (interpolate from surface/TOA)
- Handle magnetic field model (IGRF or user-specified)

## Implementation Approach

### Phase 1: Interface Design
1. Define polarized absorption structure
2. Create Zeeman RT interface module
3. Add flags to Options for Zeeman polarization

### Phase 2: LBL Integration
1. Wrap LBL code for CRTM calling convention
2. Compute off-diagonal optical depths
3. Store in extended absorption structure

### Phase 3: RT Solver Modification
1. Detect Zeeman-affected channels (50-60 GHz, 118 GHz)
2. Use matrix transmittance for those channels
3. Propagate Stokes vector through layers

### Phase 4: Fast Model (ODZeeman)
1. Extend ODZeeman to output polarization terms
2. Train coefficients for off-diagonal absorption
3. Validate against LBL reference

### Phase 5: TL/AD
1. Linearize matrix transmittance calculations
2. Implement adjoint of matrix operations
3. Validate with finite difference

## Sensors Affected

| Sensor | Channels | Zeeman Impact |
|--------|----------|---------------|
| AMSU-A | 50-60 GHz | Significant |
| ATMS | 50-60 GHz | Significant |
| SSMIS | 50-60 GHz | Significant |
| MHS | 118 GHz (if applicable) | Moderate |

## Testing

- [ ] Unit test: LBL interface returns correct structure
- [ ] Unit test: Matrix transmittance computation
- [ ] Unit test: RT with identity matrix = scalar RT
- [ ] Regression test: Zeeman off reproduces current results
- [ ] Validation: Compare to standalone Zeeman RT codes
- [ ] Validation: AMSU-A 50-60 GHz observations

## References

- Rosenkranz, P. W. (1993). Absorption of microwaves by atmospheric gases. Chapter in Atmospheric Remote Sensing by Microwave Radiometry.
- Lenoir, W. B. (1968). Microwave spectrum of molecular oxygen in the mesosphere. JGR.
- Han, Y., et al. (2007). Zeeman splitting effect on AMSU-A. IEEE TGRS.

## Acceptance Criteria

1. LBL Zeeman code callable from CRTM
2. Off-diagonal optical depths computed and stored
3. RT solver uses matrix transmittance for Zeeman channels
4. Results validated against reference calculations
5. TL/AD implemented and validated
6. No regression when Zeeman polarization disabled

## Estimated Complexity

- Interface design: Moderate
- LBL integration: Low (code exists)
- RT solver modification: High (core algorithm change)
- Fast model extension: High (coefficient training)
- TL/AD: Moderate

## Labels

`enhancement`, `MW`, `polarization`, `zeeman`, `atmospheric-absorption`
