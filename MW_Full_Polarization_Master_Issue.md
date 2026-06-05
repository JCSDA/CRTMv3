# Microwave Full Polarization (4-Stokes) Simulation Capability

## Overview

This master issue tracks the implementation of robust full polarization (4-Stokes: I, Q, U, V) simulation capability for the microwave (MW) spectral region in CRTM v3.

Currently, CRTM v3 has infrastructure for 4-Stokes simulation (`MAX_N_STOKES=4`, `Options%n_Stokes`), but MW is fundamentally implemented as a **2-component model** (V/H polarizations only). This issue coordinates the work needed to enable operational full-polarization MW simulation.

## Current State

- UV/VIS: Full Stokes capable
- IR: Full Stokes capable
- **MW: Limited to 2-component (V/H) only**

### What Works
- `n_Stokes` parameter exists and can be set to 4
- RT solvers (ADA, AMOM) have some multi-Stokes infrastructure
- Zeeman LBL code exists (but disconnected)

### What's Missing
- Surface emissivity models output only 2 components (V/H)
- No oriented particle scattering (6 elements sufficient for random orientation)
- No Faraday rotation
- Zeeman not integrated with operational RT
- No cross-polarization coupling at surface

---

## Priority Items

### High Priority

#### 1. Extend MW Surface Emissivity/Reflection Models to Full 4x4 Matrix
**Files:** `src/SfcOptics/CRTM_MW_Water_SfcOptics.f90`, `CRTM_FastemX.f90`, MW land/ice/snow modules

**Current:** Only diagonal elements populated
```fortran
SfcOptics%Emissivity(i,1)  ! V only
SfcOptics%Emissivity(i,2)  ! H only
SfcOptics%Reflectivity(i,1,i,1)  ! V-V only
SfcOptics%Reflectivity(i,2,i,2)  ! H-H only
! All off-diagonal = 0
```

**Required:**
- [ ] Compute Q, U, V Stokes emissivity components
- [ ] Populate full 4x4 reflection matrix including off-diagonal terms
- [ ] Add cross-polarization coupling from surface roughness/foam
- [ ] Update FASTEM and all MW surface models

---

#### 2. Add Oriented Particle Scattering Capability
**Files:** `src/AtmScatter/CRTM_CloudScatter.f90`, `src/RTSolution/Common_RTSolution.f90`, CloudCoeff

**Current:**
- 6 phase elements sufficient for randomly oriented symmetric particles (spheres, spheroids)
- RT solver assumes symmetry-based parameterization (α₁, α₂, α₃, α₄, β₁, β₂)
- No support for preferentially oriented particles

**Note:** The existing 6-element structure correctly handles the full 4×4 Mueller matrix for randomly oriented particles due to symmetry constraints. Expansion to 16 elements is NOT needed for this case.

**Required for oriented particles (ice crystals, canted raindrops):**
- [ ] Add orientation distribution parameters to cloud/aerosol definitions
- [ ] Implement orientation-dependent scattering calculations
- [ ] Add canting angle distribution support
- [ ] Create new coefficient data for oriented particle optical properties
- [ ] Modify RT solver to handle asymmetric scattering matrices when needed

---

#### 3. Integrate Zeeman LBL Code into Operational RT Solver
**Files:** `src/Zeeman/src_lbl/O2_zeeman_lbl.f90`, `src/AtmAbsorption/ODZeeman/`, `src/RTSolution/`

**Current:** Zeeman LBL computes off-diagonal transmittance (`tau12_Re`, `tau12_Im`) but is disconnected from main RT

**Required:**
- [ ] Create interface between Zeeman module and CRTM forward model
- [ ] Propagate off-diagonal transmission matrix through RT equations
- [ ] Add magnetic field vector as input to operational CRTM calls
- [ ] Enable Zeeman-induced cross-polarization in RTSolution

---

#### 4. Implement Ionospheric Faraday Rotation Module
**Files:** New module needed in `src/`

**Current:** Not implemented - no ionosphere module exists

**Required:**
- [ ] Create `CRTM_Faraday_Module.f90`
- [ ] Add TEC (Total Electron Content) as input parameter
- [ ] Compute rotation angle from magnetic field and electron density
- [ ] Apply polarization rotation matrix in RT solution
- [ ] Add gyrofrequency calculations

---

### Medium Priority

#### 5. Implement Rain/Ice Depolarization Physics
**Files:** `src/AtmScatter/`, RT solvers

**Current:** No differential extinction or depolarization from precipitation

**Required:**
- [ ] Add rain-induced depolarization model
- [ ] Implement cross-polarization extinction
- [ ] Add raindrop shape (oblate spheroid) effects
- [ ] Model ice particle depolarization

---

#### 6. Add Circular Polarization Generation Mechanisms
**Files:** Surface modules, RT solvers

**Current:** No mechanism to generate Stokes V

**Required:**
- [ ] Identify physical sources of circular polarization in MW
- [ ] Implement V-Stokes generation from relevant processes
- [ ] Add circular polarization coupling in reflection matrix

---

### Architecture Changes

#### 7. RT Solver Updates
**Files:** `src/RTSolution/UWisc_SOI/SOI_CRTM_Forward_Module.f90`, `src/RTSolution/CRTM_RTSolution.f90`

**Current:** SOI hardcoded to `SFCOPTICS_N_STOKES = 2`

**Required:**
- [ ] Remove hardcoded 2-Stokes limitation in SOI
- [ ] Ensure all RT solvers support n_Stokes=4 for MW
- [ ] Update surface boundary conditions for full 4x4 matrix
- [ ] Validate ADA/AMOM/SOI with 4-Stokes MW

---

#### 8. Testing and Validation
**Files:** `test/mains/`

**Current:** No regression tests for MW full-polarization

**Required:**
- [ ] Add unit tests for 4-Stokes MW surface emissivity
- [ ] Add regression tests for Zeeman polarization
- [ ] Add validation against WindSat/SSMIS polarimetric observations
- [ ] Create benchmark cases for Faraday rotation

---

## Sub-Issues

This master issue will be broken into the following sub-issues:

| Priority | Sub-Issue | Components |
|----------|-----------|------------|
| High | #xxx | MW Surface 4x4 Reflection Matrix |
| High | #xxx | Oriented Particle Scattering |
| High | #xxx | Zeeman Integration with RT Solver |
| High | #xxx | Faraday Rotation Module |
| Medium | #xxx | Rain/Ice Depolarization |
| Medium | #xxx | Circular Polarization Mechanisms |
| Arch | #xxx | RT Solver 4-Stokes Support |
| Test | #xxx | MW Full-Pol Validation Suite |

---

## References

- WindSat polarimetric radiometer documentation
- SSMIS polarimetric channels
- Meissner & Wentz ocean emissivity model
- Rosenkranz O2 Zeeman absorption model

## Related Sensors

Full-pol MW capability would benefit simulation of:
- WindSat (all 4 Stokes)
- SSMIS (V, H, plus some polarimetric)
- Future polarimetric MW radiometers

---

## Labels

`enhancement`, `MW`, `polarization`, `multi-issue`
