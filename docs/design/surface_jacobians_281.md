# Surface-Property Jacobians — Assessment & Implementation Plan (issue #281)

Status: living design doc. Tracks the analytic TL/AD effort for the surface
parameters listed in [JCSDA/CRTMv3#281](https://github.com/JCSDA/CRTMv3/issues/281).

## Progress (2026-05-26)

- **Phase 0 (done):** `test/mains/unit/Unit_Test/test_Land_Jacobian.f90` (ctest
  `test_Unit_Land_Jacobian`) validates analytic MW-land `LAI`/`Vegetation_Fraction`
  Jacobians from both the K-matrix (adjoint) and tangent-linear against central
  finite differences. Passes (agreement to ~5 sig figs on the surface-sensitive
  amsua_n19 window channels).
- **Phase 1 (done):** analytic TL/AD for MW-land `LAI` + `Vegetation_Fraction`.
  - `NESDIS_LandEM` gained optional analytic outputs `dEV_dvlai`/`dEH_dvlai`
    (via new optional derivatives in `Canopy_Optic` and `Two_Stream_Solution`);
    forward emissivity is bit-identical.
  - `CRTM_MW_Land_SfcOptics` `iVar_type` now caches the per-angle partials; the
    `_TL`/`_AD` zero-stubs are replaced with real routines.
  - The 3 dispatcher call-sites in `CRTM_SfcOptics.f90` pass `iVar%MWLSOV`,
    `Surface_TL`/`Surface_AD` (templated on MW water).
- **Reference-data consequence:** 18 MW-over-land `k_matrix`/`adjoint` reference
  tests (atms_n21, atms_npp, amsua_n19) now differ because `Surface_K`/`Surface_AD`
  `LAI`+`Vegetation_Fraction` are correctly nonzero (were 0). `Atmosphere_K` and
  `RTSolution_K` are unchanged; forward radiances are bit-identical. The `.bin`
  references are build-local (not git-tracked, regenerated on first run); fresh
  references must be accepted to green the suite.
- **Phase 2 (done):** analytic TL/AD for `Soil_Moisture_Content` (mv -> `Soil_Diel`
  esoil -> roughened r23 -> two-stream). Validated in `test_Land_Jacobian`
  (dTb/dSMC ~ -45 K/K at 23.8 GHz; AD == TL == FD).
- **Phase 3 (done, 2026-07-19):** analytic TL/AD for `Soil_Temperature` and the
  emissivity part of `Land_Temperature`. Temperature enters the emissivity by two
  paths: the soil dielectric (`Soil_Diel` `eswo`/`tauw` polynomials, reusing the
  Phase-2 r23 chain now factored into `Roughened_R23_Deriv`) and the canopy/soil
  thermal ratio `gsect0` in `Two_Stream_Solution` (the dominant term). The
  `Soil_Temperature` out-of-range aliasing (`t_soil <- t_skin`) is resolved in the
  forward: the aliased input gets a zero derivative and its sensitivity is
  re-attributed to `Land_Temperature`. The `Land_Temperature` emissivity part
  ACCUMULATES onto the existing skin-T Planck emission Jacobian
  (`CRTM_Compute_SurfaceT_AD`). Also fixes a latent bug: the forward already
  depended on skin T via `gsect0`, but that term was dropped from
  `Surface_K%Land_Temperature`, making the below-cutoff Land_Temperature Jacobian
  ~3-4x too large; it now matches FD. `Canopy_Water_Content` confirmed as a
  structural zero (never consumed by the forward; asserted in the test). Branch
  `feature/btj_landem_temperature_jacobians`, commit `aa31bb7`; suite 215/215
  after regenerating the MW-over-land Surface references.
- **Land MW analytic surface Jacobian set is now COMPLETE below 80 GHz**: LAI,
  Vegetation_Fraction, Soil_Moisture_Content, Soil_Temperature, Land_Temperature
  all analytic and FD-validated; Canopy_Water_Content structurally zero. Above
  80 GHz LandEM is gated off (constant emissivity) by design.

## TL;DR / scope correction

The issue's surface list reads like "a few missing Jacobians," but the audit shows it is
really two different problems:

1. **A plumbing problem (easy, well-templated):** the MW Land/Snow/Ice SfcOptics TL/AD
   routines are *pure zero-stubs* — they don't even take a `Surface` argument and just set
   emissivity/reflectivity to 0. The surrounding machinery (dispatcher, K-matrix, adjoint
   driver) already carries everything needed.
2. **A physics problem (the real cost):** the underlying NESDIS emissivity models
   (`NESDIS_LandEM`, the AMSU/ATMS/SSMI/SSMIS SnowEM & SeaIceEM modules) **have no
   tangent-linear/adjoint code at all** (the only existing linearization is `OceanEM_TL_SSTW`
   for water). Producing analytic Jacobians means hand-differentiating those models.

Of the 12 surface parameters listed, the honest breakdown is:

- **3 are genuinely tractable analytically and high value** (MW land: LAI, Vegetation_Fraction,
  Soil_Moisture_Content).
- **4 are tractable only in restricted regimes** (soil/land temp emissivity sensitivity;
  snow depth & ice thickness via the shared LandEM path — TBs frozen as observations).
- **5 are NOT implementable without a model change** because the forward model never consumes
  them: Canopy_Water_Content, Snow_Density, Snow_Grain_Size (MW), Ice_Density, Ice_Roughness.

## What works today (so we don't re-derive it)

| Sensitivity | Status | Mechanism |
|---|---|---|
| Land/Water/Snow/Ice **surface temperature** (emission term dI/dTs) | works | `CRTM_Compute_SurfaceT_AD` — `CRTM_SfcOptics.f90:340` |
| Water: **Salinity, Wind_Speed, Wind_Direction** | works | FASTEM / `OceanEM_TL_SSTW` — `CRTM_MW_Water_SfcOptics.f90:705` |
| IR Snow: **Snow_Grain_Size, Snow_Temperature** | works | IRSSEM — `CRTM_IR_Snow_SfcOptics.f90:606` |

Important nuance: surface *temperature* already produces a Jacobian via the Planck/emission
term. What's missing for temperatures is only the *emissivity's* dependence on temperature
(a secondary effect). The structural parameters (soil moisture, veg, LAI, snow depth, ice
thickness) produce **exactly zero** today.

## The machinery is already wired — only the leaves are stubs

The dispatcher and drivers already pass the surface state and its perturbations all the way
down; the sub-model stubs simply discard them:

- `CRTM_Compute_SfcOptics_TL(Surface, SfcOptics, Surface_TL, GeometryInfo, SensorIndex, ChannelIndex, SfcOptics_TL, iVar)` — `CRTM_SfcOptics.f90:1231`
- `CRTM_Compute_SfcOptics_AD(Surface, SfcOptics, SfcOptics_AD, GeometryInfo, SensorIndex, ChannelIndex, Surface_AD, iVar)` — `CRTM_SfcOptics.f90:1857`
- The dispatcher `iVar_type` already has dedicated per-submodel cache slots:
  `MWLSOV` (land), `MWSSOV` (snow), `MWISOV` (ice) — currently empty placeholders.
- **MW_Water is a complete working template** for forward→TL→AD (`CRTM_SfcOptics.f90:544 / 1322 / 2191`).
- **K-matrix needs no changes:** `Surface_K(ln,m)` is passed directly as the `Surface_AD`
  output of `CRTM_Compute_SfcOptics_AD`, and the `+`/`-`/zero operators in
  `CRTM_Surface_Define.f90:1685-1769` already cover all 12 fields. The moment a sub-model AD
  writes `Surface_AD%Soil_Moisture_Content`, it flows into `Surface_K` automatically.

So per surface family the wiring work is mechanical: cache `iVar` in the forward routine,
give the TL/AD routines the MW_Water signature, and update the 3 dispatcher call-sites.

## Feasibility matrix (verified against the forward code)

### Land — `NESDIS_LandEM` (smooth: dielectric mixing + Fresnel + two-stream)

| Parameter | Flows to emissivity? | Nature of dependence | Difficulty |
|---|---|---|---|
| **LAI** | Yes (`vlai = LAI*veg_frac`, `LandEM:210`) | Linear in canopy optical depth; all downstream smooth | **Easy — DONE (Phase 1)** |
| **Vegetation_Fraction** | Yes (`LandEM:205,210`) | Linear; only caveat is `MIN/MAX` clip at [0,1] | **Easy — DONE (Phase 1)** |
| **Soil_Moisture_Content** | Yes (`Soil_Diel`, `LandEM:220/429/436`) | Smooth dielectric mixing; clip at [0,1] + `vmc>0` guard | **Moderate — DONE (Phase 2)** |
| **Soil_Temperature** | Yes — soil dielectric + `gsect0` thermal ratio; aliased to skin T when out of [100,350] | Smooth in range; aliasing resolved to a zero derivative + re-attribution to Land_Temperature | **Hard — DONE (Phase 3)** |
| **Land_Temperature** (emissivity part) | Yes (`gsect0` thermal ratio in two-stream) | Smooth; accumulates onto the skin-T emission Jacobian | **Hard — DONE (Phase 3)** |
| **Canopy_Water_Content** | **No** — never passed into `Compute_MW_Land_SfcOptics`; absent from `NESDIS_LandEM` | n/a | **Structural zero (asserted in test)** |

`NESDIS_LandEM` is also the fallback for the snow and ice angle corrections, so differentiating
it once benefits multiple paths.

### Snow (MW) — sensor-specific, brightness-temperature-driven

The MW snow chain is dominated by discrete snow-type **classification from observed TBs**,
QC gates (e.g. Ts in [150,290] K), depth thresholds (0.1/10/50/100 mm), and lookup tables.
Snow properties enter mostly as perturbations within narrow regimes.

| Parameter | Flows to emissivity? | Difficulty |
|---|---|---|
| **Snow_Depth** | Yes, but through threshold branches + lookups; smooth only in LandEM shallow-snow regime (<=10 mm) and `ems_adjust` polynomials | **Hard / restricted** |
| **Snow_Temperature** (emissivity part) | Yes, partly smooth, many hard gates | **Hard / restricted** |
| **Snow_Density** | **No** — MW models never read it (IR-only) | **Not feasible** |
| **Snow_Grain_Size** | **No** in MW (IR path already has it) | **Not feasible (MW)** |

### Ice (MW) — empirical TB regression + discrete ice-type `minloc` classification

| Parameter | Flows to emissivity? | Difficulty |
|---|---|---|
| **Ice_Temperature** (emissivity part) | Yes; smooth for SSMI/SSMIS, classification-discontinuous for ATMS/AMSU/MHS | **Moderate->Hard** |
| **Ice_Thickness** | Only SSMI/SSMIS, via LandEM angle correction (clamped [0.1,10] mm); other sensors use fixed values | **Hard / restricted, small effect** |
| **Ice_Density** | **No** — never consumed anywhere | **Not feasible** |
| **Ice_Roughness** | **No** — never consumed (LandEM uses hardcoded sigma=1) | **Not feasible** |

## Recommended plan (analytic TL/AD)

### Phase 0 — Validation oracle first
Build a finite-difference Jacobian check (perturb each `Surface` field → forward → compare to
K-matrix). This is the *test* for the analytic code, not a substitute for it; it mirrors CRTM's
existing TL/AD consistency tests. Without it, analytic surface Jacobians can't be trusted.
- New test under `test/mains/regression/k_matrix/` (or `unit/`) exercising an MW land scene.

### Phase 1 — Land tier 1: LAI + Vegetation_Fraction (Easy, highest value/effort ratio)
1. Define `MWLSOVar_type` fields to cache canopy/two-stream intermediates (replace the empty
   placeholder `iVar_type` in `CRTM_MW_Land_SfcOptics.f90:87`).
2. Add an `iVar` (OUT) argument to forward `Compute_MW_Land_SfcOptics` and cache intermediates.
3. Write `NESDIS_LandEM_TL` / `NESDIS_LandEM_AD` covering the canopy-optics + two-stream path
   (`Canopy_Optic`, `Two_Stream_Solution`, `Reflectance`).
4. Replace the `Compute_MW_Land_SfcOptics_TL/_AD` zero-stubs with real routines using the
   MW_Water signature `(..., Surface_TL/Surface_AD, ..., iVar)`.
5. Update the 3 dispatcher call-sites (`CRTM_SfcOptics.f90:518 / 1298 / 2225`) to pass
   `iVar%MWLSOV`, `Surface_TL`, `Surface_AD` — copy the MW_Water blocks verbatim.
6. Validate against Phase-0 FD; confirm full ctest suite still passes (no common-suite sensor
   uses MW land emissivity Jacobians, so regressions should be inert — verify).

### Phase 2 — Land tier 2: Soil_Moisture_Content (Moderate) — DONE (2026-05-26)
Extended the LandEM TL/AD through the complex `Soil_Diel` dielectric, the Fresnel reflectances,
the roughness mixing, and the two-stream `r23` dependence:
- `Soil_Diel` returns optional `desm_dvmc` (complex `d(esm)/d(vmc)`).
- `Two_Stream_Solution` returns optional `desv_dr23v`/`desh_dr23h`.
- `NESDIS_LandEM` assembles `dEV_dmv`/`dEH_dmv` via the chain mv -> esoil -> {theta_t,
  Fresnel rv0/rh0} -> roughened r23 -> emissivity (Fresnel derivative inlined).
- `CRTM_MW_Land_SfcOptics` caches the per-angle partials and adds the soil-moisture term to
  TL/AD (with the [0,1] clip handled as a zero-derivative region).
- `test_Unit_Land_Jacobian` extended to validate `Soil_Moisture_Content` (AD + TL vs FD);
  passes (e.g. amsua_n19 ch1 dTb/dSMC ~ -45 K per unit). Full ctest 200/200 after regenerating
  the MW-over-land `Surface` references (now also carry a nonzero Soil_Moisture_Content field).

### Phase 3 (optional) — Soil/Land temperature emissivity sensitivity (Hard)
Add the emissivity-vs-temperature term to the existing temperature Jacobian. Document the
threshold caveats (100/350/280 K). Lower priority since the dominant temperature Jacobian
already exists via `CRTM_Compute_SurfaceT_AD`.

### Phase 4 (optional, restricted) — Snow_Depth & Ice_Thickness via shared LandEM path (Hard)
Reuse the Phase-1/2 `NESDIS_LandEM_TL/AD` for the LandEM-default snow/ice angle-correction
path only. **Treat observed brightness temperatures (`Surface%SensorData%Tb`) as frozen
observations — do not differentiate them.** Document validity limits (shallow snow, valid TB
ranges, specific sensors) and accept zero/undefined derivatives at classification thresholds.

### Explicitly out of analytic scope (recommend documenting in the issue)
Canopy_Water_Content, Snow_Density, Snow_Grain_Size (MW), Ice_Density, Ice_Roughness — these
have **no forward dependence**, so their analytic Jacobian is identically zero. They cannot be
made nonzero without first adding the physics to the forward model.

## Effort & risk
- **Phase 0 + Phase 1** (FD test + LAI/Veg analytic TL/AD, one surface family, full validation):
  the natural first deliverable; self-contained and templated on MW_Water.
- **Phase 2**: incremental on the same LandEM TL/AD.
- **Phases 3–4**: higher effort, restricted validity, lower marginal value — defer until there's
  a concrete DA use case.
- Main risk is correctness of hand-coded adjoints; the Phase-0 FD oracle plus CRTM's existing
  TL/AD consistency-test pattern mitigates it.

## Implementation notes (for the LAI/Veg analytic path)

In the no-snow land branch of `NESDIS_LandEM`, LAI and Vegetation_Fraction enter **only** through
`vlai = Lai*veg_frac` (`LandEM:210`), which feeds `Canopy_Optic` -> `tau_v = tau_h =
0.5*vlai*(2-tv-th)` (`LandEM:349`) -> `Two_Stream_Solution`. Everything upstream (soil/canopy
dielectric, reflectances, single-scatter albedo, `gv/gh`, `kk`, `beta`, `gamma`, `gsect0`) is
independent of LAI/veg. So the analytic derivative reduces to:

```
dEmis/dvlai = (d Emis / d tau) * (d tau / d vlai),   d tau/d vlai = 0.5*(2 - tv - th)
dvlai/dLAI  = veg_frac_clipped
dvlai/dVeg  = LAI            (zero if veg_frac is clipped at 0 or 1)
```

Forward caches per-angle `dEV_dvlai`, `dEH_dvlai`, the used `LAI`/`veg_frac`, the veg-clip flag,
and whether the low-frequency model path was taken; TL/AD just apply the chain rule. Specular
reflectivity gives `dReflectivity = -dEmissivity`. This mirrors the partial-derivative caching
idiom already used by the MW water low-frequency path.
