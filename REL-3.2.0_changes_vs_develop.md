# `feature/btj_REL-3.2.0` vs `develop`: code/data changes and the regression-test brightness-temperature differences they produce

This issue catalogs the changes on `feature/btj_REL-3.2.0` relative to `develop` and identifies, for each, the brightness-temperature (TB) differences it produces in the common `ctest` regression suite (the `forward` / `tangent_linear` / `k_matrix` / `adjoint` sensor sweeps). The goal is a reference I can point others to when explaining "why did test X change."

> **Refreshed 2026-06-11** to cover the full branch through REL-3.2.0 tag prep. The original catalog (items 1–6) predated the branch's second half; this revision adds the downwelling/upwelling radiance-profile outputs (§10), the `n_Stokes > 1` vector-RT fixes (§8, #318), the TELSEM2 MW-land atlas (§8, #314), the experimental `CRTM-Exp` cloud optics + DDA-ARTS ICE_CLOUD change (§8, #320), the grazing-angle reflectivity guards (§8), and the pre-release review fixes (§11; see the Round 1 / Round 2 comments below). **None of the additions changes any common-suite TB** — they are all opt-in, land/≥200-GHz-only, scalar-path-bit-identical, or new output fields. The md5 note at the foot is also corrected to the current tarball.

## Summary of TB-affecting changes

| # | Change | Common ctests affected | Nature of the TB difference |
|---|--------|------------------------|------------------------------|
| 1 | Coefficient I/O switched to NetCDF by default; new `fix_REL-3.2.0.0.tgz` fix tree; `.nc4` → `.nc` | All sensors load NetCDF coeffs now | For most sensors the `.nc` and `.bin` coefficients are equivalent and TB is unchanged; the exceptions are items 2–4 |
| 2 | SpcCoeff NetCDF reader now loads the `NLTECoeff` / `ACCoeff` sibling files | CrIS-FSR (`cris-fsr_n21`, `cris399_npp`) — forward/TL/K/adjoint | Non-LTE radiance correction is now applied for these sensors; previously the NetCDF SpcCoeff reader returned an empty NLTE sub-structure so the correction was a no-op. Radiance/BT shift in the NLTE-sensitive channels (shortwave CO₂ band) |
| 3 | Canonical `v.abi_g18.SpcCoeff.nc` carries per-channel `Solar_Irradiance` = 1.035050 for ABI bands 1–6 | `v.abi_g18` — forward ×8, k_matrix ×7, adjoint ×2, tangent_linear ×2 | Solar source term changes ~0.03–0.1% per channel vs prior values; cascades into `SOD`, `Layer_Optical_Depth`, `Single_Scatter_Albedo`, and reflective-band `Radiance`/`Brightness_Temperature` (BT shifts ~0.006–0.07 K) |
| 4 | WMO satellite/sensor ID values in the canonical coeff files | `abi_g18` (272/617 vs the prior −999 sentinels), `cris-fsr_n21` (WMO satellite 226 vs 224) | No TB change, but the regression comparison includes the `RTSolution` metadata, so these tests differ on that field alone unless the field is excluded |
| 5 | Analytic MW-land emissivity Jacobians (issue #281, `9358caa`+`aa31bb7`) | none in the common ocean sweep; **MW-over-land Surface reference data** (k_matrix/adjoint/TL) | The MW-land emissivity TL/AD were previously identically-zero stubs; they now return analytic sensitivities to **every physical land-state variable** — `LAI`, `Vegetation_Fraction`, `Soil_Moisture_Content`, `Soil_Temperature`, and `Land_Temperature` (the last a correctness fix, not just a new output). Forward emissivity is unchanged (bit-identical), so only the Surface-Jacobian reference fields for MW-over-land scenes move. |
| 6 | CONST_MIXED_POLARIZATION (pol type 13) Distance_Ratio fix (`bdc7fb9`) | `tms_*` / TROPICS only | Dropped the erroneous `GeometryInfo%Distance_Ratio` scaling of the fixed polarization angle in the SfcOptics FWD/TL/AD pol-13 branches. Pol type 13 is TMS/TROPICS-only, so no common-suite sensor is affected; TMS/TROPICS BT changes. |

Everything else on the branch is either inert in the common test configuration, changes the default path only for MW-water channels ≥ 200 GHz (the PARMIO backend — no common-suite sensor reaches that frequency, see §8), adds new output fields without touching the existing ones (downwelling/upwelling profiles, §10), is opt-in or land/DDA-only (TELSEM2 and `CRTM-Exp`/DDA-ICE, §8), affects only the `n_Stokes > 1` path the scalar suite never runs (§8), or only changes behavior on an error/edge path the standard scenes never exercise — see §8.

## 1. NetCDF coefficient transition

Commits: `3ac90dc` (NetCDF as the default LUT format), `69e752f` (`.nc4` → `.nc` extension), `0fcb7d4` (Zeeman SSMIS TauCoeff → NetCDF), `2a7495a` + `ODSSUBIN2NC` / `ODSSU_netCDF_IO` (ODSSU SSU TauCoeff NetCDF I/O + converter), `c88c79c` / `33a23da` / `5e0a69d` (new `fix_REL-3.2.0.0.tgz` + md5sums in `Get_CRTM_Binary_Files.sh` and `test/CMakeLists.txt`), `b6168df` / `12b0394` / `34fbbeb` (tests read coeffs from the canonical `test_data` tree).

Code changes:

* `CRTM_LifeCycle.f90`: the default coefficient format flips from `Binary` to `netCDF` for `SpcCoeff`, `TauCoeff`, `CloudCoeff`, `AerosolCoeff`, and all the IR/VIS/MW EmisCoeff files, with a `Resolve_Coeff_Format` step that falls back to `Binary` if the `.nc` file is absent but the `.bin` equivalent is present.
* `CRTM_SpcCoeff.f90` / `CRTM_TauCoeff.f90`: the readers default to NetCDF (`netCDF` argument absent ⇒ NetCDF) and probe per-sensor (SpcCoeff) / per-batch (TauCoeff: ODAS/ODPS/ODSSU/ODZeeman loaders take a single `netCDF` flag, so they switch the whole batch), falling back to the alternate format if the requested one isn't on disk. ODZeeman has its own `z<sensor>.TauCoeff.nc` probe.

  *Zeeman netCDF status (verified against `fix_REL-3.2.0.0.tgz`):* the **SSMIS** Zeeman transition is complete — the tarball ships `zssmis_f16..f20.TauCoeff.nc` and SSMIS runs load them with no fallback. **AMSU-A** ships no `zamsua*.TauCoeff` in either format; AMSU-A simply has no Zeeman coefficient and runs without one (its 60 GHz Zeeman correction is not applied — unchanged from prior releases). The probe is hardened (`ed21dc7`) so a sensor that lacks the `.nc` only forces the batch to Binary when it actually has a `.bin`; an AMSU-A with neither format no longer drags a NetCDF-only SSMIS set down to a (nonexistent) Binary and no longer emits a spurious "incomplete; falling back" message. So commit `0fcb7d4`'s "complete NetCDF transition" is accurate **for SSMIS specifically**, not for the entire Zeeman family.
* New offline `ODSSUBIN2NC` converter and `ODSSU_netCDF_IO` / `ODZeeman_netCDF_IO` modules so the SSU/Zeeman TauCoeff paths exist in NetCDF form.
* The bytes the tests load are now those in `fix_REL-3.2.0.0.tgz`.

The intent is `.bin`↔`.nc` round-trip equivalence; for the bulk of the sensor sweep that holds and TB is unchanged. The non-trivial consequences are items 2–4 below.

## 2. SpcCoeff NetCDF reader: `NLTECoeff` / `ACCoeff` sibling load

Commits: `98c6815` ("fix silent NLTECoeff sibling-load truncation; prefer canonical coeff dirs"), `3ac90dc`, `c1e00de` (stage `cris-fsr_n21.NLTECoeff.nc` for the suite).

Files: `src/Coefficients/SpcCoeff/SpcCoeff_netCDF_IO.f90`, `src/Coefficients/CRTM_SpcCoeff.f90`.

* The binary `SpcCoeff` reader streams the antenna-correction (`ACCoeff`) and non-LTE-correction (`NLTECoeff`) sub-structures inline from the same file via a `DATA_PRESENT` flag. The NetCDF layout stores them as separate sibling files: `<sensor>.ACCoeff.nc`, `<sensor>.NLTECoeff.nc`. The previous NetCDF `SpcCoeff` reader did not read them — the loaded `SpcCoeff` had empty `AC` / `NLTE` sub-structures (the path string used for the sibling-file existence check was being truncated).
* The new reader locates the siblings in the canonical REL-3.2 layout (`fix/SpcCoeff/netCDF/` ↔ `fix/ACCoeff/netCDF/` ↔ `fix/NLTECoeff/netCDF/`) or, for flat layouts, next to the `SpcCoeff` file, with an oversized path buffer so the `File_Exists` check is against the full path. It validates `Sensor_Id` / `WMO_Satellite_Id` / `WMO_Sensor_Id` / `Sensor_Channel` consistency across the trio.

**TB effect:** for sensors with an `NLTECoeff` sibling — in the regression suite that's CrIS-FSR (`cris-fsr_n21`, `cris399_npp`) — the non-LTE radiance correction is now applied (it was effectively disabled on the NetCDF path before). Radiance and BT change in the NLTE-sensitive shortwave-CO₂ channels. The `forward`, `tangent_linear`, `k_matrix`, and `adjoint` tests for those sensors all reflect this.

## 3. `v.abi_g18` reflective bands — per-channel `Solar_Irradiance`

* The reflective-band solar source is `RTSolution%Solar_Irradiance = SC%Solar_Irradiance(ch) * GeometryInfo%AU_ratio2`, with `AU_ratio2` channel-independent.
* Running current code against `v.abi_g18.SpcCoeff.nc` yields `Solar_Irradiance = 1.035050` for all six ABI reflective channels (the per-channel `SC%Solar_Irradiance` it reads × the standard AU factor).
* Prior `v.abi_g18` output had a per-channel-varying value (≈1.0353, 1.0347, 1.0348, 1.0351, 1.0354, 1.0362), i.e. a different per-channel `SC%Solar_Irradiance` set was in effect.
* The per-channel solar-source delta is ~0.03–0.1% and propagates into `SOD`, `Layer_Optical_Depth`, `Single_Scatter_Albedo`, and the reflective-band `Radiance` / `Brightness_Temperature`; the BT change is ~0.006–0.07 K.

Affected common ctests: the `v.abi_g18` `forward` (×8), `k_matrix` (×7), `adjoint` (×2), and `tangent_linear` (×2) cases. (`v.abi_gr` is unaffected.)

## 4. WMO satellite/sensor ID values in the coefficient files

Carried by the canonical coeff files in `fix_REL-3.2.0.0.tgz` (see also `98c6815`):

* `abi_g18`: `WMO_Satellite_Id` / `WMO_Sensor_Id` are now the real WMO values (`272` / `617`); previously they were `−999` placeholder values. (CRTM's canonical "no id" sentinels are `1023` for satellite, `2047` for sensor.)
* `cris-fsr_n21`: `WMO_Satellite_Id` is `226` (NOAA-21); previously `224`.

No TB change. The regression comparison includes these `RTSolution` fields, so the `abi_g18` and `cris-fsr_n21` cases differ on the metadata alone (for `abi_g18` that's the *only* difference; for `cris-fsr_n21` it's in addition to the NLTE change in §2).

## 5. Analytic MW-land emissivity Jacobians (issue #281)

Commits: `9358caa` (LAI/vegetation analytic TL/AD, Phase 1), Phase 2 (soil moisture), `aa31bb7` (soil/land temperature, Phase 3). The Jacobian test is registered as `test_Unit_Land_Jacobian`, source `test_Land_Jacobian.f90`.

Files: `src/SfcOptics/CRTM_MW_Land_SfcOptics.f90`, `src/SfcOptics/NESDIS_Emissivity/NESDIS_LandEM_Module.f90`, `src/SfcOptics/CRTM_SfcOptics.f90` (3 land dispatcher call-sites), `src/SfcOptics/CRTM_SfcOptics_Define.f90` (`iVar%MWLSOV`), `test/mains/unit/Unit_Test/test_Land_Jacobian.f90`, `docs/design/surface_jacobians_281.md`.

* Previously `Compute_MW_{Land,Snow,Ice}_SfcOptics_TL/_AD` were pure zero-stubs. This change gives the **MW-land** path (NESDIS_LandEM, < 80 GHz) analytic TL/AD for **all** of its physical state variables by hand-differentiating `NESDIS_LandEM`: the canopy `vlai = LAI*Veg_Fraction` optical-depth path, the soil-moisture and soil-temperature soil dielectric, the Fresnel/roughness chain (factored into `Roughened_R23_Deriv`, shared by the moisture and temperature paths), and the canopy/soil thermal-ratio `gsect0` (soil and land temperature). Partials are cached in `iVar%MWLSOV`. Snow/ice remain zero-stubs.
* **Land/soil temperature (Phase 3, `aa31bb7`):** soil temperature enters via the soil dielectric and `gsect0`; land (skin) temperature via `gsect0`. The soil-temperature out-of-range aliasing (`t_soil <- t_skin`) is resolved in the forward — the aliased input gets a zero derivative and its sensitivity re-attributes to land temperature. The land-temperature emissivity part **accumulates** onto the existing skin-T Planck emission Jacobian (`CRTM_Compute_SurfaceT_AD`). This corrects a latent bug: the forward already depended on skin temperature through `gsect0`, but `Surface_K%Land_Temperature` dropped it, so the below-cutoff land-temperature Jacobian was ~3-4x too large; it now matches finite differences. `Canopy_Water_Content` is never consumed by the forward → its analytic Jacobian is exactly zero (asserted in the test).
* **Forward emissivity is bit-identical** — the new derivative code is gated behind `PRESENT(...)` optional arguments that the forward never supplies. So radiances/BT do not change.

**TB effect:** none in the common ocean regression sweep. The change is to the **Surface-Jacobian reference data** for MW-over-land scenes: the `LAI` / `Vegetation_Fraction` / `Soil_Moisture_Content` / `Soil_Temperature` / `Land_Temperature` columns of `Surface_K` (and the matching TL/AD outputs) change (from zero, except land temperature which had the emission-only value). A field-level diff confirms only those columns move; `Atmosphere_K` and `RTSolution_K` are unchanged. MW-over-land Surface reference files must be regenerated after this change (build-local, self-seeding on a fresh checkout).

## 6. CONST_MIXED_POLARIZATION (polarization type 13) Distance_Ratio fix

Commits: `bdc7fb9` (fix), `57c9911` (`test_CONST_MIXED_Polarization` unit test).

File: `src/SfcOptics/CRTM_SfcOptics.f90`.

* The CONST_MIXED_POLARIZATION (pol type 13) FWD/TL/AD branches were scaling the fixed polarization angle's `SIN2_Angle` term by `GeometryInfo%Distance_Ratio`; that scaling was erroneous and is dropped. The V/H-mixed cases are untouched.

**TB effect:** pol type 13 is used only by the TMS (TROPICS / tomorrow.io) family, so no common-suite sensor is affected. TMS/TROPICS BT changes. (Gotcha noted during the fix: SfcOptics pol-mixing requires `%n_Stokes == 1` while the allocation uses `MAX_N_STOKES`.)

## 7. Where to look for a given regression difference

1. Difference is only in `WMO_*` / `Sensor_Id` → §4 (coeff-file metadata).
2. CrIS-FSR `Radiance` / `Brightness_Temperature` change → §2 (NLTECoeff sibling now loaded).
3. `v.abi_g18` `SOD` / `Layer_Optical_Depth` / `Single_Scatter_Albedo` / reflective-band `Radiance` / `Brightness_Temperature` → §3 (per-channel `Solar_Irradiance`).
4. MW-over-land `Surface_K` / Surface TL/AD change in the `LAI` / `Vegetation_Fraction` / `Soil_Moisture_Content` / `Soil_Temperature` / `Land_Temperature` columns (forward BT unchanged) → §5 (analytic MW-land Jacobians, #281).
5. `tms_*` / TROPICS pol-13 `Brightness_Temperature` change → §6 (CONST_MIXED_POLARIZATION Distance_Ratio fix).
6. Difference depends on `OMP_NUM_THREADS` → not expected; that would be a bug, not one of these changes.
7. MW-water `Radiance` / `Brightness_Temperature` / Jacobian change on a sensor with channels ≥ 200 GHz (e.g. `mwr_aws`, TROPICS/`tms_*`) → §8 (PARMIO backend — now auto-loaded and auto-dispatched at ≥ 200 GHz; no common-suite sensor reaches that frequency).
8. A new `Down_Radiance` / `Downwelling_Radiance(:)` / `Upwelling_Radiance(:)` field appears in the reference, or the `RTSolution` comparison gained columns → §10 (new downwelling/upwelling profile outputs; existing `Radiance`/BT unchanged).
9. MW-over-**land** forward emissivity/BT changed wholesale (not just the Jacobian columns) → §8 (TELSEM2 atlas auto-loaded because its file is present in the coeff path; #314).
10. MW cloud-ice (DDA-ARTS) `Radiance` / Jacobian change → §8 (`CRTM-Exp`/DDA ICE_CLOUD now scatters + habit default change; #320). Mie-TAMU default LUT is unaffected.
11. `n_Stokes > 1` (vector-RT) result change → §8 / #318. The scalar (`n_Stokes = 1`) suite is bit-identical.

## 8. Changes that do not affect the standard regression scenes

These are on the branch but produce **no TB difference** in the common `ctest` configuration. Two of them (the NESDIS guards) *do* change code behavior, but only on an error/edge path the standard scenes never hit — flagged explicitly below.

* **OpenMP thread-safety / race fixes** — `6ae8c1c`, `07e91d7`, `ec1cbb1`, `8311ba3`, `fc0a49a`, `01edea6`, `288fbdb`, `53a266d`, `b94b23f`: removed unsafe `SAVE` / implicit-`SAVE` coeff scratch (NESDIS-emissivity, ODCAPS); fixed channel-thread `!$OMP` races in `CRTM_Forward/Tangent_Linear/K_Matrix` (`Error_Status` write → `REDUCTION(MAX:...)`; unindexed `RTV%`/`RTV_Clear%` → `RTV(nt)%`; an OOB chunk-bucket write and an `end_ch` OOB read; `AAvar` privatization; per-channel NLTE/Zeeman predictor reset); hardened `CRTM_ChannelInfo_Subset`.
  - *No numerical difference, by construction:* the removed `SAVE`s are on local scratch arrays that are unconditionally re-assigned from literal `data` / array-constructor values at the top of every call (vestigial `SAVE`); the race fixes only change anything with >1 OpenMP thread, and the regression suite runs single-threaded (`CRTM_Init` coerces unset/empty `OMP_NUM_THREADS` to 1). The ctest pass/fail set was verified byte-identical before/after this work.
  - New self-consistency tests added: `test_OMP_Consistency`, `test_OMP_Speedup`, `test_ChannelSubset_OMP`, `test_OMPoverChannels` (no shared reference files). README gained an "OpenMP and thread safety" section. (`JCSDA/CRTMv3#111`, `#164`.)
* **PARMIO microwave ocean-emissivity backend (a default-path change for channels ≥ 200 GHz)** — `2c2f8d4`, `28fd40f`, `9fdf70a`, `11b9bfb`, `629b6a5`, `9cebd68`, `34fbbeb`, `7689efe`, `776aa56`, `9e8702f`, `0c6ff86`, `13c7d23`, `6df91a7`, plus `src/SfcOptics/MW_Water/PARMIO_MWSSEM/*`, `src/Coefficients/.../PARMIOCoeff/*`, `src/Coefficients/CRTM_PARMIOCoeff.f90`, `test/.../parmio_tlad/*`: a LUT-driven MW ocean SSEM. **Correction to earlier drafts of this issue: PARMIO is no longer opt-in.** The `Use_PARMIO_Model` flag (both `Options%` and `SfcOptics%`) was removed (`0c6ff86`); the backend is now auto-loaded at init and auto-dispatched by channel frequency.
  - *Auto-load (`13c7d23`):* `CRTM_Init` resolves `<File_Path>/PARMIO.MWwater.EmisCoeff.nc` and loads it whenever the loaded SpcCoeff set contains at least one microwave sensor (`CRTM_LifeCycle.f90` ~1073, ~1101–1137). The caller may override the filename via the optional `PARMIOCoeff_File` argument. Missing-LUT behavior: with no explicit file, an absent LUT is non-fatal and CRTM silently continues on the FASTEM path (drop-in); an explicitly-supplied-but-absent `PARMIOCoeff_File` is a hard `FAILURE`. The LUT *is* shipped in `fix_REL-3.2.0.0.tgz` (`fix/EmisCoeff/MW_Water/netCDF/PARMIO.MWwater.EmisCoeff.nc`), so in the default deployment the LUT loads and the routing below is active.
  - *Dispatch (`7689efe`):* `CRTM_MW_Water_SfcOptics` routes a channel to PARMIO iff `CRTM_PARMIOCoeff_IsLoaded() .AND. Frequency >= PARMIO_FREQ_THRESHOLD` (= 200 GHz) — forward `:240`, TL `:471`, AD `:692`. Below 200 GHz, or with the LUT not loaded, the code is byte-identical to the original FASTEM path. So for MW-water channels ≥ 200 GHz the default ocean emissivity — and thus `Radiance` / `Brightness_Temperature` and the MW-water Jacobians — now comes from PARMIO rather than FASTEM. This is a genuine default-behavior change, not an inert opt-in.
  - *Why the common ctest suite is unaffected:* no sensor in the regression suite reaches 200 GHz — ATMS 183.31, GMI 183.25, SSMIS 183.31, MHS 190.31, AMSU-A 89, SAPHIR 183.31, AMSR 89; the committed `Simple`/`ClearSky` sweep uses `amsua_metop-a mhs_n18 ssmis_f16 amsre_aqua` (+ `atms_npp` in `check_crtm`). The 200 GHz gate deliberately excludes the ATMS/GMI/SSMIS/MHS 183–190 GHz band (`7689efe`: PARMIO gave no skill there and degraded ~88 GHz, so it was scoped to ≥ 200 GHz). The only MW sensors that cross 200 GHz are `mwr_aws` (~325 GHz) and the TROPICS/`tms_*` family (~204 GHz); their default output now reflects PARMIO. None is in the common-suite reference data, so no `forward` / `tangent_linear` / `k_matrix` / `adjoint` reference TB changed.
  - *Test coverage:* the ≥ 200 GHz PARMIO path is exercised by self-consistency drivers — `test_PARMIO_TLAD` (two-sided finite-difference TL + adjoint dot-product) and `test_PARMIO_FASTEM_DeltaSweep[_AWS|_TMS]` (up to 325 GHz; the TROPICS 204.8 GHz window channel is the stronger PARMIO-vs-FASTEM witness, `91f5fee`) — and, since `aff03b5` (issue #311), by a stored-reference `mwr_aws` `ClearSky` regression across `forward`/`tangent_linear`/`adjoint`/`k_matrix` (19 channels, 4 at ~325 GHz), gated on `AWS_COEFFS_PRESENT AND PARMIO_LUT_PRESENT`. The gap flagged in earlier drafts of this document (no truth-file comparison for the ≥ 200 GHz default path) is closed.
  - *Thread-safety:* the PARMIO compute path (`CRTM_PARMIO.f90`, `_TL`, `_AD`) holds no writable module `SAVE` state — only the read-only LUT (`PARMIOC`), mirroring FASTEM's `MWwaterC` — so it is safe under the OpenMP-over-channels parallelism. (The bulk of the `CRTM_LifeCycle.f90` diff remains the Binary→NetCDF default-format flip — item 1, not this.)
* **nvfortran support** — `6bf5757`: new `cmake/compiler_flags_NVHPC_Fortran.cmake`; split the rank-8 PARMIOCoeff `Rdown` LUT into per-polarization rank-7 `Rdown_v` / `Rdown_h` (nvfortran caps array rank at 7).
  - *No numerical difference:* the on-disk file already stores `Rdown` per-polarization, so this is a memory-layout reshape with no value change; touches only PARMIO code plus test files; the `REAL(16)` → `REAL(fp)` edit is in two convergence *unit tests* (tolerance 0.1), not the library.
* **NESDIS ATMS snow / sea-ice emissivity guards** — `3a36f5f`, `240520b` (companion to `JCSDA/CRTMv3#192`): the diagnosis-based emissivity routine (`ATMS_SNOW_ByTBTs_D` / `ATMS_SeaICE_ByTbTs_D`) now runs only when the five window-channel TBs are `PRESENT`, `SIZE >= 5`, and all finite and within `[50, 500]` K; otherwise the default/by-type emissivity is kept. Previously it was called unconditionally (and the snow path read out of bounds when `Tbs` was absent or shorter than 5 — the #192 crash).
  - *This is a real behavior change, but confined to the error/edge path:* for the standard regression scenes (valid, in-range ATMS/AMSU window-channel Tbs) the diagnosis path runs exactly as before ⇒ identical TB. It only diverges for malformed / out-of-range / missing Tb inputs, which the standard tests don't produce. So: no observed effect in the suite, not "unconditionally identical."
* **Argument-interface / hygiene** — `b94b23f` (`FitCoeff_*_Create` assumed-shape `dimensions` arg, re-applying the #192 fix correctly) and the ODCAPS `ODCAPS_AtmAbsorption.f90` / `ODCAPS_Predictor.f90` edits: thread-safety / argument-shape cleanup, no numeric change.
* **Lifecycle wiring** — `9cebd68`, `34fbbeb`: PARMIO obs-space drivers moved onto the integrated `CRTM_Init` lifecycle; the default RT path is untouched.
* **Repo cleanup / version bump** — removed `CRTM_V30_TEST/`, `README_JEDI.md`, `Set_CRTM_Environment.sh`, `NOTES`, the deprecated `*_NC` unit-test variants; dropped dead `Zeeman_Utility.f90` from the lib build (kept for the offline `BeCoeff_ASC2NC` tool); `LICENSE` / `VERSION.cmake` / `CRTM_Version.inc` → v3.2.0; `README.md` refreshed; per-compiler flag-file updates (GNU/Intel/IntelLLVM/Cray/XL/NVHPC).
* **`n_Stokes > 1` vector-RT (polarized scattering) fixes** (`JCSDA/CRTMv3#318`) — `c95ac47`, `56037ca`, `3c65c8b`, `6ec846c`, `4c800db`, `81c4b16`, `3d3ce8f`, `5dd13bf`, `94760ef`, plus the Round-1 `Normalize_Phase` TL/AD mirror (`fe5104c`): corrected the ADA scattering-layer thermal source (intensity-slot-only guard `MOD(i-1,n_Stokes)==0`), the Kirchhoff column-sum bound (`n_Streams → n_Streams*n_Stokes`), the satellite intensity-row special case, the polarized phase-block normalization, the BT adjoint-seed routing to `Stokes(1)`, the K-matrix `SfcOptics%n_Stokes`-from-`Opt` sync, and propagated all of it to TL/AD/K. These were the ~30–44× cloudy-radiance inflation (#318) and its Jacobian consistency.
  - *No common-suite TB change, by construction:* every change reduces **exactly** to the prior scalar code at `n_Stokes = 1`, and the entire regression suite runs scalar (`Options%n_Stokes` defaults to 1). The `n_Stokes > 1` path is reachable only with a ≥ 6-phase-element cloud LUT (the `CRTM-Exp` scheme below); stock LUTs are hard-rejected by the forward guard.
  - *New coverage:* `test_VectorRT_TLADK` (Round-1, `1cc9a5b`) — TL-vs-FD on both Stokes components, full-Stokes adjoint dot-product (1e-12 tolerance, fault-injection calibrated), K-vs-AD, and an `n_Stokes=1` scalar control; gated on the `CloudCoeff_Exp_Full6.nc` LUT from #320. This is the first in-repo Jacobian coverage of the path. Remaining deferred items (fractional-cloud `n_Stokes>1` adjoint combine; surface V/H↔I/Q decoupled polarization) are tracked in #318.
* **TELSEM2 microwave land-emissivity atlas** (`JCSDA/CRTMv3#314`) — `3697a04`, `922662b`, `3d3c6d3`, `c3ac530`, plus `src/Coefficients/.../MW_Land/TELSEM2/*`, `src/Coefficients/CRTM_MWlandCoeff.f90`, integration in `CRTM_MW_Land_SfcOptics.f90` / `CRTM_LifeCycle.f90`: an optional climatological MW land-emissivity atlas (lat/lon/month), ported from RTTOV.
  - *Why the common suite is unaffected:* the standard sweep is ocean, and the ctest harness deliberately stages the atlas under a **non-default name** so the auto-load does not fire (`test/CMakeLists.txt`). The atlas-derived emissivity is treated as a constant (zero TL/AD), internally consistent since it depends only on lat/lon/month.
  - *⚠️ Behavior caveat (default-path, land users):* `CRTM_Init` **auto-loads** `<File_Path>/TELSEM2.MWland.EmisCoeff.nc` if present, and the file **ships in `fix_REL-3.2.0.0.tgz`** (`fix/EmisCoeff/MW_Land/netCDF/`). A user who flattens the fix tree into one coefficient directory will silently (a) switch all MW-land forward emissivity from `NESDIS_LandEM` to TELSEM2, and (b) lose the new #281 analytic land Jacobians (the atlas path leaves TL/AD zero). This is a drop-in policy mirroring PARMIO; it needs release-note prominence or an explicit opt-in — tracked in #314.
* **Experimental `CRTM-Exp` cloud optics + DDA-ARTS ICE_CLOUD change** (`JCSDA/CRTMv3#320`) — `7578d69`, `d04b001`, `bcb9ed4`, `7765d3d`, `f8d8dc9`, `e1eeea0`, plus `CloudCoeff_Exp_{Define,netCDF_IO}.f90`, the `CRTM_CloudScatter.f90` scheme gating, `CRTM_Parameters.f90` (`MAX_N_LEGENDRE_TERMS` 16→64): an opt-in 6-phase-element ('full-Mueller') MW cloud-optics scheme (`Cloud_Model='CRTM-Exp'`), plus a `Data_Type` discriminator (Mie-TAMU vs DDA-ARTS) on `CloudCoeff`.
  - *Default Mie-TAMU path bit-identical:* `Data_Type` is derived at load from exactly the `ALL(Reff_MW>0)` predicate `develop` evaluated per call, never read from file, so stock `.bin`/`.nc` coefficients are fully backward-compatible; the `CRTM-Exp` scheme is reachable only by the exact `Cloud_Model` string and fails loudly on a mismatched file. `MAX_N_LEGENDRE_TERMS` 16→64 is memory-only (loops bounded by the actual term count).
  - *⚠️ Behavior caveat (DDA-ARTS users):* for the DDA-ARTS cloud database, ICE_CLOUD now goes through the full scattering branch (was a non-scattering shortcut) **and** the default ICE_CLOUD habit changed `IceSphere(18) → IconCloudIce(6)` (`d04b001`). Both are intended (AWS 325 GHz O−B improvement) but silently change radiances/Jacobians for existing DDA users; no regression pins the new values. Needs release-note coverage. The common Mie-TAMU suite is unaffected.
* **Grazing-angle MW-water reflectivity guards** — `e662b02` (FastemX), `1c4d4ce` (PARMIO), test `6898e2c` (`test_Grazing_SfcOptics`): clamp the MW-water reflection-correction to a physical range at the near-grazing Gaussian quadrature angles the scattering RT uses (without the guard, the FASTEM-fit reflection correction extrapolates to ~1e35 and produced −1e15 K TBs at ≥ 200 GHz scattering channels). The clamp fires only above ~84°; standard scenes never reach those angles, so no common-suite TB change. (Round-2 `4c80915` added 85° TL/AD/K cases to `test_Downwelling_TLADK`.)

## 9. Regression baselines converted from binary to netCDF (TB-neutral)

The `ctest` regression **baselines/results** (`RTSolution{,_K,_AD,_TL}`,
`Atmosphere`, `Surface`) were switched from binary (`.bin`) to netCDF
(`.nc`). This is a **test-infrastructure change only** — it changes the
on-disk format of the self-seeded reference files in `build/test/results`,
not the radiative transfer, so it produces **no TB difference** and is
independent of items 1–6. The reference files self-seed on first run, so
nothing is committed; the conversion is a hard switch (drivers write/read
only `.nc`).

New / changed code:

* `CRTM_Surface_Define.f90`, `CRTM_Atmosphere_Define.f90`: netCDF
  read/write/inquire added for the rank-2 (`n_Channels × n_Profiles`,
  K-matrix) and rank-1 (profile-only, adjoint) objects. Each element is
  flattened into a packed `REAL(fp)` record (Surface: one var
  `Surface_Data(n_Channels,n_Profiles,n_Fields)`; Atmosphere: a
  variable-length record covering the nested `Cloud(:)`/`Aerosol(:)`),
  mirroring the existing `CRTM_RTSolution_Define` netCDF idiom. Profile-only
  files store the true `n_Channels` (`0`) as a global attribute with the
  channel dimension `MAX(n_Channels,1)`.
* `CRTM_RTSolution_Define.f90`: pre-existing netCDF gaps fixed so the
  K-matrix/adjoint objects round-trip — `n_Layers=0` (`RTSolution_K`/`_AD`
  carry only the scalar adjoint seed), the optional-argument segfault in
  `CRTM_RTSolution_InquireFile`, the missing `Reflectance`/`Reflectance_clear`
  reads, and per-element `RT_Algorithm_Name` (it varies by channel/profile
  for scattering sensors). The forward RTSolution path was already netCDF.
* New round-trip unit tests `test_Surface_netCDF_io` and
  `test_Atmosphere_netCDF_io` (ctest count 206 → 208).
* All `forward`/`k_matrix`/`adjoint`/`tangent_linear`/`Aerosol_Bypass`
  drivers flipped to `NetCDF=.TRUE.` + `.nc` baseline names.

## 10. Downwelling / upwelling radiance-profile outputs (new fields, TB-neutral for existing ones)

Commits: `64b17db`..`c3ac530` (the "downwelling"/"upwelling" series, Phases 1–5).

A new family of `RTSolution` outputs, fully differentiated (TL/AD/K) across all
three solvers (Emission/SOI/ADA):

* `Down_Radiance` — surface downwelling radiance (scalar).
* `Downwelling_Radiance(:)` — level-resolved downwelling profile (surface→TOA).
* `Upwelling_Radiance(:)` — level-resolved upwelling profile.

Opt-in via the new `Options%Compute_Down_Radiance`,
`Compute_Down_Radiance_Profile`, `Compute_Up_Radiance_Profile` flags (all default
`.FALSE.`); the clear-sky (emission) surface `Down_Radiance` is always computed.
The fractional-cloud (TCC) combine of the profiles mirrors the `Radiance` combine
in FWD/TL/AD/K. The legacy `Obs_4_downward_P` aircraft-observer hack was retired
(Phase 4) with prior aircraft behavior preserved.

**TB effect:** none on existing fields — the forward `Radiance`/`Brightness_Temperature`
are bit-identical. These are *additional* output fields; the always-on emission
`Down_Radiance` now populates in the (self-seeding) reference files, but no
committed truth file changes (the references self-seed on first run, per §9).
The scattering downwelling/profile outputs are off by default → bit-identical to
legacy when unused. Verified by `test_Downwelling_TLADK` (14 cases: TL-vs-FD,
adjoint dot-product, K-vs-AD, and the surface-profile==scalar identity, across
clear/SOI/ADA × overcast/fractional, plus Round-2 grazing-angle cases).

## 11. Pre-release review fixes

A full validity/completeness review of the branch (2026-06-11) produced two
rounds of pre-tag fixes to branch-only code — see the dedicated comments below
("Pre-release fix thread"). An earlier review pass (2026-05-29/30), two
post-doc commits, and the 2026-07-18 comprehensive review are listed here as
well. **None changes any common-suite TB** except the NLTE-staging repair in
the 2026-07-18 round (13 cris399_npp/iasi_metop-b baselines regenerated with
NLTE correctly active; suite 215/215 after). Summary:

* **Earlier review pass, 2026-05-29/30** (`bdc7fb9`+`57c9911`, `9ca8e98`,
  `aff03b5`, `a84de47`, `91f5fee`, `c89b2e4`, `8bd411c`, `ed21dc7`):
  CONST_MIXED_POLARIZATION `Distance_Ratio` fix + unit test (§6); repo hygiene
  for the tag; the `mwr_aws` ≥ 200 GHz stored-reference regression (§8, #311);
  repair of the PARMIO delta-sweep that compared PARMIO against itself; the
  TROPICS 204.8 GHz delta-sweep witness; documentation of the shared read-only
  `Predictor` OMP invariant (H1); PARMIO reflection-correction gating unified
  across FWD/TL/AD (H2); Zeeman format-probe hardening (§1).
* **Round 1** (`b19cd61`, `9eda275`, `78aac3d`, `8f75a59`, `fe5104c`, test `1cc9a5b`):
  TL `Down_Radiance` OMP race; SOI `Compute_*` gating; MW-land soil-moisture
  Jacobian `+Inf` at SMC=0 (#281); `Resolve_Coeff_Format` extension authority;
  the `Normalize_Phase` `n_Stokes>1` TL/AD mirror (#318); `test_VectorRT_TLADK`.
* **Round 2** (`f7000b6`, `8f39e93`, `4be36b1`, `6b188e7`, `4c80915`): K-matrix
  Exp-scheme Legendre hook; `Options` `Compute_*` flag plumbing
  (`SetValue`/`Inspect`/`Equal`); rank-1 netCDF `n_Channels/=0` parity guard;
  TELSEM2 longitude-wrap + `class1` validation; grazing-angle + ≥200 GHz PARMIO
  test cases + `OMP_Speedup` `RUN_SERIAL`.
* **Post-doc commits** (`518dd41`, `ef8e7c1`): `ASYMTX` balance-loop
  non-termination guard in `CRTM_Utility` (robustness only — no result change
  for previously-converging matrices); per-habit `Reff`→`Dm` conversion in the
  `CRTM-Exp` MW CloudCoeff reader (opt-in Exp LUT path only; default conversion
  factor 1.0, so existing Exp LUTs and all non-Exp paths are unaffected).
* **Comprehensive review, 2026-07-18** (13-area adversarially-verified pass over
  the full branch diff). Fixes, none of which changes any common-suite TB
  except where noted:
  - **Multi-sensor single-call `ln` offset** (Forward/TL/K): per-thread channel
    counter dropped the previous sensors' cumulative offset, so a single call
    with `ChannelInfo(1:2+)` overwrote sensor 1's outputs (preexisting on
    develop; invisible to the suite, which passes one sensor per call). Fixed
    in all three modules + new `test_Unit_MultiSensor_SingleCall` (FWD/TL/K
    combined-vs-per-sensor, bit-exact; fails on the unfixed code). Also
    guarded per-sensor `RTV_Create` re-allocation (the `error in allocate
    Pff 5014` noise on sensors ≥ 2).
  - **Split-path coefficient loading**: SpcCoeff/TauCoeff were loaded from
    `Effective_NC_Path` unconditionally, breaking `File_Path` (Binary) +
    `NC_File_Path` (netCDF) split trees; now selected per requested format.
    PARMIO/TELSEM2 drop-in auto-load additionally probes `NC_File_Path`, and
    CRTM_Init reports the ≥ 200 GHz FASTEM fallback when no PARMIO LUT loads.
  - **PARMIO physics** (≥ 200 GHz MW-water only): (a) azimuthal harmonics were
    applied to the invalid-azimuth sentinel (`Sensor_Azimuth_Angle` default
    999.9 — the common DA configuration); now gated like FastemX (azimuthal
    mean, FWD/TL/AD-consistent, sentinel case pinned in `test_PARMIO_TLAD`).
    (b) 3rd/4th-Stokes reflectivity was `1−e` (≈1) instead of the FastemX
    convention 0 — polarimetric (`n_Stokes=4`) runs would have reflected ~100%
    of downwelling U/V. Both change TBs only for ≥ 200 GHz MW-water scenes
    without a valid sensor azimuth (the `mwr_aws` references, which set a
    valid azimuth, are unchanged). Reader now also loads the group-boundary
    threshold attributes, and the forward LUT interpolation gained the
    `Is_Allocated` guard the TL/AD already had.
  - **VMOM sensor-angle normalization** (`Normalize_Phase`, `n_Stokes>1`,
    `n_Streams<nZ`): normalized column `nZ*n_Stokes` (the last polarized
    component) instead of the sensor-angle intensity element
    `(nZ-1)*n_Stokes+1`; scalar branch and all other blocks confirm the
    convention. Vector-RT (opt-in) results only.
  - **NLTE/ACCoeff staging** (test data, TB-affecting for 13 tests): the
    `crtm_test_input` entries for `cris399_npp`/`iasi_metop-b`/`airs_aqua`
    NLTECoeff and `amsua_n19` ACCoeff pointed into `SpcCoeff/netCDF/` where
    the release tree does not keep them → dangling symlinks → NLTE silently
    off in the stored baselines. Paths corrected to `NLTECoeff|ACCoeff/netCDF/`
    and the 13 affected references regenerated (CrIS 4.3 µm channels move by
    up to ~8.5 K on the daytime profile — the correct, NLTE-on radiances).
    The SpcCoeff reader now emits an INFORMATION message when a sibling
    AC/NLTE file is absent so this cannot recur silently.
  - **Exp scheme + non-MW sensors**: the CRTM-Exp branch stomped
    `n_Legendre_Terms=0` for every sensor, which also truncated the AEROSOL
    phase expansion for IR/VIS channels in cloudy profiles; the override is
    now MW-only (non-MW keeps the stream-based order; Exp cloud optics remain
    zero for non-MW as designed).
  - **Zeeman mixed-format visibility**: a Zeeman candidate whose coefficient
    exists only in the format the batch probe did not choose now gets a loud
    WARNING (was: silent loss of the Zeeman correction).
  - **Hardening/minor**: TELSEM2 `MODULO` 360.0 longitude edge + netCDF-reader
    index/geometry validation; ODSSU netCDF WMO fallbacks use the CRTM invalid
    sentinels (1023/2047) and OPTRAN dimensions are conformance-checked;
    Surface netCDF reader validates `n_Surface_Fields`; Atmosphere netCDF
    writer rejects zero-size arrays; level-profile FWD/TL/AD copies clamp to
    conformant extents; FD probes in `test_VectorRT_TLADK` /
    `test_Downwelling_TLADK` / `test_MW_O3_TLAD` now check forward
    `Error_Status`; `find_package(OpenMP REQUIRED)`; `NDEBUG` no longer leaks
    into `DEBUG` builds (case-insensitive build-type match);
    `Convert_TELSEM2_Atlas` prints usage instead of defaulting to a personal
    path; stale testinput header comment rewritten (all-netCDF staging; the
    tarball ships no `.bin` coefficients, so the binary coefficient read path
    is suite-untested for this release).
  - **Reviewed, deliberately unchanged**: `CRTM_Init` forcing 1 OpenMP thread
    when `OMP_NUM_THREADS` is unset is the documented `bb9991b` design (note
    for release notes: `OMP_SET_NUM_THREADS` is process-global and affects the
    host application); the "legacy CloudScatter OOB under vector RT" review
    finding is refuted — all four entry modules reject cloudy `n_Stokes>1`
    runs with a `<6`-element LUT before that code can run.

## Appendix: commits `develop..HEAD` (oldest → newest)

> **Note:** the list below is the original snapshot through `57c9911`. The
> branch has since added the downwelling/upwelling (§10), `n_Stokes>1` (§8/#318),
> TELSEM2 (§8/#314), and `CRTM-Exp`/DDA-ICE (§8/#320) series plus the §11
> pre-release fixes; for the full current list use `git log --oneline develop..HEAD`.

```
fa64893  updating internal versions to v3.2.0 in preparation for REL-3.2.0
69edea8  Merge remote-tracking branch 'origin/develop' into feature/btj_REL-3.2.0
f3ee1c8  minor comment change
2c2f8d4  Add PARMIO microwave ocean emissivity backend
28fd40f  Wire PARMIO opt-in selector and obs-space regression test scaffolding
9fdf70a  Refresh PARMIO TLAD reference values for Meissner-fix LUT; add diagnostic probes
11b9bfb  Extend PARMIO validation: AWS TB sweep + 4-frequency V/H emissivity sweep
a7a6114  Merge branch 'develop' into feature/btj_REL-3.2.0
629b6a5  Merge branch 'feature/btj_ML_emissivity_from_parmio' into feature/btj_REL-3.2.0
43ea155  minor update to README.md regarding version number and history
d167bfa  removing old CRTM_V30_TEST
7771e35  removing Set_CRTM_Environment.sh, no longer used in modern build systems (use cmake)
dc27b60  removing README_JEDI.md -- CRTM is JEDI ready by default
40e2052  updated LICENSE version
6872b37  removed deprecated NOTES
01edea6  OpenMP: runtime OMP_NUM_THREADS, OPENMP=OFF support, and a speedup test
3ac90dc  NetCDF default LUT format + SpcCoeff sibling-substructure read
0fcb7d4  Complete NetCDF transition for Zeeman SSMIS TauCoeff
288fbdb  OpenMP: treat empty OMP_NUM_THREADS the same as unset
613921b  Unit_AerosolScatter tests: load GOCART-GEOS5 from NetCDF
8311ba3  CRTM_Forward: index RTV by thread in Obs_4_downward warning
ec1cbb1  CRTM_Forward/TL/K: fix three OpenMP-over-channels races
fc0a49a  CRTM_K_Matrix: restore channel-thread OpenMP (JCSDA/CRTMv3#231)
12b0394  Stage PARMIO LUT + AWS coeffs in canonical test_data; demote RC residual gate
69e752f  Switch coefficient extension from .nc4 to .nc for fix_REL-3.2.0.0
2a7495a  Add ODSSU netCDF I/O and BIN2NC converter for SSU TauCoeff
b6168df  PARMIO/AWS tests: read coefficients from canonical test_data only
5e0a69d  update MD5sum hash in test/CMakeLists.txt
c88c79c  updated Get_CRTM_Binary_Files.sh with correct md5sum for current fix_REL-3.2.0.0.tgz
6bf5757  Add nvfortran support: split 8-D Rdown into V/H 7-D arrays
9cebd68  Phase 4: wire PARMIO LUT into CRTM_Init lifecycle
34fbbeb  Lift PARMIO obs-space drivers to the integrated CRTM_Init lifecycle
98c6815  SpcCoeff netCDF reader: fix silent NLTECoeff sibling-load truncation; prefer canonical coeff dirs
c1e00de  test: stage cris-fsr_n21 NLTECoeff sibling for the regression suite
6ae8c1c  Fix thread-safety issues: remove unsafe SAVE attributes and pointer initializations
07e91d7  CRTM_Forward/TL/K: aggregate channel-thread error status via reduction
33a23da  updated MD5sums for netcdf tarball fix_REL-3.2.0.0.tgz to 5777242387228359869325e1a0505f85
5d0ca34  build: drop dead Zeeman_Utility.f90 from the libcrtm sources
ce8fcf0  docs: add "OpenMP and thread safety" section to README
7c9bdd9  test: add OpenMP/serial consistency regression test (JCSDA/CRTMv3#111)
b94b23f  FitCoeff_*_Create: assumed-shape `dimensions` arg (re-applies the #192 fix correctly)
3a36f5f  NESDIS_ATMS_SnowEM: guard the diagnosis-based path against absent/short Tbs
240520b  NESDIS_ATMS_SeaICE: only run the diagnosis-based path on sane TBs (parity with #192)
53a266d  Test channel subsetting under OpenMP; harden CRTM_ChannelInfo_Subset
7689efe  PARMIO: route MW-water channels >= 200 GHz to PARMIO LUT
776aa56  PARMIO: delete ATMS-only regression tests obviated by 200 GHz gate
9e8702f  PARMIO: rewrite remaining A/B tests to two-phase CRTM_PARMIOCoeff_Load
0c6ff86  PARMIO: remove inert Use_PARMIO_Model flag
13c7d23  PARMIO: auto-load LUT from default coefficient path in CRTM_Init
6df91a7  fixing more default binary options
8dc9bc3  docs: correct Read/Write netCDF-arg doc blocks to NETCDF default
ed20667  test: extend OMP consistency check to CRTM_Tangent_Linear
9358caa  feat(SfcOptics): analytic MW land emissivity Jacobians (LAI, vegetation, soil moisture)   [§5]
43661a9  updated md5sum for tarball (supersedes 33a23da; current md5 = 056d34c0fadfd67444e69907b013a30a)
bdc7fb9  fix(SfcOptics): drop Distance_Ratio scaling in CONST_MIXED_POLARIZATION   [§6]
57c9911  test: add CONST_MIXED_POLARIZATION surface-optics unit test
```

(`69edea8`, `a7a6114` are merges pulling `develop` forward.)

> **md5sum note (updated 2026-06-11):** the tarball was re-rolled again after the
> appendix snapshot — `43661a9`'s `056d34c0fadfd67444e69907b013a30a` was superseded
> by `99a1fa8` when the TELSEM2 MW-land atlas was added to the tarball. The
> **current** `fix_REL-3.2.0.0.tgz` md5 is `3dcef94c129efb78c85cdf542fca55ae`, which
> matches both `Get_CRTM_Binary_Files.sh` and `test/CMakeLists.txt`.
