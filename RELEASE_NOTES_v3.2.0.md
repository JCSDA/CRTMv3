# CRTM v3.2.0 Release Notes

**Release date:** July 31, 2026
**Coefficient data:** `fix_REL-3.2.0.0.tgz` (~3.5 GB), auto-downloaded at build
time (or via `Get_CRTM_Binary_Files.sh`).

Developer-facing detail for every item below (commits, affected tests, TB
impact) is in `REL-3.2.0_changes_vs_develop.md`. Changes are stated relative
to CRTM v3.1.4; the v3.1.4 tag and the `develop` branch point differ only in
the default coefficient-data location (no library code differences), so the
same catalog applies against either baseline. A per-sensor inventory of the
shipped coefficient files is in `REL-3.2.0_coefficient_inventory.md`.

## Highlights

- **netCDF coefficient transition.** netCDF is now the canonical coefficient
  format: SpcCoeff, TauCoeff (ODPS/ODAS/ODSSU/Zeeman), CloudCoeff,
  AerosolCoeff, and emissivity LUTs all load from `.nc` files, and the fix
  tarball ships netCDF only. Regression baselines are netCDF as well.
- **PARMIO MW-water emissivity (new, default-on at and above 200 GHz).** A
  LUT-driven physical-reference ocean emissivity backend for sub-mm sounders
  (AWS, TROPICS class). When `PARMIO.MWwater.EmisCoeff.nc` is present
  (`CRTM_Init` auto-loads it; it ships in the fix tarball), MW-water channels
  at or above 200 GHz route to PARMIO; below 200 GHz (every legacy sensor)
  the FASTEM path is byte-identical to v3.1.x.
- **TELSEM2 MW-land emissivity atlas (new, opt-in).** Enabled by the new
  `CRTM_Init` argument `Use_MWland_Atlas=.TRUE.` (which auto-resolves the
  default-named `TELSEM2.MWland.EmisCoeff.nc` from the coefficient path) or by
  passing an explicit `MWlandCoeff_File`. Without the opt-in, MW land
  emissivity uses NESDIS_LandEM exactly as before, even if the atlas file is
  present in the coefficient directory.
- **Level-resolved downwelling/upwelling radiance profiles (new outputs).**
  Opt-in via `Options%Compute_Down_Radiance_Profile` /
  `Compute_Up_Radiance_Profile`; fully differentiated (FWD/TL/AD/K) and
  combined correctly for fractional cloud. The surface downwelling radiance
  `RTSolution%Down_Radiance` is a first-class output (always populated on the
  emission path; enabled via `Options%Compute_Down_Radiance` for scattering).
- **Vector radiative transfer (`Options%n_Stokes > 1`) usability fixes**,
  including nonzero Jacobians, polarized phase-matrix normalization, and
  SfcOptics Stokes-dimension handling.
- **Analytic MW-land surface Jacobians** (issue #281) — the NESDIS_LandEM
  microwave land path (< 80 GHz) now returns analytic TL/AD/K sensitivities for
  **every physical land-state variable**: LAI, vegetation fraction, soil
  moisture, soil temperature, and land (skin) temperature. `Canopy_Water_Content`
  has an exactly-zero Jacobian (the forward never consumes it). This also
  corrects the `Land_Temperature` Jacobian below 80 GHz, which previously omitted
  the emissivity's skin-temperature dependence (through the LandEM `gsect0`
  thermal ratio) and was therefore too large; it now matches finite differences.
  Forward radiances are unchanged.
- **MW scene-ozone transmittance component** (`GROUP_MW_O3`, Group_Index=7)
  for microwave sensors.
- **UV scene-NO2 transmittance component** (`GROUP_UV_NO2`, Group_Index=8,
  issue #340) for UV-VIS air-quality spectrometers (TEMPO, GEMS class). A
  sixth ODPS component carries scene-variable NO2 absorption: supply an NO2
  profile in the Atmosphere (HITRAN id 10, ppmv) and NO2-sensitive radiances
  respond; omit it and the coefficient file's reference climatology applies.
  Full FWD/TL/AD/K support, verified by a machine-precision predictor-level
  transpose test and an end-to-end TL/AD/K parity test. Group 1/2/3/7
  coefficient files never reach the new code.
- **UV sensors can now run the forward operator** (issue #339). The
  surface-optics dispatch previously had no UV (Sensor_Type 4) branch, so
  every UV SpcCoeff (the shipped OMPS family included) failed CRTM_Forward
  with "Unrecognised sensor type". UV channels now share the VIS Lambertian
  surface-optics path, and a UV-only sensor list loads the VIS surface
  emissivity LUTs at CRTM_Init. MW/IR/VIS behavior is unchanged.
- **CRTM-Exp cloud-optics schema (experimental, opt-in).** A new
  habit-resolved cloud LUT format selected explicitly with
  `Cloud_Model='CRTM-Exp'`; the default cloud path is unchanged.
- **SNICAR visible snow emissivity scheme (new, opt-in).** A SNICAR-based
  VIS-snow reflectance LUT (`SNICAR.VISsnow.EmisCoeff.nc`, shipped in the fix
  tarball) selectable through the VISsnowCoeff scheme machinery, alongside
  updated IR snow emissivity modules. The default (NPOESS) snow surface path
  is unchanged.
- **ODPS transmittance-algorithm modernization** (issue #343). The ODPS group
  system was rebuilt on a single group registry with load-time validation of
  `Group_Index` and the `Component_ID`/`Absorber_ID` rosters (malformed or
  mislabeled coefficient files are now rejected at load with a clear message
  instead of computing garbage), per-component predictor kernels for FWD, TL,
  and AD, and file-roster-driven dispatch. Results are bit-identical for valid
  coefficient files; Zeeman-reserved group indexes are refused (the historical
  OMPS "Group 4" failure mode).
- **Long coefficient paths** (issue #238). Coefficient file paths are now
  carried in deferred-length strings instead of fixed 80/128/256-character
  buffers, so deep installation paths no longer truncate silently;
  initialization through a ~300-character path is regression-tested.
- **Fastem1 SST Jacobian corrected.** On the legacy Fastem1 MW-water path
  (`Options%Use_Old_MWSSEM=.TRUE.`, frequency >= 20 GHz) the emissivity's
  sea-surface-temperature derivative was silently dropped, so
  `Surface_K%Water_Temperature` carried only the skin-emission term. The
  Jacobian is now complete and validated against finite differences. The
  default (FastemX) path was never affected.
- **Runtime OpenMP control.** `OMP_NUM_THREADS` is honored at run time (no
  longer captured at configure time).
- **Expanded self-checking test coverage.** New baseline-independent checks
  include general TL-vs-FD and adjoint-consistency tests across the three
  main sensor types (#280), multi-sensor single-call bit-consistency, ODPS
  group-validation and long-path initialization tests, a DDA-ARTS ICE_CLOUD
  behavior pin, and OpenMP thread-count consistency tests (#111).

## Breaking and behavior changes

1. **RTSolution file formats changed incompatibly.** The netCDF reader
   requires variables absent from files written by earlier versions
   (per-element `RT_Algorithm_Name`, `Reflectance`, `Downwelling_Radiance`,
   the `n_Layers` global attribute), and the binary record grew. RTSolution
   files written by pre-3.2.0 code cannot be read; regenerate archived files.
2. **`Options%Obs_4_downward_P` removed** (compile-breaking). Migrate to
   `Options%Compute_Down_Radiance` / `Compute_Down_Radiance_Profile` and read
   `RTSolution%Down_Radiance` / `RTSolution%Downwelling_Radiance(:)`.
3. **Coefficient wrapper I/O defaults flipped Binary → netCDF**
   (`CloudCoeff_ReadFile`/`WriteFile`/`InquireFile` and the analogous wrapper
   modules). External callers reading `.bin` files through these routines must
   now pass `netCDF=.FALSE.` explicitly.
4. **`CRTM_ChannelInfo_Subset` hard-fails on duplicate or non-member channel
   lists** (previously silent misbehavior: stalled merges and silently
   deactivated channels).
5. **DDA-ARTS cloud optics: `ICE_CLOUD` now scatters** (the legacy
   non-scattering shortcut applies only to Mie-TAMU tables), and its default
   DDA habit changed from IceSphere to IconCloudIce. Users of DDA-ARTS
   CloudCoeff tables with `ICE_CLOUD` in their profiles will see different
   brightness temperatures (validated against a sub-mm sounder; tropical O−B
   at 325 GHz moved from ~+13 K to ~−0.6 K). Default (Mie-TAMU) cloud optics
   are bit-identical.
6. **PARMIO is presence-activated; TELSEM2 is opt-in.** Placing
   `PARMIO.MWwater.EmisCoeff.nc` in the coefficient directory switches MW
   water emissivity physics at and above 200 GHz; removing it restores FASTEM.
   `CRTM_Init` prints an INFORMATION message when a sensor with ≥ 200 GHz
   channels initializes without the PARMIO LUT. The TELSEM2 MW-land atlas, by
   contrast, is **never** loaded unless requested (`Use_MWland_Atlas=.TRUE.`
   or an explicit `MWlandCoeff_File`); a present-but-not-requested atlas file
   is ignored. Note that with TELSEM2 opted in, all land-parameter Jacobians
   (LAI, vegetation, soil moisture, soil/land temperature) are zero; the
   atlas is a climatology and does not depend on them; the analytic
   NESDIS_LandEM Jacobians apply only when the atlas is not loaded.
7. **OpenMP threading.** `CRTM_Init` reads `OMP_NUM_THREADS` at run time; if
   it is **unset or empty, CRTM defaults to a single thread** via
   `OMP_SET_NUM_THREADS(1)`. Because that call is process-global, it also
   affects OpenMP regions of the host application after `CRTM_Init` — export
   `OMP_NUM_THREADS` explicitly in threaded host applications (DA systems
   embedding libcrtm). The configure-time capture of `OMP_NUM_THREADS` (and
   the per-test ENVIRONMENT overrides) are gone: the environment at run time
   is what counts.
8. **Binary coefficient files are on the way out.** The v3.2.0 fix tarball
   ships no `.bin` coefficient files and the test suite no longer exercises
   the binary coefficient read path. The binary readers remain in the library
   for users with existing binary trees, but they should be considered
   deprecated (removal expected in a later v3.2.x, per the README).

## Known issues and limitations

- **Sub-mm thin frozen cloud (≈ 325 GHz):** optically thin frozen-cloud
  layers can produce nonphysical TBs through the adding-doubling/MOM path
  (small-τ matrix conditioning with high phase-function truncation orders).
  Affects sub-mm scattering scenes only; under investigation.
- **Fractional cloud with `n_Stokes > 1`:** the K/AD Jacobian seed uses the
  scalar-Radiance approximation; a full Stokes-space fractional combine
  remains future work.
- **ODSSU netCDF supports the ODPS sub-algorithm only** (the shipped SSU
  files are ODPS-based; ODAS-based SSU coefficient files remain binary-only).
- **Options binary I/O does not persist the new fields** (`n_Stokes` and the
  radiance-profile switches); they will be added with the next format
  revision.
