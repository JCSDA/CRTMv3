# GSI interface update plan: CRTM v2.4.0.1 to REL-3.2.0

**Scope note.** The subject of this document is NOAA-EMC/GSI, not this
repository. It is pinned to external commits and will go stale on GSI's
schedule rather than on ours. Treat the pinned hashes below as the
as-of state, and re-verify against GSI `develop` before acting on it.

Drafted 2026-07-27 against NOAA-EMC/GSI `develop` @ `ec8215d` (2026-07-15, shallow
clone at `/home/ben/CRTM/GSI`) and CRTM `feature/btj_REL-3.2.0` @ `670d271`.

**Re-pinned 2026-08-04** to CRTM `feature/btj_REL-3.2.0` @ `830a280`, 53 commits
later. The GSI side is unchanged and still as of `ec8215d`. Every symbol in the
API surface listed below was re-checked against the new head and is still
present; the only `CRTM_Init`-adjacent signature change in those 53 commits is
an added argument on an unrelated internal routine. So section 0 stands as
written. What changed materially is section 2, where the coefficient-staging
risk can now be quantified rather than merely flagged.

## 0. Baseline: how GSI couples to CRTM today

- **Version pins:** modulefiles load `crtm/2.4.0.1` and set `CRTM_FIX` to a
  `v2.4.0.2` fix directory (binary, flat layout) on every platform
  (`modulefiles/gsi_common.lua`, `gsi_<machine>.lua`).
- **Build:** `src/gsi/CMakeLists.txt` does `find_package(crtm REQUIRED)` and
  links `crtm::crtm`. No version constraint is stated.
- **Primary interface:** `src/gsi/crtm_interface.f90` (3577 lines;
  `init_crtm` / `call_crtm` / `destroy_crtm`), plus
  `set_crtm_cloudmod.f90`, `set_crtm_aerosolmod.f90`, `setuprad.f90`,
  `setupaod.f90` (AOD via `crtm_aod_k`), `cads.f90` (cloud detection), and the
  obs readers `read_bufrtovs/read_iasi/read_cris/read_iasing/read_atms/
  read_mws/read_saphir`.
- **API surface used** (all verified present in REL-3.2.0):
  - `crtm_module`: types (`crtm_atmosphere_type`, `crtm_surface_type`,
    `crtm_geometry_type`, `crtm_options_type`, `crtm_rtsolution_type`,
    `crtm_channelinfo_type`), create/destroy/zero/associated routines,
    `crtm_init`, `crtm_destroy`, `crtm_forward`, `crtm_k_matrix`,
    `crtm_channelinfo_subset`, `crtm_channelinfo_n_channels`,
    `ssu_input_setvalue`, `crtm_irlandcoeff_classification`, constants
    (`success`, `fp`, `microwave_sensor`, `toa_pressure`, `max_n_layers`,
    `limit_exp`, gas ids, aerosol type ids).
  - `crtm_cloudcover_define`: overlap-method constants.
  - `crtm_aod_module`: `crtm_aod_k`.
  - **Internal-module reaches (fragile by design, but still exported by
    3.2.0):** `crtm_spccoeff` (`sc`, `crtm_spccoeff_load/destroy`) for channel
    frequency/wavenumber/polarization in `setuprad`, `setupaod`, and the
    readers; `crtm_aerosolcoeff` (`AeroC%Reff`, `AeroC%RH`) in
    `set_crtm_aerosolmod`.
  - `RTSolution` fields read: `brightness_temperature`, `radiance`,
    `surface_emissivity`, `layer_optical_depth`, `total_cloud_cover`,
    `upwelling_overcast_radiance`. `Options` field set:
    `use_antenna_correction`. All still exist in 3.2.0.
  - **v3.3.0 note (TELSEM2 DA program, `feature/btj_telsem2_da`):**
    `RTSolution_K%Surface_Emissivity` now reports the *total* derivative
    d(Tb)/d(emissivity) (the reflected-downwelling term is folded in;
    values change vs 3.2.0), and a new forward field
    `RTSolution%Surface_Emissivity_Std` carries the TELSEM2 Release-2
    per-channel emissivity error std (zero = unavailable) — relevant to
    the pending UFO '#4219 per-channel surface_emissivity' work and any
    GSI emissivity-sink/error-inflation use. The full inter-channel
    emissivity covariance is available post-Init via
    `CRTM_TELSEM2_Emissivity_Uncertainty` (exported from `CRTM_Module`).
- **Good news:** GSI never touches `Options%Obs_4_downward_P` (removed in
  3.2.0), never reads RTSolution fields that changed, and the v2.4.1
  aerosol-model init arguments GSI wants are present in 3.2.0
  (`Aerosol_Model`, `AerosolCoeff_File`, `AerosolCoeff_Format`; currently
  commented out in `crtm_interface.f90` with "for crtm2.4.1" markers).

## 1. Build and environment (coordination-heavy, start first)

1. Get CRTM v3.2.0 into spack-stack / the WCOSS2 stack (coordinate with EPIC
   and NCO). Differences from the 2.4 package: v3 requires netCDF4/HDF5 at the
   CRTM level (2.4 was self-contained binary I/O), builds a shared lib by
   default, and needs CMake 3.20+ and git-lfs for source builds.
2. Verify the exported CMake target name from the v3 `crtm-config.cmake`
   matches `crtm::crtm` (v3 installs an EXPORT set; confirm the namespace, and
   add `find_package(crtm 3.2.0 REQUIRED)` version floor in
   `src/gsi/CMakeLists.txt`).
3. Update every `modulefiles/gsi_*.lua`: `crtm_ver` 3.2.0, new
   `crtm_fix_ver` (see 2), and the `CRTM_FIX` paths.

## 2. Coefficient (fix) staging

1. Produce a GSI-consumable netCDF fix set from `fix_REL-3.2.0.0.tgz`.
   Decision: keep GSI's **flat directory** convention (the 3.2.0 SpcCoeff
   reader probes flat-layout siblings for ACCoeff/NLTECoeff, so a flat tree
   works) or adopt the canonical `fix/<type>/netCDF/` tree. Recommendation:
   flat, minimal churn; symlink from the canonical tree.
2. Cross-check the GSI active-sensor roster (global satinfo + regional) against
   the 3.2.0 inventory (`REL-3.2.0_coefficient_inventory.md`, 515 sensors).
   Verify specifically: ACCoeff siblings for AMSU-A/MHS (GSI sets
   `use_antenna_correction`), NLTECoeff for CrIS/IASI/AIRS, zssmis f16-f19
   (no f20 in 3.2.0), SSU, and the viirs-m companion files that `read_cris`
   stages for CADS.

   **This is the highest-consequence item in the migration, and it fails
   silently.** Measured 2026-08-04 by running the v3.1.4 and v3.2.0 libraries
   over 84 profiles against identical coefficients (harness in
   `release_wrap_2026-08/code_delta_314_vs_320/`):

   - A netCDF `SpcCoeff` read with the `NLTECoeff` sibling absent loses the
     non-LTE correction entirely. On AIRS that is **up to 36.4 K** on the
     4.3 um shortwave CO2 channels in daylight, and the **temperature
     Jacobians on those channels change by a median of 23 percent, up to 90
     percent**. The affected set is channels 1900-2114. At night the
     difference is exactly zero, since the correction is solar-driven.
   - With the `ACCoeff` sibling absent and `use_antenna_correction` set, AMSU-A
     loses antenna correction: **up to 1.225 K, mean 0.611 K, on every
     channel**. A systematic bias rather than scatter, which for assimilation
     is the worse shape.
   - There is **no warning and no error**. The file is simply never opened.
     Verified directly: deleting the sibling changes the answer on zero rows.

   Practical consequence for GSI: the sibling files are not optional companions
   to be staged when convenient. Omitting them produces a run that completes
   cleanly and is wrong, and the failure is invisible in anything short of a
   radiance comparison. Stage them, then confirm by deleting one in a scratch
   copy and checking that the radiances move.

   Note this bites only the netCDF path. GSI's current v2.4 binary `SpcCoeff`
   streams both substructures inline, which is why the issue does not exist
   today and appears only on migration.
3. The 2.4 Big_Endian/Little_Endian binary split disappears; netCDF is
   portable. ~~Blocker note: the 3.2.0 tarball itself must be re-rolled first
   (staging tree is 143 files ahead; tracked in the release docs).~~
   **Updated 2026-08-06: final campaign roll, published.** The archive
   was re-rolled 2026-08-06 after the coefficient regeneration campaign staged
   its 24-file replacement set (this supersedes the 2026-08-04/05 rolls). Size
   3,377,422,223 bytes, md5 `88995873986cf2b077808a75d1c56f83`, verified
   against the staging tree at 1440 members, all netCDF, zero differences.
   Both pins (`Get_CRTM_Binary_Files.sh`, `test/CMakeLists.txt`) carry the new
   value.

   For GSI: the 2026-08-06 upload is live (server `Content-Length:
   3377422223`, confirmed 2026-08-06); GSI can pull the coefficient set the
   normal way and the md5 pin will verify it.

## 3. Code changes in GSI (small; compile-driven)

1. `crtm_interface.f90` `crtm_init` call: no required change (v3 accepts the
   2.4 argument list). Optionally uncomment the `Aerosol_Model` /
   `AerosolCoeff_File` / `AerosolCoeff_Format` block to pin the aerosol
   scheme explicitly (v3 default table is the classic GOCART-heritage set the
   GSI aerosol code assumes; pinning makes that assumption loud).
2. Replace hardcoded `.SpcCoeff.bin` (and CADS `.TauCoeff.bin`) existence
   probes with `.nc` names, or drop the manual probes and rely on the load
   functions' error status: `read_bufrtovs.f90`, `read_iasi.f90`,
   `read_cris.f90` (viirs-m names), `read_iasing.f90`, `cads.f90`.
   `crtm_spccoeff_load` in 3.2.0 defaults to netCDF with binary fallback, so
   the load calls themselves are fine once the probe filenames are fixed.
3. Confirm the aerosol id constants and `AeroC` field names compile clean
   (they exist in 3.2.0; `set_crtm_aerosolmod` also has a duplicated import at
   two call sites to update if names shift).
4. No changes needed for: cloud-cover overlap imports, SSU input, AOD K-matrix
   path, channel-subset logic (note: 3.2.0 hard-fails on duplicate or
   non-member channel subset lists where 2.4 silently misbehaved; if GSI ever
   feeds a bad list this now aborts loudly, which is the desired behavior).

## 4. Behavioral deltas to expect and validate (the real work)

1. **OpenMP caution (highest operational risk):** `CRTM_Init` in 3.2.0 reads
   `OMP_NUM_THREADS` at run time and, if it is unset or empty, calls
   `OMP_SET_NUM_THREADS(1)` process-globally. GSI is an OpenMP+MPI code: any
   job that does not export `OMP_NUM_THREADS` would have its host threading
   clamped after `init_crtm`. Audit GSI job cards / ush scripts to guarantee
   the variable is always exported.
2. **Radiance-level differences vs 2.4** will appear from the physics distance
   between 2.4.0 and 3.2.0 (cloudy-sky solver work, corrected coefficient
   metadata, netCDF-canonical coefficients, NLTE siblings). Expect small
   shifts in hyperspectral shortwave channels and MW window channels;
   validate per-sensor, not in aggregate.
3. **Validation ladder:**
   a. GSI ctest/regression suite vs a 2.4 control (expect diffs; catalog them).
   b. Single-observation experiments per sensor class (MW sounder, hyperspectral
      IR, geo IR, MW imager, SSU/SSMIS special paths).
   c. Low-resolution cycled experiment (2-4 weeks) comparing O-B/O-A statistics,
      penalty, and bias-correction spin-up per sensor against the 2.4 control.
   d. Full-resolution parallel before any operational or quasi-operational use.
4. **Performance check:** wall-clock and memory of setuprad with v3 (netCDF
   coefficient load at init, larger RTSolution). Confirm no thread
   oversubscription.

   **Correction (measured 2026-08-06).** An earlier revision of this item
   claimed CRTM-internal channel OpenMP was "inert under GSI's per-profile call
   pattern". That was wrong, and backwards. GSI passes one profile per call
   (`crtm_interface.f90`, `atmosphere` is `dimension(1)`), and CRTM only reached
   its channel-threading path when threads outnumbered profiles, so a single
   profile with `OMP_NUM_THREADS` greater than 1 was the *worst* case rather
   than an inert one: every spare thread took a slice of the channels along with
   a full set of per-thread scratch structures. Measured on a 22-channel
   microwave sensor, one profile on 16 threads ran about 30 times slower than
   the same build on one thread.

   This is fixed in CRTM (see behavior-change entries 20 and 21 in the release
   notes): channels are no longer split below a minimum work threshold, and the
   single-profile case is now at worst break-even and usually a modest gain.
   Two things still follow for GSI:

   - A hybrid MPI/OpenMP GSI on **v3.1.4 or earlier** is exposed to the old
     behavior. If such a configuration is in use, `OMP_NUM_THREADS=1` for the
     CRTM portion recovers the loss with no code change.
   - Threads give GSI little either way, because one profile is the smallest
     unit CRTM can divide well. The scaling lives in profile-level parallelism,
     which GSI's call pattern does not expose. If CRTM throughput inside
     setuprad ever becomes a bottleneck, batching profiles per call is the
     change that would matter, and it is a GSI-side change, not a CRTM one.

## 5. Rollout sequencing

1. Branch `feature/crtm_rel320` in GSI; changes from section 3 plus CMake floor.
2. Stack + fix staging (sections 1-2) on one dev platform (Hera or Orion) first.
3. Validation ladder (section 4) on that platform; document per-sensor deltas
   (the CRTM-side catalog `REL-3.2.0_changes_vs_develop.md` maps most of them).
4. PR to NOAA-EMC/GSI develop with the validation evidence; coordinate the
   global-workflow fix-version bump (`crtm_fix_ver`) separately, since
   global-workflow pins its own CRTM fix path.
5. Regional (RRFS/3DRTMA) follow-up after global acceptance; they share
   crtm_interface but have their own job cards (OpenMP audit again).

## Open questions for EMC/JCSDA

- spack-stack timeline for a crtm@3.2.0 recipe (depends on the release tag).
- Whether GSI wants the new 3.2.0 capabilities exposed as namelist options in
  a follow-on PR (PARMIO is automatic at >= 200 GHz once the LUT is staged;
  TELSEM2 land atlas opt-in; downwelling radiance outputs for ground-based DA;
  MW-land analytic Jacobians arriving as nonzero Surface_K columns that the
  radiance-Jacobian QC and bias code have never seen: verify nothing assumes
  those columns are zero).
- CADS reads SpcCoeff/TauCoeff directly by filename; confirm the CADS
  team's netCDF switch or keep binary copies of just those files as a bridge.
