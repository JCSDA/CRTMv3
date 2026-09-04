# JEDI/UFO status and refresh plan for CRTM REL-3.2.0

**Scope note.** The subject of this document is the JEDI bundle and UFO, not
this repository. It is pinned to a local bundle state and will go stale on
their schedule rather than on ours. Re-verify before acting on it. The bundle
was subsequently refreshed on 2026-07-27 to 54/54, with the UFO
`Use_MWland_Atlas` opt-in added and five IR GsiHofX failures shown to be
pre-existing rather than caused by REL-3.2.0; none of that was pushed.

Drafted 2026-07-27 from the local bundle at
`~/jedi/jedi-bundle_REL-3.2.0_testing` (last built 2026-07-08). Companion to
`gsi_crtm_rel320_interface_plan.md`.

## Where things stand

The port already exists. The bundle pins four repos to `feature/btj_REL-3.2.0`
branches, and each local checkout is in sync with its origin feature branch:

| Repo | Local head | Port content | Behind its develop |
|---|---|---|---|
| crtm (JCSDA/CRTMv3) | `a8a7329` 2026-07-18 | the release branch itself | 54 commits behind origin feature branch |
| ufo | `b3130f1f2` 2026-07-02 | 12 commits ahead of develop | 18 |
| fv3-jedi | `9d214655` 2026-06-06 | 4 ahead (test-data staging only) | 8 |
| mpas-jedi | `95d13d7` 2026-06-06 | 4 ahead (test-data staging only) | 9 |

**UFO operator port is substantive and current with the v3 API:**
`CRTM_Init` is called with the full v3 keyword set (`NC_File_Path`,
`Aerosol_Model`/`AerosolCoeff_Format`/`AerosolCoeff_File`, `Cloud_Model` +
format/file, all per-surface EmisCoeff files); `Options` usage includes the
new `Compute_Down_Radiance`, `Compute_Down_Radiance_Profile`,
`Compute_Up_Radiance_Profile` (exposed as ObsDiagnostics, `f739b9888`) and
`n_Stokes` (`b3130f1f2`); geometry lat/lon/date is set for location-aware
surface models (`388187d4e`, the TELSEM2 prerequisite); netCDF coefficient
staging and `.nc` migration are done (`1508c88bf`, `ebc7eccd4`); cloud-mapping
logging and the ICE_CLOUD non-scattering warning are in (`f8f3bf89f`).
fv3-jedi and mpas-jedi carry only test-data symlink staging (netCDF
coefficients, NLTE/AC sibling directories).

## Gap analysis: what changed under the port since 2026-07-02

The bundle crtm (`a8a7329`) predates six merged CRTM PRs plus the doc work
(now at `670d271`). The ones with UFO-side consequences:

1. **TELSEM2 flipped to opt-in** (PR #334, `2ffdf22`, 2026-07-19). The UFO
   port made TELSEM2 usable (lat/lon/date) but exposes no activation control;
   it implicitly relied on presence-activation, which no longer exists. With
   current CRTM, a staged atlas file is ignored unless `CRTM_Init` gets
   `Use_MWland_Atlas=.TRUE.` or an explicit `MWlandCoeff_File`. UFO must add
   a YAML option (say `mw land emissivity atlas: telsem2`) that passes one of
   those through; without it the TELSEM2 path is unreachable from JEDI.
2. **ICE_CLOUD warning wording is now conditionally stale** (`f8f3bf89f`
   warns ICE_CLOUD is non-scattering at MW). True for Mie-TAMU tables; false
   for DDA-ARTS tables on REL-3.2.0 (ICE_CLOUD now scatters via IconCloudIce).
   Low priority: soften the message or condition it on the loaded table type.
3. **No other API breaks:** UFO does not use `Obs_4_downward_P` (removed),
   and the ODPS modernization (#343), UV NO2 (#341), path-length (#238),
   SNICAR (#324), and Fastem1 Jacobian fix are API-neutral for UFO. The UV
   forward-operator enablement (#339) is an opportunity: OMPS/TEMPO-class
   SpcCoeffs can now run, which UFO could exercise later.
4. **Baseline sensitivity:** the CRTM-side NLTE/AC staging repair changed
   CrIS 4.3 micron radiances (up to ~8.5 K daytime). UFO/fv3-jedi/mpas-jedi
   staged their NLTE/AC symlinks on 2026-06-06; after refreshing crtm,
   any ufo-data/fv3-jedi-data reference values for hyperspectral IR need
   re-verification (correct NLTE-on values may differ from stored refs).

**Merge-conflict outlook for `ufo` develop sync (18 commits):** three touch
`src/ufo/operators/crtm`: #4219 (per-channel surface_emissivity
ObsDiagnostics, overlaps the port's ObsDiagnostics additions), #4214
(zero-out Jacobian option, TLAD code), #4173 (Fortitude lint churn). Expect a
real but modest merge; everything else upstream is non-CRTM filters/operators.
fv3-jedi (8) and mpas-jedi (9) upstream commits are unrelated to the staging
commits; trivial merges expected.

## Refresh plan

1. **Update crtm in the bundle** to `670d271` (or the release tag once cut):
   `git -C crtm pull` or a fresh ecbuild configure (the bundle uses `UPDATE`).
2. **Sync each feature branch with its develop** (ufo first, then fv3-jedi,
   mpas-jedi): merge `origin/develop`, resolving the operators/crtm overlap
   (#4219/#4214/#4173). Keep the merge separate from any new feature work.
3. **UFO code follow-ups on the feature branch:**
   a. Expose the TELSEM2 opt-in (`Use_MWland_Atlas` or `MWlandCoeff_File`)
      through the operator YAML config (required for the atlas to be
      reachable at all).
   b. Refresh the ICE_CLOUD warning wording (Mie-TAMU-only claim).
   c. Optional: surface the OpenMP note (CRTM_Init clamps to 1 thread when
      OMP_NUM_THREADS is unset; process-global) in UFO docs, same caution as
      GSI.
4. **Rebuild the bundle and run the test ladder:** ufo ctests (crtm operator
   suite), then fv3-jedi and mpas-jedi radiance/hofx tests. Expect and triage
   hyperspectral-IR reference diffs from the NLTE staging repair; regenerate
   ufo-data / fv3-jedi-data references where the new values are the correct
   (NLTE-on) ones.
5. **Coefficient staging:** after the fix tarball re-roll (release blocker
   tracked in the CRTM release docs), refresh the bundle's
   `test-data-release/fix_REL-3.2.0.0` and re-verify md5-pinned downloads.
6. **PR sequencing to jcsda-internal:** once the release is tagged, the ufo
   feature branch (post-sync) goes to ufo develop referencing the crtm tag;
   fv3-jedi/mpas-jedi staging PRs follow; the bundle CMakeLists flips crtm
   from `BRANCH feature/btj_REL-3.2.0` to the `v3.2.0` tag. Coordinate with
   the JEDI infrastructure team on the spack-stack crtm recipe (same
   dependency as the GSI plan, section 1).

## Notes

- The local build tree predates all of the above (2026-07-08); do not trust
  cached ctest results.
- soca/coupling/rttov and the rest of the bundle ride on develop and need no
  CRTM-related action.
