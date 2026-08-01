# CRTM v3.2.0 Release Notes

**Status:** release candidate. The library code is frozen; the coefficient
tarball is not yet published.

**Coefficient data:** the v3.2.0 coefficient tarball has **not been rolled or
published**. Two coefficient efforts are still in progress and both change files
that would ship in it: the IASI-NG regeneration onto the 2026 gas epoch, and the
GeoXO `gxi`/`gxs` rework. The tarball will be rolled once those land, and this
section will then carry its size and checksum.

Until then the checksum pinned in `test/CMakeLists.txt` and
`Get_CRTM_Binary_Files.sh` still refers to the **June 2026** tarball, which
predates the whole of the July coefficient campaign. A default build therefore
downloads a coefficient tree that cannot reproduce this release's test suite.
Build a release candidate against the staging tree instead:

```
cmake .. -DFIX_FILE_PATH=<path-to>/fix_REL-3.2.0.0/fix
```

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
- **Vector radiative transfer (`Options%n_Stokes > 1`) substantially repaired.**
  The surface is now handed to the solver in the Stokes basis rather than
  (V,H); the third and fourth Stokes components survive the surface aggregation
  and the azimuthal accumulation instead of being zeroed at both ends; the
  clear-sky path gained a vector solver, so a cloud-free polarimetric run no
  longer returns Q = U = V = 0; the fractional-cloud clear/cloudy combine is
  differentiated across all Stokes components; and `RTSolution%Radiance` is now
  the channel-polarized measurement rather than Stokes I (see behavior change
  14). Tangent-linear, adjoint and K are verified at `n_Stokes = 4` against
  finite differences, the adjoint dot-product identity and K against AD, for
  both atmospheric and surface control variables including wind direction.
  **This path has not been validated against an external model; see the known
  issues before using U or V quantitatively.**
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
- **New sensor: `gems2_amethyst`** (Weather Stream GEMS2 24-channel
  microwave sounder, 118.75 GHz oxygen bank plus 160-183.31 GHz humidity
  bank, on the GEMS2-Amethyst smallsat). Generated with crtm-coeffgen
  (MonoRTM, ECMWF84); brightness-temperature validation against
  line-by-line truth averages 0.04 K. WMO ids are invalid-value sentinels
  until C-5/C-8 assign codes; the 118.75 GHz line-center channels carry no
  Zeeman treatment. Unrelated to the Korean GEMS UV spectrometer
  (`gems_gk2b`) despite the acronym.
- **FY-3 microwave family completed (12 sensors, FY-3C through FY-3G).**
  New coefficient pairs for MWHS-2 (FY-3C/D and the E-variant on FY-3E/F),
  MWTS-2 (FY-3C/D), MWTS-3 (FY-3E/F), MWRI (FY-3C/D), MWRI-2 (FY-3F;
  instrument failed in 2025, coefficient serves historical reprocessing),
  and MWRI-RM (FY-3G). Generated with crtm-coeffgen from the NWP-SAF
  passband definitions and validated with forward, weighting-function, and
  adjoint physics checks; the 57 GHz line-splitting channels reproduce the
  AMSU-A weighting-function progression.
- **AMSR3 humidity-channel bandwidths corrected to WMO OSCAR.** The bundled
  AMSR3 spectral response function was too wide on its three highest channels:
  165.5 GHz by a factor of 1.25, 183.31 +/- 7 by 2.35, and 183.31 +/- 3 by 2.72,
  all against the per-sideband figures published by WMO OSCAR. The other 18 of
  21 channels matched OSCAR exactly and are unchanged, and STAR independently
  revised the same three. Against 2152 collocated observations the 183.31 +/- 3
  bias improves from +4.270 K to +3.422 K and 183.31 +/- 7 from +3.118 K to
  +2.906 K, with channels 1 to 18 unmoved. The larger effect is on the
  Jacobians: the 183.31 +/- 3 water vapour Jacobian error against a
  reference-free tiled-channel truth falls from 9.468 percent to 0.338 percent,
  which retires the belief that ODPS could not represent a wide double-sideband
  channel. That was never an ODPS limitation, only a band 2.72 times too wide.
  Recorded caveat: 165.5 GHz moved the wrong way on observation-minus-background
  (+3.951 K to +4.168 K), but that channel sits inside a 3.4 to 4.2 K
  common-mode bias shared by every variant including STAR's, so
  observation-minus-background cannot arbitrate it; the change rests on OSCAR,
  on STAR's independent revision, and on the 18-of-21 exact match. Spectroscopy,
  training profiles and algorithm are unchanged, and only the TauCoeff differs
  (the regenerated SpcCoeff is identical in every data variable).
- **INSAT-3DS visible sensors completed.** `v.imgr_insat-3ds` and
  `v.sndr_insat-3ds` previously shipped a SpcCoeff with no TauCoeff and
  could not pass `CRTM_Init`; both now carry TauCoeffs generated from the
  measured ISRO SRFs (crtm-coeffgen, ECMWF84) and regenerated SpcCoeffs
  whose centroids match the previously shipped files to better than
  0.03 nm.
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
- **Intel builds fixed: TELSEM2 atlas load and the PRA polarization angle.**
  Two defects that a GNU-only test suite had been passing, both found by adding
  an `ifx` build to the release verification.

  On Intel, `CRTM_Init` segmentation faulted whenever the TELSEM2 atlas was
  requested, on any machine with the stock 8 MB Linux stack limit. The atlas is
  large: `n_data` is 2,770,889, so `cell_number` is 11 MB and `emissivity` is
  155 MB. `nf90_get_var` takes an assumed-shape dummy and passes it down to an
  F77 layer that takes an assumed-size one, and the compiler bridges the two
  with a contiguous copy-in temporary. That temporary is created inside the
  netCDF library's own compiled code, so it follows the flags netCDF was built
  with and not CRTM's, which is why no CRTM compiler flag can prevent it. The
  reader now takes the atlas in bounded slices, so every temporary stays small
  however netCDF was built. This is a read-path change only; the values loaded
  are identical.

  Separately, the `PRA_POLARIZATION` surface-optics branch divided zero by zero
  at nadir scan angle. The shared denominator reduces exactly to
  `|sin(phi)|*sqrt(1 + sin(theta_f)^2)`, and both numerators vanish with it, so
  the expression was undefined there and returned whatever the compiler folded
  it to: GNU gave a polarization weight of 1, which selects the **opposite**
  polarization to the correct limit of 0, and Intel gave a NaN that propagated
  into the radiance, the weighting functions and the adjoint. The singularity is
  removable, and passing the two numerators to `ATAN2` removes it rather than
  special-casing it. This affects `gems2_amethyst` and `gems2_beryl` only, both
  new in this release, so no previously released sensor changes behavior. The
  expression had been duplicated in the forward, tangent-linear, adjoint and
  Stokes-projection paths and is now one shared function.
- **Runtime OpenMP control.** `OMP_NUM_THREADS` is honored at run time (no
  longer captured at configure time).
- **Expanded self-checking test coverage.** New baseline-independent checks
  include general TL-vs-FD and adjoint-consistency tests across the three
  main sensor types (#280), multi-sensor single-call bit-consistency, ODPS
  group-validation and long-path initialization tests, a DDA-ARTS ICE_CLOUD
  behavior pin, multi-sensor OMPS UV and TEMPO UV+VIS physics verifications
  (each registered when its pre-release coefficient pairs are present), and
  OpenMP thread-count consistency tests (#111).

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

   The 200 GHz figure is a safety floor, not physics. It was placed where the
   traditional sounding sensors stop so that enabling PARMIO could not disturb
   operational channels. `Options%Use_PARMIO_MWSSEM` drops the floor and runs
   PARMIO everywhere its table has data. Table coverage is checked separately
   and is never relaxed, including when a caller opts in, because the
   alternative is a confident number computed at the wrong frequency: the
   interpolator otherwise clamps silently, and a 204.78 GHz channel was being
   evaluated at 229 GHz. Channels where PARMIO is wanted but uncovered fall
   back to FASTEM.
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

9. **Duplicate `_j2` sensor aliases removed from the fix tree.** Six
   sensors shipped twice under both a `_j2` and an `_n21` name for the same
   satellite (JPSS-2 = NOAA-21). Five of those pairs are identical in every
   data variable; the sixth (`atms_j2` / `atms_n21`) agrees to
   floating-point round-off only (largest difference 5.7e-14 in
   `Frequency`, with the derived Planck and band-correction coefficients
   differing in their last bits from independent computation), which is
   radiometrically identical but not literally bit-identical. Otherwise the
   pairs differ only in the internal `Sensor_Id` string and the creation
   timestamp. A seventh entry, `v.viirs-m_j2`, had no `_n21` counterpart in
   v3.1.4 and is a straight rename rather than a de-duplication. The `_j2`
   copies are gone; use the `_n21` name:

   | removed | use instead |
   |---|---|
   | `atms_j2` | `atms_n21` |
   | `atms_j2-srf` | `atms_n21-srf` |
   | `cris-fsr_j2` | `cris-fsr_n21` |
   | `viirs-i_j2` | `viirs-i_n21` |
   | `viirs-m_j2` | `viirs-m_n21` |
   | `v.viirs-i_j2` | `v.viirs-i_n21` |
   | `v.viirs-m_j2` | `v.viirs-m_n21` |

   `cris-fsr_j2.NLTECoeff.nc` went with them (identical to
   `cris-fsr_n21.NLTECoeff.nc`, and orphaned once its SpcCoeff was removed).

   `CRTM_Init` resolves coefficients by filename, so any caller configured
   with a `_j2` sensor id must be updated or initialization will fail. The
   duplication was a maintenance hazard as much as dead weight: nothing in the
   tree recorded that the two names were meant to be twins, so regenerating
   one would have silently left the other stale.

   Note this does **not** apply to the visible-channel files that share
   content across detector variants (`v.imgrD1..D8_gNN`, `v.sndrD1..D4_gNN`,
   `v.mi-l/m_coms`). Those are genuinely distinct sensors whose parent IR
   channels differ; they share a common visible channel by instrument design
   and are all retained.

10. **Twelve further sensor renames.** These are name changes only: in every
    case the replacement file is identical to the removed one in every data
    variable, and only the filename and the internal `Sensor_Id` differ.
    `CRTM_Init` resolves coefficients by filename, so a caller configured with
    an old name will fail to initialize and must be updated.

    | removed | use instead | why |
    |---|---|---|
    | `viirs-i_j1` | `viirs-i_n20` | JPSS-1 became NOAA-20 at launch. Every one of these files already carried `WMO_Satellite_Id` 225, which is NOAA-20 in WMO C-5, so the `_j1` name contradicted the file's own content. |
    | `viirs-m_j1` | `viirs-m_n20` | as above |
    | `v.viirs-i_j1` | `v.viirs-i_n20` | as above |
    | `v.viirs-m_j1` | `v.viirs-m_n20` | as above |
    | `v.viirs-dnb_j1` | `v.viirs-dnb_n20` | as above |
    | `mwi_metop-sg-a1` | `mwi_metop-sg-b1` | platform correction. MWI flies on Metop-SG-B, not Metop-SG-A. |
    | `tms_tomorrow-s01_v4` | `tms_tomorrow-s01_v4-STAR` | lineage relabel, so the v4 delivery is tagged to the organization it came from. Six sensors, `s01` through `s06`. |

    The six `tms_tomorrow-sNN_v4` entries follow the same pattern and are not
    listed individually.

11. **Five products withdrawn.** Two were duplicates under a nonsensical name
    and three were never-flown or notional instruments:

    | withdrawn | why |
    |---|---|
    | `airs_g13` | identical in every data variable to `airs281_aqua`, which still ships. Its own `WMO_Satellite_Id` is 784 (Aqua) and its sensor id is 420 (AIRS), so the `g13` suffix never described the content. |
    | `iasi_g13` | identical in every data variable to `iasi616_metop-a`, `-b` and `-c`, all of which still ship and which are distinguished from each other only by their WMO satellite ids (4, 3 and 5), as they should be for one instrument design on three platforms. `iasi_g13` carried WMO satellite 1022, which identifies no platform. `iasi_g13.NLTECoeff.nc` went with it. |
    | `ssmis_f20` | DMSP F-20 was cancelled and never launched. The file carried `WMO_Satellite_Id` 1023, the invalid-value sentinel, because no satellite id was ever assigned. |
    | `zssmis_f20` | the Zeeman companion to the above, withdrawn with it. |
    | `atms-ng_v1` | a notional next-generation ATMS with 1169 channels and placeholder WMO ids (satellite 1, sensor 1). No such instrument exists. |

    Nothing that still ships is lost by any of these: users of `airs_g13` or
    `iasi_g13` should switch to the correctly named file, which holds the same
    numbers.

12. **Twelve unnamed CloudCoeff development artifacts removed.** The v3.1.4
    tree carried `test_new.bin_type0` through `test_new.bin_type10` and
    `test_new.bin_MIESNOW` under `CloudCoeff/Little_Endian/`, and the netCDF
    transition converted them along with everything else. They carry no title,
    history or comment attribute of any kind, nothing in the library or the
    test suite references them, and their names describe a conversion run
    rather than a product. Removing them takes about 540 MB off the tarball.
    The named cloud tables are all retained, including the microphysics-scheme
    variants (`CloudCoeff.GFDLFV3`, `CloudCoeff.Thompson08`, `CloudCoeff.WSM6`)
    and the TAMU tables, none of which the test suite exercises either.

13. **OMPS replaced by per-platform NOAA-20 and NOAA-21 products.** The two
    shipped OMPS files were unusable and mislabelled, and have been retired in
    favor of four regenerated products:

    | removed | replaced by |
    |---|---|
    | `u.omps-npAllFOV_j2` | `u.omps-np_n20` (151 ch), `u.omps-np_n21` (158 ch) |
    | `u.omps-tcAllFOV_j2` | `u.omps-tc_n20` (196 ch), `u.omps-tc_n21` (198 ch) |

    Three separate defects motivated this:

    - **The files could not be loaded at all.** Both carried
      `Group_Index=4`, which has been Zeeman-reserved with zero components
      since 2008, so `CRTM_Predictor_Create` failed outright. The replacements
      are `Group_Index=8` (`GROUP_UV_NO2`), the UV variant carrying a scene-NO2
      component, and all four now pass `CRTM_Init`.
    - **The platform labels were wrong, in opposite directions.**
      `u.omps-npAllFOV_j2` was labelled NOAA-21 (WMO 226) but its channel set
      matches the NOAA-20 grid (rms 0.04 nm, max 0.07 nm; every other
      platform's grid is at least 7 times farther); NOAA-21's nadir profiler
      natively has 158 channels reaching 245 nm, not 151. `u.omps-tcAllFOV_j2`
      really was NOAA-21 content, but indexed with a three-channel offset.
    - **Only one platform was represented** where two instruments exist.

    Channel numbering now follows each platform's own SRF. For total column,
    old channel *N* corresponds to `u.omps-tc_n21` channel *N+3*; the old file
    omitted n21 channels 1 to 3 (298.1 to 298.9 nm) and extended three channels
    past its red end. Channel selections carried over from the old files must
    be re-mapped, not reused.

    The per-platform channel sets are verified against primary sources: the
    JPSS NOAA-21 OMPS SDR validated-maturity record (nadir profiler 158
    channels, nadir mapper 198) and the published instrument table in Yan et
    al. 2024, doi:10.3390/rs16234488 (nadir profiler SNPP 147 / NOAA-20 151 /
    NOAA-21 158; nadir mapper 196 / 196 / 198). Note the same record flags
    NOAA-21 nadir mapper radiances below 302 nm (roughly `u.omps-tc_n21`
    channels 1 to 10) as not validated for operational use.

14. **`RTSolution%Radiance` is now the channel-polarized measurement when
    `n_Stokes > 1`, not Stokes I.** The emergent Stokes vector is projected onto
    the channel's polarization, so a vertically polarized channel reports I+Q
    where it previously reported I. `Brightness_Temperature`, which is derived
    from it, moves with it. `RTSolution%Stokes` is unchanged and still holds the
    physical (I, Q, U, V). Anyone already running `n_Stokes > 1` on a polarized
    channel will see the reported radiance and brightness temperature change by
    the polarization difference, which over ocean is order 20 percent of the
    signal. The projection weights are taken from the scalar path's own
    polarization handling, so the two now agree: a run with the channel forced
    to pure vertical or pure horizontal reproduces the corresponding scalar
    radiance to machine precision. Nothing at `n_Stokes = 1` changes.

15. **`RTSolution_AD%Radiance` and `RTSolution_K%Radiance` are now honoured as
    input seeds on the vector path.** Previously `%Radiance` was an output alias
    for `Stokes(1)` but not an input one, and seeding it for an `n_Stokes > 1`
    run silently produced a zero Jacobian; only `%Stokes` or
    `%Brightness_Temperature` worked. Both now work. Code that seeded
    `%Radiance` and `%Stokes(1)` together to work around this will now double
    count.

16. **`RT_Algorithm_Id = RT_SOI` with `n_Stokes > 1` is now an error.**
    Previously the vector branch was taken before the algorithm selector was
    consulted, so the caller silently received ADA results labelled as SOI. SOI
    has no vector solver, so this is now rejected with a message naming RT_ADA.

17. **`CRTM_MWwaterCoeff_Load_FASTEM` discards the previously loaded scheme and
    reports failure.** Switching scheme, for example FASTEM6 to FASTEM4, used to
    leave the shared coefficient structure deallocated while the function
    returned SUCCESS, because the underlying setter rejects a shape mismatch by
    destroying the target. Scheme switching now works, and a failed load returns
    FAILURE instead of SUCCESS. Related: a coefficient dimension mismatch used to
    terminate the program with a Fortran runtime formatting error rather than
    reporting cleanly.

## Known issues and limitations

- **Sub-mm thin frozen cloud (≈ 325 GHz):** optically thin frozen-cloud
  layers can produce nonphysical TBs through the adding-doubling/MOM path
  (small-τ matrix conditioning with high phase-function truncation orders).
  Affects sub-mm scattering scenes only; under investigation.
- **The FASTEM to PARMIO handover at 200 GHz is a step, not a blend.** The two
  models are independent and are not reconciled at the boundary, so a sensor
  with channels either side of 200 GHz sees a discontinuity in ocean
  emissivity there. Measured at -1.42 K mean in brightness temperature for
  TROPICS channel 12. This is a consequence of switching models at a frequency
  rather than a defect in either one, and it is why the floor sits where no
  operational sounding channel crosses it. PARMIO's own published validation
  stops at 165.5 GHz, so its whole default dispatch range is above the range
  its authors validated; the table labels that band
  `extrapolated-experimental` in its `confidence_label` variable, and callers
  should read that label rather than assume the table is uniform.
- **ODSSU netCDF supports the ODPS sub-algorithm only** (the shipped SSU
  files are ODPS-based; ODAS-based SSU coefficient files remain binary-only).
- **Options binary I/O does not persist the new fields** (`n_Stokes` and the
  radiance-profile switches); they will be added with the next format
  revision.

### Polarimetric (`Options%n_Stokes > 1`) radiative transfer

The vector path is functional and internally verified, but it has **not been
validated against an independent radiative transfer model or against
observations**. Treat it as a capability under development rather than a
production-ready product. The specific limitations follow.

- **No external reference.** Every check on this path compares CRTM against
  itself: tangent-linear against finite differences, the adjoint dot-product
  identity, K against AD, physical invariants such as the polarization bound,
  and agreement with the scalar path in the limits where the two must agree.
  None of that can fix a convention that is consistently wrong. In particular
  the **sign convention of the third Stokes component** between the surface
  models and the solver is unverified: every internal test passes unchanged if
  the sign of U is flipped. Anyone using U or V quantitatively should establish
  the sign against an external reference first.
- **The default microwave water surface has no polarimetric model.** The
  `CRTM_Init` default is FASTEM6, whose azimuth model parameterises the
  vertical and horizontal components only and returns the third and fourth
  Stokes components as identically zero. A polarimetric surface requires
  `MWwaterCoeff_Scheme = 'FASTEM4'` (or FASTEM5), or PARMIO, which is used
  automatically at and above 200 GHz when its lookup table is present. Note
  that `MWwaterCoeff_File` does **not** select the model; the scheme argument
  does.
- **Microwave only.** The coupled (V,H) to Stokes surface branch exists only in
  the microwave section of the surface optics. Infrared and visible sensors
  populate the first Stokes component alone, so a vector run there returns
  Q = U = V = 0 from the surface regardless of geometry.
- **Aerosols contribute no polarization.** The shipped `AerosolCoeff` carries a
  single phase element, so a vector run mixes polarized cloud scattering with
  unpolarized aerosol scattering. This is deliberate and is not blocked, since
  a scalar aerosol table must not prevent a polarized cloud run.
- **The surface reflects no U or V.** Both FASTEM and PARMIO set the third and
  fourth Stokes reflectivity components to zero, so U and V are emitted by the
  surface but never reflected. This is exact for clear sky, where the
  downwelling reaching a specular surface is unpolarized and symmetry about the
  meridional plane forbids a reflected U. It is an approximation once the
  downwelling is itself polarized by scattering.
- **Polarimetric cloud lookup tables are not validated.** The `n_Stokes > 1`
  scattering path requires a six-phase-element table (the experimental
  `CRTM-Exp` scheme). Those tables have not been validated for full-Stokes
  work, so polarized scattering results carry that uncertainty independently of
  the code.
- **Phase-matrix normalization is inconsistent for below-diagonal polarized
  blocks.** `Normalize_Phase` scales each row's intensity and polarized
  elements together and then applies an intensity-only symmetry copy, so a
  block below the diagonal ends up with its (1,1) element carrying one row's
  normalization and its polarized elements another. This does not produce a
  polarization-bound violation with the shipped coefficients, where the
  measured worst ratio is 0.65, but it is not the correct polarized symmetry
  treatment.
- **Mixed-polarization channels are projected at the sensor angle.** For the
  V/H-mixed, constant-mixed and PRA polarizations the vector path applies the
  polarization mixing once, to the emergent radiance at the sensor angle, which
  is where a receiver projects. The scalar path instead applies it to the
  surface emissivity at every quadrature angle. The two coincide when there is
  a single angle, and differ slightly for a scattering mixed-polarization
  channel.
- **`plus45L`, `minus45L`, `RC` and `LC` polarizations are treated as
  vertical**, inherited unchanged from the scalar path so that the two agree.
  These are placeholders, not the true projections.
- **`RT_Algorithm_Id = RT_SOI` is not supported with `n_Stokes > 1`** and is now
  rejected; see the behavior changes above. SOI has no vector solver.
