# Handoff: convert ctest results/baselines from binary to netCDF

**Status: COMPLETE. Retained as a record, not as pending work.** The
conversion landed on 2026-05-31 in six commits ending at `ee6b908`, "test:
flip regression baselines from binary to netCDF (hard switch)". Everything
below is written in the forward tense because it was a working handoff
between sessions; read it as a description of what was done and why, not as
a task list. The suite counts quoted in the ground state are long superseded,
206 then 208 at the time, 239 as of 2026-08-01.

Worth keeping for two things the conversion turned up: four latent K-matrix
and adjoint bugs in the RTSolution netCDF path, and the convention that
adjoint output is rank-1 while k_matrix is rank-2.

## Goal (user request)
Store the ctest regression **baselines/results** in **netCDF** instead of binary. Decided scope:
**full conversion** (RTSolution + Atmosphere + Surface), **hard switch** (drivers write/read
only `.nc`; delete stale `.bin` baselines and re-seed). Confirmed by the user.

## Ground state (verified this session)
- Branch `feature/btj_REL-3.2.0`, HEAD `ed21dc7`, **8 commits ahead of origin**, working tree clean.
- Clean from-scratch build verified: `rm -rf build && cmake .. -DCMAKE_BUILD_TYPE=RELEASE && make -j8 && ctest` ⇒ **206/206**.
- GNU gfortran 13.3 in PATH builds the repo. The GNU build dir is `build/` (Intel `build-ifx`, nv `build-nv` exist but nv ctest is environment-blocked; ignore for this task).
- **WARNING for the implementer:** in the prior session, multi-line file *reads* were intermittently
  corrupted (Read tool returned plausible-but-fabricated lines). `grep`/`sed`/`md5sum` were reliable.
  **Trust the compiler and round-trip tests, not file reads.** After every edit: (1) `grep` to confirm
  it landed, (2) build crtm from scratch and READ the exit code, (3) run the round-trip test. Never
  claim a build/test passed without an observed exit code. (A green incremental ctest can run against a
  STALE object — always force a clean compile of the edited module.)

## Why this splits into easy + hard
| Baseline | Driver calls | netCDF I/O today | Work |
|---|---|---|---|
| RTSolution (`.RTSolution{,_K,_AD,_TL}.bin`) | ~27 drivers (fwd/TL/AD/K) | **EXISTS** — `CRTM_RTSolution_WriteFile(...,NetCDF=.TRUE.)` already works (test_Simple, test_Aerosol_Bypass already use it with `.nc`) | flip the driver calls only |
| Atmosphere (`.Atmosphere.bin`) | K-matrix drivers | **NONE** | write netCDF I/O from scratch |
| Surface (`.Surface.bin`) | K-matrix drivers | **NONE** | write netCDF I/O from scratch |

## Authoritative facts gathered (grep-verified)

### The template to mirror: `src/RTSolution/CRTM_RTSolution_Define.f90`
- Public `CRTM_RTSolution_WriteFile` (~L2194) routes:
  `binary=.TRUE.; IF(PRESENT(NetCDF)) binary=.NOT.NetCDF; IF(binary) ..._Binary ELSE ..._NetCDF`.
- Private workers: `CRTM_RTSolution_WriteFile_NetCDF` (~L2441), `..._ReadFile_NetCDF` (~L1587),
  `..._InquireFile_NetCDF` (~L1106), plus a `CreateFile_netCDF` helper (defines dims+vars).
- netCDF idiom: `USE netcdf`; pattern is `NF90_INQ_VARID` then `NF90_PUT_VAR`/`NF90_GET_VAR`,
  each followed by an `IF (NF90_Status /= NF90_NOERR) THEN msg=...; CALL Write_Cleanup(); RETURN`.
  Cleanup subroutine closes the file if `Close_File` and sets `err_stat=FAILURE`.
- Dimension/variable names are module CHARACTER PARAMETERs (e.g. `PROFILE_DIMNAME`,`LAYER_DIMNAME`,
  `CHANNEL_DIMNAME`,`STOKE_DIMNAME`, `CHANNEL_VARNAME`, etc.). Define analogous ones for Atm/Sfc.

### Dispatch surface for Surface (`src/Surface/CRTM_Surface_Define.f90`, md5 b1876c2..., committed blob 999eb2f)
- `CRTM_Surface_WriteFile` and `CRTM_Surface_ReadFile` are **generic INTERFACEs**:
  - `INTERFACE CRTM_Surface_WriteFile` ⇒ `Write_Surface_Rank1` (L1374), `Write_Surface_Rank2` (L1468)
  - `INTERFACE CRTM_Surface_ReadFile`  ⇒ `Read_Surface_Rank1` (L1061),  `Read_Surface_Rank2` (L1183)
  - `CRTM_Surface_InquireFile` is a single FUNCTION (L910).
- Binary record workers: `Write_Record`/`Read_Record` (near end; `END MODULE` at L2508).
- Module header USEs (L21-44): `Type_Kinds:fp`, `Message_Handler:SUCCESS,FAILURE,WARNING,INFORMATION,Display_Message`,
  `File_Utility:File_Open,File_Exists`, `CRTM_SensorData_Define` (SensorData type + Create/Destroy/Associated/Zero).
  **Must ADD `USE netcdf`.** `INTEGER,PARAMETER :: ML=256` exists (L135). No SL param — SensorData uses `STRLEN`.
- The test baselines call **rank-2** only: `Surface_K(n_Channels, N_PROFILES)`.
- `OPERATOR(==)`/`CRTM_Surface_Compare` exist and are PUBLIC — USE them in the round-trip test.

### `CRTM_Surface_type` fields to serialize (all SCALAR per element + nested SensorData)
REAL(fp): Land_Coverage, Water_Coverage, Snow_Coverage, Ice_Coverage, Land_Temperature,
Soil_Moisture_Content, Canopy_Water_Content, Vegetation_Fraction, Soil_Temperature, LAI,
Water_Temperature, Wind_Speed, Wind_Direction, Salinity, Snow_Temperature, Snow_Depth,
Snow_Density, Snow_Grain_Size, Ice_Temperature, Ice_Thickness, Ice_Density, Ice_Roughness.
INTEGER: Land_Type, Soil_Type, Vegetation_Type, Water_Type, Snow_Type, Ice_Type.
LOGICAL: Is_Allocated (always .TRUE.; can skip or store as int).
Nested `TYPE(CRTM_SensorData_type) :: SensorData`: n_Channels(INT), Sensor_Id(CHAR STRLEN),
WMO_Satellite_ID(INT), WMO_Sensor_ID(INT), and IF n_Channels>0: Sensor_Channel(:)(INT), Tb(:)(fp).
**In ALL test baselines SensorData%n_Channels==0** (no driver sets it) — so SensorData arrays need
not round-trip for the tests, but store n_Channels (and handle >0 for correctness if cheap).

### `CRTM_Atmosphere_type` fields (`src/Atmosphere/CRTM_Atmosphere_Define.f90`, END MODULE; binary Read_Record L2474 / Write_Record L2607; InquireFile L1508)
Scalars/dims: n_Layers, n_Absorbers, n_Clouds, n_Aerosols, n_Added_Layers, Climatology(INT), Add_Extra_Layers(LOG).
Arrays: Absorber_ID(:)(INT J), Absorber_Units(:)(INT J), Level_Pressure(0:K), Height(0:K),
Pressure(K), Temperature(K), Absorber(K,J), Cloud_Fraction(K), Relative_Humidity(K).
Nested **populated** in baselines (so MUST round-trip):
- `Cloud(:)` — test_Simple N_CLOUDS=1, test_ScatteringSwitch=2, test_AOD=0. `CRTM_Cloud_type`
  (`src/Atmosphere/Cloud/CRTM_Cloud_Define.f90`): Type(INT), n_Layers(INT), plus 4 allocatable
  REAL(fp)(K): Effective_Radius, Effective_Variance, Water_Content, Water_Density.
- `Aerosol(:)` — test_AOD N_AEROSOLS=3, test_Simple=1, test_ScatteringSwitch=2. `CRTM_Aerosol_type`
  (`src/Atmosphere/Aerosol/CRTM_Aerosol_Define.f90`): Type(INT), n_Layers(INT), 3 allocatable
  REAL(fp)(K): Effective_Radius, Effective_Variance, Concentration.
Use the existing `CRTM_Cloud_Create`/`CRTM_Aerosol_Create`/`CRTM_Atmosphere_Create` to allocate on read.

### Driver files to convert (the hard half — Atmosphere+Surface, all rank-2, all K-matrix/adjoint)
`test/mains/regression/k_matrix/{test_Simple,test_ClearSky,test_SOI,test_SSU,test_AOD,test_Zeeman,
test_User_Emissivity,test_ChannelSubset,test_ScatteringSwitch,test_VerticalCoordinates}/...f90`,
`test/mains/regression/adjoint/{test_Simple,test_ClearSky}/...f90`,
`test/mains/unit/Unit_Test/{test_Aerosol_Bypass_k_matrix,test_Aerosol_Bypass_adjoint}.f90`.
RTSolution-only drivers (easy half) = the rest of the ~27 that call `CRTM_RTSolution_WriteFile`.
Baseline filename literals in drivers: `'.Atmosphere.bin'`, `'.Surface.bin'`, `'.RTSolution.bin'`,
`'.RTSolution_K.bin'`, `'.RTSolution_AD.bin'`, `'.RTSolution_TL.bin'` (+ `.GFS./.NAM.` variants in
test_VerticalCoordinates). The drivers also `InquireFile` then `ReadFile` the same path.

## Recommended implementation order (each phase gated on clean build + test)
1. **Surface netCDF I/O** (simplest: all-scalar). Add `_NetCDF` workers for rank-2 (and rank-1 for
   interface symmetry; rank-1 unused by tests but the generic interface members should share the
   optional `NetCDF` arg signature), add optional `NetCDF` to InquireFile, route. Suggested schema:
   dims `n_Channels`,`n_Profiles`; each scalar field a 2-D var (n_Channels×n_Profiles); SensorData
   n_Channels as a 2-D int var (skip the arrays while max==0, or implement padded 3-D if doing it fully).
2. **Round-trip unit test** `test/mains/unit/Unit_Test/test_Surface_netCDF_io.f90`: build Surface(3,2)
   with varied nonzero values, write `NetCDF=.TRUE.`, read back, compare with `==`/`CRTM_Surface_Compare`.
   Print `SURFACE NETCDF ROUNDTRIP PASS`/STOP 0. Register in `test/CMakeLists.txt`.
3. **Atmosphere netCDF I/O** (incl. Cloud/Aerosol). Schema needs an `n_Layers` dim, `n_Absorbers` dim,
   and per-(channel,profile) cloud/aerosol counts; the variable-length Cloud/Aerosol arrays are the
   tricky part — store Max_Clouds/Max_Aerosols dims and the per-cell n_Clouds/n_Aerosols, pad arrays.
   (Because results are written+read by the same build, you have total schema freedom; only requirement
   is exact round-trip.) Round-trip test likewise.
4. **Flip all baseline drivers** to `NetCDF=.TRUE.` + `.nc` extensions (Write, Read, Inquire calls).
   Mirror exactly how `regression/forward/test_Simple` already does RTSolution netCDF (lines ~215-290).
5. **Hard switch:** delete stale `build/test/results/**/*.bin`, reseed (run ctest once to create `.nc`
   references), run ctest again to confirm compare-mode passes. Then **clean from-scratch** build+ctest
   to certify: must be 206/206 (or 206/206 + any new round-trip tests registered).
6. Update `REL-3.2.0_changes_vs_develop.md` (untracked working-tree doc) with a bullet on the baseline
   format change.

## Acceptance criteria
- `rm -rf build && cmake .. -DCMAKE_BUILD_TYPE=RELEASE && make -j8 && ctest` ⇒ all pass, 0 build errors,
  observed via exit codes (not assumed).
- No `.bin` baseline files remain under `build/test/results/`.
- New Surface/Atmosphere round-trip unit tests registered and passing.
- The two new netCDF I/O implementations mirror the RTSolution template's error-cleanup idiom.

## Commit guidance
Separate commits: (a) Surface netCDF I/O + test, (b) Atmosphere netCDF I/O + test, (c) driver flip +
reseed. Commit each only after an observed clean-build + passing test. Branch policy: commit on
`feature/btj_REL-3.2.0`; user pushes manually.
