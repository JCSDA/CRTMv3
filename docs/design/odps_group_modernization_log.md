# ODPS group modernization: implementation log

Branch: `feature/btj_odps_group_modernization` (off `feature/btj_uv_no2_group8`).
Purpose: running record of exactly what changed and why, kept current as Tiers
0 through 2 land. This file is the source for the developer-facing GitHub issue
to be posted on JCSDA/CRTMv3 after Tier 2 is complete. Design rationale lives in
`odps_group_modernization.md`; forensic background in
`crtm-coeffgen/docs/omps_group4_odps_group_assessment.md`.

## Tier 0: validation and guardrails

Status: implemented; full-suite verification in progress.

What changed, file by file:

1. `src/AtmAbsorption/ODPS/ODPS_Predictor.f90`
   - New PUBLIC parameters `RESERVED_ZSSMIS_GROUP = 4` and
     `RESERVED_ZAMSUA_GROUP = 5`, declared beside the group table with a
     comment making this module the single owner of the ODPS group index
     space (previously the reserved indexes existed only as magic numbers in
     `ODZeeman_Predictor.f90`).
   - New PUBLIC parameter `VALID_GROUPS = (/1, 2, 3, 7, 8/)`.
   - New PUBLIC function `ODPS_Validate_Group(Group_Index, Component_ID,
     Absorber_ID, Message)`: verifies the group is supported and that the
     file's component and absorber rosters match the group's compiled maps
     exactly, in content and order (the predictor code addresses components
     positionally). Returns a specific, self-explanatory message for:
     Zeeman-reserved indexes, unknown indexes, roster size mismatches, and
     roster content/order mismatches.
   - Component-ID registry documented in place: the IDs are heritage
     transmittance-production "molecule set" numbers (12 wet, 13 dry, 14
     ozone, 15 wco, 101 effective-mol1, 113 effdry, 114 effozo, 118 to 121
     effective trace gases, 122 NO2 CRTM extension). New IDs must be
     coordinated with crtm-coeffgen.
   - The `HUGE()` poison values returned by `ODPS_Get_Component_ID` and
     `ODPS_Get_Absorber_ID` for out-of-range queries are replaced by a named
     PUBLIC sentinel `ODPS_INVALID_ID = -1` (also now used by
     `ODPS_Get_Ozone_Component_ID`).

2. `src/AtmAbsorption/ODZeeman/ODZeeman_Predictor.f90`
   - `ODPS_gINDEX_ZSSMIS` and `ODPS_gINDEX_ZAMSUA` are now derived from the
     `RESERVED_*` constants via `USE ODPS_Predictor` instead of hardcoding 4
     and 5 independently. Values unchanged; single source of truth.

3. `src/Coefficients/CRTM_TauCoeff.f90`
   - After the ODPS container load, every `TC%ODPS(i)` is validated with
     `ODPS_Validate_Group`; failure aborts `CRTM_Load_TauCoeff` with the
     sensor name and the specific reason. Before this change, a bad file
     surfaced later as a generic "Error allocating predictor structure" from
     `CRTM_Predictor_Create` with no sensor context (this is exactly how the
     OMPS Group-4 files failed).
   - The same validation runs over every nested ODPS sub-structure of each
     loaded ODSSU sensor (guarded on the ODPS pointer being associated, since
     ODSSU can instead carry ODAS sub-structures).
   - Zeeman companion files are exempt by construction: they load into the
     separate ODZeeman container, which this check does not touch.

4. `test/mains/unit/Unit_Test/test_ODPS_Group_Validation.f90` (new) plus
   CMake registration as `test_Unit_ODPS_Group_Validation`: 17 cases covering
   all five canonical rosters (accepted), the three Zeeman-reserved indexes,
   unknown/fill/negative indexes, and size, content, and order mismatches
   (all rejected with non-blank messages). Includes the literal OMPS roster
   (Group 4, Component_ID [13, 14], Absorber_ID [3]) as a must-reject case.

Compatibility statement: no coefficient file changes. A pre-change sweep of
all 568 ODPS structures (including ODSSU sub-blocks) in fix_REL-3.2.0.0
confirmed every group 1/2/3/7/8 roster matches its compiled map exactly, so
the new validation rejects only the nine Group-4 files, which never loaded
successfully anyway (four OMPS: clean failure replaces cryptic failure; five
zssmis: unaffected, they use the ODZeeman path).

## Tier 1: group registry and per-component predictor kernels

Status: in design; key verified facts recorded here for the implementation.

Verified structure of the current monoliths (needed to preserve bit-identical
results):

- FWD `ODPS_Compute_Predictor` (driver at 739; contained IR at 869, MW at
  1127 by pre-change numbering): driver computes integrated variables (Tz,
  Tzp, GAz, GAzp, GATzp) in a layer loop sized by `N_ABSORBERS_G(Group_ID)`,
  then dispatches IR (groups 1, 2, 8) or MW (groups 3, 7). Within the IR
  layer loop, component blocks execute in order: DRY, WCO, OZO, CO2, WLO,
  then the group-1 extension block (CO2 predictor 11, WLO predictors 16-18,
  CO, CH4, N2O), then NO2 (group 8). Forward X assignments are independent
  (no cross-accumulation), so kernel execution order is bit-irrelevant in
  FWD and TL.
- TL mirrors FWD exactly (value + derivative pairs); NO2 TL block sits after
  the group-1 block. MW TL has two quirks to preserve verbatim: it zeroes
  `X(:, 8:, COMP_DRY_MW)` inside the layer loop, and one H2OdH2OTzp_TL term
  indexes `GATzp_TL(k, ABS_H2O_IR)` (harmless: both constants are 1).
- AD is NOT a strict reverse of FWD and its accumulation order determines
  bit-level results. Verified per-layer order in IR_AD: recompute forward
  scalars; then adjoint blocks DRY, WCO, NO2 (group 8, self-contained
  through to Absorber_AD), OZO, CO2 (1-10), WLO (1-15), group-1 extension
  (WLO 16-18 AND CO2 predictor 11 together), CO, CH4, N2O; then the group-1
  variable-chain epilogue (CH4, CO, N2O chains to Absorber_AD); then the
  common epilogue (O3 chain, H2O chain, DT2/T2 chains, absorber and
  temperature propagation). MW_AD: DRY, WET, then a fully self-contained
  OZO block (group 7), then the common epilogue. Any kernel decomposition
  must call AD kernels in exactly this order.
- Group variants are encoded by per-component predictor counts (WLO 18 in
  group 1 vs 15 in group 2/8; CO2 11 vs 10), so kernels can key their
  extensions on the roster's count instead of the group index.

Planned decomposition (bit-identity preserving): keep the single layer loop
per basis; replace the interleaved IF(Group_ID) blocks with per-component
kernel calls in the monolith's exact statement order; kernels are contained
subroutines host-associating the shared per-layer scalars.

## Tier 2: file-driven allocation and dispatch, Zeeman dual-acceptance

Status: pending.
