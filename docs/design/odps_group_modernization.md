# Evaluation: modernizing the ODPS transmittance group system

Date: 2026-07-26
Status: evaluation / proposal, no code changes.
Companion document: `crtm-coeffgen/docs/omps_group4_odps_group_assessment.md`
(the Group_Index=4 forensic assessment; this document builds on its verified
findings and does not repeat the history).

## 1. Purpose

The ODPS "group" system works, but it has accreted since 2008 and its recent
growth (GROUP_MW_O3, GROUP_UV_NO2) plus the OMPS Group-4 incident exposed how
cumbersome and ad hoc it has become. This document inventories the actual design
debt, separates it from the parts that are sound, and proposes a tiered
modernization path with effort and risk called out per tier.

## 2. What a "group" actually is today

A Group_Index bundles four orthogonal things into one hardcoded integer:

1. A shared predictor basis class: the IR-style basis (groups 1, 2, 8) or the
   MW-style basis (groups 3, 7), selected in `ODPS_Compute_Predictor`
   (`ODPS_Predictor.f90:823-828` and the TL/AD copies at 1378 and 2013).
2. A component roster: which effective-transmittance components exist and in
   what order (`COMPONENT_ID_MAP_G*`, `N_COMPONENTS_G`).
3. An absorber roster: which variable gases the coordinate mapping must supply
   (`ABSORBER_ID_MAP_G*`, `N_ABSORBERS_G`).
4. A per-component predictor recipe: which formulas, and how many, feed each
   component (`N_PREDICTORS_G*` plus roughly 2,000 lines of positional
   assignments, three times over for FWD, TL, and AD).

Because all four axes are welded to one index, adding one scene gas to an
existing sensor class has each time required a whole new group: MW plus ozone
became group 7; UV/VIS (group-2 basis) plus NO2 became group 8. The next such
need (IR plus NO2? MW plus N2O?) means another group, another set of parallel
arrays, and another pass through three copies of the predictor code.

## 3. What is sound and should be kept

- The optical-depth application side is already data-driven. `ODPS_AtmAbsorption.f90`
  loops components using the file's own per-channel `n_Predictors(j, channel)`
  and `Pos_Index(j, channel)` (lines 172-198, 359-386, 563-588). Nothing there
  cares about the group.
- Absorber gathering is already data-driven. `ODPS_CoordinateMapping.f90:246`
  matches `Atm%Absorber_ID` against the file's `Absorber_ID` list, not against
  the group table.
- The physics must be compiled code. Predictor functional forms cannot be
  meaningfully "data-driven" without an expression interpreter, which is the
  wrong trade for an operational Fortran library. Any modernization keeps
  compiled predictor kernels; the question is how they are organized and
  dispatched.
- The ODPS file format itself is adequate. It already stores `n_Components`,
  `Component_ID(:)`, `Absorber_ID(:)`, and per-channel predictor counts, which
  is precisely the metadata a modern dispatch would want. No format change is
  required for most of what follows.

## 4. Design debt inventory

Each item is verified against source on branch `feature/btj_uv_no2_group8`.

D1. Two sources of truth, never cross-checked. The file carries the component
    and absorber rosters; the library carries compile-time tables keyed by
    Group_Index; nothing ever verifies they agree. `Component_ID` from the file
    is read, stored, and never consulted by the runtime (the predictor code is
    purely positional). A file whose component count happens to match a group
    but whose components mean something else would compute silently wrong
    radiances. The OMPS Group-4 files failed loudly only because the group-4
    table entry is zero; with a nonzero entry they would have "worked".

D2. No load-time validation surface. An unsupported Group_Index is discovered
    as a generic "Error allocating predictor structure" deep inside
    `CRTM_Predictor_Create`, with no mention of the sensor, the group, or what
    would be supported. The OMPS diagnosis took a forensic session largely
    because of this.

D3. Magic reserved indices and split ownership. Groups 4 to 6 are zero-padded
    holes in the ODPS table (`ODPS_Predictor.f90:93-95`) whose meaning
    (`ODPS_gINDEX_ZSSMIS = 4`, `ODPS_gINDEX_ZAMSUA = 5`) lives in a different
    module (`ODZeeman_Predictor.f90`), while the actual Zeeman routing decision
    is made in a third place by WMO sensor ID (`CRTM_TauCoeff.f90:633,655`).
    The Group_Index in a Zeeman file is decorative; the Group_Index in an ODPS
    file is load-bearing. Same field, two unrelated meanings.

D4. Parallel-array table growth. A new group touches `N_G`, three dimension
    arrays, a `GROUP_*` constant, a `COMPONENT_ID_MAP_G*`, an
    `ABSORBER_ID_MAP_G*`, an `N_PREDICTORS_G*`, the accessor SELECT CASE blocks
    (`ODPS_Predictor.f90:3506-3536`), and the compute dispatch in three places.
    That is roughly ten scattered edit sites for what is conceptually one row of
    a table.

D5. FWD/TL/AD triplication at the wrong granularity. The three copies of the
    predictor computation are monolithic per basis class (one giant contained
    subroutine each), so every roster change re-touches all three monoliths.
    Triplication of formulas is inherent to hand-written TL/AD; triplication of
    roster and dispatch logic is not.

D6. Undocumented identifier conventions. Component IDs are heritage
    molecule-set numbers (12 = wet, 13 = dry, 14 = ozone, 101/113/114 =
    effective variants) documented only in a generation-side parameter file in
    another repository. This directly enabled the OMPS confusion, and the
    library exports accessors (`ODPS_Get_Component_ID`, `ODPS_Get_Absorber_ID`,
    `ODPS_Get_Ozone_Component_ID`) that no code outside the module calls, with a
    CASE DEFAULT that returns `HUGE()` in the hope of inducing a downstream
    failure (`ODPS_Predictor.f90:3517-3518`).

D7. Generator mirroring. `crtm-coeffgen` re-implements the group tables in
    Python (`odps/predictors.py`, `taucoeff/netcdf_io.py`) and must be kept in
    lockstep with the Fortran by hand. There is no shared, machine-readable
    definition of what a group means.

One honest complication to record: components are not fully independent within
a group. Group 1's water-line component includes CH4 and CO cross-terms as
predictors 16 to 18, and its CO2 component includes a CO term as predictor 11
(`ODPS_Predictor.f90:1013-1020`). The same Component_ID (101, WLO) means an
18-predictor recipe in group 1 and a 15-predictor recipe in group 2. So
predictor kernels cannot be keyed on Component_ID alone; they key on
(Component_ID, basis class), and a few kernels consume other components' gas
variables. Any refactor must model this honestly rather than pretending
components are independent.

## 5. Constraints

- Bit-for-bit reproducibility of existing sensors is the acceptance criterion
  for any refactor tier below Tier 3. The 200+ ctest baselines plus the
  machine-precision TL/AD parity tests (for example test_UV_NO2_TLAD) are the
  safety net.
- TL and AD stay hand-written (project convention); the goal is to shrink what
  must be written by hand, not to introduce automatic differentiation.
- The binary and netCDF TauCoeff formats keep working unchanged through Tier 2.
- Fortran 2008, no external dependencies.

## 6. Modernization tiers

### Tier 0: hygiene and guardrails (small; no numerical change possible)

1. Load-time validation in the TauCoeff loader: when an ODPS file is loaded,
   verify Group_Index is a supported, non-reserved group AND that the file's
   n_Components, Component_ID list, and Absorber_ID list exactly match the
   group's compiled maps. On mismatch, fail CRTM_Init with the sensor name, the
   offending values, and the supported set. This closes D1 and D2 outright: the
   OMPS files would have produced a one-line, self-explanatory error.
2. Reserved-index ownership: move the Zeeman gINDEX constants (or at minimum a
   mirror of them, compile-checked) into ODPS_Predictor next to the table, so
   one module states the whole index space. Replace the zero-padding comment
   with named RESERVED entries.
3. Document the component-ID registry (molecule-set heritage, the effective +100
   convention, who assigns new IDs) in ODPS_Predictor.f90 and in crtm-coeffgen
   docs. NO2 = 122 already follows the pattern; write the pattern down.
4. Delete or wire up the dead accessors; replace the HUGE() sentinel with an
   explicit error path.

This tier is shippable in a 3.2.x timeframe and is worth doing regardless of
any later tier.

### Tier 1: registry plus component kernels (medium; the substantive cleanup)

1. Replace the parallel arrays with a single derived-type registry:

   ```fortran
   TYPE :: ODPS_Group_Spec
     CHARACTER(16) :: Name              ! 'IR_HIRES', 'UV_NO2', ...
     INTEGER       :: Basis             ! BASIS_IR or BASIS_MW (or RESERVED_ZEEMAN)
     INTEGER       :: n_Components
     INTEGER       :: Component_ID(MAX_COMPONENTS)
     INTEGER       :: n_Absorbers
     INTEGER       :: Absorber_ID(MAX_ABSORBERS)
     INTEGER       :: n_Predictors(MAX_COMPONENTS)
   END TYPE
   TYPE(ODPS_Group_Spec), PARAMETER :: GROUP_REGISTRY(N_G) = [ ... ]
   ```

   All accessors become one bounds-checked lookup. Adding a group becomes one
   registry row. D4 closed.
2. Break the per-basis monoliths into per-component kernel subroutines (FWD, TL,
   AD triples): dry_ir, wlo, wco, ozo, co2, n2o, co, ch4, no2, edry_mw, wet_mw,
   ozo_mw. The shared layer quantities (T, DT, gas ratios, secant products) stay
   in a basis-level "context" derived type computed once per layer sweep and
   passed to kernels, which resolves the cross-component coupling cleanly (the
   WLO group-1 kernel variant reads CH4_A and CO_A from the context). Kernel
   selection is a loop over the registry roster dispatching on
   (Component_ID, Basis), replacing the interleaved IF(Group_ID == ...) blocks.
   D5 reduced to its irreducible core (formula triples).
3. Positional layout stays exactly as is (component k of the file is roster
   entry k, enforced by the Tier 0 validation), so results are bit-identical and
   the file format is untouched.

Payoff test: after Tier 1, "add scene gas X to class Y" is one registry row,
one kernel triple (if the gas is new), and the generator change. Today the same
change is a new group threaded through roughly ten sites and three monoliths
(this is precisely what the GROUP_MW_O3 and GROUP_UV_NO2 efforts had to do).

### Tier 2: allocate and dispatch from the file (medium-large; optional beyond Tier 1)

1. Allocate the predictor structure from the file's own dimensions
   (n_Components, MAXVAL of per-channel n_Predictors, n_Absorbers) instead of
   the registry, keeping the registry as the validator and kernel-selector.
2. Dispatch kernels by the file's Component_ID sequence rather than by roster
   position. At that point Group_Index degrades to a basis tag plus a validation
   convenience, and a well-formed file whose components are all known kernels
   loads even if its exact roster combination was never enshrined as a group.
   A future dry+ozone-only UV file (the physics the OMPS files actually
   contain) would load without anyone defining a group for it.
3. This is the right moment to also normalize the Zeeman side-load: select the
   ODZeeman path from file metadata (Algorithm or a sub-algorithm field)
   instead of hardwired WMO IDs, removing the last out-of-band dispatch (D3
   fully closed). Requires a coordinated coefficient re-release or a
   dual-acceptance window, so schedule it with a release, not a patch.

### Tier 3: schema and cross-algorithm unification (large; v4 material)

Per-component recipe version tags in the TauCoeff schema; one strategy-style
interface over ODAS, ODPS, ODSSU, and ODZeeman; shared machine-readable group
and component registry consumed by both the library and crtm-coeffgen (D7
closed at the root, for example a small JSON or Fortran namelist shipped with
the fix set and code-generated into both languages). Worth designing only as
part of a deliberate v4 / TauCoeff-Release-bump effort.

## 7. Recommendation

1. Do Tier 0 now. It is small, cannot change numbers, and converts the entire
   OMPS class of failure into a one-line diagnostic. It also documents the
   conventions while the forensic context is fresh.
2. Adopt Tier 1 as the modernization target for the next minor release cycle.
   It is a pure refactor with bit-identical acceptance against the existing
   baselines, and it is where the "cumbersome and ad hoc" feeling actually
   lives: ten edit sites and three monoliths per group become one registry row
   and a kernel set.
3. Treat Tier 2 as opportunistic follow-on once Tier 1 has soaked, with the
   Zeeman de-gating explicitly scheduled alongside a coefficient release.
4. Defer Tier 3 to a v4 conversation.

## 8. Validation strategy for Tiers 1 and 2

- Bit-identical forward, TL, AD, and K baselines across all ODPS sensors in the
  suite (groups 1, 2, 3, 7, 8 are all exercised by current ctests).
- The existing machine-precision TL/AD/K parity tests re-run unchanged.
- New negative tests: a synthetic TauCoeff with a reserved Group_Index, one
  with a Component_ID/roster mismatch, and one with an unknown Component_ID,
  each asserting the specific CRTM_Init error message (Tier 0 deliverable).
- OpenMP consistency test unchanged (registry is PARAMETER data; kernels are
  pure; no new shared state).
