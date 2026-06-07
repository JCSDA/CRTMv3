# Design Plan: First-Class Downwelling Radiance in CRTM

**Baseline:** `feature/btj_REL-3.2.0` @ 6898e2c
**Branch:** `feature/btj_downwelling_radiance`
**Goal:** Replace the `Obs_4_downward_P` option hack with an always-on, fully-differentiated
(Forward / Tangent-Linear / Adjoint / K-Matrix) downwelling-radiance output.

**Locked scope (per review):**
- Surface scalar `Down_Radiance` across **all** solvers (Emission, ADA, SOI), always-on, fully differentiated.
- **Retire** `Obs_4_downward_P` once the formal output covers the need.
- **No** level-resolved profile field for now (deferred — see §10 R4).
- Phases delivered: **0 → 1 → 2 → 4** (Phase 3 profile dropped).

> ⚠️ Consequence of "surface only + retire option": `Down_Radiance` is strictly the surface
> value (`e_Level_Rad_DOWN(n_Layers)`). The retired `Obs_4_downward_P` gave downwelling at an
> arbitrary pressure level. Arbitrary-level downwelling returns only with the deferred profile
> (Phase 3). Retired because the option is forward-only and inert in TL/AD/K — nothing in a DA
> pipeline can depend on it. (CORRECTION: an earlier draft called it "not publicly settable" — it
> was in fact a *public* field of the public `CRTM_Options_type`, set directly by the old
> `test_Downwelling_Radiance`, so its removal is a breaking API change. See Phase 4.)

---

## 1. Recommendation

Make downwelling radiance a first-class RTSolution output that **mirrors the existing upwelling
outputs** (`Up_Radiance` scalar, `Upwelling_Radiance(:)` profile). Complete the symmetry with
`Down_Radiance` (scalar, at surface — *field already exists and is fully plumbed*). The surface
downwelling value is **already computed and already differentiated** inside the surface-reflection
term every solver evaluates (`R·I_down`); Phase 1 simply routes that existing sensitivity out
through the existing field. No new physics, no API/IO churn for the scalar.

The AD/K **seed** for the downwelling Jacobian is `RTSolution_{AD,K}%Down_Radiance` (uniformly).
`%Radiance` always means TOA upwelling.

## 2. Why `Obs_4_downward_P` is wrong (not just inelegant)

(file:line on baseline)
- Overwrites `RTSolution%Radiance`/`%Stokes`/`%Brightness_Temperature` with the downwelling value
  (`Common_RTSolution.f90:1186,1201,1222-1244`) — destroys TOA output for that channel.
- Forward-only; **inert in TL/AD/K** (TL/AD modules have zero references; K reads
  `RTV%obs_4_downward%rt` at `CRTM_K_Matrix_Module.f90:1417` but nothing sets it). Proven
  numerically: FD-vs-TL ratio ≈ 0.002, K Jacobian ≈ 0.
- SOI silently broken (`s_Level_Rad_DOWN` filled only by ADA, gated).
- Magic sentinel (`>0`, default `-1`), nearest-level pressure snap, half-plumbed Option
  (missing from SetValue/Inspect/Equal/IO), duplicates the aircraft machinery.

## 3. Assets already in place

- `RTSolution%Down_Radiance` field (`CRTM_RTSolution_Define.f90:240`), fully plumbed
  (IO/Inspect/Compare/Equal/arithmetic, `DWNR_VARNAME` at :150), and **already set every forward
  emission call** to surface downwelling (`Common_RTSolution.f90:1206`). → scalar output needs
  **no new field, no netCDF changes** — only differentiation.
- Emission computes the full downwelling profile unconditionally
  (`Emission_Module.f90:116-135`); its TL/AD already compute the surface downwelling sensitivity
  as a running scalar for the surface-reflection term (TL `:229-249`, AD `:368-390`).
- Verification harness already drafted: `test/mains/unit/Unit_Test/test_Downwelling_TLADK.f90`
  (FD-vs-TL, adjoint dot-product, K-vs-AD; TOA control + downwelling).

## 4. Design decisions

- **D1** Output via RTSolution fields, not an option, not a second structure. Symmetric with upwelling.
- **D2** Emission (clear-sky) `Down_Radiance` is always-on (cheap). Scattering (ADA/SOI)
  `Down_Radiance` is **opt-in** via `Options%Compute_Down_Radiance` (default .FALSE.) because the
  adding-doubling downward sweep adds real per-call cost. (Revised after the perf finding; the
  earlier "always-on for all solvers" intent is kept only for the cheap emission path.)
- **D3** AD/K seed via `%Down_Radiance` uniformly (resolves today's K-vs-AD seed inconsistency).
  `%Radiance` stays TOA upwelling.
- **D4** Differentiate by exposing existing intermediates, not new physics.
- **D5** Retire `Obs_4_downward_P`, `RTV%obs_4_downward`, the `%Radiance` overwrite, the ADA gate.

## 5. Coupling with the surface-reflection-correction work

The reflected-downwelling surface term is literally `R·I_down`:
- clear sky:  `reflectivity(1,1)·e_Level_Rad_DOWN(nL)` (`Emission_Module.f90:142`)
- scattering: `matmul(s_Level_Refl_UP(:,:,nL)=reflectivity, …)` (`ADA_Module.f90:152,216`)

The grazing-angle guards (`1c4d4ce`, `e662b02`, `6898e2c`) clamp `R` to a physical band and split
its TL/AD per-Stokes via `iVar%Reflectivity_Clamped` (`CRTM_FastemX.f90`, `CRTM_PARMIO*`).

- **C1** Exposing/differentiating `I_down` means the TL/AD of `R·I_down` must use the same clamp
  flags; a clamped Stokes component uses `d(1-emissivity)` not `dR`. The downwelling-output
  Jacobian must follow that branch or TL≠ADᵀ at grazing.
- **C2** Add a grazing-angle case to the downwelling verification (reuse `test_Grazing_SfcOptics`
  angles za=82/84/86°).
- **C3** Angular/profile downwelling (if ever added) runs on the same Gauss streams (~86°) — inherit
  the guard's conditioning protection.

## 6. Per-solver cost (drives phasing)

| Solver | FWD surface | TL/AD surface |
|---|---|---|
| Emission | free (have it) | exists as a running scalar — route it out (LOW) |
| ADA | extract surface value (currently behind the gated `DO 20`) | **new** TL/AD (HIGH — dominant effort) |
| SOI | accumulate `s_Level_IterRad_DOWN` over orders | **new** TL/AD accumulation (MODERATE) |

## 7. Phased implementation

### Phase 0 — Acceptance harness first
`test_Downwelling_TLADK.f90`: FD-vs-TL (central diff), adjoint dot-product
`⟨TL·x,TL·x⟩=⟨x,AD·TL·x⟩`, K-vs-AD equality; run for a TOA control + downwelling, seeding
downwelling via `%Down_Radiance`. Clear-sky + (later) scattering + grazing cases. Defines
"correct" before source changes; fails on baseline, passes after each phase.

### Phase 1 — Emission (clear-sky) surface `Down_Radiance`, differentiated  *(low risk)*
- FWD: already populated (`Common_RTSolution.f90:1206`).
- TL: expose `e_Level_Rad_DOWN_TL(nL)` from `CRTM_Emission_TL`; set `RTSolution_TL%Down_Radiance`
  in `Assign_Common_Output_TL`.
- AD: seed `down_rad_AD` from `RTSolution_AD%Down_Radiance` in `Assign_Common_Input_AD`; feed
  `CRTM_Emission_AD` surface adjoint.
- K: inherits AD path; seed `RTSolution_K%Down_Radiance`. Extend `Pre_Process_RTSolution_AD/_K`
  with a `%Down_Radiance` seed path that does not collide with the `%Radiance`/BT seed.
- Gate on Phase-0 clear-sky test; regenerate clear-sky TL/AD baselines.

### Phase 1 — OUTCOME (implemented)
- Emission surface `Down_Radiance` now differentiated FWD/TL/AD/K. Verified by
  `test_Unit_Downwelling_TLADK` at machine precision (clear-sky atms_npp): FD-vs-TL 5.7e-11,
  adjoint dot-product 2.4e-16, K-vs-AD exact; TOA control likewise exact.
- Touch points: `Emission_Module.f90` (optional `down_rad_TL_out` / `down_rad_AD_in`),
  `CRTM_RTSolution.f90` (`CRTM_Compute_RTSolution_TL` sets `RTSolution_TL%Down_Radiance`;
  `_AD` extracts/zeros `RTSolution_AD%Down_Radiance` seed and feeds `CRTM_Emission_AD`).
- Known, expected fallout: `test_tangent_linear_Simple_{atms_n21,amsua_n19}` now fail ONLY on
  `Down_Radiance` (previously TL=0; now the correct nonzero TL). These are CLOUDY scenes, so the
  baselines are intentionally NOT regenerated here — see the clear/cloudy combine in Phase 2.
  Forward / adjoint / k_matrix Simple all still pass.

### Phase 2 — ADA + SOI surface `Down_Radiance`, differentiated  *(high effort — DONE & VERIFIED)*
- ADA FWD: **DONE** (commit 64b17db). Un-gated the downward sweep so `s_Level_Rad_DOWN` is always
  computed; scattering branch of `Assign_Common_Output` sets `Down_Radiance` from the finalized
  surface value; clear/cloudy `Total_Cloud_Cover` combine added. Verified TOA upwelling unchanged
  (BT delta = 0). FWD scattering `Down_Radiance` not yet FD-verified (needs ADA TL).
  - PERF: un-gating runs the adding-doubling downward sweep (2 matinv/layer: Inv_Gamma2 +
    Inv_Gamma3) on EVERY scattering RT call. Acceptable per "always-on", but a surface-only
    optimization (full per-layer Inv_Gamma2 recursion, Inv_Gamma3 finalization only at the surface
    level) would roughly halve the added cost. Deferred; correctness first.
- Opt-in flag: **DONE**. `Options%Compute_Down_Radiance` + `RTV%Compute_Down_Radiance`, plumbed in
  all four drivers (FWD/TL/AD/K), gates the ADA downward sweep and the scattering `Down_Radiance`
  assignment. Default off → existing scattering results unchanged. Full suite green (212/212);
  emission TL baselines regenerated (clear-sky, FD-verified). (A stale aerosol baseline from a
  cross-branch shared build/ was also regenerated — not a regression; this branch has no
  scatter/aerosol source changes.)
- SOI FWD/TL/AD: **DONE & VERIFIED.** FWD accumulates the per-order surface downwelling into
  `s_Level_Rad_DOWN(Index_Sat_Angle,n_Layers)` (gated on `Compute_Down_Radiance`); TL adds optional
  `down_rad_TL_out` (sum over orders of `s_IterRad_DOWN_TL` at the surface); AD adds optional
  `down_rad_AD_in` (per-order seed into `s_IterRad_DOWN_AD`, exact transpose). Both SOI call sites in
  `CRTM_RTSolution.f90` pass the new args. Verified by `test_Unit_Downwelling_TLADK` (cloudy ICE_CLOUD,
  RT_SOI, overcast TCC=1): FD-vs-TL 6.6e-11, adjoint dot-product 2.4e-16, K-vs-AD exact; scattering
  values differ from clear-sky (cloud effect confirmed, 22/22 channels scattering). Full suite 212/212.
  - Harness note: scattering scenes use Cloud_Fraction=1 (overcast, TCC=1) so the missing fractional
    TL/AD `Down_Radiance` combine doesn't bite; TL/AD/K cloud structures made congruent with FWD.
- ADA TL/AD: **DONE & VERIFIED.** Dense TL/AD of the downward adding recursion (matmul + matinv)
  + Inv_Gamma3 finalization. Verified by `test_Unit_Downwelling_TLADK` on a strongly-scattering
  overcast SNOW scene (Reff=500um, TCC=1, RT_ADA): FD-vs-TL best |FD/TL-1| = 3.2e-10, adjoint
  dot-product = 0 (exact transpose), K-vs-AD = 0 (exact). All six harness cases green; full suite
  212/212.
  - Interface mirrors Emission/SOI: optional `down_rad_TL_out` (CRTM_ADA_TL) / `down_rad_AD_in`
    (CRTM_ADA_AD) + optional `Index_Sat_Angle`, wired at the four ADA call sites in
    CRTM_RTSolution.f90 (TL n_Stokes>1 & ==1; AD n_Stokes>1 & ==1). Both gated on
    `RTV%Compute_Down_Radiance .AND. PRESENT(...)` (ADA's downward RTV intermediates are only
    populated when the flag is on, unlike SOI which always iterates downward).
  - **FWD bug found & fixed (was latent):** `CRTM_ADA` copied the finalized profiles back into the
    working arrays (`s_Level_Rad_DOWN = s_Level_Rad_DOWNT; s_Level_Rad_UP = s_Level_Rad_UPT`)
    unconditionally inside the downward-sweep `IF`. With `Compute_Down_Radiance` on this CLOBBERED
    the INTERMEDIATE `s_Level_Rad_UP` that the existing upward TL/AD read → corrupted the **TOA**
    Jacobian (FD/TL≈0.345), not just the downwelling one. The adjoint dot-product still passed (TL &
    AD consistently wrong) so only the FD check caught it. Fix: gate the copy-back on
    `aircraft%rt .or. obs_4_downward%rt` (the only forward-only observers that read the finalized
    full profiles); `Compute_Down_Radiance` now leaves `s_Level_Rad_DOWN`/`_UP` as intermediates for
    TL/AD. Output extraction reads `s_Level_Rad_DOWNT` directly for ADA/VMOM — but **solver-specific**:
    SOI stores its final surface value in `s_Level_Rad_DOWN` (no DOWNT), so Common_RTSolution branches
    on `RTV%RT_Algorithm_Id == RT_SOI`.
  - Implementation facts (mirror the FWD DO 20 + finalization): output is the FINALIZED
    `s_Level_Rad_DOWNT(n1,n_Layers)`; the TL re-runs the downward adding recursion using the saved
    RTV intermediates (Inv_Gamma2/2T, Refl_Trans_DOWN, s_Level_Refl_DOWN, Inv_Gamma3, layer matrices)
    and re-calls CRTM_AMOM_layer_TL/_AD per scattering layer (layer-indexed, so reproduces the FWD
    exactly; redundant w/ the upward call but correctness-first); the surface-boundary `s_rad_up_TL` /
    `s_refl_up_TL` are CAPTURED before DO 10 (which overwrites them); the AD surface contributions are
    added AFTER the upward DO 10 "=" surface block (so they `+=`) and before the total_opt_AD→T_OD_AD
    propagation (the downward AMOM_layer_AD calls accumulate total_opt_AD too).
- Clear/cloudy `Down_Radiance` combine in FWD/TL/AD/K: **DONE & VERIFIED.** Mirrors the existing
  Radiance/Stokes combine including the cloud-fraction (TCC) sensitivity term; all gated on
  `Opt%Compute_Down_Radiance` (so flag-off scattering scenes keep the original Down_Radiance=0,
  matching baseline). Verified with a fractional (Cloud_Fraction=0.5) SOI scene in
  test_Unit_Downwelling_TLADK: FD-vs-TL 6.3e-11, adjoint 2.4e-16, K-vs-AD exact. Full suite 212/212.
  KEY ADA TAKEAWAY (why this was done first): the fractional/TCC coupling lives entirely in the
  driver combine, and both clear & cloudy sub-solves receive their Down_Radiance seed via the same
  Assign_Common_Input_AD extraction — so the ADA solver only needs the SAME optional
  down_rad_TL_out / down_rad_AD_in interface as Emission/SOI; nothing ADA-combine-specific.
- **Clear/cloudy combine:** for fractional-cloud scenes the driver combines `%Radiance`/`%Stokes`
  from the clear (emission) and cloudy (scattering) sub-calls but NOT `%Down_Radiance`. Add the
  same cloud-fraction combine for `Down_Radiance` (FWD/TL/AD/K). This is what makes cloudy
  `Down_Radiance` physically complete and resolves the 2 deferred TL Simple failures.
- Gate on Phase-0 scattering + grazing tests; regenerate scattering/cloudy TL/AD baselines.

### Phase 4 — Retire `Obs_4_downward_P`  *(DONE)*
- Removed: the Options field (`CRTM_Options_Define`), `obs_4_downward_type` + `RTV%obs_4_downward`
  (`RTV_Define`), the FWD parse block (`CRTM_Forward_Module`), both `%Radiance`/downward-overwrite
  branches (`Common_RTSolution`, scattering + emission), and the ADA gate refs (now keyed on
  `aircraft%rt .or. Compute_Down_Radiance`). `CRTM_Get_PressureLevelIdx` /
  `AIRCRAFT_PRESSURE_THRESHOLD` kept (still used by the aircraft observer). No TL/AD/K dead imports
  existed on this baseline (the channel-threaded refactor had already dropped them).
- `test_Downwelling_Radiance` migrated: sets `Options%Compute_Down_Radiance = .TRUE.` and regresses
  the full RTSolution (now `%Radiance` = TOA upwelling, `%Down_Radiance` = surface downwelling) across
  all 8 registered sensors; baselines regenerated (self-healing — not git-tracked). TL/AD/K coverage
  is provided rigorously by `test_Unit_Downwelling_TLADK` (FD-vs-TL, adjoint dot-product, K-vs-AD for
  clear/SOI/ADA) rather than duplicated into the multi-sensor forward regression.
- **API note:** `Obs_4_downward_P` was a *public* field of the public `CRTM_Options_type` (the old
  test set it directly), so its removal is a breaking change for any external code that referenced it.
  Capability change: arbitrary-pressure-level downwelling is gone; `%Down_Radiance` is surface-only
  (the deferred `Downwelling_Radiance(:)` profile, §10 R4, would restore arbitrary levels).
- **Release note:** a distributed `fix_*` tarball / stale `build/test/results` carrying the old
  `test_Downwelling_Radiance` baselines (downwelling-at-320-hPa overwriting `%Radiance`) must be
  regenerated; clean builds self-heal.
- Full suite 212/212 after retirement.

## 8. Verification

Primary gate: `test_Downwelling_TLADK` — FD-vs-TL central diff (ratio→1 ~1e-6), adjoint
dot-product (machine precision), K-vs-AD (exact), per solver + TOA control + grazing. TOA control
must stay green throughout (regression guard on the primary radiance path).

## 9. Baseline / IO impact

- Scalar `Down_Radiance` is already serialized and already nonzero in forward; making TL/AD
  RTSolution carry `Down_Radiance_TL/_AD` changes TL/AD baselines → regenerate. Baselines are not
  in git; clean builds self-heal, but the distributed `fix_REL-3.2.0.0` tarball / stale
  `build/test/results` must be regenerated.

## 10. Risks & open questions

- **R1** ADA TL/AD of the downward path (Phase 2) is the dominant effort and main correctness risk.
- **R2** Grazing consistency (C1): downwelling Jacobian must match the clamp-flag branch or the
  adjoint test fails at grazing — explicit test required.
- **R3** Aircraft (upwelling-at-altitude) shares the same hack pattern and is forward-only for
  TL/AD/K. Out of scope; candidate for the same first-class treatment later.
- **R4** Arbitrary-level downwelling returns only with the deferred profile field (Phase 3). Flag if
  sub-mm/aircraft downwelling-at-level becomes a near-term requirement.
- **R5** n_Stokes>1: emission branch currently broadcasts a scalar downwelling into all Stokes
  (`Common_RTSolution.f90:1201`) — write only `Stokes(1)` for emission.
