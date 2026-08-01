# Design brief: realistic MW land-surface *state* Jacobians below ~80 GHz

**Goal (BTJ, revised 2026-07-19 after the architecture decision):** realistic
land-surface Jacobians with respect to **land geophysical state** (soil moisture,
LAI, vegetation), using the TELSEM2 atlas as the forward emissivity and
`NESDIS_LandEM` as the sensitivity source. Enable any remaining missing Jacobians
in the process.

> **Scope change — read this.** The original goal said "realistic Jacobians at
> **all** MW frequencies." **That is no longer the target.** BTJ chose
> Architecture 1 (physical state) alone; see the decision record below. Arch 1 is
> physically capped at ~80 GHz, so above that band this program delivers a
> *correctly near-zero* soil-moisture sensitivity and **nothing** for the
> variables that actually drive high-frequency land emission. Achieving
> "realistic at all frequencies" would have required Architecture 2 (emissivity
> as control variable), which was considered and **not** selected.

**Status (corrected 2026-08-01):** scoping complete, architecture decided, and
the implementation has since landed. The header previously read "implementation
not started", which its own inline update blocks below already contradicted.
Verified by `git branch --contains`:

- `6c8f420` TELSEM2/LandEM additive-anomaly hybrid MW-land Jacobians, and
  `58c7179` the `Use_MWland_Atlas` opt-in gate, both 2026-07-19. These exist
  **only** on `feature/btj_telsem2_landem_hybrid` and have **not** been merged
  into `feature/btj_REL-3.2.0`.
- `aa31bb7` analytic MW-land Soil/Land temperature emissivity Jacobians,
  the #281 Phase 3 item. This one **has** propagated to
  `feature/btj_REL-3.2.0`.

**Repo/branch:** `~/CRTM/CRTMv3`. The former `feature/btj_lifecycle_header_docs`
was merged as PR #332 and deleted.

**Release context:** the v3.2.0 date of 2026-07-31 has passed. The original
claim that none of this work fits that window is now only partly true, since
`aa31bb7` did land on the release branch. The TELSEM2 hybrid itself remains
unmerged and is a 3.2.x/3.3 program. Of the punch-list bug fixes, only item 1
(Fastem1 SST) was ever realistically 3.2.0-eligible; see the verified
punch-list, which downgrades the other two.

> ## ✅ DECISION RECORD (BTJ, 2026-07-19)
>
> **Architecture 1 only — physical state as the control vector.** Emissivity-space
> (Arch 2) is not being pursued.
>
> **Accepted consequences:**
> - "Realistic surface Jacobians at all MW frequencies" is **formally dropped** as
>   a deliverable. Above ~80 GHz the honest answer stays ~zero/absent.
> - `RTSolution_K%Surface_Emissivity` stays unpopulated (verified: still zero
>   occurrences in `CRTM_K_Matrix_Module.f90` and `CRTM_Adjoint_Module.f90`).
>   The Arch-2 write-up below is retained as the record of the road not taken —
>   it remains the only route to high-frequency surface sensitivity if the
>   requirement ever returns.
> - TELSEM2's dropped std/covariance stays dropped; no background-error work.
>
> **First implementation step:** un-shadow the `9358caa` LandEM Jacobians, which
> the TELSEM2 early `RETURN` at `CRTM_MW_Land_SfcOptics.f90:254` currently zeroes
> (verified in tree; the inline comment there confirms `iVar%Compute` stays
> `.FALSE.`). Reconciling TELSEM2-forward with LandEM-sensitivity is the core of
> the program. The TELSEM2 opt-in gate (§ Implementation surface) is a
> prerequisite, not an optional extra.

> ## ✅ IMPLEMENTED (BTJ, 2026-07-19) — Arch-1 hybrid, first step done
>
> The atlas early `RETURN` is removed. `Compute_MW_Land_SfcOptics`
> (`src/SfcOptics/CRTM_MW_Land_SfcOptics.f90`) now runs an **additive
> emissivity-anomaly hybrid**: below the 80 GHz cutoff, when the atlas is valid
> it keeps the atlas emissivity as the forward value **and** calls
> `NESDIS_LandEM` into scratch purely to harvest `dEV/dEH_dvlai`, `dEV/dEH_dmv`,
> which flow through the existing TL/AD unchanged. The LandEM derivative is used
> **unscaled** (chosen over the ratio `e_atlas/e_landem` form: minimal assumption,
> ~2% numerical difference here, and the reported sensitivity is exactly LandEM's
> already-validated derivative).
>
> **Accepted trade (Option A, chosen 2026-07-19):** the forward emissivity over
> atlas-covered land is the state-independent atlas value, so the reported
> Jacobian is **NOT** the finite-difference derivative of the forward (FD of the
> actual forward w.r.t. soil state = 0). TL↔AD/K stay mutually exact (verified to
> 1.4e-14). This is fine for emissivity-sink-style DA; a strict incremental
> variational check will flag the forward/Jacobian mismatch. The consistent
> alternative (Option C — a true anomaly forward carrying a background state
> reference through SfcOptics) was deferred as architecturally invasive.
>
> **Verification:** forward is bit-identical on every non-atlas path (215/215
> ctest; ocean-fallback diff exactly 0). `test_TELSEM2_MWland` was rewritten from
> "assert zero Jacobian" to the hybrid contract: non-zero atlas-path Jacobian
> (|K|≈44.8 K/unit SMC), TL==K to 1.4e-14, atlas still driving the forward
> (spatial/seasonal/ocean-fallback checks retained).
>
> **NOW ALSO DONE — TELSEM2 opt-in gate (BTJ, 2026-07-19).** The atlas no longer
> activates from mere file presence. `CRTM_Init` gained an optional
> `Use_MWland_Atlas` (default `.FALSE.`); an explicit `MWlandCoeff_File` also
> counts as opt-in. Without opt-in the atlas is skipped even when the
> default-named file is on the path, so MW-land optics use `NESDIS_LandEM` (with
> its #281 Jacobians). `Use_MWland_Atlas=.TRUE.` auto-resolves the default-named
> atlas; if absent it warns and falls back (non-fatal). Missing explicit file =
> hard error (unchanged). `test_TELSEM2_MWland` stages the atlas under the
> default name too and proves the gate both ways; the default-named symlink is
> inert for the rest of the suite (215/215). Commit `58c7179`.
>
> **STILL not done:** reconciling with the invalid-cell fall-through for
> snow/ice; anything above the 80 GHz cap (out of scope by the Arch-1 decision).
> The line references in the sections below predate these changes and may be a
> few lines off.

---

## The central finding — read this first

**There is no physical model in this tree — and effectively none in the standard
MW-RT world — that produces a *realistic* land emissivity Jacobian from
geophysical state above ~80 GHz.** This is not a wiring gap that effort closes.
It is physics plus a missing control vector:

1. Above 80 GHz every operational MW land/snow/ice path returns a **hardcoded
   constant** (0.95 / 0.90 / 0.92) with zero derivatives. `NESDIS_LandEM` — the
   only model with real physics *and* analytic Jacobians — is hard-gated to
   < 80 GHz at every call site (`CRTM_MW_Land_SfcOptics.f90:213`,
   `CRTM_MW_Snow_SfcOptics.f90:183`).
   **Precision:** "zero derivatives" here means the *emissivity-state* Jacobians
   (∂Tb/∂LAI, ∂Vegetation_Fraction, ∂Soil_Moisture, ∂Soil_Temperature) — the
   sensitivities that route through the emissivity model. The land skin-**temperature**
   Jacobian ∂Tb/∂`Land_Temperature` is *not* zero and is frequency-independent:
   it comes through the surface Planck emission term (`CRTM_SfcOptics.f90`
   `SurfaceT_AD`, weighted by `Land_Coverage`), so a constant ε = 0.95 still emits
   ε·B(T) and a skin-T perturbation still moves Tb at 183 GHz. `Soil_Temperature`
   is *sub-surface* (feeds only LandEM's emissivity physics, not the skin-T
   emission), so it behaves like the emissivity-state terms, not like
   `Land_Temperature`.
2. At 89–183 GHz the soil is radiatively invisible (sub-mm penetration). The
   physically correct ∂e/∂(soil moisture) there is **~zero**. High-frequency land
   emissivity is driven by vegetation water content and snow volume scattering.
3. CRTM's surface type does not carry the variables that drive that physics as
   *differentiable state*: `Canopy_Water_Content` is never read by any SfcOptics
   module; `Snow_Density` and `Snow_Grain_Size` are never read by any MW model
   (grain size is used by IR only); snow density/grain size are hardcoded inside
   LandEM (`va = 0.4 + 0.0004*depth`, `rad = 0.5 + 0.005*depth`).

So "realistic land Jacobian at 183 GHz w.r.t. soil moisture" is a category error —
the honest value is ~zero. And "realistic Jacobian w.r.t. the variables that *do*
matter at 183 GHz" cannot be produced because neither the forward physics nor the
control variables exist in the tree. **No amount of LandEM extension fixes this.**

TELSEM2 is the only thing giving realistic, spatially/seasonally varying
emissivity up to ~190 GHz — but it responds to **geography and season, not
surface state**, carries **no uncertainty** (the reference atlas's std/covariance
is explicitly dropped in the CRTM port, `TELSEM2_Atlas_Module.f90:12-14`), and has
**identically-zero TL/AD** by construction. Its high-frequency values at 166/183
are a *twice-extrapolated* projection off the 85 GHz anchor (native channels are
only 19/22/37/85 GHz), and only for water-like surface classes; other classes are
held flat at the 85 GHz value. There is **no frequency-range guard** — it
silently extrapolates (`interp_freq2`, `TELSEM2_Atlas_Module.f90:444-484`).

---

## Two architectures can deliver "realistic Jacobians at all frequencies"

They answer **different scientific questions**. The choice is the whole decision.

### Architecture 1 — physical hybrid (TELSEM2 forward + LandEM sensitivities)

TELSEM2 supplies absolute emissivity; LandEM supplies normalized ∂e/∂state.
Makes **land geophysical state** (soil moisture, LAI, vegetation) observable.

- **Works, and is the right tool, below ~80 GHz** — where LandEM is valid and
  where soil/vegetation genuinely modulate emission. Coincidentally exactly the
  band LandEM already covers.
- **Cannot deliver realistic *state* sensitivity above 80 GHz.** Best case it
  returns a correctly-near-zero soil-moisture Jacobian and *nothing* for the
  variables that actually matter there (vegetation water, snow) because that
  physics and those control variables are absent. So it does not, by itself,
  satisfy "realistic at all frequencies."
- Acceptance criterion must reward *correctly small where the surface is opaque*,
  not magnitude — a validation that chases non-zero produces a confidently wrong
  operator.

### Architecture 2 — emissivity as the control variable (retrieve/adjust e)

Treat surface emissivity itself as the control/sink variable; expose ∂Tb/∂e.
Makes **emissivity** observable/adjustable at **every** frequency. This is the
operational mainstream (ECMWF, NOAA dynamic-emissivity / emissivity-sink DA).

- **Delivers a realistic, complete surface Jacobian at all frequencies** — because
  it never asks a physical model to know the geophysics. The emissivity increment
  absorbs whatever the true surface does.
- CRTM **already computes ∂Tb/∂e internally** in the SfcOptics_K machinery; it is
  simply **not mapped out**. `RTSolution_K%Surface_Emissivity` exists as a field
  (`CRTM_RTSolution_Define.f90:237`) but is **never populated** in
  `CRTM_K_Matrix_Module.f90` or `CRTM_Adjoint_Module.f90` (verified: zero
  occurrences). Wiring it out is tractable, contained, and independently useful.
- TELSEM2 is used in its *natural* role here: the emissivity **background/first
  guess**. This is where its dropped **uncertainty** matters — an emissivity-sink
  scheme needs a background error, so reinstating TELSEM2's std/covariance (or
  supplying a spectral-persistence assumption) becomes real work.
- Changes the control vector: emissivity per channel, not soil moisture / LAI.

### ~~Recommendation: build both, split by what is being retrieved~~ — NOT ADOPTED

> **Superseded by the decision record above.** BTJ selected Arch 1 alone. The
> split-by-frequency recommendation below is retained for context and as the
> re-entry point if the "all frequencies" requirement ever returns.

- **< 80 GHz:** Architecture 1 — physical Jacobians so land *geophysical state* is
  observable where the physics supports it (reconcile TELSEM2 forward with the
  `9358caa` LandEM sensitivities; today TELSEM2's early `RETURN` at
  `CRTM_MW_Land_SfcOptics.f90:254` *shadows* those Jacobians to zero).
- **All frequencies, incl. > 80 GHz:** Architecture 2 — wire out ∂Tb/∂e and adopt
  TELSEM2 as the emissivity background. This is the only route that literally
  satisfies "realistic surface Jacobians at all frequencies," and it is the
  operationally proven one.

The two compose cleanly: below 80 GHz a user can choose state-space (Arch 1) or
emissivity-space (Arch 2); above 80 GHz only Arch 2 is physically honest.

**Open decision for BTJ — this determines everything downstream:** is the DA goal
to make *emissivity* observable at all frequencies (Arch 2 — achievable), or to
make *land geophysical state* observable (Arch 1 — physically capped at ~80 GHz
no matter the effort)? "Realistic Jacobians at all frequencies" is reachable for
the former and physically impossible for the latter.

---

## Full missing-Jacobian punch-list (from the #281 audit)

Beyond the TELSEM2 question, "enable any remaining missing Jacobians" resolves to:

### Verified against the tree, 2026-07-19 (scope check before 3.2.0)

All three were confirmed to exist. Two are **larger than originally written** —
neither is a one-line fix, and only item 1 is 3.2.0-eligible.

**1. Water SST Jacobian silently zero (Fastem1) — REAL, smallest, 3.2.0-eligible.**
`iVar%dEH_dTs`/`dEV_dTs` are declared and zero-initialised
(`CRTM_MW_Water_SfcOptics.f90:99,101`) and read by TL (`:531,533`) and AD
(`:761,767`) — but **never assigned anywhere in `src/`** (tree-wide grep: zero
assignments). `Fastem1` emits only `dEH_dWindSpeed`/`dEV_dWindSpeed` (`:313-314`).

*Blast radius is narrower than feared:* the legacy branch splits by frequency, and
the low-frequency leg goes to `LowFrequency_MWSSEM` (`:294`), which carries a
proper internal-variable TL/AD and handles Ts correctly. Only the **`Fastem1` leg
of the legacy `.NOT. Use_New_MWSSEM` path** is affected. Default (`Use_New_MWSSEM`
→ `Compute_FastemX`) and PARMIO are unaffected — which is why no ctest catches it.

*Fix options:* (a) differentiate Fastem1's permittivity/Fresnel — correct, real
work; (b) **finite-difference `dE/dTs` around the existing Fastem1 call and cache
it in `iVar`** — ~15 lines, contained, correct to FD accuracy, the sane 3.2.0 move;
(c) document + warn that the legacy path has no SST Jacobian. Recommend (b).

**2. RH K-matrix leg — REAL, but the original framing is WRONG. Not 3.2.0.**
Confirmed zero occurrences of `Relative_Humidity` in `CRTM_K_Matrix_Module.f90`
and `CRTM_Adjoint_Module.f90`. **But mapping it into `Atmosphere_K` would be a
mistake.** `Atm%Relative_Humidity` is a *derived diagnostic*, computed from T/P/H2O
(`CRTM_Atmosphere_Define.f90:130` uses `PPMV_to_MR, MR_to_RH` — the only public
entries of `CRTM_Relative_Humidity`, and **no TL/AD forms are exported**). It is
not an independent control variable; exposing it alongside q would let a DA system
double-count.

*The real defect:* `CRTM_AerosolScatter.f90:786` accumulates
`Atm_AD%Relative_Humidity`, and that sensitivity is then **silently dropped** — the
chain back to `Atm_AD%Temperature` / `Atm_AD%Absorber(H2O)` does not exist. So
aerosol hygroscopic-growth sensitivity is missing from the **T and q** Jacobians.
Fixing it means writing/wiring the AD of the `MR_to_RH`/`PPMV_to_MR` chain, not
adding a field. The "#285 resolves" annotation is **not** reflected on this branch.

**3. Cloud `Water_Density` adjoint — REAL, and it affects passive MW, not just radar.**
Confirmed commented out at `CRTM_K_Matrix_Module.f90:1131` and
`CRTM_Adjoint_Module.f90:633`. Two corrections to the original note:

- **Uncommenting will not compile.** The commented calls use a two-argument form
  `Calculate_Cloud_Water_Density(Atm, Atm_AD)`, but only a **one-argument forward
  routine exists** (`CRTM_Active_Sensor.f90:127-160`). The adjoint has to be
  *written*; forward is the trivial `Water_Density = Water_Content / dZ_m` (`:157`).
- **This is not active-sensor-only.** `Water_Density` is consumed by the passive
  cloud optics — `CRTM_CloudScatter.f90:343` (FWD), `:614` (TL), `:941` (AD out) —
  where it indexes the **MW** cloud LUT (`CloudC%Water_Density_MW`, `:1390-1392`).
  That is the actual source of the "MW/IR inconsistency" in the issue: the IR path
  indexes on effective radius, the MW path on water density.

*Consequence:* `Atm_AD%Cloud%Water_Density` accumulates and is then thrown away, so
the **MW cloud `Water_Content` Jacobian is missing its density-indexing term.**
Also note `dZ_m` depends on `Height`, so a complete adjoint has a height leg too.

### Real, but need forward-model work (differentiable state exists, no derivative wired)
| Item | Where | Effort |
|---|---|---|
| Soil/Land temperature emissivity sensitivity | `CRTM_MW_Land_SfcOptics.f90:295-296` | Read forward, no derivative. #281 design "Phase 3" (Hard; lower value — bulk T Jacobian already exists via `CRTM_Compute_SurfaceT_AD`). |
| Snow `Snow_Depth`/`Snow_Temperature` | `CRTM_MW_Snow_SfcOptics.f90` | TL/AD are **zero-stubs that don't even take `Surface_TL/_AD`**; `iVar` is `INTEGER::Dummy`. Needs signature change + real iVar, following `9358caa`. |
| Ice `Ice_Thickness`/`Ice_Temperature` | `CRTM_MW_Ice_SfcOptics.f90` | Same zero-stub situation. `NESDIS_SIce_Phy_EM` (Fresnel, responds to T+salinity) exists but is **commented out** at the call site. |
| Cloud `Effective_Variance` | `CRTM_CloudScatter.f90` | Zero occurrences; required only if scattering LUTs depend on size-distribution variance. |

### Structurally impossible without new forward physics (identically zero today)
Per `surface_jacobians_281.md` (a **local, unposted** adjudication): these have
**no forward dependence**, so their analytic Jacobian is identically zero and
cannot be made non-zero without first adding the physics —
`Canopy_Water_Content`, `Snow_Density`, `Snow_Grain_Size` (MW), `Ice_Density`,
`Ice_Roughness`. Note these are exactly the variables that would matter for
high-frequency land/snow — reinforcing that Architecture 2 is the realistic path
to high-frequency surface sensitivity.

**Process note:** the "not implementable" adjudication and the proposal to close
these five out of #281 live only in `docs/design/issue-281-comment.md`, a **draft
never posted to GitHub**. #281's three actual comments do not contain it. So this
is our analysis of the tree, not a JCSDA-accepted decision — it needs posting and
agreement.

---

## Implementation surface (`src/SfcOptics/CRTM_MW_Land_SfcOptics.f90`)

| What | Where |
|---|---|
| `iVar_type` (cached derivatives + clip flags) | ~line 150 |
| `DEFAULT_EMISSIVITY = 0.95`, `FREQUENCY_CUTOFF = 80.0` | 213–214 |
| TELSEM2 atlas path | 227–254 |
| **`IF (atlas_valid) RETURN` — shadows the Land Jacobian** | **254** |
| LandEM path, populates derivative cache | 278–306 |
| High-frequency constant fallback | 312–316 |
| `Compute_MW_Land_SfcOptics_TL` (guard `IF (.NOT. iVar%Compute) RETURN` @392) | 366 |
| `Compute_MW_Land_SfcOptics_AD` | 467 |

For Architecture 2, the mapping to add is in `CRTM_K_Matrix_Module.f90` around the
`SfcOptics_K` solve (~line 1586) → `RTSolution_K(...)%Surface_Emissivity`, mirrored
in `CRTM_Adjoint_Module.f90`.

**Also recommended (both agents):** gate TELSEM2 behind an explicit option, not
mere file-presence. Today dropping `TELSEM2.MWland.EmisCoeff.nc` into the coeff
tree silently switches all land forward emissivity to the atlas *and* zeroes the
land Jacobian (`CRTM_LifeCycle.f90:1250-1263`, no opt-in). That is a foot-gun
independent of this work and a prerequisite for offering Arch 1/Arch 2 as a choice.

---

## Test surface

`test/mains/unit/Unit_Test/test_TELSEM2_MWland.f90` **asserts zero surface
Jacobians as a pass criterion**. It goes red the moment any of this works —
rewrite it, do not delete it (it is the only atlas-path regression guard).
`test_Land_Jacobian.f90` (ctest `test_Unit_Land_Jacobian`) validates the < 80 GHz
LandEM Jacobians via AD/TL-vs-FD on amsua_n19; extend it with the SMC = 0 boundary
case flagged in #281 comment 3. New TL/AD/K parity tests should copy
`test_MW_O3_TLAD.f90` (commit `b9c525a`): self-adapting FD/adjoint/K-≡-AD checks.

---

## What was searched (so it isn't re-searched)

No prior design, session, or issue covered the hybrid. Checked: all session
transcripts (TELSEM2 × hybrid/blend/scale/anomaly/ratio/Jacobian), CRTMv3 `src/`,
`docs/`, root `*.md`, open issues incl. #281 body+comments. Related-but-not-on-topic
sessions: `79967ce7` (polsir_experiment, Jun 8–10, TELSEM2 integration + ATMS obs),
`7e6054e8` (Jun 8 precursor), `14a59ce2` (Jul 16, `~/jedi/parmio_experiment/
telsem2_land/`, atlas maps + coastline masking; a ready-made Arch-1/Arch-2
validation harness with ATMS obs/geovals).

## Governing issues / docs
- **#314** TELSEM2 land atlas — scoped TL/AD as identically zero; the caveat notes
  at `REL-3.2.0_changes_vs_develop.md:124` and `RELEASE_NOTES_v3.2.0.md:73-75`
  route here.
- **#281** missing Jacobians — the `9358caa` LandEM analytic Jacobians +
  `docs/design/surface_jacobians_281.md`.
- User Guide upgrade notes (`sec:v320_upgrade_notes`, crtm-documentation PR #26)
  document the zero-Jacobian consequence; update on any behaviour change.
