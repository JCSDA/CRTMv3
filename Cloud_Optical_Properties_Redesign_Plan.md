# Cloud Physical & Optical Properties — Redesign Plan

**Status:** Draft / analysis pinned 2026-06-01
**Owner:** B. Johnson
**Scope:** Redesign how CRTMv3 represents and consumes cloud (hydrometeor) optical
properties, built from in-house DDSCAT (DDA) single-particle data.
**Relationship:** This is the cloud-optics foundation for the work staged in
`MW_Full_Polarization_Master_Issue.md` and `issues/` (esp. `02_Oriented_Particle_Scattering.md`,
`05_Rain_Ice_Depolarization.md`).

---

## 1. Data inventory (in-house DDA archive)

Archive root (slow DrvFs mount): `D:\WSL_temp_storage` → `/mnt/d/WSL_temp_storage`.

Per **particle** (`geom_param.dat`, `shape.dat`):
- Mass = N_dipoles · d³ · ρ_ice (d = 50 µm dipole spacing, ρ_ice = 917 kg/m³)
- D_max (maximum diameter), D_eq (volume-equivalent)
- Projected area, surface area, volume
- Circumscribing-ellipsoid axis lengths → **aspect ratio**
- Multiple density definitions (circumscribing-ellipsoid based, D_max based)
- Box-counting fractal dimension
- Morphological provenance, e.g. `cmb:blr-003->scl-050` (combine → blur → scale):
  a **controlled densification continuum** at fixed D_max.

Per **(habit, density-class, size, frequency)** (`ddscat07/<freq>GHz/w000r000.avg`):
- Orientation-averaged over **900 orientations × 2 incident polarizations**
- Scalar efficiencies: Qext, Qabs, Qsca, g=<cosΘ>, <cos²Θ>, Qbk, Qpha, Qpol, dQpha
- **Full angular Mueller matrix** S₁₁, S₁₂, S₂₁, S₂₂, S₃₁, S₄₁ … at 1° resolution, 0–180°

Habits present: `pristine`, `aggregate`, `snowfake`. Frequencies (25):
3, 5, 11, 14, 19, 24, 36, 89, 94, 150, 166, 176, 180, 186, 190, 205, 240, 325, 380, 425, 462, 500, 640, 683, 874 GHz
(spans GMI / ATMS / MHS / AMSR2 / ICI / sub-mm).

Assembled DB: `openSSP.nc` (5.3 GB single-scattering-property database; ncdump
currently core-dumps reading its header over the DrvFs mount — read locally).

Current build pipeline:
- `parse_avg_dat_files.pl` — **extracts scalar efficiencies only** (discards the angular
  Mueller matrix and the Qpol/dQpha polarization info).
- `mass-dimension_eval.py` / `create_csv_15.py` — fit M = a·Dᵇ, benchmark vs.
  Thompson / Heymsfield / Brandes.

---

## 2. Current paradigm in CRTM (baseline)

Shipped LUT: `build/test_data/3.2.0/.../CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2024.nc`
(Release 3, Version 4; Moradi 2024; MW solid from Eriksson et al. 2018 ARTS-SSDB, MW
liquid + IR/VIS from Mie/Simmer + Yang).

- Bulk, **PSD-pre-integrated scalars**: ke, w, g, kb + **scalar** phase function
  (`n_Phase_Elements = 1`, `n_Legendre_Terms = 39`).
- MW **liquid** indexed by (Temperature, R_eff, Frequency): 8 T, 100 R_eff, 886 freq.
- MW **solid** indexed by (**Density**, R_eff, Frequency): 18 densities (density used as a
  **habit / m–D proxy**, `DDA_Habits_S_MW`), **no temperature axis**.
- IR/VIS indexed by (Density, R_eff, Frequency).

Runtime (`src/AtmScatter/CRTM_CloudScatter.f90`, `CloudCoeff_Define.f90`):
- Input per layer per cloud: **Water_Content + Effective_Radius** (+ Water_Density fallback).
- Lagrange interpolation (1D/2D/3D) over freq / R_eff / T / density; TL/AD/K versions exist.
- Cloud types: WATER(1) ICE(2) RAIN(3) SNOW(4) GRAUPEL(5) HAIL(6) + extended DDA habits 7–25.
  Type→LUT mapping is coarse (e.g. ICE → 1-D in frequency only).
- τ_ext = ke·WC; SSA = τ_scat/τ_ext; P₁₁ Legendre coeffs accumulated by backscatter weighting.
- RT is scalar by default; vector infra exists (`n_Stokes` up to 4) but MW is effectively V/H.

---

## 3. Diagnosis — four structural limitations

1. **R_eff + one fixed PSD is the wrong coordinate for frozen hydrometeors.** Cannot
   distinguish populations with equal R_eff but different PSD width/shape and m–D; cannot
   respond to the host model's actual PSD.
2. **"Density as habit proxy" conflates two independent axes** — intrinsic structure (habit)
   vs. the mass–size relation (m–D). `geom_param.dat` shows ρ varies *with* D within a habit.
3. **Scalar phase function** (`n_Phase_Elements=1`) discards the 6 independent elements
   (α₁–α₄, β₁, β₂) vector RT needs — the very parameterization `Common_RTSolution.f90` names.
4. **No consistency contract with host microphysics** (Thompson/Morrison/P3/GFDL); the
   forward-operator m–D/PSD mismatch is the dominant all-sky MW assimilation bias.

---

## 4. Proposed method — habit-resolved, two-moment, full-Stokes cloud optics

Four layers, ordered by impact-per-effort. Layers 1–2 capture most of the benefit at low
architectural risk; Layers 3–4 are the novel, publishable contribution.

**Layer 1 — Normalized two-moment PSD indexing (replace R_eff).**
Normalized gamma: N(D) = N_w·f(μ)·(D/D_m)^μ·exp(−(4+μ)·D/D_m). Index bulk properties by
**(D_m, μ)** per (habit, freq, T). Maps cleanly to 1-/2-/3-moment microphysics (mass+number→D_m;
reflectivity→μ). Keep **pre-integrated** first so the interpolation machinery and TL/AD/K are
near-mechanical. (Conceptually RTTOV-SCATT hydrotables; differentiator is Layers 2–4.)

**Layer 2 — Store the full 6-element phase matrix.**
Re-parse `.avg` to keep S₁₁, S₁₂, S₂₂, S₃₃, S₃₄, S₄₄ and expand in **generalized spherical
functions / Wigner-d** (not plain Legendre — only P₁₁ expands in ordinary Legendre) → Greek
constants (α₁–α₄, β₁, β₂). Bump `n_Phase_Elements` 1→6 in `CloudCoeff_Define.f90` + netCDF IO.
Prerequisite for the full-pol roadmap. Plan **δ-M / δ-fit truncation + TMS correction** for the
forward peak at high frequency (39 terms is too few ≥325 GHz).

**Layer 3 — m–D as a first-class, host-aligned axis (kills limitation #2).**
Tabulate kernels in (habit, m–D) / mass-space; bulk density becomes diagnostic, not an index.
Operator uses the *same* m–D as the coupled microphysics scheme.

**Layer 4 — Oriented particles → extinction matrix (ambitious / novel).**
Average over azimuth (β,φ) but keep tilt distributed about a canting mean (current `ddscat.par`
averages all of β,θ,φ). Yields the 4×4 **extinction matrix** (κ_I, κ_Q) + azimuth-dependent
phase matrix for `issues/02` & `issues/05`; gives ZDR/KDP/A_DP consistency with polarimetric radar.
New axes: canting-mean, canting-width.

**Novel contribution unique to this archive:** a **continuous density / rime-fraction axis**
(from `cmb:blr->scl` morphological densification) — the P3 "predicted particle properties"
philosophy applied to radiative transfer, full-Stokes, on explicit-geometry DDA (not
soft-sphere / Maxwell-Garnett EMA). With kb already stored, make it **self-consistent for
active + passive** (GPM-DPR / EarthCARE + GMI / ICI).

---

## Design contract (owner-confirmed 2026-06-01)

1. **Legacy is the untouched default.** The existing MIE_TAMU / DDA_ARTS CloudCoeff path and its
   results are byte-for-byte unchanged unless the user *explicitly* opts in. CRTM-as-a-black-box
   library users see no behavior change. The new scheme is **experimental / opt-in**.
2. **Opt-in via `Cloud_Model`.** Reuse the existing `Cloud_Model` argument (threaded
   `CRTM_Init` → `CRTM_CloudCoeff_Load`; default `'CRTM'`). A new value (e.g. `'CRTM-Exp'`) selects
   the experimental loader + scattering path. The scheme is **explicitly selected, not auto-detected
   from file contents** (today MIE vs DDA is inferred via `ALL(Reff_MW>0)` in `CRTM_CloudScatter`) —
   so swapping a coefficient file can never silently change a black-box user's results.
3. **Format is free to change** (owner: not tied to the current structure). The experimental
   CloudCoeff drops the parts that are accuracy/solver-hostile:
   - **Remove the `{4,6,8,16}` `lOffset` block packing** — it bakes the §7 Legendre↔stream coupling
     into the file (costs up to ~2.7 K on DDA habits, §7.4). Store **one full/variable-length GSF
     (Greek-constant) expansion + its own truncation order**; the solver reads that order.
     Decoupling is the *default* for the new scheme — no `Options` flag needed.
   - **Replace Reff-only / density-as-habit-proxy indexing** with explicit **habit axis + PSD-moment
     axes (D_m, μ)** (Layer 1), **`n_Phase_Elements` up to 6** (Layer 2), and **temperature for solid**.
4. **Self-contained.** The experimental path has its own Define/IO/scatter read code; legacy
   modules and regression baselines are not modified. Collapsing the new scheme to 1 phase element +
   legacy truncation must still reproduce legacy results (sanity bridge, not a format constraint).

## 5. Migration path (respects CRTM architecture)

Infra is more ready than it looks: `n_Phase_Elements` already a dimension; solid MW already
uses a non-R_eff index; TL/AD/K interpolation already exists in `CRTM_CloudScatter.f90`.

1. **Parser first, no CRTM change** — emit 6 angular elements + GSF coefficients; validate
   energy/symmetry (e.g. F₂₂ ≤ F₁₁).
2. **Prototype LUT** keyed by (habit, D_m, μ, T), `n_Phase_Elements=6`, offline normalized-gamma
   integration of single-particle kernels. Keep old LUT for regression.
3. **Runtime, scalar-RT first** — confirm 6-element LUT collapsed reproduces today's intensity
   results (no-regression safety net), then enable vector consumption.
4. **TL/AD/K** for new (D_m, μ) axes (mechanical given existing Lagrange derivatives).
5. **Oriented / extinction-matrix** as an optional path gated by an orientation flag
   (random → existing fast path; never regress the symmetric case).

---

## 6. Open decisions

- **Ice refractive index is single-temperature** (`ior_266K.dat`); solid MW has no T axis.
  Add one? (matters most for low-freq emission/absorption).
- **Online PSD integration vs. pre-integrated tables** — start pre-integrated.
- **GSF vs. Legendre + truncation scheme (δ-M / δ-fit)** — decide before high-freq tables.
- **Habit taxonomy** — discrete list vs. continuous-density axis.
- **Target** — operational CRTM (assimilation, TL/AD) vs. research/publication prototype
  (changes differentiability-vs-flexibility emphasis).

---

## 7. Active thread — Legendre terms vs. solver streams

**Question:** Is the number of phase-function Legendre terms genuinely *pinned* to the number
of RT solver streams, or is that a soft accuracy convention? Other models (RT4, Evans &
Stephens) impose no such hard coupling.

**Answer (verified in source): it is coupled by an explicit code identity, not by mathematics.**
The chain is:

```
maxReff ──> CRTM_Compute_nStreams (Mie parameter heuristic) ──> n_Full_Streams (∈ {4,6,8})
        ──> AtmOptics%n_Legendre_Terms = n_Full_Streams           [the hard identity]
        ──> lOffset selects the {4,6,8,16}-term LUT block
        ──> CRTM_Phase_Matrix sums l = mth_Azi … n_Legendre_Terms-1
        ──> RTV%n_Streams = n_Full_Streams/2  ──> Double_Gauss_Quadrature
```

Verified call sites:
- **Streams from Mie parameter** — `src/RTSolution/CRTM_RTSolution.f90:1242-1254`. `MieParameter = 2π·maxReff·ν`; then `<0.01 → 4`, `<1 → 6`, `else → 8`. **Auto path caps n_Full_Streams at 8 ⇒ only 4 hemispheric Gaussian streams**, regardless of how forward-peaked the particle is.
- **The hard identity** — `src/CRTM_Forward_Module.f90:961`: `AtmOptics(nt)%n_Legendre_Terms = n_Full_Streams` (mirrored in TL:1023, AD:892, K:1120).
- **LUT packaging constraint** — `src/AtmScatter/CRTM_CloudScatter.f90:260-266` (and `CRTM_AerosolScatter.f90:208-213`): the coefficient file only stores Legendre blocks for {4,6,8,16} terms (`lOffset` 0/5/12/21); the `CASE DEFAULT` falls back to HG/Rayleigh and carries the comment `! Is this correct?`.
- **Solver truncation** — `src/RTSolution/Common_RTSolution.f90:1930`: `DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1`, using `Phase_Coefficient(l,1,k)` (element index **1** = P₁₁ only ⇒ scalar).
- **delta-M scaling present, TMS/single-scatter correction absent on main path** — `src/AtmOptics/CRTM_AtmOptics.f90:414-421` truncates the forward peak using the moment at index `n_Legendre_Terms` and rescales τ, ω. No Nakajima-Tanaka/TMS single-scatter restoration on the operational chain (only in legacy `RTSolution/UWisc_SOI/`).
- **User override exists** — `src/CRTM_Forward_Module.f90:946-948`: `Options%Use_N_Streams` / `Options%n_Streams` lets you force the stream count (→ up to the 16-term LUT block), the lever for sensitivity tests.

**Why it half-makes-sense and where it doesn't:**
- The identity `n_Legendre_Terms = n_Full_Streams` *is* the classical discrete-ordinate aliasing
  bound (with N hemispheric streams, 2N-point Gauss quadrature resolves ≤ 2N moments; here
  n_Full_Streams = 2N moments). So as a *resolution ceiling for the multiple-scattering source
  integral* it is legitimate.
- It is **not** a mandate to globally truncate the phase function there. RT4 (Evans & Stephens
  1991) / RT3, DISORT, LIDORT decouple the expansion order (NLEGEN) from the quadrature order
  (NUMMU/NSTR): they keep the full expansion for the **exact single-scattering** term and use
  delta-M **plus** a TMS/Nakajima-Tanaka single-scatter correction so few streams stay accurate
  for forward-peaked particles. CRTM instead (a) hard-wires moments ≡ streams, (b) caps the auto
  stream count at 4 hemispheric, (c) ships only {4,6,8,16}-term LUT blocks, and (d) has no TMS
  correction on the main path. Those four are CRTM implementation choices, not physics.

**Why this matters for the redesign (Layers 1–2):** a richer, high-order phase matrix for
high-frequency forward-peaked ice/snow will be **clipped to ≤8 terms and 4 streams** by the auto
path, and the off-diagonal phase elements are not even read (`Phase_Coefficient(l,1,k)`). To
actually exploit better cloud optics we must:
1. **Decouple** `n_Legendre_Terms` from `n_Full_Streams` (store a variable/full expansion length
   in CloudCoeff; stop the verbatim assignment; choose streams for accuracy independently).
2. **Re-pack the LUT** away from fixed {4,6,8,16} `lOffset` blocks to a full expansion (+ stored
   truncation order per entry).
3. **Add a TMS/Nakajima-Tanaka single-scatter correction** (or raise the stream cap) so accuracy
   survives at few streams.
4. Read phase elements **2–6** in `CRTM_Phase_Matrix` for vector RT (ties to Layer 2).

**Quick experiment available now:** set `Options%Use_N_Streams=.TRUE.`, `Options%n_Streams=14`
(→16-term block) for a high-freq scattering case and compare TB/reflectance vs. the auto path to
quantify the truncation error the identity is hiding.

### 7.1 Measured stream/Legendre sensitivity (GMI, 2026-06-01)

Driver: `build/stream_test/stream_sensitivity.f90` (standalone, links `build/lib/libcrtm.so`).
GMI, US-standard profile, fully-overcast frozen column over layers 50–75 (T≈228–250 K),
default `CloudCoeff.nc` (Mie-sphere, 38 Legendre terms, 1 phase element). Forward run at forced
`n_Streams ∈ {4,8,16}`. Because of the §7 identity, each run simultaneously sets
n_Legendre_Terms = n_Full_Streams, so the sweep measures the **coupled** (streams+truncation) error.

Two gotchas found while building the test (both now in the driver):
- **`Cloud_Fraction` defaults to 0** ⇒ column treated as clear, cloud stripped (`maxReff→0`,
  auto streams=4, TB=clear). Must set `Atm%Cloud_Fraction>0`.
- **`Options%n_Legendre_Terms` is vestigial** — never read in the forward path; line 961
  overwrites it. The only public dial is `n_Streams`. (Direct confirmation of the hard coupling.)

ΔTB = TB(n16) − TB(n) brightness-temperature error vs. the 16-stream result:

| GMI band | clear→cloudy (snow) | snow n16−n4 | snow n16−n8 | graupel n16−n4 | graupel n16−n8 |
|----------|--------------------:|------------:|------------:|---------------:|---------------:|
| ≤23.8 GHz | ~3–10 K | <0.06 K | <0.01 K | <0.06 K | <0.01 K |
| 36.5 GHz | ~15 K | ~0.1–0.14 K | ~0.01 K | ~0.08 K | <0.01 K |
| 89 GHz | ~50 K | ~0.3 K | <0.03 K | ~0.15 K | ~0.02 K |
| 166 GHz | ~90 K | **~2.0 K** | ~0.09 K | ~0.19 K | <0.01 K |
| 183 GHz | ~85 K | **~2.0–2.5 K** | ~0.10 K | ~0.20 K | <0.01 K |

Findings:
- The truncation error is large where it matters: **~2.0–2.5 K at 166–183 GHz for snow** across
  the available range (4→16). It grows monotonically with frequency (size parameter).
- The **auto heuristic lands at 8 streams** for ≥36.5 GHz here, which is within ~0.1 K of the
  16-stream result — so for these strong-scattering cases auto is *adequate but not converged*,
  and the residual is non-monotonic in channel (e.g. 18.7H, 36.5H flip sign), i.e. genuine
  phase-matrix representation error, not a smooth bias.
- Snow (low density) shows ~10× the sensitivity of graupel (higher density, less forward-peaked).
- **You cannot buy the missing ~0.1–2.5 K of phase fidelity without also adding streams** — the
  identity forbids it. This is the concrete radiative cost of the §7 coupling, and the headroom
  the Layer-2 redesign (full phase matrix, decoupled truncation) is competing for.

Caveat: this uses the **Mie-sphere** default LUT. Real DDA aggregates (the archive) have less
forward-peaked but more structured phase functions; the sensitivity will differ and should be
re-measured against the DDA CloudCoeff once a habit/Reff cloud setup is wired up.

### 7.2 Where the auto heuristic under-resolves — size/WC/sensor map (2026-06-01)

Driver: `build/stream_test/stream_map.f90`. Snow column, Reff swept 50–1500 µm, WC = 0.1 and
0.5 kg/m²/layer, sensors GMI / ATMS / MHS. Per channel: `max|TB_auto − TB_16|` over the size
sweep, the stream count the heuristic chose there, and the full-range `max|TB_n4 − TB_16|`.

Headline numbers (snow, WC=0.5):

| Channel | scat. depression | max\|auto−16\| | auto NFS there | max\|n4−16\| |
|---------|-----------------:|--------------:|---------------:|-------------:|
| GMI 166 V/H | ~130 K | ~0.51 K | 8 | ~2.3 K |
| GMI 183±3/±7 | ~90–110 K | ~0.23 K | 6–8 | ~2.6–2.7 K |
| ATMS 183 sounders (17–22) | 70–130 K | 0.20–0.52 K | 6–8 | 2.2–2.7 K |
| MHS 89/157/183/190 (1–5) | 70–134 K | 0.20–0.56 K | 6–8 | 2.0–2.9 K |
| ATMS 50–60 GHz O₂ (10–15) | 0 K | 0.000 K | — | 0.000 K |

Three structural conclusions:

1. **The Mie-parameter heuristic is largely self-protecting.** Because it scales streams with
   size parameter, `max|auto−16|` stays ≲0.5 K even though the full 4→16 range costs up to
   **~2.9 K**. Window/sounding channels with no scattering show exactly 0 (good sanity check).

2. **A systematic mid-size-parameter blind spot.** At **Reff ≈ 200–250 µm, 166–190 GHz, optically
   thick (WC=0.5)** the Mie parameter sits just below 1, so the heuristic picks **6 streams** and
   loses a consistent **~0.20–0.24 K** (same sign). This lands squarely on the operational
   183 GHz humidity sounders (ATMS 18–22, MHS 3–5, GMI 12–13) that dominate all-sky MW
   assimilation. It appears only at high WC — i.e. deep convection / heavy snow.

3. **A hard ceiling at 8 streams.** The auto selector never exceeds `n_Full_Streams=8`, so for the
   most strongly scattering cases (large Reff, 166 GHz) even auto's best leaves **~0.5 K** vs. 16
   streams. The heuristic *cannot* close this — only forcing `Options%n_Streams` can.

Implications for the redesign:
- **Cheap, high-value fix:** raise the auto stream cap above 8 and/or make the Mie-parameter
  thresholds depend on scattering optical depth (closes both #2 and #3 without new LUTs).
- These are **Mie spheres** — a smooth phase function that *under*-populates high Legendre orders.
  Real DDA aggregates put energy into higher orders, so the truncation error here is a **floor**;
  the DDA case (the archive) is expected to be worse and is the one that matters. Re-measure
  against the DDA CloudCoeff is the natural next experiment.
- The hard coupling (§7) means any extra phase-matrix fidelity from the redesign is **gated by the
  stream count** — so decoupling + raising the cap is a prerequisite for the optics work to pay off.

### 7.3 DDA habits vs. Mie — the truncation error is 5–6× larger (2026-06-01)

Driver: `build/stream_test/dda_stream_map.f90`, GMI, `CloudCoeff_DDA_Moradi_2024.nc`. DDA mode
indexes solid hydrometeors by **habit** (cloud type) and **`Water_Density`** (size axis), not Reff —
confirmed: the file has no `Reff_MW` (only `PSD_Computed_Reff_MW`), so `ALL(Reff_MW>0)` is false and
the DDA branch (`CRTM_CloudScatter.f90:1314`) fires. Truncation measured with **forced** streams
(Reff-independent). WC=0.5 kg/m²/layer, `Water_Density` swept 1e-6–1e-2 (LUT range).

Worst-case TB error over the size sweep at 166/183 GHz (GMI ch 10–13):

| Habit (LUT index) | scat. depression | max\|n4−16\| | max\|n8−16\| | max\|auto−16\| |
|-------------------|-----------------:|------------:|------------:|---------------:|
| Mie sphere (§7.2 ref) | ~90–130 K | ~2.5 K | ~0.1 K | ~0.5 K |
| SNOW default (→SectorSnowflake) | ~150–170 K | ~1.3 K | ~0.01 K | ~0.1 K |
| LargePlateAggregate | ~130–150 K | ~8 K | ~0.2 K | ~1.2 K |
| EvansSnowAggregate | ~70–85 K | ~13 K | ~0.4 K | **~2.3 K** |
| SixBulletRosette | ~115–135 K | **~16 K** | ~0.4 K | **~2.9 K** |

Findings:
- **Realistic DDA aggregates/rosettes are far more forward-/structure-peaked than Mie spheres.**
  The full 4→16 truncation range reaches **~16 K** (SixBulletRosette) vs. ~2.5 K for Mie, and the
  operational **auto path leaves ~2–3 K** at 166–183 GHz vs. ~0.5 K for Mie. These are exactly the
  channels + particle types used in all-sky precipitation assimilation.
- **Habit dominates.** The default `SNOW_CLOUD` mapping (→SectorSnowflake) is benign (~0.1 K auto),
  but the realistic aggregate/rosette habits are **20–30× more sensitive**. Whatever habit the
  redesign exposes directly sets how badly the §7 coupling bites.
- **8 terms is the practical knee** for these habits (n8→n16 ≈ 0.2–0.4 K), but the auto heuristic
  selects fewer (it keys off the nominal Reff and the Mie-parameter table), so it under-resolves by
  ~2–3 K. Forcing `Options%n_Streams=8` already recovers most of it — cheap validation of the fix.
- Confirms the §7.2 caveat: the Mie numbers were a **floor**. For the in-house DDA archive (mostly
  aggregates) the truncation/stream resolution controls **multi-Kelvin** radiance error, so
  decoupling Legendre order from streams + raising the cap is not optional polish — it is required
  for the DDA optics to be usable at high frequency.

### 7.4 Prototype: decouple Legendre truncation from streams — the fix is nearly free (2026-06-01)

**Source change (working tree, gated, NOT committed):** `src/CRTM_Forward_Module.f90` after line 961
now honors the previously-vestigial `Options%n_Legendre_Terms`:
```fortran
AtmOptics(nt)%n_Legendre_Terms = n_Full_Streams
IF ( Opt%n_Legendre_Terms > 0 ) AtmOptics(nt)%n_Legendre_Terms = Opt%n_Legendre_Terms
```
Default (`==0`) is byte-for-byte unchanged (config A below reproduces §7.3). Driver:
`build/stream_test/decouple_test.f90`. Three configs per DDA habit at the worst Water_Density cell:
- **A** = default auto (coupled: ~6 Legendre / 4 hemispheric streams)
- **B** = auto streams, **Legendre forced to 16** (DECOUPLED — 16 Legendre, *same* 4 streams)
- **C** = forced 16 streams (reference: 16 Legendre / 9 streams)

Decomposition of the baseline error A−C into a **Legendre part (A−B)** and a **stream part (B−C)**,
GMI 183±7 GHz (ch 13), WC=0.5:

| Habit | A−C (baseline) | A−B (Legendre) | B−C (streams) |
|-------|---------------:|---------------:|--------------:|
| SNOW default | 0.08 K | 0.08 K | 0.003 K |
| LargePlateAggregate | 1.13 K | 1.13 K | 0.004 K |
| EvansSnowAggregate | 2.45 K | 2.44 K | 0.003 K |
| SixBulletRosette | 2.67 K | 2.67 K | 0.004 K |

**Conclusion — the entire multi-Kelvin error is Legendre truncation, not stream count.**
- Config B (16 Legendre, **unchanged** 4-stream auto quadrature) recovers the reference to **within
  ~0.004 K**. The quadrature stream count (4→9 hemispheric) contributes essentially nothing (B−C).
- So the cheap fix wins: **use the richest available Legendre block (16) regardless of the stream
  count.** This is near-free — no extra streams, so no extra doubling-adding cost (which scales
  ~streams²–³); only the phase-matrix sum runs 16 vs 6 terms. "Raising the stream cap" is **not**
  needed for accuracy here and would just burn compute.
- This is exactly the RT4/DISORT separation of expansion order from quadrature order, confirmed
  numerically in CRTM.

Caveats / remaining work:
- Forward-only prototype. TL/AD/K carry the identical assignment (`CRTM_Tangent_Linear_Module.f90`
  ~1023, `CRTM_Adjoint_Module.f90` ~892, `CRTM_K_Matrix_Module.f90` ~1120) and need the same gating.
- The LUT only stores up to a 16-term block (`lOffset` {4,6,8,16}); 16 is the ceiling available
  today. Very peaked habits may want >16 terms — that needs a re-packed CloudCoeff (Layer 2), which
  also removes the {4,6,8,16} constraint so the truncation order can be chosen freely.
- For a permanent fix, drive the truncation from the CloudCoeff (store the needed order per entry,
  or always use the richest block) rather than from a user Option.

**Update (2026-06-01):** prototype edit **reverted** from `CRTM_Forward_Module.f90` and library
rebuilt (legacy default fully restored, per Design contract #1). The permanent home for "use the
full expansion, truncation LUT-driven" is the experimental scheme — see
**`CloudCoeff_Experimental_Schema_v1.md`** for the v1 format + dispatch design.
