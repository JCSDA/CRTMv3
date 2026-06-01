# Experimental CloudCoeff — Schema v1 (design spec)

**Status:** Design draft 2026-06-01. Companion to `Cloud_Optical_Properties_Redesign_Plan.md`
(see its Design contract + §7 findings). Implements Layers 1–2 as an **opt-in, experimental**
scheme; the legacy MIE_TAMU / DDA_ARTS path is untouched and remains the default.

## 0. Why this format (carry-over from findings)

- §7.4 proved the multi-Kelvin DDA error is **Legendre truncation, not stream count**, and that the
  fix is nearly free. So the format stores **one full / variable-length expansion with its own
  truncation order**, read by the solver independently of the stream count.
- The legacy file bakes the coupling in: `CloudCoeff_type` carries `n_Stream_Sets`, `n_Streams()`,
  `Legendre_Offset()` and packs Legendre into `{4,6,8,16}` blocks (`lOffset`). **All of that is
  removed here.**
- Legacy indexes solid by density-as-habit-proxy or Reff-only, scalar phase function
  (`n_Phase_Elements=1`), no temperature for solid. v1 replaces these with explicit **habit** +
  **PSD-moment (Dm, μ)** axes, **6 phase elements**, and **temperature for all phases**.

## 1. Scope of v1

- **Microwave only** (extend to IR/VIS later). Spectrum chosen by the active scheme.
- **Opt-in** via `Cloud_Model='CRTM-Exp'` (see §5). Default `'CRTM'` ⇒ legacy, byte-for-byte.
- **Host-transparent size input**: v1 keeps the existing cloud inputs (`Water_Content`, `Effective_Radius`,
  `Type`); the runtime maps `Effective_Radius → Dm` internally for the assumed PSD, and uses a single
  shape (`n_Mu=1`). The format already carries a `μ` axis so 2-moment (host-supplied `Dm,μ`) is a
  data-only upgrade later (needs a `CRTM_Cloud_type` field — deferred).
- **Phase matrix**: store all 6 Greek constants; scalar RT uses element 1 only, vector RT (the
  polarization roadmap) uses all 6. No format change needed to turn polarization on.

## 2. Physical model

Bulk properties are **pre-integrated** over a normalized-gamma PSD
`N(D) = N_w · f(μ) · (D/Dm)^μ · exp(−(4+μ)·D/Dm)` for each (habit, Dm, μ, T, frequency), from the
in-house DDA single-particle data. The scattering matrix is expanded in **generalized spherical
functions** (Wigner-d); stored coefficients are the standard 6 "Greek constants" for randomly
oriented particles with a symmetry plane:

| elem | constant | phase-matrix meaning |
|-----:|----------|----------------------|
| 1 | α₁(l) | P₁₁ (intensity) — the only one scalar RT needs |
| 2 | α₂(l) | P₂₂+P₃₃ |
| 3 | α₃(l) | P₂₂−P₃₃ |
| 4 | α₄(l) | P₄₄ |
| 5 | β₁(l) | P₁₂ = P₂₁ |
| 6 | β₂(l) | P₃₄ = −P₄₃ |

(Matches `Common_RTSolution.f90`'s α₁–α₄,β₁,β₂ convention, so the solver consumes them directly.)

## 3. netCDF dimensions (v1 sizes are placeholders)

```
n_Frequency      = 30      // MW LUT frequency nodes (interpolated to channels at runtime)
n_Dm             = 40      // mass-weighted mean diameter grid
n_Mu             = 1       // PSD shape; 1 in v1, grows for 2-moment
n_Temperature    = 5       // applies to ALL phases (incl. solid — new)
n_Habit          = 20      // explicit habit axis (DDA archive habits)
n_Legendre       = 64      // L_max: full expansion length (was effectively ≤16)
n_Phase_Elements = 6       // full Mueller (random orientation); 1 collapses to scalar
```

## 4. netCDF variables (CDL sketch)

```
// --- scheme identity ---
:Scheme        = "CRTM-Exp" ;     // explicit; loader matches Cloud_Model
:Release = 1 ; :Version = 1 ;
:PSD           = "normalized_gamma" ;
:Orientation   = "random" ;        // future: "azimuth_avg_canted" for Layer 4

// --- coordinate axes ---
double Frequency(n_Frequency) ;            // GHz
double Dm(n_Dm) ;                          // mass-weighted mean diameter, microns
double Mu(n_Mu) ;                          // gamma shape parameter (dimensionless)
double Temperature(n_Temperature) ;        // K

// --- habit axis + provenance/microphysics (Layer 3 hooks) ---
string Habit_Name(n_Habit) ;
int    Habit_Phase(n_Habit) ;              // 0=liquid, 1=frozen
double mD_a(n_Habit) ;                     // mass-dimension prefactor  m = a*D^b
double mD_b(n_Habit) ;                     // mass-dimension exponent
double Habit_Density(n_Habit) ;            // nominal/effective bulk density (diagnostic)

// --- bulk optical properties (PSD-integrated) ---
//     dims fastest->slowest as in CRTM CDL convention
double ke(n_Frequency, n_Temperature, n_Mu, n_Dm, n_Habit) ;   // mass extinction, m^2/kg
double ka(n_Frequency, n_Temperature, n_Mu, n_Dm, n_Habit) ;   // mass absorption, m^2/kg
double g (n_Frequency, n_Temperature, n_Mu, n_Dm, n_Habit) ;   // asymmetry (convenience)
double kb(n_Frequency, n_Temperature, n_Mu, n_Dm, n_Habit) ;   // mass backscatter (active)
// derived at runtime: ks = ke-ka ;  w = ks/ke   (storing ka avoids precision loss when w->1)

// --- phase-matrix expansion + its OWN truncation order (the decoupling) ---
int    n_Legendre_Eff(n_Frequency, n_Temperature, n_Mu, n_Dm, n_Habit) ; // effective order per entry
double pcoeff(n_Phase_Elements, n_Legendre, n_Frequency, n_Temperature, n_Mu, n_Dm, n_Habit) ;
```

Notes:
- `pcoeff(1, 0:L-1, ...)` is α₁ — collapsing to element 1 + a low `n_Legendre_Eff` reproduces the
  legacy scalar behavior (the contract's sanity bridge).
- `n_Legendre_Eff` is the order beyond which the δ-fit/GSF coefficients are negligible; the solver
  uses `min(n_Legendre_Eff, MAX_N_LEGENDRE_TERMS)` and is **never tied to the stream count**.
- v1 size estimate at the placeholder dims: `pcoeff` ≈ 6·64·30·5·1·40·20·8 B ≈ **370 MB**
  (cost drivers are `n_Mu` and `n_Phase_Elements`; scalar v1 ≈ 60 MB, `n_Mu=3` ≈ 1.1 GB).

## 5. Fortran derived type (sketch)

```fortran
TYPE :: CloudCoeff_Exp_type
  INTEGER(Long) :: Release = 1, Version = 1
  LOGICAL       :: Is_Allocated = .FALSE.
  CHARACTER(SL) :: Scheme = 'CRTM-Exp'
  ! dimensions
  INTEGER(Long) :: n_Frequency=0, n_Dm=0, n_Mu=0, n_Temperature=0, &
                   n_Habit=0, n_Legendre=0, n_Phase_Elements=0
  ! axes
  REAL(Double), ALLOCATABLE :: Frequency(:), Dm(:), Mu(:), Temperature(:)
  ! habit metadata
  CHARACTER(SL), ALLOCATABLE :: Habit_Name(:)
  INTEGER(Long), ALLOCATABLE :: Habit_Phase(:)
  REAL(Double),  ALLOCATABLE :: mD_a(:), mD_b(:), Habit_Density(:)
  ! bulk optics : (Frequency, Temperature, Mu, Dm, Habit)  [store ka not w; w=ks/ke derived]
  REAL(Double), ALLOCATABLE :: ke(:,:,:,:,:), ka(:,:,:,:,:), g(:,:,:,:,:), kb(:,:,:,:,:)
  ! per-entry truncation + expansion
  INTEGER(Long), ALLOCATABLE :: n_Legendre_Eff(:,:,:,:,:)
  REAL(Double),  ALLOCATABLE :: pcoeff(:,:,:,:,:,:,:)   ! (Phase,Legendre,Freq,T,Mu,Dm,Habit)
END TYPE
```

New modules (legacy `CloudCoeff_*` untouched): `CloudCoeff_Exp_Define.f90`,
`CloudCoeff_Exp_netCDF_IO.f90`.

## 6. Dispatch flow (opt-in, self-contained)

```
CRTM_Init(..., Cloud_Model='CRTM-Exp', CloudCoeff_File='CloudCoeff_Exp_*.nc')
   └─ CRTM_CloudCoeff_Load(Cloud_Model, file, ...)
        IF Cloud_Model=='CRTM-Exp' : load CloudCoeff_Exp_netCDF_IO -> module CloudC_Exp
                                     set module flag  Active_Cloud_Scheme = CRTM_EXP
        ELSE                       : legacy load (unchanged)

CRTM_Compute_CloudScatter(...)            ! per channel
   IF Active_Cloud_Scheme==CRTM_EXP :
        Get_Cloud_Opt_MW_Exp(habit=Type, Dm=f(Effective_Radius), mu=Mu_default, T, freq)
            -> ke, w, g, kb, pcoeff(1:6, 0:Le-1), Le=n_Legendre_Eff
        CScat%n_Legendre_Terms = MAX(CScat%n_Legendre_Terms, Le)     ! LUT-driven, NOT streams
        accumulate pcoeff(1:6,...) into AtmOptics%Phase_Coefficient
   ELSE : legacy Get_Cloud_Opt_MW (unchanged)
```

The stream count (`RTV%n_Streams`) still comes from the Mie heuristic / `Options%n_Streams` — §7.4
showed that's adequate; only the truncation order is freed.

## 7. Enabling changes outside the LUT (small, additive)

1. **`CRTM_Parameters.f90`**: raise `MAX_N_LEGENDRE_TERMS` 16 → 64 (sizes `RTV%Pleg`,
   `AtmOptics%Phase_Coefficient`; memory grows ∝ this — verify footprint). `MAX_N_PHASE_ELEMENTS`
   is already 6. **Legacy results unchanged** (legacy still requests ≤16).
2. **`CRTM_AtmOptics`**: ensure the combined `n_Legendre_Terms = MAX` over constituents
   (cloud-exp may exceed aerosol/molecular) and the phase-matrix sum runs to that order.
3. **`CRTM_CloudScatter`**: add the `Active_Cloud_Scheme` branch + `Get_Cloud_Opt_MW_Exp`
   (forward first; TL/AD/K mirror later).
4. **Cloud input**: v1 needs none. 2-moment (v1.1) adds `Dm`/`μ` (or number conc.) to
   `CRTM_Cloud_type` — `Effective_Variance` exists but is unused and could carry μ.

## 8. Decisions — LOCKED 2026-06-01

1. **Size coordinate = `Dm`** (mass-weighted mean diameter). Runtime maps the host `Effective_Radius`
   → `Dm` internally (per habit/μ), so v1 is host-transparent and the file is 2-moment-ready. Reff
   as the axis was rejected (conflates PSD scale with moment definition).
2. **`MAX_N_LEGENDRE_TERMS = 64`.** The archive spans to 874 GHz (ICI/sub-mm), where size parameters
   give convergence orders ~40–60; 64 is safe headroom. Runtime cost is negligible (`Pleg`,
   `Phase_Coefficient` working sets are ~100s of KB); only the LUT grows, and the near-zero
   high-order tails compress under netCDF deflate. Per-entry `n_Legendre_Eff` means low-freq/small
   entries still use few terms. Legacy unaffected (requests ≤16). *Veto only if file size is a hard
   constraint → fall back to 32 (MW ≤183 GHz).*
3. **Store `ke` + `ka` + `kb` (+ `g`); derive `ks=ke−ka`, `w=ks/ke`.** Revised from "store `w`":
   storing absorption directly avoids precision loss in the emission term when `w→1`
   (matches legacy, which already carries `ka_*`/`kb_*`).
4. **Frequency grid = native DDA nodes (~25–30), log interpolation** to channel frequencies — matches
   how the data was computed; no resampling artifacts.
5. **Habit key = existing extended cloud-type integers** (`EvansSnowAggregate`, …) for dispatch;
   carry `Habit_Name` string for provenance/validation.

Additional locked points:
- **`n_Phase_Elements = 6` stored and populated** from the DDA Mueller matrices (the archive has
  them). v1 scalar-RT validation consumes element 1 only; elements 2–6 light up with the vector-RT
  polarization work — no format change to flip on.
- **`n_Mu = 1`** in v1 (single assumed shape, host-transparent). 2-moment (`n_Mu>1`) is a data-only
  upgrade that additionally needs a `CRTM_Cloud_type` shape field — deferred to v1.1.
- **netCDF deflate enabled** on `pcoeff` (crushes the near-zero high-order tails).
```
