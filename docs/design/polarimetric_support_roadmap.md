# Polarimetric (n_Stokes > 1) support: rigorous path forward

Written 2026-07-30, after the surface Stokes-basis fix (`4fe0191`, merged
`8bf010e`). Every factual claim below carries a file and line citation and was
verified against the code rather than against comments or prior documentation.
Claims that were *not* verified are listed separately in "Open questions", and
should be treated as unknown rather than as likely-fine.

## Why this document exists

A defect in the surface-to-solver handoff survived from the first v3.0
instantiation until 2026-07-30. On the `n_Stokes > 1` path the microwave surface
optics were handed to the radiative transfer solver in the (V,H) basis the
surface models produce, while the solver consumes a Stokes source vector. The
vertical emissivity was read as Stokes I and the horizontal emissivity as
Stokes Q. At 23.8 GHz and 45 degrees over ocean the solver received I = 0.547
and Q = 0.342 where the correct values are I = 0.444 and Q = 0.103, so the
polarization difference entered 3.3 times too large.

It survived that long for a specific and correctable reason. Every test of the
vector path is a self-consistency test: tangent-linear against finite
differences, the adjoint dot-product identity, and K against AD. All three ask
whether the derivative code is the true derivative of the forward code. None
asks whether the forward code is physically correct. A forward model with the
wrong surface basis is smooth and perfectly differentiable, so it passes all
three. `test_VectorRT_TLADK` passed both before and after the fix.

The organizing principle of everything below follows from that:

> **You cannot validate the middle of a chain without trusted ends.** Establish
> external ground truth first, pin the convention at every interface second,
> then fix inward from the boundaries. Every phase exits on a test that
> demonstrably fails against the unfixed code.

That last clause is not a formality. The fix committed on 2026-07-30 was
accepted only after reverting it, rebuilding, and confirming the new test failed
all four of its assertions by 0.10 to 0.24 in emissivity units.

## Established facts

These are settled and should not be re-litigated.

**The solver state vector is standard Stokes (I, Q, U, V), with slot 1 as total
intensity.** Two independent lines of evidence agree.

1. In the v3.0 alpha `ADA_Module.f90`, unpolarized layer thermal emission is
   written with a stride, `DO i = 1, nZ, RTV%n_Stokes`, placing the Planck source
   in the first slot of each angle and leaving the rest zero. The reflected
   cosmic background uses the same stride. Unpolarized radiation in standard
   Stokes is (B, 0, 0, 0); in modified Stokes (Iv, Ih, U, V) it would be
   (B/2, B/2, 0, 0) and would require filling two slots. It does not.
2. The scalar branch of `CRTM_SfcOptics.f90` states the mapping outright:
   `SECOND_STOKES_COMPONENT` is `0.5*(Emissivity(:,1) - Emissivity(:,2))`,
   unpolarized is `0.5*(Emissivity(:,1) + Emissivity(:,2))`, `VL_POLARIZATION`
   is `Emissivity(:,1)` and `HL_POLARIZATION` is `Emissivity(:,2)`. Component 1
   is therefore eV and component 2 is eH, and I and Q are their half-sum and
   half-difference.

**The handoff performs no conversion of its own.** `Reshape_Surf_Opt`
(`src/CRTM_Utility/CRTM_Utility.f90:281`) flattens `emissivity(angle, m)`
directly into the source vector consumed by `CRTM_ADA` and `CRTM_Emission`
(`src/RTSolution/CRTM_RTSolution.f90:249`). Whatever sits in
`SfcOptics%Emissivity(:,1:n_Stokes)` is what the solver integrates.

**The surface conversion is now correct** in the forward, tangent-linear and
adjoint routines, and is pinned by `test_VectorRT_SurfaceBasis`, which recovers
eV and eH by temporarily forcing the channel to pure V and pure H polarization
on the scalar path and asserts the vector path returns their half-sum and
half-difference, together with the reflectivity identities. It requires no cloud
lookup table and no reference radiances.

**Nothing at `n_Stokes = 1` is affected.** All executable lines of the fix are
inside the `ELSE` of `IF (SfcOptics%n_Stokes == 1)`. `CRTM_Options_type%n_Stokes`
and `RTV_type%n_Stokes` both default to 1, the propagation is guarded by
`IF (Opt%n_Stokes > 0)`, and all four entry modules set `SfcOptics%n_Stokes` from
`RTV%n_Stokes`. `SfcOptics%n_Stokes` does default to 0 in its own type
definition, which would route into the vector branch, but the only library
caller of `CRTM_Compute_SfcOptics` is `Common_RTSolution.f90:370`, reached solely
through those entry modules. The second caller,
`src/RTSolution/UWisc_SOI/SOI_CRTM_Forward_Module.f90:664`, never sets
`n_Stokes`, but it is not compiled into the library: it is absent from
`src/CMakeLists.txt` and contributes zero symbols to `libcrtm.so`.

## Known defects and gaps

Each of these was verified on 2026-07-30. None is speculative.

| # | Gap | Evidence |
|---|-----|----------|
| 1 | Reported radiance is total intensity, not the channel's polarized measurement. There is no projection of the Stokes vector onto the channel polarization. | `src/RTSolution/Common_RTSolution.f90:1279`, `RTSolution%Radiance = RTSolution%Stokes(1)` inside the `n_Stokes > 1` block |
| 2 | U and V are computed by the surface model and then discarded. FastemX receives the full four-component slice, but the coverage aggregation copies only components 1 and 2 into the local array, which is zero-initialised. The solver always sees U = V = 0 at the surface. | Full slice passed at `CRTM_MW_Water_SfcOptics.f90:276`; aggregation at `CRTM_SfcOptics.f90:546, 575, 604, 633`; zero-init at `:509` |
| 3 | The fractional-cloud K seed is defective for vector runs. Inside a `DO ks = 1, n_Stokes` loop it assigns to a scalar with no `ks` index, so each iteration overwrites the last and only the final Stokes component survives. The correct line is present, commented out, directly above. | `src/CRTM_K_Matrix_Module.f90:1335` (loop), `:1339` (defect), `:1336` (commented correct form) |
| 4 | Selecting SOI with `n_Stokes > 1` silently yields ADA. The vector branch precedes the algorithm dispatch entirely, so `RT_Algorithm_Id` is never consulted. | `src/RTSolution/CRTM_RTSolution.f90:248` versus the `ELSE IF (RTV%Scattering_RT)` dispatch that follows |
| 5 | Aerosols contribute no polarization. The shipped `AerosolCoeff.nc` carries `n_Phase_Elements = 1`, and the scatter routine fills `MIN(n_Phase_Elements, AeroC%N_PHASE_ELEMENTS)`. A vector run therefore mixes polarized cloud scattering with unpolarized aerosol scattering. | Coefficient file dimension; `src/AtmScatter/CRTM_AerosolScatter.f90:316` |
| 6 | Polarimetric support is microwave-only. The coupled-polarization branch exists solely in the microwave section; infrared and visible set component 1 only. | Section banners at `CRTM_SfcOptics.f90:517`, `:852`, `:992`; branch spans `:655` to `:845` |
| 7 | U and V are never exercised by any test. `test_VectorRT_TLADK` runs at `n_Stokes = 2` with a scalar control. | test header and setup |
| 8 | No polarized radiance has ever been compared against a reference outside CRTM. | Whole-repository survey of tests touching `n_Stokes` |

## Open questions

Unknown, and to be resolved in Phase 1 rather than assumed.

- **Frame convention at the surface.** Whether the surface V/H frame coincides
  with the radiative transfer meridional plane at non-zero azimuth. If it does
  not, a rotation is required, and its absence is invisible at nadir, which is
  where most casual testing happens.
- **Phase matrix ordering.** `CloudCoeff_Exp_Define.f90:12` documents six phase
  elements as "alpha1..alpha4, beta1, beta2", the Hovenier and Mishchenko
  expansion convention for randomly oriented particles with a plane of symmetry.
  That convention pairs correctly with standard Stokes in the meridional frame,
  so it is consistent with the solver. But it is one comment line and nothing
  verifies the code honours that ordering.
- **Reflectivity structure.** The surface models fill only diagonal angle terms,
  `Reflectivity(i,j,i,j)`. Whether a polarimetric surface requires off-diagonal
  angle coupling is unresolved.
- **Azimuthal Fourier decomposition.** Components 1 and 2 are accumulated with a
  cosine weighting and components 3 and 4 with a sine weighting
  (`Common_RTSolution.f90` around `:1276`). This is the standard convention, but
  it has not been verified against the solver's internal expansion.

## Phased plan

### Phase 0. External ground truth

Nothing downstream is meaningful without a reference that is not CRTM. Begin
with closed-form cases, which are unambiguous and cheap: an isothermal
non-scattering slab over a specular Fresnel surface has an analytic (I, Q), and
single-scatter Rayleigh has known polarization. Follow with an independent code,
either RTTOV-SCATT's polarimetric path or a PolRadtran or VDISORT reference, for
scattering configurations.

*Exit:* a registered test comparing CRTM vector output against a closed-form
solution for at least one clear-sky and one single-scatter configuration.

### Phase 1. Convention audit, written down and pinned

Document, then verify against code, the convention at every interface: surface
model output, solver state vector, phase matrix expansion ordering and reference
frame, azimuthal Fourier assignment, and the surface-to-meridional frame
relationship. Resolve all four open questions above.

*Exit:* a design note plus assertion tests pinning each interface independently,
so that a later change violating one fails immediately rather than silently.

### Phase 2. Complete the surface

Carry U and V through the coverage aggregation. Note that the naive change from
`1:2` to `1:nL` breaks the scalar path, which requires both V and H even at
`nL = 1` for its polarization mixing; the correct form is `1:MAX(2,nL)`, applied
at twelve sites across the forward, tangent-linear and adjoint routines. Add the
frame rotation if Phase 1 shows one is needed, and resolve the reflectivity
structure question.

*Exit:* `test_VectorRT_SurfaceBasis` extended to four Stokes components, plus a
test proving FastemX's U and V survive to the solver input.

### Phase 3. Fix the output side

Project the emergent Stokes vector onto the channel polarization for `%Radiance`
and brightness temperature, and rework the adjoint seeding, which currently
assumes `Radiance` is `Stokes(1)` and seeds from `Stokes(1:2)`. This is the
largest single item, because it reaches the brightness-temperature adjoint and
the K seed. A forward-only change here would recreate exactly the class of
forward-versus-Jacobian inconsistency this document exists to prevent.

*Exit:* tangent-linear against finite differences, adjoint dot-product, and K
against AD at `n_Stokes = 4` on a polarized channel, plus a cross-check that a
channel with pure V or pure H polarization agrees between the scalar and vector
paths.

### Phase 4. Interior and dispatch

Fix the fractional-cloud K seed (gap 3). Decide the aerosol story (gap 5), which
is a coefficient-generation question as much as a code question. Make SOI
combined with `n_Stokes > 1` an explicit error rather than a silent substitution
(gap 4).

*Exit:* fractional-cloud K against AD at `n_Stokes > 1`; explicit failure tests
for the unsupported combinations.

### Phase 5. Scope decision on infrared and visible

Either extend the coupled branch beyond microwave or state plainly in the
documentation that polarimetric support means microwave. The present situation,
where the capability is described generally but implemented for one sensor
class, is the part most likely to mislead.

## Immediate actions, ahead of the phases

**Convert silent wrongness into loud refusal.** Have initialization or the
forward entry point reject, or at minimum warn on, the combinations now known to
be unsupported: `n_Stokes > 1` together with a single-phase-element aerosol
coefficient, with SOI selected, with an infrared or visible sensor, or with
fractional cloud on the K path. This is contained, it protects users throughout
the work above, and it costs a fraction of any single phase.

**Correct the release documentation.** The v3.2.0 notes list vector radiative
transfer among the highlights in terms that read as readiness. The known-issues
section should carry gaps 1 through 6.

## Sequencing

Phases 0 and 1 are prerequisites and should not be compressed under schedule
pressure, because skipping exactly that groundwork is what allowed the original
defect to persist. Phases 2 and 3 are the substance and are independent enough
to proceed in parallel. Phase 4 is contained and can be done at any point after
Phase 1. Phase 5 is a product decision rather than an engineering one.

A closing note on test design. For every phase, the acceptance criterion is a
test that fails against the unfixed code. Self-consistency checks remain
valuable and should be kept, but they must never again be the only coverage of a
physics path, because they are structurally incapable of detecting the class of
defect described at the top of this document.
