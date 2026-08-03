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

**The surface Stokes frame is the meridional frame, and no rotation is
required.** Resolved 2026-07-31; this was the largest open question and the
answer is negative, so no code change follows from it. Three independent lines
of evidence agree.

1. The solver's Stokes basis is the per-direction meridional frame. The
   polarized phase matrix is assembled from the generalized spherical functions
   `Pplus` (R_l^m) and `Pminus` (T_l^m) at `Common_RTSolution.f90:2017` and
   `:2019`, which is the standard meridional-frame azimuthal Fourier expansion
   with the scattering-plane rotations folded in analytically. Consistent with
   that, a whole-tree search of `src/RTSolution` and `src/SfcOptics` finds no
   rotation matrix of any kind.
2. The surface models refer their vector to the plane of incidence, which for a
   plane-parallel atmosphere is the same plane: both are spanned by the local
   zenith and the propagation direction, so the angle between them is zero by
   construction at every azimuth. The code shows this directly in the parity of
   the azimuthal model. `Azimuth_Emissivity_Module.f90:139-142` builds the
   vertical and horizontal components from `cos(m*phi)` and the third and
   fourth from `sin(m*phi)`, where `phi` is the relative azimuth; PARMIO does
   the same at `PARMIO_Azimuth_Module.f90:79-91`. Even V and H with odd U and V
   is exactly the signature of a reference frame lying in the view plane, since
   reflecting the scene through that plane maps `phi` to `-phi` and must leave
   I and Q alone while flipping the handedness of U and V.
3. Measured, not argued. `test_VectorRT_SurfaceFrame` mirrors an ocean scene
   through the view plane at 45 degrees and recovers I and Q bit-identical and
   U and V exactly negated. `test_VectorRT_StokesOutput` repeats the
   measurement through the full scattering solver at `n_Stokes = 4` and gets I
   and Q bit-identical with U and V odd to 3.4e-13 against an intensity scale
   of 2.2e-1.

The residual risk in this area was not the frame but the *sign* convention of
the third Stokes component between FASTEM and the solver. **RESOLVED 2026-08-02**:
This was tested against RTTOV's FASTEM5 implementation and the U and V signs
match perfectly at all azimuths.

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

Each of these was verified on 2026-07-30, and the 2a to 2d entries on
2026-07-31. None is speculative.

| # | Gap | Evidence |
|---|-----|----------|
| 1 | ~~Reported radiance is total intensity, not the channel's polarized measurement.~~ **FIXED 2026-07-31.** `RTSolution%Radiance` is now the emergent Stokes vector projected onto the channel polarization, in forward, tangent-linear and adjoint, and the brightness temperature follows from it. The weights are derived from the scalar branch of `CRTM_SfcOptics` rather than from first principles, so the two paths agree by construction: every case there is `a*eV + b*eH (+ c*e3 + d*e4)`, and with `eV = I+Q`, `eH = I-Q` that is `w = (a+b, a-b, c, d)`. `%Stokes` is untouched and remains the physical (I,Q,U,V). Three things had to move with it: `Pre_Process_RTSolution_AD` mirrored `%Radiance` into `%Stokes(1)` as a workaround for the seed being ignored, which now double counts and is removed; all four entry modules re-asserted `Radiance = Stokes(1)` after their fractional-cloud combine, silently undoing the projection, and now combine the projected radiance linearly instead; and that combine needed its own adjoint split, without which the cloudy column took the whole seed. Verified against the unfixed code by reverting: the V and H checks fail by 2.25e-3, the full polarization difference. | Projection and its transpose in `Common_RTSolution.f90` (`Channel_Polarization_Weights`); fractional combine in all four entry modules. Pinned by `test_VectorRT_ScalarLimit` (vector `%Radiance` equals the pure-V and pure-H scalar runs at 1.4e-16) and by four checks in `test_VectorRT_TLADK`: TL against finite difference on `%Radiance`, and the adjoint dot product seeded through `%Radiance`, each overcast and fractional |
| 2 | ~~U and V are computed by the surface model and then discarded by the coverage aggregation.~~ **FIXED 2026-07-31.** Corrected only at the microwave *water* sites in the forward, tangent-linear and adjoint routines, not at twelve sites as originally prescribed. The land, snow and ice models write components 1 and 2 and never define 3 and 4, and nothing zeroes `SfcOptics%Emissivity` between calls, so extending their aggregation would have read stale values from a previous water channel. Water is the only microwave surface with a polarimetric model, and the accumulator is already zeroed, so leaving the other three at `1:2` is also the physically correct contribution. | Aggregation now `1:nS` with `nS = MAX(2,nL)` at `CRTM_SfcOptics.f90` water blocks; land/snow/ice fill only 1:2, e.g. `CRTM_MW_Land_SfcOptics.f90:257`; no entry zeroing in any of the four models. Pinned by `test_VectorRT_SurfaceFrame` |
| 2a | U and V were then annihilated a second time, on output. `RTSolution%Stokes(3:4)` accumulated with a weight of `SIN(mth_Azi*dphi)`, and `n_Azi` is set above zero only for visible channels while the polarimetric surface branch exists only for microwave, so `mth_Azi` is always 0 on every path where `n_Stokes > 1` is meaningful and the weight was always `SIN(0)`. **FIXED 2026-07-31** in the forward, tangent-linear and adjoint accumulations. | `Common_RTSolution.f90` accumulation blocks; `n_Azi` at `CRTM_Forward_Module.f90:993` versus `:1011`. Pinned by `test_VectorRT_StokesOutput` |
| 2b | The shipped default microwave water model has no third or fourth Stokes azimuth model at all. `CRTM_Init` defaults to FASTEM6, whose azimuth routine parameterises the vertical and horizontal components only and returns components 3 and 4 as identically zero. A polarimetric run has a real surface U and V only on FASTEM4/5 or PARMIO. Not a defect in itself, but it means gap 2 was invisible in the default configuration, and it constrains what a polarimetric user can actually run. | `CRTM_LifeCycle.f90:802`; `Azimuth_Emissivity_F6_Module.f90:187-188` versus `Azimuth_Emissivity_Module.f90:139-142` and `PARMIO_Azimuth_Module.f90:79-91` |
| 2c | ~~`MWwaterCoeff_File` does not select the microwave water model. The file-based load is commented out and selection is by the `MWwaterCoeff_Scheme` string; the filename argument survives only in a diagnostic message. Passing `MWwaterCoeff_File='FASTEM4...'` silently leaves FASTEM6 loaded.~~ **FIXED 2026-07-31.** The argument now selects the model. The shipped names are exactly `<SCHEME>.MWwater.EmisCoeff.<ext>`, so the scheme is recovered as the leading component of the base name, derived before the `File_Path` join so it works with or without a path. An explicit `MWwaterCoeff_Scheme` remains the direct selector and still wins; a disagreement between the two is now reported rather than resolved silently. The load message named the file it never read and now names the model actually loaded. The dead file-based loader was left commented rather than resurrected: it was only ever needed for the FASTEM5 binary tables, and reviving it would change the scalar path for no benefit here. Backward compatible, because the default name derives to FASTEM6 and no in-repo caller passes the argument at all. This mattered beyond CRTM: JEDI/UFO selects the model through exactly this argument, building `TRIM(MWwaterCoeff)//".MWwater.EmisCoeff.nc"` from its own yaml key (`ufo_crtm_utils_mod.F90:741`), so a JEDI user writing `MWwaterCoeff: FASTEM4` to obtain a polarimetric surface was given FASTEM6 and U = 0. Verified by disabling the fix and rebuilding: the FASTEM4 request returns U = V = 0 exactly and the test fails. | Derivation and helper in `CRTM_LifeCycle.f90`; pinned by `test_MWwaterCoeff_FileSelects`, which asserts both directions so it cannot pass by always loading FASTEM4 |
| 2d | On the **scalar** path, a channel declared `THIRD_STOKES_COMPONENT` or `FOURTH_STOKES_COMPONENT` receives an identically zero surface emissivity. Those cases read `Emissivity(:,3)` and `(:,4)` from the same accumulator, which at `nL = 1` still stops at component 2. Latent only: no sensor in the shipped test suite uses polarization 3 or 4 (a scan of all `testinput/*.SpcCoeff.nc` finds only 1, 5, 6, 9, 10, 11, 13, 14). Left unfixed deliberately, being a scalar-path behaviour change outside this work's scope; the one-line fix is to raise the water aggregation floor from `MAX(2,nL)` to 4. | `CRTM_SfcOptics.f90:683-691` |
| 3 | **FIXED 2026-07-31**, after being mis-diagnosed twice. Two real defects sat behind it, both on the microwave path and both invisible to every self-consistency test. First, the fractional-cloud adjoint and K seeds handled `%Radiance` alone while the forward combined every Stokes component, so the clear column was never seeded (the vector path reads `%Stokes`, never `%Radiance`) and the cloudy column was never scaled by the cloud cover. Second, `RTV_Clear%n_Stokes` was set only in `CRTM_Forward_Module`, so the tangent-linear, adjoint and K-matrix ran their clear column **scalar** while the forward ran it vector, and the two disagreed on the forward radiance itself. The adjoint dot-product identity measured the first as a factor of two (0.99 relative) and the second as 5.2e-3; TL-vs-FD and K-vs-AD passed throughout both. Now 1.7e-15. Original mis-diagnosis follows for the record. | Seeds at `CRTM_Adjoint_Module.f90:1262`, `CRTM_K_Matrix_Module.f90:1492`; `RTV_Clear%n_Stokes` now set in all four entry modules. Pinned by the fractional-cloud block of `test_VectorRT_TLADK` |
| 3-orig | ~~The fractional-cloud K seed overwrites a scalar inside a `DO ks` loop.~~ **MIS-DIAGNOSED; corrected 2026-07-31.** There is no overwrite: an `IF( ks == 1 )` guard sits on that assignment (`:1338`) and has since April 2025, and the whole block is inside a visible/ultraviolet sensor test so it is unreachable on the microwave-only polarimetric path. The real gap is elsewhere and is a forward-versus-adjoint asymmetry: the K-matrix **forward** fractional-cloud combine handles every Stokes component (`:1413-1419`) but the **adjoint seed** on the microwave path handles `%Radiance` alone (`:1496-1500`), so the polarized clear/cloudy sensitivities are never propagated. It cannot be fixed on its own, because the clear-sky half has no vector solver to seed from; see gap 9. | `src/CRTM_K_Matrix_Module.f90:1331-1342` (VIS/UV block, guard at `:1338`), `:1413-1419` (forward Stokes combine), `:1492-1500` (adjoint seed, scalar only) |
| 9 | ~~**The clear-sky / non-scattering path has no vector solver.**~~ **FIXED 2026-07-31** by `CRTM_Emission_Stokes` and its tangent-linear and adjoint siblings. In the absence of scattering the atmosphere is polarization neutral, so for k >= 2 the whole solution is a surface boundary value transported upward with no source of its own: `S_k(sfc) = e_k*B + R_k1*D`, then multiplied by the layer transmittances. `test_VectorRT_ScalarLimit` now passes at 1.1e-16 in I and 8.2e-17 in Q and is registered; `test_VectorRT_TLADK` gained a clear-sky block covering TL against finite difference, the adjoint dot product and K against AD on the new path. One trap worth recording: `CRTM_Emission_AD` **zeroes** `T_OD_AD`, `Planck_Surface_AD`, `emissivity_AD` and `reflectivity_AD` on entry, so polarized adjoint contributions accumulated before that call are erased. Accumulating them directly left the adjoint non-transpose at 6.4e-5 while TL-vs-FD and K-vs-AD both still passed; only the dot-product identity caught it. They are held in locals and added back afterwards. Original description follows. | |
| 9-orig | **The clear-sky / non-scattering path has no vector solver.** With `n_Stokes > 1` and no significant scattering, the dispatch calls the scalar `CRTM_Emission` with the flattened (angle, Stokes) arrays. With `n_Angles = 1` the flattening happens to place Stokes I first, so `emissivity(n_Angles)` is e_I and `reflectivity(1,1)` is R_II and the **intensity is correct to 1.1e-16**, measured. But no Q, U or V is produced, and `Assign_Common_Output` fills `Radiance(1)` only, leaving `Radiance(2:n_Stokes)` read before assignment. A clear-sky polarimetric run, which is exactly what ocean wind-vector retrieval needs, returns Q = 0. An earlier claim that the intensity was about 9 percent high was an artefact of `test_VectorRT_ScalarLimit` never copying `Atm(2) = Atm(1)`, not a CRTM defect; that test bug is fixed. Fixing this needs a vector emission path in forward, tangent-linear and adjoint, and it is the prerequisite for gap 3. | `src/RTSolution/CRTM_RTSolution.f90:275-292` (dispatch), `Emission_Module.f90:139-141` (scalar surface boundary), `Common_RTSolution.f90` `Radiance(1)` assignment in the emission branch. Measured by `test_VectorRT_ScalarLimit` |
| 4 | ~~Selecting SOI with `n_Stokes > 1` silently yields ADA.~~ **FIXED 2026-07-31.** Now an explicit failure with a message naming the alternative, rather than a substitution the caller cannot see. Refusal was chosen over a warning because the caller receives a different algorithm's numbers under the name of the one they asked for. Verified by reverting the guard and rebuilding: without it the SOI call returns SUCCESS. | Guard at the head of the `n_Stokes > 1` branch, `src/RTSolution/CRTM_RTSolution.f90`. Pinned by `test_VectorRT_Unsupported`, which also runs an RT_ADA vector control so the test cannot pass by rejecting everything |
| 5 | Aerosols contribute no polarization. The shipped `AerosolCoeff.nc` carries `n_Phase_Elements = 1`, and the scatter routine fills `MIN(n_Phase_Elements, AeroC%N_PHASE_ELEMENTS)`. A vector run therefore mixes polarized cloud scattering with unpolarized aerosol scattering. | Coefficient file dimension; `src/AtmScatter/CRTM_AerosolScatter.f90:316` |
| 6 | ~~Polarimetric support is microwave-only.~~ **WRONG IN BOTH DIRECTIONS; corrected 2026-08-03 by measurement.** The claim was inferred from the surface code alone and never tested. Two things are true instead. First, the **visible already produces real polarized radiance** and always did: `CRTM_MoleculeScatter.f90:161-172` fills the polarized Rayleigh expansion coefficients under an `n_Stokes` gate, with matching tangent linear and adjoint, and visible/ultraviolet is the only sensor class where the azimuthal Fourier machinery is switched on (`CRTM_Forward_Module.f90:1065`). Measured on `v.abi_g18` clear-sky: Q/I of 47.8 percent at 0.47 um falling to 37.1 percent at 0.86 um, U exactly zero at relative azimuth 0 and 6.4e-17 at 180 while reaching 0.437 at 90, V below 1.9e-14, and a vector-versus-scalar intensity difference of up to 0.34 percent. So warning that the visible is unsupported would have been wrong. What the visible lacks is a polarized *surface*, not polarized radiative transfer. Second, the **infrared** returned exactly zero on all channels at all azimuths, confirming it had no polarization source at all; that half is now addressed, see gap 10. | Measured by probe against a clean build; the visible source is `CRTM_MoleculeScatter.f90:161-172` and its TL at `:279-290` |
| 10 | ~~The infrared has no polarized surface.~~ **FIXED 2026-08-03** on `feature/btj_polarimetric_ir_surface`, deliberately held out of v3.2.0. `CRTM_IR_Water_SfcOptics.f90` already computed `rv` and `rh` separately from the Wieliczka and Friedman seawater refractive index and then averaged them away on the next line, the same shape as the microwave defect that started this work. The V/H *shape* now comes from Fresnel and the *magnitude* from the validated WuSmith or Nalli lookup table, so `e_I = (e_V+e_H)/2` is the table value by construction and reduces to `e_Q = e*(rh-rv)/(2-rv-rh)`, exactly zero at nadir. Note the convention differs from microwave on purpose: infrared models write **(I,Q)** rather than (V,H), because an infrared channel has no polarization mixing and component 1 is already the intensity emissivity the scalar path consumes. Infrared land, snow and ice stay at component 1: they have no Fresnel basis, and "unpolarized rough emitter" is the honest statement. Polarized infrared cloud scattering remains refused, see gap 11. Measured on `abi_g18` at 60 degrees: window-channel Q/I = 2.2 percent, about 1.4 K of V minus H brightness temperature, Q strictly increasing with view angle, and the opaque 6.2 um channel acquiring 1.1e-16. | `CRTM_IR_Water_SfcOptics.f90` forward, TL and AD; water-only aggregation in `CRTM_SfcOptics.f90`. Pinned by `test_VectorRT_IR_Surface`, verified to fail against the unfixed code |
| 11 | **A scalar phase lookup table cannot stand in for a polarized one, and the difference is not conservative.** With one phase element, `Phase_Coefficient(l,2:6,k)` are all zero and `Pff` is zeroed before assembly, so the assembled 4x4 block has only its (1,1) entry nonzero. That is a *perfect depolarizer*: every scattering event annihilates Q, U and V while scattering I normally. It is not a neutral placeholder but a strong and false physical claim, which is why the existing guard refusing `n_Stokes > 1` with a sub-six-element cloud table is correct and must stay. Recorded because "just use the scalar tables for now" is the obvious suggestion and it is the wrong one. | Assembly at `Common_RTSolution.f90:2319` and `:2341-2400`; guard at `CRTM_Forward_Module.f90:705`, mirrored at TL `:693`, AD `:642`, K `:739` |
| 7 | ~~U and V are never exercised by any test. `test_VectorRT_TLADK` runs at `n_Stokes = 2` with a scalar control.~~ **NO LONGER TRUE, struck 2026-07-31.** `test_VectorRT_TLADK` now runs at `n_Stokes = 4` and differentiates dU and dV, including dU/d(wind direction). `test_VectorRT_Physics` checks U and V vanish at relative azimuth 0 and 180, `test_VectorRT_SurfaceFrame` pins their odd parity on both surface backends, `test_VectorRT_StokesOutput` proves they survive the Fourier accumulation, and `test_VectorRT_StokesSign` pins their sign. Retained here only so the claim is not re-derived from an older revision. | superseded by the verification table below |
| 8 | No polarized **radiance** has been compared against a reference outside CRTM. Narrowed 2026-08-02: the polarimetric *surface* now has an external check, since the third and fourth Stokes signs were compared against RTTOV's FASTEM5 and match at every relative azimuth. That closes the sign question but validates one interface, not the emergent radiance, and it does not exercise the solver, the phase matrix or the transport at all. | Whole-repository survey of tests touching `n_Stokes`; the RTTOV comparison is external to the repository and is recorded in `docs/design/polarimetric_conventions.md` rather than as a registered test |

## Verification status, 2026-07-31

What has been measured, as distinct from argued. Each entry names the
instrument so it can be re-run.

| Claim | Evidence |
|-------|----------|
| The default scalar path is unchanged, bit for bit | `dump_scalar_fullprec` built in this tree and at the branch point, 90 values across clear, overcast and fractional-cloud scenes, identical to all 17 digits. Necessary because the regression suite compares at `DEFAULT_N_SIGFIG` (= `SP_N_SIGFIG`, about six figures) and cannot see a change in the last bits |
| The full four-component Jacobian chain is correct | `test_VectorRT_TLADK` at `n_Stokes = 4`: TL against finite difference on dI, dQ, dU and dV, all non-degenerate (dU/dWC = 2.7e-6 agreeing to 6.1e-12), the four-component adjoint dot product exact, and K against AD. Before this, every Jacobian check ran at `n_Stokes = 2`, so the (1,3) (3,1) (2,3) (3,2) (3,3) (2,4) (4,2) (3,4) (4,3) (4,4) blocks and the whole U/V chain had never been differentiated |
| The emergent Stokes vector is physically admissible | `test_VectorRT_Physics`, clear sky so no cloud lookup table is involved and a failure could not be blamed on coefficient quality. Polarization bound `I^2 >= Q^2+U^2+V^2` holds with margin -7.1e-7, I positive, and the polarized part is genuinely non-zero (1.7e-3) so the bound is not vacuous |
| `n_Stokes = 3` works, and truncation is consistent | Same test. Running one scene at `n_Stokes` 2, 3 and 4 returns bit-identical I and Q, and bit-identical U between 3 and 4. This is the only exercise of the three-component truncation anywhere; the solver guards U with `n_Stokes > 2` and V with `n_Stokes == 4`, so it is a distinct path |
| U and V are the azimuthal signal, not something leaking into those slots | Same test: both vanish to 2.7e-21 at relative azimuth 0 and 180, where every odd harmonic must |
| The `(1,1)` positivity clamp does not fire in practice | Instrumented count over a realistic `n_Stokes = 4` scattering run on the shipped CRTM-Exp table: **0 clamps in 7092 assembled elements**. The clamp is a dormant hazard, not an active defect |
| The surface Jacobian is correct on the vector path, including the observable the capability exists for | `test_VectorRT_TLADK`: dU and dV against wind direction (2.56e-7, agreeing with finite differences to 1.9e-10) and dQ against wind speed, plus wind speed, wind direction, sea surface temperature and salinity in the adjoint dot product and in K against AD. Before this, `Sfc_TL` and `Sfc_AD` were zeroed in every vector test and `Sfc_AD` and `Sfc_K` were never read, so a completely broken surface adjoint satisfied the whole suite. The dot product now reports the surface share of the inner product, 0.3 percent in the overcast scattering blocks and 54 percent clear-sky, so it cannot silently become vacuous |
| PARMIO's polarimetric surface uses the same frame convention as FASTEM, and its Jacobians are correct | `test_VectorRT_SurfaceFrame` runs both backends at the surface interface: amsua_n19 channel 1 (23.8 GHz, FASTEM) and mwr_aws channel 16 (325 GHz, above the PARMIO gate). Both give U and V reaching the solver exactly and the exact even/odd mirror signature, so two independent surface models agree on the convention. `test_VectorRT_PARMIO_TLAD` then covers PARMIO through the full radiative transfer chain on TROPICS channel 12: dU/d(wind direction) against finite differences at 1.96e-10, the adjoint dot product at 1.6e-16 with a 40.8 percent surface share, and K against AD exact. Note the FASTEM channels in that test carry no U at all, because it takes the FASTEM6 default, so the signal is unambiguously PARMIO's |
| Which channels can see the PARMIO surface at all | Only two shipped sensors exceed the 200 GHz gate. mwr_aws is at 325.15 GHz, on a water-vapour line, and is opaque: measured Stokes U there is 1e-16 to 1e-13, against 3.7e-5 on the FASTEM channels. TROPICS is at 204.783 GHz, between the 183 and 325 GHz lines, and is **not** opaque: U/I is 1.6e-3, comparable to the best FASTEM window channels. An earlier revision of this document asserted that both sat on lines and that PARMIO could therefore never influence a top-of-atmosphere radiance. The 325 GHz half was measured; the 204.78 GHz half was inferred and is wrong. The gate is a compile-time parameter and a policy choice, not a data limit: the PARMIO table itself spans 1.4 to 700 GHz |
| The new RTV state is thread safe | `RTV` is a per-thread array (`RTV(n_channel_threads)`), so `e_Rad_UP_Stokes` is per thread by construction |
| **The polarized cloud phase matrix drives a real observable, and it is the only thing that can produce the sub-millimetre signal** | Measured in the jedi bundle on PolSIR, not in this suite, because it needs the ablated CRTM-Exp tables. `CloudCoeff_Exp_Full6_noAllPol.nc` zeroes phase elements 2 to 6 exactly and leaves element 1 untouched (verified directly against the full table), so the ablation is real rather than a no-op. At 40 degrees the 683 GHz polarization difference between the two identically-polarized channels, ch5 (VL_MIXED) and ch6 (HL_MIXED), which share frequency, wavenumber, Planck and band-correction coefficients and differ only in polarization, is: cloud polarization ON, mean -0.010974 K, range -0.218 to +0.121; cloud polarization OFF, mean +0.000004 K, max +0.007263 K. The OFF case reproduces the **scalar** path for the same scene to 0.0009 K, and the scalar path structurally cannot represent cloud polarization at all. So the entire 683 GHz signal is cloud-scattering induced. At nadir the same ablation changes the split by ~1e-5 K, which is correct physics rather than a null result: 180 degree backscatter suppresses the polarization, though the polarized elements still couple into multiply-scattered **intensity** there by up to 0.128 K, almost entirely from the frozen habit (liquid contributes 3e-5 K). An earlier reading of the nadir ablation alone as "the polarized elements do nothing" was wrong, and wrong because nadir is the one geometry that cannot see the effect |
| The polarized infrared surface is correct, and adding it did not move the scalar path | `test_VectorRT_IR_Surface` on `abi_g18`, clear sky and night-time so neither a cloud table nor the solar BRDF is involved. Q exactly zero at nadir where rv equals rh; strictly increasing across 0, 15, 30, 45, 60 degrees, reaching Q/I = 2.24e-2 on the 11.2 um window channel; intensity **bit-identical** between `n_Stokes` 1 and 4 (difference exactly 0.0, which is the construction check on `e_I = (e_V+e_H)/2`); U and V exactly zero; the opaque 6.2 um channel at 1.07e-16 against the window channel's 2.12, so no Q is being injected anywhere but the surface; the polarization bound holding with margin 0.221. Jacobians: TL against central finite differences at 4.3e-8 for wind speed and 5.4e-11 for water temperature, the four-component adjoint dot product at 2.8e-16, K against AD at 9.8e-16. Verified to FAIL against the unfixed code by ablating the forward block and rebuilding: five assertions fail, including Q growing with view angle. Note which two did **not** fail under ablation, the adjoint dot product and K-versus-AD, since a transpose of nothing is still a valid transpose. That is the structural blindness this document opens with, reproduced deliberately |
| Adding the infrared surface left every scalar number untouched | `dump_scalar_fullprec` before and after, 1128 full-precision values, **zero differing lines**. Coverage spans forward, tangent linear, adjoint and K, a microwave sounder, a visible imager and an infrared imager, ocean and land, and clear, overcast and fractional cloud. The instrument was extended to carry the infrared sensor first: without that it would have certified bit-identity while never executing the modified path, which is a certificate of nothing. One finite-difference trap worth keeping: the emissivity table is tabulated on a wind grid of 0, 0.5, ... 10, 12.5, 15 m/s, and a central difference taken exactly at a knot straddles two interpolation cells and disagrees with the analytic tangent linear by 1.8e-5. Mid-cell it is 4.3e-8. That is the differencing, not the derivative |
| The convention is written down and cannot drift silently | `docs/design/polarimetric_conventions.md` states the adopted azimuth and Stokes-sign convention; `test_VectorRT_StokesSign` pins the sign per backend. Flipping U's sign consistently across the forward, tangent-linear and adjoint routines was shown by measurement to leave `test_VectorRT_SurfaceFrame`, `test_VectorRT_StokesOutput` and `test_VectorRT_TLADK` all passing, so before this test nothing in the suite could detect a convention change. FASTEM and PARMIO were also measured to agree in the sign of U, which establishes the two independently fitted coefficient sets share a convention, but not on its own that CRTM's `phi` origin matches it. That last step was closed externally on 2026-08-02 against RTTOV's FASTEM5, which gave matching U and V4 signs at every relative wind azimuth |

## Open questions

- **Normalize_Phase applies inconsistent normalization to below-diagonal
  polarized blocks.** Found 2026-07-31 while bounding the clamp. For each row
  `i` the routine scales that row's intensity and polarized elements by the
  same factor, then performs an intensity-ONLY symmetry copy,
  `Pff(j1,i1) = Pff(i1,j1)`. The `(1,1)` loop covers columns `>= i` while the
  polarized loop covers all columns, so a block below the diagonal ends up with
  its intensity element carrying row `i`'s normalization and its polarized
  elements carrying row `j`'s. Reconciling them needs the polarized symmetry
  relations rather than a guess, so nothing was changed. This is not
  hypothetical: it is why the stress case of `test_PhaseMatrix_Invariants`
  still shows a ratio near 2 after `Bound_Phase_Block` reduced it from 5.6e6.
  It does not produce a bound violation on realistic coefficients, where the
  measured ratio is 0.65.


Unknown, and to be resolved in Phase 1 rather than assumed.

- ~~**Frame convention at the surface.**~~ **RESOLVED 2026-07-31**, negatively:
  the two frames coincide identically at every azimuth and no rotation is
  required. Promoted to "Established facts" above, with the evidence and the
  two tests that measure it. The one thing this did *not* settle was the sign
  convention of the third Stokes component between FASTEM and the solver, which
  could not be decided against CRTM alone; it moved to Phase 0 and was closed
  there on 2026-08-02 against RTTOV.
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

  **Refined 2026-08-03 by measurement, on the visible path where `n_Azi > 0` and
  the Fourier sum is actually exercised.** The adjoint dot-product residual on a
  fixed clear-sky visible scene, ozone control vector, is:

  | `n_Stokes` | residual |
  |---|---|
  | 1 (seeded through `%Radiance`) | 2.41e-13 |
  | 2 | 1.94e-13 |
  | 3 | 9.71e-11 |
  | 4 | 1.87e-10 |

  The degradation appears **exactly when U enters**, and I and Q close to the
  same precision as the scalar path. The structural reading is that U and V ride
  sine harmonics, so their m = 0 term vanishes identically and they are built
  entirely from the higher Fourier components, whose adding-doubling solves are
  more poorly conditioned; I and Q take most of their value from the
  well-conditioned m = 0 solve. Two alternatives were ruled out by measurement
  rather than argument: cancellation in the test's own inner product (the
  factor `sum|t|/|sum t|` is exactly 1.0), and the phase-matrix positivity
  clamp (instrumented on this exact scene, it fires **0 times in 19200**
  assembled blocks, which also extends the earlier CRTM-Exp-only result to
  Rayleigh).

  What is still not established is whether the remaining 1.9e-10 is purely
  conditioning or hides a small genuine non-transpose in the sine-weighted
  accumulation. `test_VectorRT_VIS_Rayleigh` therefore asserts the two tiers
  separately, tight on the cosine components and loose on the sine ones, so a
  real break in I or Q cannot hide behind the looser bound.

  Note also that `Bound_Phase_Block` is called from the **forward** phase-matrix
  assembly only, and has no tangent-linear or adjoint counterpart. That is
  harmless while it never fires, which is now measured on both a microwave
  scattering scene and a visible Rayleigh scene, but it is a latent
  forward-versus-Jacobian inconsistency and should be treated as one if any
  future coefficient set makes the clamp active.

## Phased plan

### Phase 0. External ground truth

Nothing downstream is meaningful without a reference that is not CRTM. Begin
with closed-form cases, which are unambiguous and cheap: an isothermal
non-scattering slab over a specular Fresnel surface has an analytic (I, Q), and
single-scatter Rayleigh has known polarization. Follow with an independent code,
either RTTOV-SCATT's polarimetric path or a PolRadtran or VDISORT reference, for
scattering configurations.

This phase also owned the **third Stokes sign question**, which Phase 1 pinned
but could not validate. **Closed 2026-08-02.** CRTM's adopted convention is
self-consistent and cannot now drift silently, and its relative-azimuth origin
was confirmed to match the one the FASTEM azimuth coefficients were regressed
under, by comparison against RTTOV's FASTEM5 at nonzero wind direction: the U
and V4 signs match at every relative azimuth. No internal instrument could have
determined this, because V and H ride cosine harmonics and are even in the
azimuth, so they are structurally blind to a global sign error in U. It mattered
in practice because such an error stays invisible until it reaches O minus B,
where it doubles the innovation on the one observable a polarimetric instrument
exists to measure.

*Exit:* external confirmation of the third Stokes sign convention, **met**. A
registered test comparing CRTM vector output against a closed-form solution for
at least one clear-sky and one single-scatter configuration, **still open**.

### Phase 1. Convention audit, written down and pinned

Document, then verify against code, the convention at every interface: surface
model output, solver state vector, phase matrix expansion ordering and reference
frame, azimuthal Fourier assignment, and the surface-to-meridional frame
relationship.

**Mostly done, 2026-07-31.** The surface-to-meridional frame relationship is
resolved and pinned by two tests, and the azimuthal Fourier assignment is
resolved as far as the accumulation weights go (gap 2a, including the proof
that the m = 0 phase matrix is block diagonal in {I,Q} and {U,V} because
`Pminus` vanishes there).

The surface convention is now written down and pinned:
`docs/design/polarimetric_conventions.md` is the design note this phase asked
for, stating the relative azimuth definition in terms of CRTM's own
`Wind_Direction` (direction-toward, opposite the meteorological convention) and
`Sensor_Azimuth_Angle` (satellite-to-FOV, clockwise from north), the harmonic
form against its primary source, and `U = T(+45) - T(-45)` per WindSat and
RTTOV. The authoritative statement sits at the single point the angle is formed
(`CRTM_MW_Water_SfcOptics.f90`), with pointers from all three azimuth backends.
`test_VectorRT_StokesSign` is the assertion test, and it was verified to fail
against a deliberately sign-flipped build while the rest of the suite passed.

Note what this does and does not settle. The convention is now explicit,
self-consistent across FASTEM4/5, FASTEM6 and PARMIO, and protected against
silent drift. Whether CRTM's `phi` origin matches the one the FASTEM
coefficients were regressed under was tested against RTTOV (Phase 0) on
2026-08-02, and **the sign conventions match exactly**. This Phase 0 open
question is officially closed.

Phase matrix expansion ordering and the reflectivity structure remain open.

*Exit:* a design note plus assertion tests pinning each interface independently,
so that a later change violating one fails immediately rather than silently.

### Phase 2. Complete the surface

**Largely done, 2026-07-31.** U and V now travel from the surface model to the
solver and out to `RTSolution%Stokes` (gaps 2 and 2a). No frame rotation is
needed, per the Phase 1 resolution. The correction is `1:MAX(2,nL)` at the
microwave *water* aggregation sites only, in the forward, tangent-linear and
adjoint routines; see gap 2 for why the other three surface types must stay at
`1:2` rather than the twelve sites originally prescribed.

Still open in this phase: the reflectivity structure question (whether a
polarimetric surface needs off-diagonal angle coupling), and the fact that the
water reflectivity's third and fourth components are hard-zeroed by FASTEM
itself (`CRTM_FastemX.f90:469`), so U and V are emitted but never reflected.

*Exit:* met for the pass-through. `test_VectorRT_SurfaceFrame` proves the
surface model's U and V reach the solver input; `test_VectorRT_StokesOutput`
proves they reach the reported Stokes vector and that the solver preserves the
frame. Both fail against the unfixed code, the first at exactly zero against a
surface U of 1.2e-2 in emissivity units, the second at exactly zero on all 19
channels.

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

**Started 2026-08-03, on `feature/btj_polarimetric_ir_surface` and deliberately
held out of v3.2.0.** The premise of this phase turned out to be wrong, and it
was wrong because it had been reasoned from the surface code rather than
measured. See gap 6: the visible was already producing real polarized radiance
through Rayleigh scattering, so "state plainly that polarimetric support means
microwave" would have documented something untrue.

What is now done: the infrared water surface is polarized, from Fresnel shape
on lookup-table magnitude (gap 10). The scalar path is unchanged bit for bit,
verified by `dump_scalar_fullprec` over 1128 full-precision values across
forward, tangent linear, adjoint and K, three sensor classes, two surfaces and
three cloud states, with zero differing digits. That instrument had to be
extended to carry an infrared sensor first, because as it stood it would have
certified bit-identity without ever executing the modified path.

What remains:

- The **visible surface** is still unpolarized, so a visible run has correct
  molecular polarization and no sun-glint polarization. That is the single
  largest remaining physical gap, because glint is strongly polarizing and is
  what visible and shortwave polarimetry is largely about.
- **Infrared land, snow and ice** stay at component 1 by choice, not oversight.
- **Polarized infrared cloud scattering** stays refused, correctly, per gap 11.
- The **visible vector path has no registered test at all**. Every
  `test_VectorRT_*` runs a microwave sensor, and the infrared one added here is
  the first outside it. The visible numbers above are the first time anyone has
  looked.

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

**Make the backend selection honest. Half done, 2026-07-31.** Two gaps
compounded into a silent-zeros trap for exactly the users this capability
targets: `MWwaterCoeff_File` selected nothing (gap 2c), while the FASTEM6
default has no third or fourth Stokes azimuth model and returns both as
identically zero (gap 2b). Together, a user could select what looks like a
polarimetric backend, get a successful run, and receive U = 0 that is
indistinguishable from a scene with no polarimetric signal.

**Both halves are now done.** Gap 2c is fixed and the filename argument is
honoured, so the JEDI path works as written. The second half is closed too:
`CRTM_Forward` now warns when `n_Stokes > 1` is requested while the loaded
microwave water backend has no third or fourth Stokes azimuth model, naming
FASTEM4 and PARMIO as the alternatives. It warns rather than fails, because
the configuration is legitimate for the intensity and refusing it would break
callers who set `n_Stokes` globally but only want I.

Two details worth keeping. The check sits before the profile loop, so it is
evaluated once per call and outside the parallel region, and it is latched to
once per loaded scheme on top of that: unlatched it produced 168 copies in
`test_VectorRT_TLADK` alone, which buries the message rather than delivering
it. And the conditions are nested rather than combined with `.AND.`, because
Fortran does not guarantee short-circuit evaluation and querying the latch
consumes it.

The discriminator is `CRTM_MWwaterCoeff_HasPolarimetric`, tied to the measured
surface by `test_MWwaterCoeff_FileSelects` rather than left to agree only with
itself: it must be true exactly when the surface actually produces a
polarimetric signal.

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
