## Status update — MW-land surface Jacobians complete below 80 GHz (for v3.2.0)

The microwave **land** surface portion of this issue is now fully implemented and
validated. The NESDIS_LandEM path (< 80 GHz) returns analytic TL/AD/K
sensitivities for **every physical land-state variable**:

| Variable | Status |
|---|---|
| `LAI` | ✅ analytic (Phase 1) |
| `Vegetation_Fraction` | ✅ analytic (Phase 1) |
| `Soil_Moisture_Content` | ✅ analytic (Phase 2) |
| `Soil_Temperature` | ✅ analytic (Phase 3) |
| `Land_Temperature` | ✅ analytic emissivity part + existing skin-T emission (Phase 3) |
| `Canopy_Water_Content` | ✅ exactly zero — never consumed by the forward (asserted in the test) |

**How the temperatures enter (Phase 3).** Temperature reaches the LandEM
emissivity by two paths, both differentiated exactly: the soil dielectric
(`Soil_Diel`, temperature-dependent water permittivity) → roughened Fresnel `r23`
→ two-stream, and the canopy/soil thermal ratio `gsect0` inside the two-stream
(the dominant term). The soil dielectric → `r23` chain is shared with the
soil-moisture Jacobian. When `Soil_Temperature` is out of the LandEM valid range
it is aliased to the skin temperature in the forward; the analytic path mirrors
that — the aliased input gets a zero derivative and its sensitivity is
re-attributed to `Land_Temperature`.

**Latent bug fixed.** The forward emissivity already depended on skin temperature
through `gsect0`, but that term was dropped from `Surface_K%Land_Temperature`. So
the below-cutoff land-temperature Jacobian carried only the Planck emission term
and was ~3–4× too large (the emissivity part is sizable and *negative*, partially
cancelling the emission term). It now matches finite differences.

**Validation.** `test_Unit_Land_Jacobian` (source `test_Land_Jacobian.f90`)
checks AD ≡ TL ≡ central finite-difference of the forward for all five nonzero
variables, plus a guarded `Canopy_Water_Content == 0`. Forward radiances are
bit-identical; a field-level diff confirms only the `Land_Temperature` and
`Soil_Temperature` Surface-Jacobian columns move (`Atmosphere_K` and
`RTSolution_K` unchanged). Full regression suite passes (215/215).

### Disposition of the remaining listed parameters

**Structurally zero (no forward dependence — a valid Jacobian is zero).** These
cannot be made nonzero without adding forward physics; proposing to close them
out of this issue (or split to a forward-enhancement issue):
`Canopy_Water_Content`, `Snow_Density`, `Snow_Grain_Size` (MW), `Ice_Density`,
`Ice_Roughness`.

**IR/VIS land.** Out of scope by construction: IR/VIS land emissivity is a
lookup table indexed by `Land_Type` (a classifier), with no continuous
physical-state dependence, so the physical-state Jacobians are identically zero
there (the modules' TL/AD set emissivity sensitivity to zero). The skin-temperature
Jacobian is frequency-independent and already provided by
`CRTM_Compute_SurfaceT_AD`.

**Snow / ice (`Snow_Depth`, `Snow_Temperature`, `Ice_Thickness`,
`Ice_Temperature`).** Still open. These use separate models (NESDIS SnowEM /
SeaIceEM) whose TL/AD are zero-stubs that don't take `Surface_TL/_AD`; a real
implementation needs signature changes and per-model differentiation. Tracked as
follow-on (Phase 4), not part of the v3.2.0 land completion.

Design notes: `docs/design/surface_jacobians_281.md`.
