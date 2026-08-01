# CRTM polarimetric conventions

Status: adopted 2026-07-31. Applies to vector radiative transfer, `Options%n_Stokes > 1`.

Companion to `polarimetric_support_roadmap.md`, which tracks the capability as
a whole. This document is the Phase 1 design note; the sign question it leaves
open belongs to Phase 0 there.

This document fixes the sign and angle conventions for the polarimetric
(third and fourth Stokes) components in CRTM, states where each one is
implemented, and records what is verified against what is still open.

Nothing in CRTM consumed U or V before the vector radiative transfer work,
so adopting an explicit convention now costs nothing and prevents the
components from drifting silently later.

## 1. The Stokes vector

CRTM carries `(I, Q, U, V)` in the standard Stokes basis, not modified
Stokes. The surface models produce `(e_V, e_H, e_U, e_V4)` and the
conversion to the solver basis is

    e_I = (e_V + e_H) / 2
    e_Q = (e_V - e_H) / 2

with U and V4 passing through unchanged. This is applied in
`CRTM_SfcOptics.f90` on the coupled-polarization branch.

The third Stokes component follows the standard radiometric definition

    U = T(+45 degrees) - T(-45 degrees)

This is the definition used by WindSat, whose measurements the FASTEM
azimuth coefficients were fitted to, and by RTTOV. It is not a free choice
for CRTM: we evaluate coefficients that were regressed under it, so
adopting the opposite sign would put us at odds with both the coefficients
and with any instrument reporting U.

## 2. The relative azimuth angle

There is exactly one place where the relative azimuth is formed, and all
three microwave water backends read it from there
(`CRTM_MW_Water_SfcOptics.f90`):

    phi = Surface%Wind_Direction - Sensor_Azimuth_Angle        [degrees]

Each backend then applies `phi_radians = phi * DEGREES_TO_RADIANS` with no
further reflection, offset, or sign change.

The two inputs are defined by CRTM as:

| quantity | definition | source |
|---|---|---|
| `Wind_Direction` | direction the wind blows **toward**, clockwise from North. Zero is a wind blowing toward the north, i.e. a southerly. Deliberately opposite to the meteorological convention. | `CRTM_Surface_Define.f90`, `DEFAULT_WIND_DIRECTION` |
| `Sensor_Azimuth_Angle` | azimuth of the horizontal projection of the line **from the satellite to the FOV**, clockwise from North, 0 to 360 | `CRTM_Geometry_Define.f90:312-316` |

So `phi = 0` means the wind blows toward the same compass azimuth as the
satellite-to-FOV horizontal projection.

Stating this plainly matters because CRTM's own history contains a 180
degree change. The legacy Fastem3 path carries both forms, one commented
out (`CRTM_Fastem3.f90:599-600`):

```fortran
! version 8_5      phi = (wind10_direction-Sat_Azimuth_Angle)*pi/180.0_fp
      phi = PI - (wind10_direction-Sat_Azimuth_Angle)*pi/180.0_fp    ! version 8_7
```

Under `phi -> PI - phi` the identities `cos(m phi) -> (-1)^m cos(m phi)` and
`sin(m phi) -> -(-1)^m sin(m phi)` mean that reflection flips the **odd**
harmonics of V and H and the **even** harmonics of U and V4. Different
components, different harmonics, no error message. The modern path
(FastemX, FASTEM6, PARMIO) uses the un-reflected form above.

## 3. The harmonic expansion

From the primary source, Liu et al., *FASTEM-4 Validation*,
NWPSAF-MO-VS-045 (2011), equations 2a-2d:

    E_v = ... + SUM_m c_m cos(m phi_R)
    E_h = ... + SUM_m d_m cos(m phi_R)
    E_3 =       SUM_m e_m sin(m phi_R)
    E_4 =       SUM_m g_m sin(m phi_R)

Cosine for V and H, sine for the third and fourth Stokes components.

Note that the report never defines the origin or sense of `phi_R` itself.
It says only "a relative azimuth angle". That omission is the reason
section 5 below exists.

## 4. Implementation status per backend

| backend | file | phi | V, H | U (3rd) | V4 (4th) |
|---|---|---|---|---|---|
| FASTEM4 / FASTEM5 | `Azimuth_Emissivity_Module.f90:139-142` | `Az * D2R` | `cos(m phi)` | `sin(m phi)` | `sin(m phi)` |
| FASTEM6 (**CRTM default**) | `Azimuth_Emissivity_F6_Module.f90:144,187-188` | `Az * D2R` | `cos phi, cos 2phi` | **identically zero** | **identically zero** |
| PARMIO | `PARMIO_Azimuth_Module.f90:72,88-90` | `Az * D2R` | `cos phi, cos 2phi` | `sin phi, sin 2phi` | `sin phi, sin 2phi` |

All three agree on the angle convention and on the cosine/sine parity.

**FASTEM6 is the default and has no third or fourth Stokes model.** A
polarimetric run over water therefore has a real surface U and V4 only on
FASTEM4 or PARMIO. Note that `CRTM_MWwaterCoeff_Load_FASTEM` accepts only
`FASTEM4` and `FASTEM6`; any other scheme string, FASTEM5 included, is a hard
error rather than a fallback.

Either `MWwaterCoeff_Scheme` or `MWwaterCoeff_File` selects the model. The
filename form used to select nothing at all and was fixed on 2026-07-31: the
shipped names are `<SCHEME>.MWwater.EmisCoeff.<ext>` and the scheme is now
recovered from the leading component. If both arguments are given and
disagree, the scheme wins and the disagreement is reported. See roadmap gap
2c, pinned by `test_MWwaterCoeff_FileSelects`.

The default itself used to be a silent trap: requesting `n_Stokes > 1` while
FASTEM6 is loaded succeeds and returns U = 0, indistinguishable from a scene
with no polarimetric signal. `CRTM_Forward` now warns on that combination,
naming FASTEM4 and PARMIO as the alternatives. It warns rather than fails,
since the configuration is legitimate for the intensity. The warning is
latched to once per loaded scheme, so a finite-difference driver calling the
forward model hundreds of times reports it once.

`CRTM_MWwaterCoeff_HasPolarimetric` is the public query behind it, and
answers whether the loaded scheme carries an azimuth model for the third and
fourth Stokes components at all.

## 5. PARMIO dispatch and its coverage

PARMIO serves a microwave water channel when three things hold, and
`PARMIO_Is_Active_At` is the single predicate that answers it. Everything
that needs to know asks that rather than re-deriving the rule, so the call
sites cannot drift apart:

1. the LUT is loaded;
2. the frequency is at or above the default dispatch floor of 200 GHz, **or**
   the caller set `Options%Use_PARMIO_MWSSEM`. The floor is a safety gate
   rather than physics: it sits where the traditional sounding sensors stop,
   so that enabling PARMIO could not disturb operational channels while the
   implementation was still being shaken out. Nothing at or above 200 GHz was
   exercised operationally, so nothing could regress. Obs-space validation
   against ATMS-NPP supports the same placement from the other side: FASTEM6
   is competitive through the whole ATMS band (top 183.31 GHz) and PARMIO's
   advantage appears where FASTEM6 extrapolates;
3. **the table actually holds data at that frequency.**

`Use_PARMIO_MWSSEM` is a logical in `CRTM_Options_type`, sitting beside
`Use_Old_MWSSEM`, which is the existing switch for selecting the microwave
water surface model. That placement is deliberate. `CRTM_Init` is about what
to load; which model runs is a runtime choice, and `Use_Old_MWSSEM` is the
precedent. An earlier revision put a real-valued frequency threshold on
`CRTM_Init` instead, which matched nothing there: of that routine's optional
arguments, 30 are CHARACTER (paths, filenames, formats, scheme names), 4 are
LOGICAL feature toggles and 2 are MPI process IDs. A boolean also states the
actual intent better than a magic number, since the floor is a safety gate
rather than a physical boundary.

Following the precedent set by `Compute_Up_Radiance_Profile` and the other
newer Options components, `Use_PARMIO_MWSSEM` is deliberately excluded from
the Options binary record and takes its type default on read, so the on-disk
format is unchanged.

The third condition is a hard requirement and is not relaxed by opting in.
The coefficient groups are gridded separately either side of the permittivity
switch and their grids do not meet it: `sss_nominal_m` stops at 183.31 GHz and
`sss_nominal_h` starts at 229 GHz, so 183.31 to 229 selects a group with
nothing in it. The interpolator clamps an out-of-range query to the nearest
grid node and returns a confident number computed somewhere else, so before
this check a 204.78 GHz channel was being evaluated at **229 GHz**, about
24 GHz away, with nothing in the result to say so.

With the shipped LUT that makes the effective floor 229 GHz rather than 200.
It becomes 200 by itself if the 183 to 229 band is ever filled in, with no
constant to update. The gap is a grid-spacing choice in the offline generator
(`parmio/scripts/parmio_lut_grid.py`, `PRODUCTION_FREQS` jumps 183.31 to 229),
not a limitation of the model, so it is fillable by regenerating the table.

Clamping on the remaining axes is still possible and is now reported once per
run: the table spans zenith 0 to 65 degrees, wind 1 to 25 m/s and SST -2 to
30 C, all of which real scenes exceed.

One consequence worth recording. `test_VectorRT_PARMIO_TLAD` used TROPICS
channel 12 at 204.78 GHz, which is exactly the band with no data, so it was
demonstrating PARMIO Jacobians on clamped coefficients. It now lowers the
floor to reach 91.319 GHz, well inside `sss_nominal_m`, where the signal is
real: |U/I| is 8.7e-3 against 1.6e-3 before.

## 5. Verified, and what is not

Verified by measurement:

- U and V4 reach the vector solver from the surface model unchanged, on
  both the FASTEM and PARMIO backends (`test_VectorRT_SurfaceFrame`).
- I and Q are even under `phi -> -phi`; U and V4 are odd
  (`test_VectorRT_SurfaceFrame`, with a non-degeneracy floor so the
  assertions cannot pass on zeros).
- U and V4 survive the azimuthal Fourier accumulation into
  `RTSolution%Stokes` (`test_VectorRT_StokesOutput`).
- The surface `(V,H)` to Stokes `(I,Q)` conversion
  (`test_VectorRT_SurfaceBasis`).
- The adopted sign itself, per backend, at `phi = +90` over ocean at 45
  degrees zenith and 12 m/s wind (`test_VectorRT_StokesSign`):

  | backend | U(+90) | V4(+90) | U(0) | U(180) |
  |---|---|---|---|---|
  | FASTEM4, amsua_n19 ch1, 23.8 GHz | -2.162863e-03 | -9.539306e-05 | 0.0 | -2.3e-18 |
  | PARMIO, mwr_aws ch16, 325 GHz | -2.091229e-03 | -3.286526e-05 | 0.0 | -6.2e-19 |

  FASTEM and PARMIO agree in sign. Since the two coefficient sets were
  fitted independently, that agreement establishes they were regressed
  under a common convention. It does not establish that CRTM's `phi`
  origin matches that convention, because an error there flips both
  together.

The blindness claimed above was measured, not argued. Flipping the sign of
the third Stokes component consistently across the forward, tangent-linear
and adjoint routines of `Azimuth_Emissivity_Module` produces:

| test | result under a consistent sign flip |
|---|---|
| `test_VectorRT_StokesSign` | **fails** |
| `test_VectorRT_SurfaceFrame` | passes |
| `test_VectorRT_StokesOutput` | passes |
| `test_VectorRT_TLADK` | passes |

Flipping only the forward routine additionally breaks `test_VectorRT_TLADK`,
but that is the tangent-linear disagreeing with the forward, not the suite
detecting the convention change.

**Open.** Whether the `phi` origin defined in section 2 matches the origin
the FASTEM azimuth coefficients were fitted under is *not* established.
The FASTEM-4 report does not define `phi_R`, and no accessible RTTOV or
NWP SAF document states it either. Everything CRTM can check internally is
even in `phi` or checks parity only, and is structurally blind to a global
sign error in U.

Closing it requires an external reference. Cheapest first:

1. One RTTOV run at a nonzero wind direction, same geometry and
   coefficients, comparing the sign of the third Stokes emissivity. RTTOV
   consumes the same model from the same coefficient lineage.
2. Failing that, WindSat published upwind/downwind harmonic amplitudes;
   comparing against real data pins it directly.

Until then `test_VectorRT_StokesSign` pins the convention *as adopted here*
so that any change is deliberate and reviewable. It is not evidence that
the sign is correct against nature, and the test says so.

## References

- Liu, Q., et al. (2011). *FASTEM-4 Validation.* NWP SAF, NWPSAF-MO-VS-045.
- Gaiser, P. W., et al. (2004). The WindSat spaceborne polarimetric
  microwave radiometer. *IEEE TGRS* 42(11).
- Saunders, R., et al. (2018). An update on the RTTOV fast radiative
  transfer model (currently at version 12). *GMD* 11, 2717-2737.
- Kilic, L., et al. (2023). Development of SURFEM-Ocean based on the
  PARMIO radiative transfer model. *Earth and Space Science.*
