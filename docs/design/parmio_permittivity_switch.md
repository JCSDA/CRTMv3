# The PARMIO permittivity switch

Status: decided on evidence, not yet implemented. The switch should be removed
and Meissner used across the whole table, because that is PARMIO's own reference
configuration and the one behind SURFEM-Ocean in RTTOV. Sections 1 to 6 are the
measurement, section 7a is the citation that settles it, section 7b records an
argument that pointed the wrong way and why.

Implementing it means regenerating the LUT, reverting the installed table first,
and changing emissivity at and above 200 GHz by up to 0.056 in e_v. Nothing here
does any of that.

Measured 2026-08-01 against PARMIO at `/home/ben/CRTM/parmio` (NWP SAF clone,
master + 2 local commits). Suite green at 239/239 before and after; nothing in
this note changes code or coefficients.

Model provenance throughout is taken from the PARMIO reference paper rather than
inferred from the code: E. Dinnat and coauthors, "PARMIO: A Reference Quality
Model for Ocean Surface Emissivity and Backscatter from the Microwave to the
Infrared", Bull. Amer. Meteor. Soc., 104 (4), E742-E748, 2023,
doi:10.1175/BAMS-D-23-0023.1.

## 1. The switch is ours, not PARMIO's

The LUT is built in three groups and the group selection switches the sea water
dielectric model at `PERMITTIVITY_SWITCH_GHZ`:

    sss_dependent    f <= 10.65 GHz     Meissner, SSS axis active
    sss_nominal_m    10.65 < f < 200    Meissner, SSS = 35 only
    sss_nominal_h    f >= 200 GHz       high-frequency tabulated dielectric

That 200 GHz is a local choice. Upstream PARMIO switched at 28.837 GHz, which is
simply the bottom edge of the tabulated model's range: use the table wherever it
is defined, and Meissner only below it. The value was moved to 200 GHz in parmio
commit `5c4579c` (2026-05-07), whose message states the reason, that the Meissner
branch should cover the ATMS and AMSU window channels. `parmio_pool.py` still
carries the original rule as `--perm-mode production`.

VERIFIED: the shipped table was built with the 200 GHz rule, not the 28.8 GHz
one. The installed LUT's groups meet exactly at 199.90 and 200.00 GHz, and its
`permittivity_policy` attribute reads "Meissner below 200.0 GHz; high-frequency
tabulated dielectric at and above 200.0 GHz".

So the description of the switch in `PARMIOCoeff_Define` is correct as a
description. What was wrong was calling the resulting step "PARMIO's
two-permittivity construction" and "a question for the model rather than for
CRTM". PARMIO is continuous in frequency within either model. The step is
manufactured by our own table, and it is ours to place.

The move was also in the right direction against the model's own documentation.
The ISSI team's stated default is Meissner "in the microwaves" (section 2), and
switching away from it at 28.8 GHz does not implement that, 28.8 GHz being
nowhere near the top of the microwave by any convention. Moving the switch to
200 GHz brings the table closer to the documented default, not further from it.

## 2. What the two models are

`m`, `epsilon_MW`: Meissner and Wentz 2004, with the 2012 and 2014 corrections.
A double-Debye form fitted to satellite and laboratory data. Dinnat et al. 2023
describes it as "an empirical parameterization adjusted and validated using
remote sensing observations from 1.4 to 89 GHz". At 200 GHz it is an
extrapolation of a bit more than a factor of two.

**This is PARMIO's documented default in the microwave.** Dinnat et al. 2023,
page E744: "The team selected the Meissner and Wentz parameterization as the
default configuration for the reference model in the microwaves."

`h`, `epsilon_hifreq`: a tabulated complex refractive index, converted by
eps = (n + ik)^2. Three properties matter and none of them are obvious from the
name:

- Its lowest node is 0.962 cm^-1, which is 28.837 GHz. That is where the table
  stops, not where the physics changes. Dinnat et al. 2023 describes it as "a
  high-frequency model developed for the project at frequencies from 28.8 GHz
  and up to 449,677 GHz", so the floor is a deliberate model bound and not an
  accident of tabulation.
- Below 800 GHz the nodes are spaced 28.9 GHz apart. The two bracketing 200 GHz
  are 173.372 and 202.279. A query at 200 GHz is a straight line drawn between
  two points 29 GHz apart.
- Temperature enters as a two-point linear fit, n(T) = n_T0 + T*coef_n, built
  from data at 273 K and 298 K. Salinity enters as a correction that the table
  header describes as reliable only over 500 to 5000 cm^-1, and which is set to
  zero everywhere else. Through the entire microwave the table is pure water.

The underlying measurements (Rowe et al. 2020, Pinkley and Williams 1976,
Newman et al. 2005) are far infrared and optical work. The microwave end is the
tail of that table, not its subject.

VERIFIED by independent reimplementation: reading `refindex_hifreq.dat` in
Python and repeating the interpolation reproduces PARMIO's own reported
permittivity at 200 GHz, 5.5442 + 7.6820i against the 5.54 + 7.68i in its
output.

## 3. The models cross, but not at a fixed frequency

Sweep: 300 frequencies from 29 to 700 GHz, both models, zenith 45 degrees,
U10 10 m/s, SSS 35, foam off, at three sea surface temperatures. 600 PARMIO
jobs, all clean. Driver `parmio_pool.py --perm-mode both`, analysis
`parmio/scripts/perm_crossover_analysis.py`, output under
`parmio/Outputs/sweep/perm_crossover/`.

Each pair crosses exactly once. The crossing moves with temperature:

| SST | e_v crossing | e_h crossing |
|---|---|---|
| 0 C | 29.9 GHz | 30.0 GHz |
| 15 C | 43.0 GHz | 43.0 GHz |
| 30 C | 155.1 GHz | 153.9 GHz |

That is a 125 GHz spread over the ocean temperature range. **There is no fixed
frequency at which the two models agree.** Placing the switch where they cross
is available only for one temperature at a time.

Below the crossing the tabulated model gives the higher emissivity, above it the
lower. Away from the crossing the disagreement is large and it does not close
again anywhere in the band:

    SST = 15 C, h minus m in e_v
      30 GHz   +0.028
      43 GHz    0.000   (crossing)
      90 GHz   -0.029
      200 GHz  -0.032
      325 GHz  -0.022
      683 GHz  -0.020

## 4. 200 GHz is close to the worst available choice

At 15 C the disagreement at 200 GHz is 99.5 percent of its maximum over the
whole 29 to 700 GHz band. The switch sits almost exactly on the peak.

Step in e_v at each candidate switch frequency, by SST, with the rough
brightness temperature equivalent:

| switch | 0 C | 15 C | 30 C | worst |
|---|---|---|---|---|
| 28.84 GHz (upstream) | +0.003 (1.0 K) | +0.031 (8.9 K) | +0.037 (10.8 K) | 10.8 K |
| 43.0 GHz (15 C crossing) | -0.039 (11.1 K) | +0.001 (0.3 K) | +0.022 (6.3 K) | 11.1 K |
| 155 GHz (30 C crossing) | -0.059 (16.9 K) | -0.028 (7.9 K) | 0.000 (0.0 K) | 16.9 K |
| 200 GHz (current) | -0.056 (16.1 K) | -0.032 (9.1 K) | -0.009 (2.7 K) | 16.1 K |

No fixed switch does better than about 11 K worst case over SST. Moving the
switch trades which temperatures are penalised; it does not remove the step.
Blending across a band converts the step into a ramp of the same size and adds a
second arbitrary parameter, the bandwidth.

## 5. What this does not turn on

Salinity is not a discriminator, contrary to how the tabulated model's missing
salinity correction first appears. Over the realistic open ocean range, Meissner
gives (SST 15 C, zenith 45, U10 10):

    200 GHz   30 to 37 psu spans 0.00044 in e_v   (0.13 K)
    325 GHz   30 to 37 psu spans 0.00149 in e_v   (0.43 K)
    683 GHz   30 to 37 psu spans 0.00082 in e_v   (0.23 K)

So the table having no salinity above 10.65 GHz costs a few tenths of a kelvin,
and the LUT's own decision to drop the SSS axis above 10.65 GHz is defensible on
the same numbers. Comparing 0 psu against 35 psu makes the effect look ten times
larger and is not a relevant comparison for ocean.

## 6. What it does turn on

Below 200 GHz the evidence is one-sided. At 0 C the two models differ by 0.078
in e_v at 77.5 GHz, about 21 K, in the middle of the operational window
channels. This is the region where Meissner is fitted and validated and where
the tabulated model is a two-point linear temperature extrapolation of far
infrared data sampled every 29 GHz. Meissner should be used there. The move from
28.8 to 200 GHz in `5c4579c` was the right direction, whatever the rationale
recorded at the time.

Above 200 GHz the question is settled too, and it is settled by the reference
literature rather than by anything we can measure here. Kilic et al. 2023
answers it directly: the PARMIO configuration adopted for the reference model
and used to build SURFEM-Ocean uses Meissner and Wentz across the whole range,
500 MHz to 700 GHz. Section 7a has the wording.

So both sides of our switch should be Meissner, and the tabulated model has no
role anywhere in the microwave. It is the infrared dielectric. The step is not
an expression of a real scientific uncertainty; it is a configuration error that
applies an infrared model to sub-millimetre channels.

## 7. Options, with the measurement attached to each

Recorded as they stood before Kilic et al. 2023 was read, because the reasoning
for discarding three of them is still the reasoning, and because option 3 was
withdrawn on this list for a bad reason that is worth not repeating.

1. **Move the switch.** Ruled out as a fix by section 3. There is no frequency
   where the models agree across temperature, and every alternative is worse
   than 200 GHz for some part of the ocean.

2. **Keep the switch at 200 GHz and record the step as a known limitation.**
   Changes nothing, costs nothing, and leaves a manufactured discontinuity at
   the frequency where the two models disagree most at mid-latitude SST. Now
   also known to leave an infrared dielectric applied to sub-millimetre
   channels, so this is no longer a neutral do-nothing.

3. **One model throughout, and drop the third group.** Meissner from 1.4 to
   700 GHz is continuous by construction and has no switch to place. It does
   not leave the table alone: it changes the emissivity at every frequency at
   and above 200 GHz, which is exactly the range served by default, while
   changing nothing below, which is the range that is currently unreachable
   without opting in. See section 7b for the size of that change. **This is the
   answer.** It is what PARMIO's own reference configuration does.

4. **Settle it externally.** Done, by reading rather than by measuring. See
   section 7a.

5. **Blend across a transition band.** Moot. There are not two locally valid
   models to blend between.

## 7a. Settled by SURFEM-Ocean: Meissner throughout

Dinnat et al. 2023 records that PARMIO was used as the reference to train
SURFEM-Ocean, a neural network fast emissivity model that "extends the frequency
coverage of the previous fast ocean surface emissivity model for microwave
frequencies FASTEM to 0.5-700 GHz", and that it has shipped in RTTOV since
version 13.2 in December 2022 and targets ECMWF Cycle 49r1.

That range is exactly ours, so the configuration PARMIO was run in to generate
SURFEM-Ocean's training set is the community's de facto answer, and it is
already operational. Kilic et al. 2023 gives it plainly.

Section 2.1, Dielectric Constants: "The dielectric constants from Meissner and
Wentz (2012) used in Remote Sensing Systems (RSS) ocean emissivity model have
been chosen for PARMIO. [...] We perform a comparison of the flat ocean
emissivity to evaluate the extrapolation of the dielectric constant model for
the low frequencies down to 500 MHz and for the high frequencies up to 700 GHz."

Section 3: "In the following, PARMIO is used with the configuration described
above, that is, with the dielectric constants from Meissner and Wentz (2004,
2012), the wave spectrum from Durden and Vesecky (1985) with the amplitude
coefficient multiplied by 1.25, and the new foam coverage [...]".

And for the fast model itself: "It is estimated in SURFEM with the dielectric
constant module from Meissner and Wentz (2004, 2012) and the Fresnel equations."

So Meissner is used across the entire 500 MHz to 700 GHz range, the
extrapolation to 700 GHz was looked at deliberately and accepted, and its
Figure 1 compares Klein and Swift, Ellison, and Meissner over that range. The
high-frequency tabulated dielectric appears nowhere.

VERIFIED: the strings "Rowe", "Pinkley", "hifreq" and "tabulated dielectric"
do not occur anywhere in the 22 pages of Kilic et al. 2023. The tabulated model
is not part of the microwave reference configuration at all.

That answers section 6 without a measurement campaign. It also means the useful
follow-on is no longer "which model is right" but "does our LUT reproduce
SURFEM-Ocean", which is a direct comparison against RTTOV 13.2 over 200 to
700 GHz and checks our interpolation at the same time.

VERIFIED: PARMIO's own validation against observations stops at 165.5 GHz. The
sensors listed in Dinnat et al. 2023 are SMAP at 1.4 GHz, AMSR2 from 6.9 to
89 GHz, GMI from 10.6 to 166 GHz, and ATMS between 23.8 and 165.5 GHz. Nothing
above 166 GHz was compared to observations. So the entire default PARMIO
dispatch range in CRTM, which is 200 GHz and above, lies outside anything PARMIO
itself was validated against. That is a sharper statement than the
`confidence_label` in `parmio_lut_grid.py`, which calls 24 to 225 GHz
"extrapolated-defensible" and only above 225 GHz "extrapolated-experimental".

## 7b. What option 3 changes, and the argument that briefly blocked it

Replacing the tabulated group with Meissner raises the emissivity everywhere it
applies. m minus h in e_v, with the rough surface Tb equivalent:

| freq | SST 0 C | SST 15 C | SST 30 C |
|---|---|---|---|
| 200.00 GHz | +0.056 (15.3 K) | +0.032 (9.1 K) | +0.009 (2.9 K) |
| 204.78 GHz | +0.056 (15.2 K) | +0.032 (9.1 K) | +0.010 (3.0 K) |
| 229.00 GHz | +0.051 (14.0 K) | +0.030 (8.5 K) | +0.009 (2.8 K) |
| 325.15 GHz | +0.031 (8.6 K) | +0.022 (6.3 K) | +0.013 (3.9 K) |
| 448.00 GHz | +0.025 (6.9 K) | +0.021 (6.2 K) | +0.017 (5.2 K) |
| 683.00 GHz | +0.025 (6.9 K) | +0.020 (5.7 K) | +0.014 (4.2 K) |

e_h moves further, to +0.070 at 200 GHz and 0 C.

### The FASTEM argument, and why it was wrong

Option 3 was briefly withdrawn on the following comparison. At 325 GHz, PARMIO
minus FASTEM in e_v, from `parmio_fastem_vh_sweep.csv` with the m minus h shift
applied at matched state:

| SST | as shipped (h) | under one model (m) |
|---|---|---|
| 0 C | +0.0015 | +0.0329 |
| 15 C | +0.0124 | +0.0343 |

Read as "FASTEM is closer to the tabulated model, so the tabulated model is
better above 200 GHz". That reading is wrong, for three reasons, all of which
were available before the comparison was made.

- FASTEM's own dielectric is Ellison et al. 1998, and Kilic et al. 2023 reports
  that Meissner and Ellison differ by only 0.009 in flat-ocean emissivity at
  200 GHz. So the dielectric cannot account for a 0.034 gap. That gap is
  roughness and foam, which PARMIO and FASTEM treat differently. The caveat was
  written down at the time and then not applied.
- Kilic et al. 2023 states that FASTEM "produces unrealistic emissivity
  calculations at frequencies above 200 GHz (below 0 or higher than 1)".
  A model documented as unphysical in a band cannot referee that band.
- FASTEM is the model PARMIO and SURFEM-Ocean exist to replace above 200 GHz.
  Agreement with it is not evidence of correctness there.

The general lesson is the one in the working rules: agreement with an existing
implementation is an inference, not verification. The primary source settled in
one reading what the proxy had pointed the wrong way on.

Recommendation: option 3. Meissner throughout, and delete the third group. This
matches PARMIO's own reference configuration, matches SURFEM-Ocean and therefore
RTTOV, removes the switch rather than relocating it, and stops applying an
infrared dielectric to sub-millimetre channels.

## 7c. Full Stokes is unaffected, and is the point of the exercise

Moving to Meissner throughout does not cost polarimetric capability. The
azimuthal harmonics are non-zero across the whole band in both configurations,
and the dielectric choice barely touches them. At SST 15 C, zenith 45 degrees,
U10 10 m/s, in kelvin:

| freq | U1 (m) | U1 (h) | U2 (m) | U2 (h) | V2 (m) | V2 (h) |
|---|---|---|---|---|---|---|
| 90 GHz | -0.721 | -0.693 | -1.710 | -1.836 | +0.121 | +0.147 |
| 200 GHz | -0.536 | -0.539 | -1.034 | -1.169 | +0.049 | +0.056 |
| 325 GHz | -0.372 | -0.385 | -0.804 | -0.882 | +0.038 | +0.040 |
| 683 GHz | -0.181 | -0.195 | -0.696 | -0.745 | +0.028 | +0.026 |

At the switch the dielectric change moves U1 by 0.6 percent and U2 by about
12 percent, against 4 percent in e_v. The step discussed in this note is
concentrated in the V and H pair, not in the third and fourth Stokes terms, so
the polarimetric signal is the part of the table least disturbed by the
decision.

SURFEM-Ocean is full Stokes on the same harmonic structure we use:

    e_p = e_p0 + e_p1 cos(phi) + e_p2 cos(2 phi)      V and H
    e_q =       e_q1 sin(phi) + e_q2 sin(2 phi)       S3 and S4

which is the 14-term layout of our own groups.

Worth recording because it validates the premise of this branch: Kilic et al.
2023 says FASTEM "was developed for use only in the frequency range of
1-200 GHz and the viewing angle range of 0-60 degrees and without a full
treatment of polarization (full Stokes vector)", and that FASTEM-6 "is the only
version recommended for most sensors, apart from those with polarimetric
channels and those beyond 200 GHz". Those two exclusions, polarimetric channels
and above 200 GHz, are exactly the two reasons CRTM reaches for PARMIO. The
independent conclusion matches ours, including our finding that FASTEM6 returns
U = 0 by default.

## 8. Scope note

By default `PARMIO_Is_Active_At` serves only f >= 200 GHz, which is exactly the
`sss_nominal_h` group. The `sss_nominal_m` group is reachable only when a caller
sets `Options%Use_PARMIO_MWSSEM`. So in the default configuration the step
described here is not reachable inside PARMIO at all: what a default user meets
at 200 GHz is the FASTEM to PARMIO handover, a separate discontinuity measured
at -1.42 K mean for TROPICS channel 12. The permittivity step becomes reachable
the moment a caller opts in, and the 21 K disagreement at 77.5 GHz and 0 C is in
the opted-in range.
