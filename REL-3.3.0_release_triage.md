# CRTM release triage: the 3.2.x patch gate and the 3.3.x punch-list

**Adopted 2026-08-13 (BTJ) on `feature/btj_REL-3.3.0`.** This file states the
rule for deciding whether work lands in a 3.2.x patch release or queues for
3.3.x, and triages the known backlog against it. It formalizes the gate that
`REL-3.2.0_changes_vs_develop.md` applied informally ("none of the additions
changes any common-suite TB"). The 3.3.x queue is tracked as epic
[JCSDA/CRTMv3#357](https://github.com/JCSDA/CRTMv3/issues/357).

## The patch gate

A change is **3.2.x-eligible** only if ALL of the following hold:

1. **Default-path forward bit-identity.** The common regression suite's
   forward results (TB/radiance) are bit-identical on every default code
   path. Opt-in paths (a non-default `RT_Algorithm_Id`, an explicitly
   requested model or option) may change only when the change is a defect
   repair, and the defect and repair are documented in the release notes.
2. **No public interface changes.** No new or altered structure components,
   argument lists, module exports, or file formats that a downstream build
   (GSI, UFO, pycrtm) links against.
3. **No new or changed default coefficient files.** The pinned tarball and
   its checksums are untouched.
4. **Jacobian changes are zero-to-correct only.** A sensitivity that is
   identically zero (or nondeterministic garbage) today may become correct in
   a patch; a Jacobian that today returns a *different nonzero value* on a
   default path moves numerical results consumers may have tuned against,
   and queues for 3.3.x.

Everything failing any clause is **3.3.x (or later)** by definition:
"substantial numerical default values have changed" (BTJ's criterion).

Test-data staging, build hygiene, dead code, and documentation are always
patch-eligible.

---

## Already in 3.3.0 (context)

The TELSEM2 DA program (merge `067c111`): Release-2 uncertainty content,
total-derivative d(Tb)/d(emissivity) (G1/G2/G3/G6 — G1 and G2 change
K-matrix values, which is what makes this 3.3.0 material), the covariance
query API, `RTSolution%Surface_Emissivity_Std` (interface addition), the
class-gated TELSEM2+LandEM hybrid state Jacobians, and the class-consistent
snow/ice dispatch. See `docs/design/telsem2_landem_jacobians_plan.md`
(2026-08-13 decision record).

---

## Queued for 3.3.x (fails the patch gate)

| # | Item | Why it fails the gate | Where (as last verified) | Effort |
|---|------|----------------------|--------------------------|--------|
| 1 | **G4: fractional-cloud d(Tb)/d(emissivity)** — the reported Jacobian is cloudy-column-only; fold in the clear-column term as `(1-cc)*clear + cc*cloudy` | Changes K values for every `Cloud_Fraction < 1` scene | Clear-column capture at the clear AD call, combine in the four driver modules (`CRTM_K_Matrix_Module.f90` ~:1630-1643 and mirrors) | Medium; completes the 3.3.0 emissivity-Jacobian story |
| 2 | **Aerosol hygroscopic-growth adjoint** — `Atm_AD%Relative_Humidity` accumulates and is silently dropped; the AD of the `MR_to_RH`/`PPMV_to_MR` chain back to T and q does not exist | Changes default-path T and q Jacobians for aerosol-laden scenes | `CRTM_AerosolScatter.f90:786`; chain routines in `CRTM_Relative_Humidity` (no TL/AD forms exported) | Medium-high: new AD code. Do NOT map RH into `Atmosphere_K` (derived diagnostic; double-counting risk) |
| 3 | **Cloud `Water_Density` adjoint** — commented-out two-argument call has no implementation; the MW cloudy `Water_Content` Jacobian is missing its density-indexing term, plus a height leg via `dZ_m` | Changes default-path MW cloudy Jacobians | `CRTM_K_Matrix_Module.f90` ~:1131, `CRTM_Adjoint_Module.f90` ~:633 (commented), forward `CRTM_Active_Sensor.f90:127-160`, consumption `CRTM_CloudScatter.f90` :343/:614/:941 | Medium: adjoint must be written, not uncommented |
| 4 | **Snow/ice physical state Jacobians** — TL/AD are zero-stubs whose signatures do not accept `Surface_TL/AD`; the physical sea-ice model `NESDIS_SIce_Phy_EM` is commented out at its call site | Signature changes; enabling the physical ice model changes the ice forward | `CRTM_MW_Snow_SfcOptics.f90`, `CRTM_MW_Ice_SfcOptics.f90` (iVar is `INTEGER :: Dummy`) | High. Composes with Phase G: physical Jacobians where NESDIS runs, atlas+uncertainty where TELSEM2 runs |
| 5 | **PARMIO activation-policy consistency** — PARMIO is presence-activated while TELSEM2 is opt-in; the same foot-gun the TELSEM2 gate closed | Making PARMIO opt-in changes MW-water defaults >= 200 GHz for anyone with the file staged | `CRTM_LifeCycle.f90` (PARMIO load block); policy note at `RELEASE_NOTES_v3.2.0.md` §"PARMIO is presence-activated; TELSEM2 is opt-in" | Low code / high communication. Decide before 3.3.0 ships a second precedent |
| 6 | **Release-2 coefficient tarball roll** — ship the staged `TELSEM2.MWland.EmisCoeff.nc` (R2, 344 MB; emissivity payload bit-identical to R1) and update both md5 pins | New default coefficient file (gate clause 3) | `Get_CRTM_Binary_Files.sh`, `test/CMakeLists.txt` (precedent commit `5659f9d`); inventory action item in `REL-3.2.0_coefficient_inventory.md` | Deferred by BTJ until REL-3.3.0 matures. Until rolled, the TELSEM2 covariance/snow-ice ctests require a locally staged R2 file |
| 7 | **New high-frequency surface physics** — canopy water content, MW snow grain size/density, ice density/roughness as differentiable state | New forward physics changes forward wherever added | Adjudicated "structurally impossible without new physics" in `docs/design/surface_jacobians_281.md`; the closure proposal in `docs/design/issue-281-comment.md` was never posted to #281 | Program-scale (3.3+/3.4). First step is administrative: post the adjudication |
| 8 | **Vector-RT (`n_Stokes > 1`) energy-consistency completion** — the (V,H)→Stokes reflectivity conversion covers the intensity block; U/V and off-diagonal consistency remain | Opt-in path only, so technically patch-eligible under clause 1 — queued to 3.3.x by blast radius, not by the letter of the gate | `CRTM_SfcOptics.f90` coupled-polarization branch (in-code follow-up note) | Medium |

## 3.2.x-eligible (passes the patch gate)

| # | Item | Why it passes | Where | Effort |
|---|------|---------------|-------|--------|
| 1 | **Fastem1 legacy SST Jacobian** — `dEV/dEH_dTs` declared, read by TL/AD, never assigned; K w.r.t. SST is silently zero on the legacy `.NOT. Use_New_MWSSEM` leg | Zero-to-correct on an opt-out path; forward untouched. Already adjudicated 3.2.0-eligible in the #281 audit | `CRTM_MW_Water_SfcOptics.f90` (~:99/:101 declarations; TL ~:531/:533; AD ~:761/:767) | Small: FD `dE/dTs` around the existing Fastem1 call (~15 lines), per audit option (b) |
| 2 | **Staged test data gaps** — `SNICAR.VISsnow.EmisCoeff.nc` and `CloudCoeff_Exp_Full6.nc` missing; the standing 2 ctest failures (`test_EmisCoeff_IO` #11, `test_VectorRT_ScalarLimit`) | Test infrastructure only | Staged fix tree / testinput symlinks | Small (locate or regenerate the files, stage, optionally gate the tests on presence) |
| 3 | **Dead-code hygiene** — `Generate_CRTM_Stats.f90` references the renamed `Compute_Switch` (does not compile; unbuilt); unbuilt UWIREMIS `mod_iratlas.f90` needs keep-or-delete | Not built; zero numerical surface | `src/Statistics/Generate_CRTM_Stats/Generate_CRTM_Stats.f90:734,811`; `src/SfcOptics/IR_Land/UWIREMIS/mod_iratlas.f90` | Small |
| 4 | **Documentation debt** — post the #281 adjudication comment (draft exists in-repo, never posted); release-notes corrections; keep this triage file current | Docs only | `docs/design/issue-281-comment.md`; GitHub #281 | Small |

## Borderline calls (decide per consumer impact)

- **SOI-class adjoint hygiene backports.** The G6 fix (SOI AD zeroing its
  surface duals on entry) repairs nondeterministic channel-carryover on an
  opt-in solver — allowed in a patch under gate clause 4 if any 3.2.x
  consumer runs SOI K-matrix multi-channel. Shipped in 3.3.0; backport only
  on demand.
- **G2/G3-style Jacobian corrections** change nonzero `RTSolution_K%
  Surface_Emissivity` values and therefore stay 3.3.x, even though they are
  defect repairs: DA consumers may have tuned against the old (biased)
  values.
