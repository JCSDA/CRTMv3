#!/usr/bin/env python3
"""
build_cloudcoeff_exp.py  --  Builder for the experimental CRTM CloudCoeff (v1).

Implements the schema in CloudCoeff_Experimental_Schema_v1.md:
  - explicit habit axis + (Dm, mu) PSD-moment axes + temperature (all phases)
  - one full / variable-length GSF (Greek-constant) expansion per entry, with a
    per-entry effective truncation order (n_Legendre_Eff) -- DECOUPLED from streams
  - 6 phase elements (alpha1..alpha4, beta1, beta2); scalar RT uses element 1 only
  - bulk optics stored as ke, ka, kb (+ g); w = (ke-ka)/ke derived at runtime

Two modes:
  --synthetic : emit a spec-compliant netCDF with smooth, plausible *synthetic*
                values (NOT physical). Phase-function peakedness rises with size
                parameter so the high-Legendre / decoupling regime is exercised.
                Lets the CRTM-side reader be developed/tested before the real LUT.
  --from-dda  : (scaffold) integrate the in-house DDSCAT archive over a normalized
                gamma PSD and project the Mueller matrices onto GSF coefficients.

ARRAY ORDERING (canonical, see schema note):
  Fortran reader declares fastest->slowest:
     bulk    : (Frequency, Temperature, Mu, Dm, Habit)
     pcoeff  : (Phase_Elements, Legendre, Frequency, Temperature, Mu, Dm, Habit)
  The netCDF file therefore stores dims in the REVERSE (C/CDL) order; this builder
  creates variables in that reversed order so the Fortran netCDF IO reads the
  canonical Fortran order directly.
"""

import argparse
import numpy as np
from netCDF4 import Dataset, stringtochar

SCHEME = "CRTM-Exp"
RELEASE, VERSION = 1, 1

# Habit key = existing CRTM extended cloud-type integers (decision #5).
# (id, name, phase[0=liquid,1=frozen], m-D a, m-D b)  -- a,b are placeholders in synthetic mode.
HABITS = [
    (4,  "SnowDefault",          1, 0.069, 2.0),
    (16, "LargePlateAggregate",  1, 0.040, 1.9),
    (14, "EvansSnowAggregate",   1, 0.020, 1.8),
    (9,  "SixBulletRosette",     1, 0.060, 2.1),
    (25, "LiquidSphere",         0, 523.6, 3.0),
]

# Native-ish MW/sub-mm frequency nodes (GHz), matching the DDA archive set (decision #4).
FREQ_GHZ = np.array([3., 5., 11., 14., 19., 24., 36., 89., 94., 150., 166., 176., 180.,
                     186., 190., 205., 240., 325., 380., 425., 462., 500., 640., 683., 874.])
C_LIGHT = 2.99792458e8  # m/s


def make_grids(n_dm=40, n_mu=1, n_t=5):
    dm_um = np.geomspace(50.0, 8000.0, n_dm)             # mass-weighted mean diameter, microns
    mu    = np.array([0.0]) if n_mu == 1 else np.linspace(0.0, 6.0, n_mu)
    temp  = np.linspace(233.0, 273.0, n_t)               # K, applies to all phases
    return dm_um, mu, temp


def synthetic_entry(freq_ghz, dm_um, mu, temp_k, phase, L_max):
    """Smooth, plausible SYNTHETIC single-entry optics. Not physical."""
    lam_um = (C_LIGHT / (freq_ghz * 1e9)) * 1e6
    x = np.pi * dm_um / lam_um                            # size parameter
    # asymmetry grows with x (more forward-peaked for larger/higher-freq particles)
    g = 0.95 * x / (1.0 + x)
    # single-scatter albedo: more scattering at larger x; liquid a bit more absorbing
    w = (0.2 if phase == 0 else 0.05) + 0.9 * x / (3.0 + x)
    w = min(w, 0.999)
    # mass extinction (m^2/kg): smooth, O(0.01-10)
    ke = 0.05 * (x ** 2) / (1.0 + 0.3 * x) * (917.0 / max(dm_um, 1.0)) * 1.0e-3 + 1.0e-3
    ka = ke * (1.0 - w)
    kb = ke * w * (1.0 - g) * 0.5
    # GSF expansion: element 1 (alpha1) = Henyey-Greenstein chi_l = g^l (dimensionless)
    l = np.arange(L_max)
    a1 = g ** l
    # other elements: synthetic placeholders (real builder projects DDA Mueller matrices)
    a2 = 0.9 * a1
    a3 = 0.1 * a1
    a4 = (0.6 + 0.4 * (phase == 1)) * a1
    b1 = np.zeros(L_max); b1[2] = -0.04 * (1.0 - g)       # small Rayleigh-like P12 bump
    b2 = 0.2 * b1
    pc = np.stack([a1, a2, a3, a4, b1, b2], axis=0)        # (6, L_max)
    # effective truncation: where alpha1 drops below tol (capped to L_max)
    tol = 1.0e-7
    sig = np.where(np.abs(a1) >= tol)[0]
    n_eff = int(sig[-1] + 1) if sig.size else 1
    return ke, ka, kb, g, pc, n_eff


def build_synthetic(path, L_max=64, n_dm=40, n_mu=1, n_t=5):
    dm, mu, temp = make_grids(n_dm, n_mu, n_t)
    NF, NT, NMU, NDM, NH, NP = len(FREQ_GHZ), len(temp), len(mu), len(dm), len(HABITS), 6

    ke = np.zeros((NH, NDM, NMU, NT, NF))
    ka = np.zeros_like(ke); kb = np.zeros_like(ke); gg = np.zeros_like(ke)
    neff = np.zeros((NH, NDM, NMU, NT, NF), dtype="i4")
    pcoeff = np.zeros((NH, NDM, NMU, NT, NF, L_max, NP))   # netCDF (reversed) order

    for ih, (hid, hname, hph, a, b) in enumerate(HABITS):
        for idm in range(NDM):
            for imu in range(NMU):
                for it in range(NT):
                    for jf in range(NF):
                        e_ke, e_ka, e_kb, e_g, pc, ne = synthetic_entry(
                            FREQ_GHZ[jf], dm[idm], mu[imu], temp[it], hph, L_max)
                        ke[ih, idm, imu, it, jf] = e_ke
                        ka[ih, idm, imu, it, jf] = e_ka
                        kb[ih, idm, imu, it, jf] = e_kb
                        gg[ih, idm, imu, it, jf] = e_g
                        neff[ih, idm, imu, it, jf] = ne
                        pcoeff[ih, idm, imu, it, jf, :, :] = pc.T   # (L_max, 6)
    write_netcdf(path, dm, mu, temp, ke, ka, kb, gg, neff, pcoeff, L_max, synthetic=True)


def write_netcdf(path, dm, mu, temp, ke, ka, kb, gg, neff, pcoeff, L_max, synthetic):
    NH = len(HABITS)
    nc = Dataset(path, "w", format="NETCDF4")
    # dimensions
    nc.createDimension("n_Frequency", len(FREQ_GHZ))
    nc.createDimension("n_Dm", len(dm))
    nc.createDimension("n_Mu", len(mu))
    nc.createDimension("n_Temperature", len(temp))
    nc.createDimension("n_Habit", NH)
    nc.createDimension("n_Legendre", L_max)
    nc.createDimension("n_Phase_Elements", 6)
    nc.createDimension("nchar", 32)
    # global attributes
    nc.Scheme = SCHEME
    nc.Release = np.int32(RELEASE); nc.Version = np.int32(VERSION)
    nc.PSD = "normalized_gamma"
    nc.Orientation = "random"
    nc.pcoeff_convention = ("element1=alpha1(P11) as HG chi_l=g**l; elements 2-6 = "
                            "alpha2,alpha3,alpha4,beta1,beta2 (GSF Greek constants)")
    nc.Title = "Experimental CRTM Cloud Optical Properties (v1)"
    nc.Note = "SYNTHETIC, non-physical values" if synthetic else "Built from DDSCAT/DDA archive"

    def cvar(name, dims, data, **kw):
        v = nc.createVariable(name, "f8", dims, **kw); v[:] = data; return v

    # axes
    cvar("Frequency", ("n_Frequency",), FREQ_GHZ).units = "GHz"
    cvar("Dm", ("n_Dm",), dm).units = "microns"
    cvar("Mu", ("n_Mu",), mu).units = "1"
    cvar("Temperature", ("n_Temperature",), temp).units = "K"
    # habit metadata
    ids   = np.array([h[0] for h in HABITS], dtype="i4")
    phase = np.array([h[2] for h in HABITS], dtype="i4")
    mda   = np.array([h[3] for h in HABITS]); mdb = np.array([h[4] for h in HABITS])
    nc.createVariable("Habit_Id", "i4", ("n_Habit",))[:] = ids
    nc.createVariable("Habit_Phase", "i4", ("n_Habit",))[:] = phase
    cvar("mD_a", ("n_Habit",), mda); cvar("mD_b", ("n_Habit",), mdb)
    names = nc.createVariable("Habit_Name", "S1", ("n_Habit", "nchar"))
    names[:] = stringtochar(np.array([h[1].ljust(32) for h in HABITS], dtype="S32"))
    # bulk optics  (netCDF dims reversed from Fortran (Frequency,Temperature,Mu,Dm,Habit))
    bdims = ("n_Habit", "n_Dm", "n_Mu", "n_Temperature", "n_Frequency")
    cvar("ke", bdims, ke, zlib=True, complevel=4).units = "m^2/kg"
    cvar("ka", bdims, ka, zlib=True, complevel=4).units = "m^2/kg"
    cvar("kb", bdims, kb, zlib=True, complevel=4).units = "m^2/kg"
    cvar("g",  bdims, gg, zlib=True, complevel=4).units = "1"
    nc.createVariable("n_Legendre_Eff", "i4", bdims, zlib=True, complevel=4)[:] = neff
    # phase expansion (netCDF dims reversed from Fortran (Phase,Legendre,Freq,T,Mu,Dm,Habit))
    pdims = ("n_Habit", "n_Dm", "n_Mu", "n_Temperature", "n_Frequency", "n_Legendre", "n_Phase_Elements")
    cvar("pcoeff", pdims, pcoeff, zlib=True, complevel=4)
    nc.close()
    print(f"wrote {path}")


def build_from_dda(path, archive_root, **kw):
    raise NotImplementedError(
        "Real DDA ingest scaffold. Steps:\n"
        "  1. Walk archive (habit/density/size/freq); read geom_param.dat (mass,Dmax,area,aspect)\n"
        "  2. Read .avg Mueller matrices S11,S12,S22,S33,S34,S44 over scattering angle\n"
        "  3. Integrate single-particle Qext/Qabs/Qbk + Mueller over normalized-gamma PSD(Dm,mu)\n"
        "     to bulk ke,ka,kb and the angular phase matrix\n"
        "  4. Project phase matrix onto GSF -> alpha1..alpha4,beta1,beta2; set n_Legendre_Eff\n"
        "  5. write_netcdf(...) with synthetic=False")


def main():
    ap = argparse.ArgumentParser(description="Build experimental CloudCoeff (v1) netCDF")
    ap.add_argument("-o", "--output", default="CloudCoeff_Exp_synthetic_v1.nc")
    ap.add_argument("--synthetic", action="store_true", help="emit synthetic spec-compliant LUT")
    ap.add_argument("--from-dda", metavar="ARCHIVE_ROOT", help="(scaffold) build from DDA archive")
    ap.add_argument("--L-max", type=int, default=64)
    ap.add_argument("--n-dm", type=int, default=40)
    ap.add_argument("--n-mu", type=int, default=1)
    ap.add_argument("--n-t", type=int, default=5)
    args = ap.parse_args()
    if args.from_dda:
        build_from_dda(args.output, args.from_dda, L_max=args.L_max)
    else:
        build_synthetic(args.output, L_max=args.L_max, n_dm=args.n_dm, n_mu=args.n_mu, n_t=args.n_t)


if __name__ == "__main__":
    main()
