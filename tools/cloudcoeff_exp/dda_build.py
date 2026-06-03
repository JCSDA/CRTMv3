#!/usr/bin/env python3
"""
dda_build.py -- build a PHYSICAL experimental CloudCoeff (v1) from the DDSCAT archive.

Per habit (one netCDF4 file):
  1. sample shape dirs from the shape_mass_dimension CSV (paths only; CSV masses are
     unreliable -- mass/size come from each .avg's aeff).
  2. parse .avg over all freqs -> single-particle cross-sections sigma_{ext,abs,sca,bk}
     (= Q * pi*aeff^2) and S_11(theta); D_eq = 2*aeff, mass = (4/3)pi aeff^3 rho_ice.
  3. PSD-integrate over a normalized gamma in D_eq, N(D)=(D/Dm)^mu exp(-(4+mu)D/Dm):
       ke=Bext/M, ka=Babs/M, kb=Bbk/M, w=Bsca/Bext ;  bulk P11 = int N*S11 dD.
  4. Legendre-project bulk P11 -> chi_l (alpha1); n_Legendre_Eff = #(|chi_l|>=tol).
  5. temperature: scale ka (=> ke, w) by Matzler ice eps''(T,f)/eps''(266,f) [option 2];
     kb, g and the phase shape are ~T-independent (replicated across the T axis).
  6. write per-habit netCDF4 via build_cloudcoeff_exp.write_netcdf (size axis = D_eq).

NETCDF4 only. Size coordinate = D_eq (mass proxy; m-D verifiable).
"""
import os, glob, csv, math, argparse
import numpy as np
from dda_parse import parse_avg, RHO_ICE
from build_cloudcoeff_exp import write_netcdf

T_REF = 266.0  # K, the temperature of the DDA refractive index (ior_266K)


def matzler_eps_ice(T, f_ghz):
    """Maetzler (2006) pure-ice complex permittivity. T[K], f[GHz]."""
    eps_re = 3.1884 + 9.1e-4 * (T - 273.0)
    theta = 300.0 / T - 1.0
    alpha = (0.00504 + 0.0062 * theta) * math.exp(-22.1 * theta)
    beta = ((0.0207 / T) * math.exp(335.0 / T) / (math.exp(335.0 / T) - 1.0) ** 2
            + 1.16e-11 * f_ghz ** 2 + math.exp(-9.963 + 0.0372 * (T - 273.16)))
    return eps_re, alpha / f_ghz + beta * f_ghz


def sample_shape_leaves(csv_path, habit, n):
    paths = []
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            p = row["path"].lstrip("./")
            if p.split("/")[0] == habit and p.endswith("/ddscat07/shape.dat"):
                paths.append(p[:-len("/ddscat07/shape.dat")])
    if n and n < len(paths):
        idx = np.linspace(0, len(paths) - 1, n).astype(int)
        paths = [paths[i] for i in idx]
    return paths


def extract_particle(root, leaf):
    avgs = sorted(glob.glob(os.path.join(root, leaf, "ddscat07", "*GHz", "w000r000.avg")))
    if not avgs:
        return None
    f, sext, sabs, ssca, sbk, s11 = [], [], [], [], [], []
    theta_ref, aeff = None, None
    for ap in avgs:
        a = parse_avg(ap)
        if "qext" not in a or not a["theta"] or "freq_ghz" not in a:
            continue
        aeff = a["aeff_um"]
        area = math.pi * (aeff * 1e-6) ** 2                       # m^2
        th, s = np.array(a["theta"]), np.array(a["s11"])
        if theta_ref is None:
            theta_ref = th
        if s.shape != theta_ref.shape:
            s = np.interp(theta_ref, th, s)
        f.append(round(a["freq_ghz"]))
        sext.append(a["qext"] * area); sabs.append(a["qabs"] * area)
        ssca.append(a["qsca"] * area); sbk.append(a["qbk"] * area)
        s11.append(s)
    if aeff is None or not f:
        return None
    return dict(D_eq=2.0 * aeff * 1e-6,
                mass=(4.0 / 3.0) * math.pi * (aeff * 1e-6) ** 3 * RHO_ICE,
                f=np.array(f), sext=np.array(sext), sabs=np.array(sabs),
                ssca=np.array(ssca), sbk=np.array(sbk), theta=theta_ref,
                s11=np.array(s11))                                # (nfreq, ntheta)


def legendre_coeffs(theta_deg, P_theta, L_max):
    """chi_l = 0.5 * int_{-1}^{1} Pn(mu) P_l(mu) dmu, Pn normalized so 0.5*int Pn dmu = 1."""
    mu = np.cos(np.deg2rad(theta_deg))
    order = np.argsort(mu)
    mu, P = mu[order], np.maximum(P_theta[order], 0.0)
    norm = 0.5 * np.trapz(P, mu)
    if norm <= 0:
        return np.zeros(L_max)
    Pn = P / norm
    chi = np.empty(L_max)
    for l in range(L_max):
        Pl = np.polynomial.legendre.legval(mu, [0] * l + [1])
        chi[l] = 0.5 * np.trapz(Pn * Pl, mu)
    return chi


def build_habit(root, csv_path, habit, habit_id, habit_name, n_shapes,
                out_path, mu_val=0.0, L_max=64, n_dm=24, n_t=5, tol=1.0e-3):
    leaves = sample_shape_leaves(csv_path, habit, n_shapes)
    parts = [p for p in (extract_particle(root, lf) for lf in leaves) if p]
    if not parts:
        raise RuntimeError("no particles parsed")
    print(f"  parsed {len(parts)}/{len(leaves)} shapes")

    # Keep only frequencies with adequate shape coverage (the archive's DDA runs
    # are complete to ~190 GHz; >=205 GHz is only partially computed). Building a
    # frequency from a handful of shapes gives a biased/extrapolated PSD integral.
    all_freqs = sorted(set(int(x) for p in parts for x in p["f"]))
    fcount = {f: sum(1 for p in parts if f in set(int(x) for x in p["f"])) for f in all_freqs}
    min_keep = max(8, int(0.7 * len(parts)))
    freqs = np.array([f for f in all_freqs if fcount[f] >= min_keep])
    dropped = [(f, fcount[f]) for f in all_freqs if fcount[f] < min_keep]
    if dropped:
        print("  WARNING undersampled freqs dropped (need >=%d shapes): %s"
              % (min_keep, ", ".join("%dGHz(%d)" % (f, c) for f, c in dropped)))
    print("  kept %d freqs: %d-%d GHz" % (len(freqs), freqs[0], freqs[-1]))
    NF = len(freqs)
    ntheta = len(parts[0]["theta"]); theta = parts[0]["theta"]
    f_index = {f: i for i, f in enumerate(freqs)}

    # gather single-particle samples per frequency
    D = np.array([p["D_eq"] for p in parts])
    order = np.argsort(D)
    Dmin, Dmax = D.min(), D.max()
    Dg = np.geomspace(Dmin, Dmax, 200)                            # fine integration grid (m)

    # interpolate sigma(D_eq) and S11(D_eq,theta) onto Dg, per frequency
    sext_g = np.zeros((NF, len(Dg))); sabs_g = np.zeros_like(sext_g)
    ssca_g = np.zeros_like(sext_g);   sbk_g = np.zeros_like(sext_g)
    s11_g = np.zeros((NF, len(Dg), ntheta))
    for jf, f in enumerate(freqs):
        rows = [(p["D_eq"], p) for p in parts if f in set(int(x) for x in p["f"])]
        rows.sort(key=lambda r: r[0])
        Dr = np.array([r[0] for r in rows])
        def col(key):
            return np.array([r[1][key][list(r[1]["f"]).index(f)] for r in rows])
        sext_g[jf] = np.interp(Dg, Dr, col("sext"))
        sabs_g[jf] = np.interp(Dg, Dr, col("sabs"))
        ssca_g[jf] = np.interp(Dg, Dr, col("ssca"))
        sbk_g[jf]  = np.interp(Dg, Dr, col("sbk"))
        s11r = np.array([r[1]["s11"][list(r[1]["f"]).index(f)] for r in rows])   # (nrow,ntheta)
        for it in range(ntheta):
            s11_g[jf, :, it] = np.interp(Dg, Dr, s11r[:, it])

    mass_g = RHO_ICE * (math.pi / 6.0) * Dg ** 3                  # kg, exact for D_eq

    # PSD-integrate over (Dm) at fixed mu -> base (T_REF) bulk properties
    Dm_grid = np.geomspace(2.0 * Dmin, 0.7 * Dmax, n_dm)          # mass-weighted mean diameter (m)
    ke0 = np.zeros((n_dm, NF)); ka0 = np.zeros_like(ke0); kb0 = np.zeros_like(ke0)
    chi = np.zeros((n_dm, NF, L_max)); neff = np.zeros((n_dm, NF), dtype=int)
    for idm, Dm in enumerate(Dm_grid):
        Nd = (Dg / Dm) ** mu_val * np.exp(-(4.0 + mu_val) * Dg / Dm)
        M = np.trapz(Nd * mass_g, Dg)
        for jf in range(NF):
            Bext = np.trapz(Nd * sext_g[jf], Dg)
            ke0[idm, jf] = Bext / M
            ka0[idm, jf] = np.trapz(Nd * sabs_g[jf], Dg) / M
            kb0[idm, jf] = np.trapz(Nd * sbk_g[jf], Dg) / M
            Pb = np.trapz(Nd[:, None] * s11_g[jf], Dg, axis=0)    # bulk P11(theta)
            c = legendre_coeffs(theta, Pb, L_max)
            chi[idm, jf] = c
            sig = np.where(np.abs(c) >= tol)[0]
            neff[idm, jf] = int(sig[-1] + 1) if sig.size else 1

    # temperature axis via Matzler absorption rescaling (option 2)
    temp = np.linspace(233.0, 273.0, n_t)
    NH, NMu = 1, 1
    KE = np.zeros((NH, n_dm, NMu, n_t, NF)); KA = np.zeros_like(KE)
    KB = np.zeros_like(KE); GG = np.zeros_like(KE)
    NEFF = np.zeros((NH, n_dm, NMu, n_t, NF), dtype="i4")
    PCO = np.zeros((NH, n_dm, NMu, n_t, NF, L_max, 6))
    ks0 = ke0 - ka0                                               # scattering part (T-independent)
    for it, T in enumerate(temp):
        for jf, f in enumerate(freqs):
            r = matzler_eps_ice(T, f)[1] / matzler_eps_ice(T_REF, f)[1]
            kaT = ka0[:, jf] * r
            keT = ks0[:, jf] + kaT
            KA[0, :, 0, it, jf] = kaT
            KE[0, :, 0, it, jf] = keT
            KB[0, :, 0, it, jf] = kb0[:, jf]
            GG[0, :, 0, it, jf] = chi[:, jf, 1]                   # asymmetry ~ chi_1
            NEFF[0, :, 0, it, jf] = neff[:, jf]
            PCO[0, :, 0, it, jf, :, 0] = chi[:, jf, :]            # element 1 = alpha1

    habits = [(habit_id, habit_name, 1, RHO_ICE * math.pi / 6.0, 3.0)]   # m=a*D_eq^b, b=3
    # write_netcdf expects Frequency from its module global; set it to our freqs
    import build_cloudcoeff_exp as B
    B.FREQ_GHZ = freqs.astype(float)
    write_netcdf(out_path, Dm_grid * 1e6, np.array([mu_val]), temp,
                 KE, KA, KB, GG, NEFF, PCO, L_max, synthetic=False, habits=habits)
    print(f"  Dm range {Dm_grid[0]*1e6:.0f}-{Dm_grid[-1]*1e6:.0f} um, "
          f"freqs {freqs[0]}-{freqs[-1]} GHz, n_Legendre_Eff {NEFF.min()}-{NEFF.max()}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/mnt/d/WSL_temp_storage/full_db")
    ap.add_argument("--csv", default="/mnt/d/WSL_temp_storage/full_db/shape_mass_dimension_a_compare.csv")
    ap.add_argument("--habit", default="aggregate")
    ap.add_argument("--habit-id", type=int, default=4)            # SNOW_CLOUD by default
    ap.add_argument("--habit-name", default="AggregateSnow")
    ap.add_argument("--n-shapes", type=int, default=25)
    ap.add_argument("--mu", type=float, default=0.0)
    ap.add_argument("-o", "--output", default="CloudCoeff_Exp_aggregate_v1.nc")
    a = ap.parse_args()
    print(f"Building experimental LUT: habit={a.habit} id={a.habit_id} n_shapes={a.n_shapes}")
    build_habit(a.root, a.csv, a.habit, a.habit_id, a.habit_name, a.n_shapes,
                a.output, mu_val=a.mu)
    print(f"wrote {a.output}")


if __name__ == "__main__":
    main()
