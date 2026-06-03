#!/usr/bin/env python3
"""
ren_build.py -- build an experimental CloudCoeff (v1) from the Ren et al. (2022)
"Snow single-scattering database" (II-TM/IGOM), MW part.  doi:10.18738/T8/LGJ9SA

Self-contained (only numpy + netCDF4) so it can run on the HPC where the ~5.8 GB
data live, emitting a small netCDF to copy back.

Layout expected (from README_Snowflake.txt):
  <root>/SnowflakeModel_MW/Scheme{1..4}/T{190,210,230,250,270}/
        isca.dat   : 70 freq x 131 size rows (size = inner loop), 7 cols:
                     wl_um, Dmax_um, vol_um3, projarea_um2, Qext, ssa, g
        P11.dat    : line 1 = 498 scattering angles (deg); then 70*131 rows x 498
  <root>/SnowflakeModel_MW/MW_freq.txt  : 70 frequencies (GHz)   [optional]

v1 scope: scalar phase function (alpha1 from P11); ke/ka/kb/g + per-entry truncation;
REAL 5-temperature optics (no Maetzler rescaling). Full 6-element GSF (alpha2-4,beta1-2
from P12/P22/P33/P43/P44) is a documented follow-on. Size coordinate: Dmax.
"""
import os, glob, argparse
import numpy as np
from numpy.polynomial.legendre import legval
from netCDF4 import Dataset, stringtochar

RHO_ICE = 917.0
C_LIGHT = 2.99792458e8
N_SIZE  = 131
TEMPS   = [190, 210, 230, 250, 270]
SCHEME_NAME = {1: "ConstDensity0.1", 2: "Thompson08", 3: "Heymsfield04", 4: "Brandes07"}
SCHEME_CLOUD_ID = {1: 4, 2: 4, 3: 4, 4: 4}     # map to SNOW_CLOUD(4); adjust per use


def read_scheme_T(root, scheme, T):
    """Return per-(freq,size) arrays for one scheme/temperature."""
    hits = glob.glob(os.path.join(root, "**", "Scheme%d" % scheme, "T%d" % T, "isca.dat"),
                     recursive=True)
    if not hits:
        raise FileNotFoundError("Scheme%d/T%d/isca.dat not found under %s" % (scheme, T, root))
    base = os.path.dirname(hits[0])
    isca = np.loadtxt(os.path.join(base, "isca.dat"))          # (70*131, 7)
    nf = isca.shape[0] // N_SIZE
    isca = isca.reshape(nf, N_SIZE, 7)
    wl, Dmax, vol, area, qext, ssa, g = (isca[..., k] for k in range(7))
    with open(os.path.join(base, "P11.dat")) as fh:
        ang = np.array(fh.readline().split(), dtype=float)     # 498 angles (deg)
    p11 = np.loadtxt(os.path.join(base, "P11.dat"), skiprows=1).reshape(nf, N_SIZE, ang.size)
    return dict(nf=nf, wl=wl, Dmax=Dmax, vol=vol, area=area, qext=qext, ssa=ssa, g=g,
                ang=ang, p11=p11)


def legendre_alpha1(theta_deg, P11_bulk, L_max):
    """alpha1(l) = 0.5*int P11(mu)P_l(mu)dmu; Ren P11 already has 0.5*int P11 dmu = 1."""
    mu = np.cos(np.deg2rad(theta_deg)); o = np.argsort(mu)
    mu, P = mu[o], np.maximum(P11_bulk[o], 0.0)
    norm = 0.5 * np.trapz(P, mu)
    if norm <= 0:
        return np.zeros(L_max)
    Pn = P / norm
    return np.array([0.5 * np.trapz(Pn * legval(mu, [0]*l + [1]), mu) for l in range(L_max)])


def build(root, scheme, out, mu_val=0.0, n_dm=30, L_max=64, tol=1e-3, freqs=None):
    # frequency axis (prefer MW_freq.txt; else derive from wavelength column)
    s0 = read_scheme_T(root, scheme, TEMPS[0])
    if freqs is None:
        ff = os.path.join(root, "SnowflakeModel_MW", "MW_freq.txt")
        freqs = (np.loadtxt(ff) if os.path.exists(ff)
                 else C_LIGHT / (s0["wl"][:, 0] * 1e-6) / 1e9)
    freqs = np.round(np.asarray(freqs, float), 3)
    NF = s0["nf"]; NT = len(TEMPS)

    # single-particle size axis (Dmax, microns) -- assume same grid across freq/T
    D = s0["Dmax"][0, :] * 1e-6                       # m
    mass = s0["vol"][0, :] * 1e-18 * RHO_ICE          # kg (ice volume * rho)
    order = np.argsort(D); D = D[order]; mass = mass[order]
    Dg = np.geomspace(D.min(), D.max(), 400)
    massg = np.interp(Dg, D, mass)
    Dm_grid = np.geomspace(2*D.min(), 0.7*D.max(), n_dm)

    KE = np.zeros((1, n_dm, 1, NT, NF)); KA = np.zeros_like(KE); KB = np.zeros_like(KE)
    GG = np.zeros_like(KE); NEFF = np.zeros((1, n_dm, 1, NT, NF), dtype="i4")
    PCO = np.zeros((1, n_dm, 1, NT, NF, L_max, 6))

    for it, T in enumerate(TEMPS):
        s = read_scheme_T(root, scheme, T)
        sext = (s["qext"] * s["area"] * 1e-12)[:, order]        # m^2  (area um^2 -> m^2)
        ssca = sext * s["ssa"][:, order]
        sabs = sext - ssca
        p11  = s["p11"][:, order, :]                            # (NF, NSIZE, nang)
        # backscatter: sigma_bk = sigma_sca * P11(180deg)
        i180 = int(np.argmax(s["ang"]))
        sbk  = ssca * p11[:, :, i180]
        nang = s["ang"].size
        for jf in range(NF):
            se = np.interp(Dg, D, sext[jf]); sa = np.interp(Dg, D, sabs[jf])
            sc = np.interp(Dg, D, ssca[jf]); sb = np.interp(Dg, D, sbk[jf])
            # interpolate P11(size, angle) onto the integration grid ONCE per frequency
            p11g = np.empty((Dg.size, nang))
            for ia in range(nang):
                p11g[:, ia] = np.interp(Dg, D, p11[jf, :, ia])
            for idm, Dm in enumerate(Dm_grid):
                Nd = (Dg/Dm)**mu_val * np.exp(-(4.0+mu_val)*Dg/Dm)
                M = np.trapz(Nd*massg, Dg)
                KE[0, idm, 0, it, jf] = np.trapz(Nd*se, Dg)/M
                KA[0, idm, 0, it, jf] = np.trapz(Nd*sa, Dg)/M
                KB[0, idm, 0, it, jf] = np.trapz(Nd*sb, Dg)/M
                w = Nd*sc                                  # sca-weighted bulk P11 -> alpha1
                P11b = np.trapz(w[:, None]*p11g, Dg, axis=0) / max(np.trapz(w, Dg), 1e-30)
                a1 = legendre_alpha1(s["ang"], P11b, L_max)
                GG[0, idm, 0, it, jf] = a1[1]
                sig = np.where(np.abs(a1) >= tol)[0]
                NEFF[0, idm, 0, it, jf] = int(sig[-1]+1) if sig.size else 1
                PCO[0, idm, 0, it, jf, :, 0] = a1

    write_netcdf(out, freqs, Dm_grid*1e6, np.array([mu_val]), np.array(TEMPS, float),
                 KE, KA, KB, GG, NEFF, PCO, L_max, scheme)
    print("wrote %s  (scheme=%s, %d freq %g-%g GHz, %d T, %d Dm)"
          % (out, SCHEME_NAME[scheme], NF, freqs.min(), freqs.max(), NT, n_dm))


def write_netcdf(path, freq, dm, mu, temp, ke, ka, kb, g, neff, pcoeff, L_max, scheme):
    nc = Dataset(path, "w", format="NETCDF4")
    for n, v in [("n_Frequency", len(freq)), ("n_Dm", len(dm)), ("n_Mu", len(mu)),
                 ("n_Temperature", len(temp)), ("n_Habit", 1), ("n_Legendre", L_max),
                 ("n_Phase_Elements", 6), ("nchar", 32)]:
        nc.createDimension(n, v)
    nc.Scheme = "CRTM-Exp"; nc.Release = np.int32(1); nc.Version = np.int32(1)
    nc.PSD = "normalized_gamma"; nc.Orientation = "random"
    nc.Source = "Ren et al. 2022 snow II-TM/IGOM, doi:10.18738/T8/LGJ9SA, %s" % SCHEME_NAME[scheme]
    nc.Note = "v1: alpha1/scalar only; full 6-element GSF pending; real 5-T optics"

    def cv(name, dims, data, **k):
        x = nc.createVariable(name, "f8", dims, **k); x[:] = data; return x
    cv("Frequency", ("n_Frequency",), freq).units = "GHz"
    cv("Dm", ("n_Dm",), dm).units = "microns"
    cv("Mu", ("n_Mu",), mu); cv("Temperature", ("n_Temperature",), temp).units = "K"
    nc.createVariable("Habit_Id", "i4", ("n_Habit",))[:] = [SCHEME_CLOUD_ID[scheme]]
    nc.createVariable("Habit_Phase", "i4", ("n_Habit",))[:] = [1]
    cv("mD_a", ("n_Habit",), [0.0]); cv("mD_b", ("n_Habit",), [0.0])
    nc.createVariable("Habit_Name", "S1", ("n_Habit", "nchar"))[:] = \
        stringtochar(np.array([("Snow_%s" % SCHEME_NAME[scheme]).ljust(32)], dtype="S32"))
    bd = ("n_Habit", "n_Dm", "n_Mu", "n_Temperature", "n_Frequency")
    for nm, arr in (("ke", ke), ("ka", ka), ("kb", kb), ("g", g)):
        cv(nm, bd, arr, zlib=True, complevel=4)
    nc.createVariable("n_Legendre_Eff", "i4", bd, zlib=True, complevel=4)[:] = neff
    cv("pcoeff", ("n_Habit","n_Dm","n_Mu","n_Temperature","n_Frequency","n_Legendre","n_Phase_Elements"),
       pcoeff, zlib=True, complevel=4)
    nc.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="dir containing SnowflakeModel_MW/")
    ap.add_argument("--scheme", type=int, default=2, choices=[1, 2, 3, 4])
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--mu", type=float, default=0.0)
    ap.add_argument("--n-dm", type=int, default=30)
    a = ap.parse_args()
    out = a.output or "CloudCoeff_Exp_RenSnow_%s.nc" % SCHEME_NAME[a.scheme]
    build(a.root, a.scheme, out, mu_val=a.mu, n_dm=a.n_dm)
