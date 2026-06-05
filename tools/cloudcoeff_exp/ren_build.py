#!/usr/bin/env python3
"""
ren_build.py (v2) -- physical experimental CloudCoeff from the Ren et al. (2022)
snow single-scattering database (II-TM/IGOM), microwave part. doi:10.18738/T8/LGJ9SA

v2 = FULL 6-element phase matrix (alpha1..alpha4, beta1, beta2) projected onto CRTM's
generalized-spherical-function basis (a port of CRTM_Utility Gl2n / Legendre_M), with
the verified normalization: stored coeff(l) = (2l+1)*chi_l, reader multiplies by 0.5
=> CRTM's (2l+1)/2*chi_l (pcoeff(0)=0.5, pcoeff(l=1)/g=1.5). Real 5-temperature optics.

Self-contained (numpy + netCDF4). Runs on the HPC where the data live.

Layout: <root>/**/Scheme{1..4}/T{190,210,230,250,270}/{isca.dat,P11..P44.dat}
  isca.dat : 70 freq x 131 size rows (size inner); cols wl_um,Dmax,vol,projarea,Qext,ssa,g
  Pxx.dat  : line1 = 498 scattering angles (deg); then 70*131 rows x 498 values
             P11 absolute (int P11 dOmega = 4pi); P12,P22,P33,P43,P44 are ratios /P11.
"""
import os, glob, argparse
import numpy as np
from netCDF4 import Dataset, stringtochar

RHO_ICE = 917.0; C_LIGHT = 2.99792458e8; N_SIZE = 131
TEMPS = [190, 210, 230, 250, 270]
SCHEME_NAME = {1: "ConstDensity0.1", 2: "Thompson08", 3: "Heymsfield04", 4: "Brandes07"}
PHASE_FILES = ["P11", "P12", "P22", "P33", "P43", "P44"]


# ---------------------------------------------------------------- GSF basis (CRTM port)
def _fact(n):
    f = 1.0
    for i in range(2, n + 1):
        f *= i
    return f

def legendre_basis(u, Lmax):                       # P^l_{0,0} = ordinary Legendre
    P = np.zeros((Lmax, u.size)); P[0] = 1.0
    if Lmax > 1: P[1] = u
    for l in range(1, Lmax - 1):
        P[l + 1] = ((2*l + 1)*u*P[l] - l*P[l - 1]) / (l + 1)
    return P

def gl2n_basis(MF, n, u, Lmax):                    # CRTM Gl2n: P^l_{MF,n}(u), n=+-2
    G = np.zeros((Lmax, u.size))
    if Lmax < 3: return G
    if MF == 0:
        G[2] = -0.25*np.sqrt(6.0)*(1.0 - u*u)
    else:
        fac = -np.sqrt(_fact(2*MF)/_fact(MF + 2)/_fact(MF - 2)) / (2**MF)
        if n == 2: G[MF] = (1 - u)**(MF/2.0 - 1)*(1 + u)**(MF/2.0 + 1)
        else:      G[MF] = (1 - u)**(MF/2.0 + 1)*(1 + u)**(MF/2.0 - 1)
        G[MF] *= fac
    if MF < 2:
        if Lmax > 3:
            f = 2.0*np.sqrt(((2 + 1)**2 - 4)*((2 + 1)**2 - MF*MF))
            G[3] = (2*3*u - MF*n)*G[2]*(2*2 + 1)/f
        ks = 3
    else:
        if Lmax > MF + 1:
            f = MF*np.sqrt(((MF + 1)**2 - 4)*(2*MF + 1.0))
            G[MF + 1] = (2*MF + 1)*(MF*(MF + 1)*u - n*MF)*G[MF]/f
        ks = MF + 1
    for k in range(ks, Lmax - 1):
        f = k*np.sqrt(((k + 1)**2 - 4)*((k + 1)**2 - MF*MF))
        G[k + 1] = ((2*k + 1)*(k*(k + 1)*u - n*MF)*G[k]
                    - (k + 1)*np.sqrt((k*k - 4)*(k*k - MF*MF))*G[k - 1]) / f
    return G


class Projector:
    """Precompute GSF basis on the scattering-angle grid; project F-matrix elements."""
    def __init__(self, theta_deg, Lmax):
        x = np.cos(np.deg2rad(theta_deg)); o = np.argsort(x)
        self.x = x[o]; self.o = o; self.L = Lmax
        self.P00 = legendre_basis(self.x, Lmax)             # F11, F44
        self.P02 = gl2n_basis(0, 2, self.x, Lmax)           # F12, F34
        self.P22 = gl2n_basis(2, 2, self.x, Lmax)           # F22+F33
        self.P2m = gl2n_basis(2, -2, self.x, Lmax)          # F22-F33

    def _coef(self, F, G):                                   # stored coeff = (2l+1)*(1/2) int F G dx
        Fs = F[self.o]
        return (2*np.arange(self.L) + 1) * 0.5 * np.trapz(Fs[None, :]*G, self.x, axis=1)

    def project(self, F11, F12, F22, F33, F34, F44):
        a1 = self._coef(F11, self.P00)
        a4 = self._coef(F44, self.P00)
        s  = self._coef(F22 + F33, self.P22)
        d  = self._coef(F22 - F33, self.P2m)
        a2 = (s + d)/2.0; a3 = (s - d)/2.0
        b1 = self._coef(F12, self.P02)
        b2 = self._coef(F34, self.P02)
        return np.array([a1, a2, a3, a4, b1, b2])           # (6, Lmax), CRTM element order


def self_test(Lmax=64, nang=498):
    """Orthonormality of the GSF + round-trip reconstruction of a Rayleigh matrix."""
    th = np.linspace(0, 180, nang); pr = Projector(th, Lmax); x = pr.x
    # (a) orthonormality: int (P^l_{mn})^2 dx ~ 2/(2l+1)
    ok = True
    for name, G in (("P00", pr.P00), ("P02", pr.P02), ("P22", pr.P22), ("P2m", pr.P2m)):
        for l in (2, 5, 10, 20):
            val = np.trapz(G[l]*G[l], x); ref = 2.0/(2*l + 1)
            if abs(val - ref) > 0.02*ref + 1e-6: ok = False; print(f"  ORTHO FAIL {name} l={l}: {val:.4e} vs {ref:.4e}")
    # (b) Rayleigh (delta=0) analytic coefficients in CRTM convention (2l+1)*chi_l:
    #     a1={1,0,1/2}, a2_2=3, a3_2=0, a4_1=3/2, |b1_2|=sqrt(6)/2, b2=0  (Hovenier/de Rooij).
    #     Checking known low-l coefficients is robust to the (2l+1)-amplified high-l
    #     quadrature noise that the fixed 498-angle grid produces.
    c = np.cos(np.deg2rad(th))
    F11 = 0.75*(1 + c*c); F22 = F11.copy(); F33 = 1.5*c; F44 = F33.copy()
    F12 = -0.75*(1 - c*c); F34 = np.zeros_like(c)
    A = pr.project(F11, F12, F22, F33, F34, F44)
    # magnitudes match literature (a2/b1 sign follows CRTM's Gl2n basis, hence abs())
    checks = [("a1[0]", A[0, 0], 1.0), ("a1[2]", A[0, 2], 0.5), ("|a2[2]|", abs(A[1, 2]), 3.0),
              ("a3[2]", A[2, 2], 0.0), ("a4[1]", A[3, 1], 1.5),
              ("|b1[2]|", abs(A[4, 2]), np.sqrt(6)/2), ("b2[2]", A[5, 2], 0.0)]
    bad = []
    for nm, v, r in checks:
        flag = abs(v - r) > 0.02*abs(r) + 0.02
        print("   %-8s = %8.4f  (ref %8.4f)%s" % (nm, v, r, "  <-- FAIL" if flag else ""))
        if flag: bad.append(nm)
    # convention-internal round-trip, truncated to the real content (avoids high-l noise)
    Lt = 8
    e11 = np.max(np.abs(A[0, :Lt] @ pr.P00[:Lt] - F11[pr.o]))
    e22 = np.max(np.abs((A[1, :Lt] + A[2, :Lt]) @ pr.P22[:Lt] - (F22 + F33)[pr.o]))
    e12 = np.max(np.abs(A[4, :Lt] @ pr.P02[:Lt] - F12[pr.o]))
    print("   round-trip(L<=%d): |dF11|=%.1e |dF22+33|=%.1e |dF12|=%.1e" % (Lt, e11, e22, e12))
    if max(e11, e22, e12) > 1e-2: bad.append("roundtrip")
    if not ok or bad:
        raise RuntimeError("GSF self-test FAILED: " + (",".join(bad) if bad else "orthonormality"))
    print("  self-test PASSED (Rayleigh magnitudes + round-trip + orthonormality)")


# ---------------------------------------------------------------- data reading
def read_node_T(root, node, T):
    """Read one (habit-node, temperature) folder. `node` is the directory name that varies by
    database: 'Scheme{N}' for the Ren snow DB (doi:10.18738/T8/LGJ9SA), 'massR_{R}' for the
    conical-graupel DB (doi:10.18738/T8/DWZXZX). Both share the isca.dat (wl,Dmax,vol,area,
    Qext,ssa,g) + 6 Pxx.dat (498-angle) format. N_SIZE is auto-detected (Ren=131, graupel=70;
    particle size is the inner loop, so it = #leading rows sharing the first wavelength)."""
    hits = glob.glob(os.path.join(root, "**", node, "T%d" % T, "isca.dat"), recursive=True)
    if not hits:
        raise FileNotFoundError("%s/T%d/isca.dat not found under %s" % (node, T, root))
    base = os.path.dirname(hits[0])
    isca = np.loadtxt(os.path.join(base, "isca.dat"))
    wl0 = isca[:, 0]
    nsize = int(np.argmax(wl0 != wl0[0])) or isca.shape[0]      # inner-loop particle-size count
    nf = isca.shape[0] // nsize
    isca = isca.reshape(nf, nsize, 7)
    P = {}
    for nm in PHASE_FILES:
        with open(os.path.join(base, nm + ".dat")) as fh:
            ang = np.array(fh.readline().split(), dtype=float)
        P[nm] = np.loadtxt(os.path.join(base, nm + ".dat"), skiprows=1).reshape(nf, nsize, ang.size)
    return dict(nf=nf, ang=ang, wl=isca[..., 0], Dmax=isca[..., 1], vol=isca[..., 2],
                area=isca[..., 3], qext=isca[..., 4], ssa=isca[..., 5], g=isca[..., 6], P=P)


def build(root, node, out, mu_val=0.0, n_dm=30, L_max=64, tol=1e-3, dm_max_um=10000.0,
          habit_id=4, habit_name=None, prov="Snow", source=""):
    print("self-test:"); self_test(L_max)
    s0 = read_node_T(root, node, TEMPS[0])
    NF = s0["nf"]; nang = s0["ang"].size; theta = s0["ang"]; NT = len(TEMPS)
    freqs = np.round(C_LIGHT/(s0["wl"][:, 0]*1e-6)/1e9, 3)
    pr = Projector(theta, L_max)

    D = s0["Dmax"][0, :]*1e-6
    mass = s0["vol"][0, :]*1e-18*RHO_ICE
    order = np.argsort(D); D = D[order]; mass = mass[order]
    Dg = np.geomspace(D.min(), D.max(), 400); massg = np.interp(Dg, D, mass)
    Dm_grid = np.geomspace(2*D.min(), min(0.7*D.max(), dm_max_um*1e-6), n_dm)

    KE = np.zeros((1, n_dm, 1, NT, NF)); KA = np.zeros_like(KE); KB = np.zeros_like(KE)
    GG = np.zeros_like(KE); NEFF = np.zeros((1, n_dm, 1, NT, NF), dtype="i4")
    PCO = np.zeros((1, n_dm, 1, NT, NF, L_max, 6))

    for it, T in enumerate(TEMPS):
        s = read_node_T(root, node, T)
        sext = (s["qext"]*s["area"]*1e-12)[:, order]
        ssca = sext*s["ssa"][:, order]; sabs = sext - ssca
        # per-particle F-matrix elements (P12.. are ratios; F34 = -F43)
        p = s["P"]; p11 = p["P11"][:, order, :]
        F = dict(F11=p11, F12=p["P12"][:, order, :]*p11, F22=p["P22"][:, order, :]*p11,
                 F33=p["P33"][:, order, :]*p11, F34=-p["P43"][:, order, :]*p11,
                 F44=p["P44"][:, order, :]*p11)
        i180 = int(np.argmax(theta)); sbk = ssca*p11[:, :, i180]
        for jf in range(NF):
            se = np.interp(Dg, D, sext[jf]); sa = np.interp(Dg, D, sabs[jf])
            sc = np.interp(Dg, D, ssca[jf]); sb = np.interp(Dg, D, sbk[jf])
            Fg = {k: np.array([np.interp(Dg, D, F[k][jf, :, ia]) for ia in range(nang)]).T
                  for k in F}                                # each (len(Dg), nang)
            for idm, Dm in enumerate(Dm_grid):
                Nd = (Dg/Dm)**mu_val*np.exp(-(4.0 + mu_val)*Dg/Dm)
                M = np.trapz(Nd*massg, Dg)
                KE[0, idm, 0, it, jf] = np.trapz(Nd*se, Dg)/M
                KA[0, idm, 0, it, jf] = np.trapz(Nd*sa, Dg)/M
                KB[0, idm, 0, it, jf] = np.trapz(Nd*sb, Dg)/M
                w = Nd*sc; wn = max(np.trapz(w, Dg), 1e-30)
                Fb = {k: np.trapz(w[:, None]*Fg[k], Dg, axis=0)/wn for k in Fg}   # bulk F(theta)
                A = pr.project(Fb["F11"], Fb["F12"], Fb["F22"], Fb["F33"], Fb["F34"], Fb["F44"])
                GG[0, idm, 0, it, jf] = A[0, 1]/3.0          # g = chi_1 = coeff(1)/(2*1+1)
                # truncation from real content (bare chi over all elements)
                lvec = (2*np.arange(L_max) + 1.0)
                chi = np.abs(A)/lvec[None, :]
                sig = np.where(chi.max(axis=0) >= tol)[0]
                NEFF[0, idm, 0, it, jf] = int(sig[-1] + 1) if sig.size else 1
                PCO[0, idm, 0, it, jf, :, :] = A.T            # (L_max, 6)

    write_netcdf(out, freqs, Dm_grid*1e6, np.array([mu_val]), np.array(TEMPS, float),
                 KE, KA, KB, GG, NEFF, PCO, L_max, prov=prov, source=source,
                 habit_id=habit_id, habit_name=habit_name)
    print("wrote %s  (material=%s, habit_id=%d, %d freq %g-%g GHz, %d T, %d Dm, n_Legendre_Eff %d-%d, 6 phase elements)"
          % (out, prov, habit_id, NF, freqs.min(), freqs.max(), NT, n_dm, NEFF.min(), NEFF.max()))


def write_netcdf(path, freq, dm, mu, temp, ke, ka, kb, g, neff, pcoeff, L_max, prov="Snow",
                 source="", habit_id=4, habit_name=None):
    nc = Dataset(path, "w", format="NETCDF4")
    for n, v in [("n_Frequency", len(freq)), ("n_Dm", len(dm)), ("n_Mu", len(mu)),
                 ("n_Temperature", len(temp)), ("n_Habit", 1), ("n_Legendre", L_max),
                 ("n_Phase_Elements", 6), ("nchar", 32)]:
        nc.createDimension(n, v)
    nc.Scheme = "CRTM-Exp"; nc.Release = np.int32(1); nc.Version = np.int32(2)
    nc.PSD = "normalized_gamma"; nc.Orientation = "random"
    nc.Source = source or ("Ren et al. 2022 snow II-TM/IGOM, doi:10.18738/T8/LGJ9SA, %s" % prov)
    nc.Note = "v2: full 6-element GSF (a1,a2,a3,a4,b1,b2); coeff=(2l+1)*chi_l (reader x0.5); real 5-T"

    def cv(nm, dims, data, **k):
        x = nc.createVariable(nm, "f8", dims, **k); x[:] = data; return x
    cv("Frequency", ("n_Frequency",), freq).units = "GHz"
    cv("Dm", ("n_Dm",), dm).units = "microns"
    cv("Mu", ("n_Mu",), mu); cv("Temperature", ("n_Temperature",), temp).units = "K"
    if habit_name is None:
        habit_name = prov
    nc.createVariable("Habit_Id", "i4", ("n_Habit",))[:] = [habit_id]  # CRTM cloud-type integer
    nc.createVariable("Habit_Phase", "i4", ("n_Habit",))[:] = [1]      # frozen
    cv("mD_a", ("n_Habit",), [0.0]); cv("mD_b", ("n_Habit",), [0.0])
    nc.createVariable("Habit_Name", "S1", ("n_Habit", "nchar"))[:] = \
        stringtochar(np.array([habit_name.ljust(32)], dtype="S32"))
    bd = ("n_Habit", "n_Dm", "n_Mu", "n_Temperature", "n_Frequency")
    for nm, arr in (("ke", ke), ("ka", ka), ("kb", kb), ("g", g)):
        cv(nm, bd, arr, zlib=True, complevel=4)
    nc.createVariable("n_Legendre_Eff", "i4", bd, zlib=True, complevel=4)[:] = neff
    cv("pcoeff", ("n_Habit","n_Dm","n_Mu","n_Temperature","n_Frequency","n_Legendre","n_Phase_Elements"),
       pcoeff, zlib=True, complevel=4)
    nc.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--scheme", type=int, default=2, choices=[1, 2, 3, 4],
                    help="Ren snow m-D scheme (ignored if --massr is given)")
    ap.add_argument("--massr", default=None,
                    help="conical-graupel mass ratio, e.g. 0.500 (doi:10.18738/T8/DWZXZX); "
                         "selects folder massR_<massr> and switches DB from Ren-snow to graupel")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--mu", type=float, default=0.0)
    ap.add_argument("--n-dm", type=int, default=30)
    ap.add_argument("--habit-id", type=int, default=None,
                    help="CRTM cloud-type integer (2=ICE,4=SNOW,5=GRAUPEL); default 4 snow / 5 graupel")
    ap.add_argument("--habit-name", default=None, help="provenance string")
    ap.add_argument("--selftest-only", action="store_true")
    a = ap.parse_args()
    if a.selftest_only:
        self_test(); raise SystemExit
    if a.massr is not None:
        node   = "massR_%s" % a.massr
        prov   = "ConicalGraupel_massR%s" % a.massr
        source = "Tang et al. 2017 conical graupel, doi:10.18738/T8/DWZXZX, massR %s" % a.massr
        hid    = a.habit_id if a.habit_id is not None else 5
        out    = a.output or "CloudCoeff_Exp_Graupel_massR%s.nc" % a.massr
    else:
        node   = "Scheme%d" % a.scheme
        prov   = "Snow_%s" % SCHEME_NAME[a.scheme]
        source = "Ren et al. 2022 snow II-TM/IGOM, doi:10.18738/T8/LGJ9SA, %s" % SCHEME_NAME[a.scheme]
        hid    = a.habit_id if a.habit_id is not None else 4
        out    = a.output or "CloudCoeff_Exp_RenSnow_%s.nc" % SCHEME_NAME[a.scheme]
    build(a.root, node, out, mu_val=a.mu, n_dm=a.n_dm,
          habit_id=hid, habit_name=a.habit_name, prov=prov, source=source)
