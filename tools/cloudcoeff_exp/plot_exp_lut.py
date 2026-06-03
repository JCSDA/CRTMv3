#!/usr/bin/env python3
"""Overview plots of an experimental CloudCoeff (v1) LUT. Usage: plot_exp_lut.py LUT.nc OUT.png"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
from netCDF4 import Dataset

lut = sys.argv[1] if len(sys.argv) > 1 else "CloudCoeff_Exp_aggregate_v1.nc"
out = sys.argv[2] if len(sys.argv) > 2 else "cloudcoeff_exp_overview.png"

d = Dataset(lut)
freq = d["Frequency"][:]; dm = d["Dm"][:]; T = d["Temperature"][:]
ke = d["ke"][:]; ka = d["ka"][:]; g = d["g"][:]; ne = d["n_Legendre_Eff"][:]; pc = d["pcoeff"][:]
hname = b"".join(d["Habit_Name"][0].astype("S1")).decode().strip()
NH, NDm, NMu, NT, NF = ke.shape
itc = NT // 2                                   # mid temperature slice
w = 1.0 - np.divide(ka, ke, out=np.zeros_like(ka), where=ke > 0)

# Dm indices to show (small / mid / large)
dmsel = [1, NDm // 2, NDm - 2]
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle(f"Experimental CloudCoeff v1 — habit '{hname}'  ({NDm} Dm × {NF} freq × {NT} T)", fontsize=14)

# (1) mass extinction vs frequency
for i in dmsel:
    ax[0,0].loglog(freq, ke[0,i,0,itc,:], "-o", ms=3, label=f"Dm={dm[i]:.0f} µm")
ax[0,0].set(xlabel="Frequency (GHz)", ylabel="ke (m²/kg)", title="Mass extinction"); ax[0,0].grid(True, which="both", alpha=.3); ax[0,0].legend(fontsize=8)

# (2) single-scatter albedo vs frequency
for i in dmsel:
    ax[0,1].semilogx(freq, w[0,i,0,itc,:], "-o", ms=3, label=f"Dm={dm[i]:.0f} µm")
ax[0,1].set(xlabel="Frequency (GHz)", ylabel="ω", title="Single-scatter albedo", ylim=(0,1)); ax[0,1].grid(True, which="both", alpha=.3); ax[0,1].legend(fontsize=8)

# (3) asymmetry vs frequency
for i in dmsel:
    ax[0,2].semilogx(freq, g[0,i,0,itc,:], "-o", ms=3, label=f"Dm={dm[i]:.0f} µm")
ax[0,2].set(xlabel="Frequency (GHz)", ylabel="g", title="Asymmetry parameter", ylim=(0,1)); ax[0,2].grid(True, which="both", alpha=.3); ax[0,2].legend(fontsize=8)

# (4) reconstructed bulk phase function P11(theta) at a mid Dm, several freqs
th = np.linspace(0, 180, 361); mu = np.cos(np.deg2rad(th)); idm = NDm // 2
for f_t in (89, 166, 874):
    jf = int(np.argmin(np.abs(freq - f_t)))
    L = int(ne[0,idm,0,itc,jf]); chi = pc[0,idm,0,itc,jf,:L,0]
    P = legval(mu, (2*np.arange(L)+1)*chi)
    ax[1,0].semilogy(th, np.maximum(P,1e-3), label=f"{round(freq[jf])} GHz (L={L})")
ax[1,0].set(xlabel="scattering angle θ (deg)", ylabel="P₁₁ (norm.)", title=f"Bulk phase function  (Dm={dm[idm]:.0f} µm)", xlim=(0,180)); ax[1,0].grid(True, which="both", alpha=.3); ax[1,0].legend(fontsize=8)

# (5) effective Legendre truncation order vs frequency
for i in dmsel:
    ax[1,1].semilogx(freq, ne[0,i,0,itc,:], "-o", ms=3, label=f"Dm={dm[i]:.0f} µm")
ax[1,1].axhline(8, color="k", ls="--", lw=1, label="legacy auto (~8)")
ax[1,1].set(xlabel="Frequency (GHz)", ylabel="n_Legendre_Eff", title="Phase-fn truncation order (LUT-driven)"); ax[1,1].grid(True, which="both", alpha=.3); ax[1,1].legend(fontsize=8)

# (6) temperature sensitivity of absorption (Mätzler) at a few freqs, mid Dm
for f_t in (89, 166, 183):
    jf = int(np.argmin(np.abs(freq - f_t)))
    ax[1,2].plot(T, ka[0,idm,0,:,jf]/ka[0,idm,0,itc,jf], "-o", ms=4, label=f"{round(freq[jf])} GHz")
ax[1,2].set(xlabel="Temperature (K)", ylabel=f"ka(T) / ka({int(round(T[itc]))} K)", title="Absorption temperature sensitivity"); ax[1,2].grid(True, alpha=.3); ax[1,2].legend(fontsize=8)

fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(out, dpi=110)
print(f"wrote {out}")
