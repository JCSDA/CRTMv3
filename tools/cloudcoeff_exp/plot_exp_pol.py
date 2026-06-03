#!/usr/bin/env python3
"""Polarized phase-matrix view of an experimental CloudCoeff (v2, 6-element).
Reconstructs the bulk scattering matrix F_ij(theta) from the stored GSF coefficients
(coeff_l = (2l+1)*chi_l) using the same basis as the builder, and plots the phase
function and the polarized ratios. Usage: plot_exp_pol.py LUT.nc OUT.png [Dm_um]"""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from ren_build import legendre_basis, gl2n_basis      # exact builder basis

lut = sys.argv[1]; out = sys.argv[2] if len(sys.argv) > 2 else "pol.png"
dm_want = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0

d = Dataset(lut); d.set_auto_mask(False)        # plain ndarrays (masked @ trips a numpy bug)
freq = d["Frequency"][:]; T = d["Temperature"][:]; dm = d["Dm"][:]
pc = d["pcoeff"][:]; ne = d["n_Legendre_Eff"][:]; L = pc.shape[5]
idm = int(np.argmin(np.abs(dm - dm_want))); itc = len(T)//2

th = np.linspace(0, 180, 721); u = np.cos(np.deg2rad(th))
P00 = legendre_basis(u, L); P02 = gl2n_basis(0, 2, u, L)
P22 = gl2n_basis(2, 2, u, L); P2m = gl2n_basis(2, -2, u, L)

def recon(c, Le):                                       # c=(L,6) stored coeffs; F = sum coeff_l * basis_l
    Le = max(int(Le), 2)
    F11 = c[:Le, 0] @ P00[:Le]; F44 = c[:Le, 3] @ P00[:Le]
    s = (c[:Le, 1] + c[:Le, 2]) @ P22[:Le]; dd = (c[:Le, 1] - c[:Le, 2]) @ P2m[:Le]
    F22 = (s + dd)/2; F33 = (s - dd)/2
    F12 = c[:Le, 4] @ P02[:Le]; F34 = c[:Le, 5] @ P02[:Le]
    return F11, F12, F22, F33, F34, F44

fig, ax = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("Ren-Thompson snow: reconstructed bulk scattering matrix  "
             f"(Dm={dm[idm]:.0f} um, T={T[itc]:.0f} K)", fontsize=13)
for f_t in (37, 89, 166, 325, 874):
    jf = int(np.argmin(np.abs(freq - f_t)))
    c = pc[0, idm, 0, itc, jf, :, :]
    F11, F12, F22, F33, F34, F44 = recon(c, ne[0, idm, 0, itc, jf])
    lab = f"{round(freq[jf])} GHz"
    ax[0, 0].semilogy(th, np.maximum(F11, 1e-3), label=lab)
    ax[0, 1].plot(th, -100.0*F12/np.maximum(F11, 1e-9), label=lab)   # degree of linear pol (%)
    ax[1, 0].plot(th, F22/np.maximum(F11, 1e-9), label=lab)
    ax[1, 1].plot(th, F44/np.maximum(F11, 1e-9), label=lab)
ax[0, 0].set(xlabel="scattering angle (deg)", ylabel="F11 (norm.)", title="Phase function F11", xlim=(0, 180))
ax[0, 1].set(xlabel="scattering angle (deg)", ylabel="-F12/F11  (%)", title="Degree of linear polarization", xlim=(0, 180))
ax[1, 0].set(xlabel="scattering angle (deg)", ylabel="F22/F11", title="F22/F11 (depolarization)", xlim=(0, 180))
ax[1, 1].set(xlabel="scattering angle (deg)", ylabel="F44/F11", title="F44/F11 (circular)", xlim=(0, 180))
for a in ax.ravel(): a.grid(True, which="both", alpha=.3); a.legend(fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(out, dpi=110)
print("wrote", out)
