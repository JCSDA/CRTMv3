#!/usr/bin/env python3
"""
dda_parse.py -- robust parsers for the in-house DDSCAT archive.

Two readers:
  parse_geom(geom_param.dat) -> physical mass / Dmax / area / density per shape
  parse_avg(w000r000.avg)    -> orientation-averaged Qext/Qabs/Qsca/g/Qbk and the
                                angular S_11(theta) phase function for one frequency.

Notes:
  - Physical mass = (DDA target Volume in um^3) * RHO_ICE, converted to kg.
    (The shape_mass_dimension*.csv masses assume a fixed 50 um dipole spacing and
     are WRONG for shapes with a different spacing -- always use geom_param.dat.)
  - AEFF/WAVE in the .avg are in physical units (um); frequency = c/lambda.
  - This archive's .avg Mueller columns are S_11 S_12 S_21 S_22 S_31 S_41 (set by
    the ddscat.par NSMELTS line); S_33/S_34/S_44 are NOT present -> only the scalar
    phase function P11 (and P12 via S_12) can be built from it. Sufficient for the
    v1 scalar LUT (alpha1); full 6-element GSF needs re-runs with more S_ij.
"""

import math

RHO_ICE = 917.0          # kg/m^3
C_LIGHT = 2.99792458e8   # m/s


def parse_geom(path):
    """Parse a geom_param.dat -> dict with physical mass_kg, dmax_um, proj_area_um2,
    vol_um3, aspect_ratio, snow_density_gcc."""
    vals = {}
    lines = [l.rstrip("\n") for l in open(path)]
    def numbers(s):
        out = []
        for tok in s.replace(",", " ").split():
            try: out.append(float(tok))
            except ValueError: pass
        return out
    for i, line in enumerate(lines):
        t = line.strip()
        if t.startswith("#"):
            key = t.lstrip("#").strip().lower()
            # the data line(s) follow a '#' label line
            # first Volume/Area block is the true geometry; a later identical
            # sub-header belongs to the "Sphere-equivalent Radii" block -> skip it.
            if key.startswith("volume") and "sfc" in key and "vol_um3" not in vals:
                v = numbers(lines[i+1]);  vals.update(vol_um3=v[0], sfc_area_um2=v[1], proj_area_um2=v[2])
            elif key.startswith("axis lengths"):                    # circumscribing ellipsoid axes
                v = numbers(lines[i+1]); vals["axes_um"] = v
            elif key.startswith("maximum diameter"):
                v = numbers(lines[i+1]); vals["dmax_um"] = v[0]
    if "vol_um3" in vals:
        vals["mass_kg"] = vals["vol_um3"] * 1.0e-18 * RHO_ICE       # um^3 -> m^3 -> kg
    if "axes_um" in vals and len(vals["axes_um"]) == 3:
        a = sorted(vals["axes_um"]); vals["aspect_ratio"] = a[0] / a[2] if a[2] else 0.0
    return vals


def parse_avg(path):
    """Parse a DDSCAT .avg -> dict(aeff_um, wave_um, freq_ghz, qext, qabs, qsca, g,
    cos2, qbk, theta[deg], s11[], s12[])."""
    d = {}
    theta, s11, s12 = [], [], []
    in_mueller = False
    col = {}
    for line in open(path):
        s = line.strip()
        if "AEFF=" in line and "aeff_um" not in d:
            d["aeff_um"] = float(line.split("AEFF=")[1].split("=")[0].split()[0])
        elif "WAVE=" in line and "wave_um" not in d:
            d["wave_um"] = float(line.split("WAVE=")[1].split("=")[0].split()[0])
        elif s.startswith("mean:") and "qext" not in d:
            v = [float(x) for x in s.split()[1:]]
            d.update(qext=v[0], qabs=v[1], qsca=v[2], g=v[3], cos2=v[4], qbk=v[5],
                     qpha=v[6] if len(v) > 6 else float("nan"))
        elif "theta" in s and "S_11" in s:           # Mueller table header
            toks = s.split()
            for j, name in enumerate(toks):
                col[name] = j
            in_mueller = True
        elif in_mueller and s and s[0].isdigit():
            toks = s.split()
            try:
                theta.append(float(toks[col["theta"]]))
                s11.append(float(toks[col["S_11"]]))
                if "S_12" in col: s12.append(float(toks[col["S_12"]]))
            except (ValueError, IndexError):
                pass
    if "wave_um" in d and d["wave_um"] > 0:
        d["freq_ghz"] = C_LIGHT / (d["wave_um"] * 1.0e-6) / 1.0e9
    d["theta"] = theta; d["s11"] = s11; d["s12"] = s12
    return d


if __name__ == "__main__":
    import sys
    base = "/mnt/d/WSL_temp_storage/full_db/aggregate/d165/cmb/d050/01_0042_002"
    g = parse_geom(base + "/geom_param.dat")
    a = parse_avg(base + "/ddscat07/089GHz/w000r000.avg")
    print("GEOM: vol_um3=%.3e mass_kg=%.3e dmax_um=%.1f proj_area_um2=%.3e aspect=%.3f"
          % (g["vol_um3"], g["mass_kg"], g["dmax_um"], g["proj_area_um2"], g["aspect_ratio"]))
    print("AVG : freq=%.1f GHz aeff=%.2f um Qext=%.4e Qabs=%.4e Qsca=%.4e g=%.4f Qbk=%.4e"
          % (a["freq_ghz"], a["aeff_um"], a["qext"], a["qabs"], a["qsca"], a["g"], a["qbk"]))
    print("AVG : n_theta=%d  theta[0,-1]=%.1f,%.1f  S_11[0]=%.4e  S_11[180deg]=%.4e"
          % (len(a["theta"]), a["theta"][0], a["theta"][-1], a["s11"][0], a["s11"][-1]))
