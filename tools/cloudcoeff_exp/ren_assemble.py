#!/usr/bin/env python3
"""
ren_assemble.py -- stack per-scheme single-habit experimental CloudCoeff files
(produced by ren_build.py) into ONE multi-habit netCDF that the CRTM 'CRTM-Exp'
reader loads in a single call. The Fortran reader (Get_Cloud_Opt_MW_Exp) selects a
habit slice by matching Habit_Id to the CRTM cloud-type integer:
    2 = ICE_CLOUD,  4 = SNOW_CLOUD,  5 = GRAUPEL_CLOUD   (CRTM_Cloud_Define.f90)

All inputs MUST share identical Frequency / Dm / Mu / Temperature grids and the same
n_Legendre / n_Phase_Elements (they do when built by ren_build.py with the same
--n-dm / L_max on the same Ren size set -- the Dmax size grid is fixed across schemes).

The same source file may back several habit IDs (e.g. ICE and SNOW both from the
Thompson08 Scheme2 optics -- below the column/aggregate threshold the Ren particle is a
single hexagonal column, a sound small-Dm cloud-ice proxy; the habits then differ only by
the host's Water_Content and Dm).

Usage:
  ren_assemble.py -o OUT.nc  ID:NAME:FILE  [ID:NAME:FILE ...]
Example (ice/snow/graupel):
  ren_assemble.py -o CloudCoeff_Exp_RenSnow_isg.nc \
    2:Ice_Thompson08_hexcol:CloudCoeff_Exp_RenSnow_snow_sch2.nc \
    4:Snow_Thompson08:CloudCoeff_Exp_RenSnow_snow_sch2.nc \
    5:Graupel_ConstDensity0.1:CloudCoeff_Exp_RenSnow_graupel_sch1.nc
"""
import argparse
import numpy as np
from netCDF4 import Dataset, stringtochar

AXIS_VARS = ["Frequency", "Dm", "Mu", "Temperature"]
BULK = ["ke", "ka", "kb", "g"]            # file dims (Habit, Dm, Mu, Temperature, Frequency)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("habits", nargs="+", metavar="ID:NAME:FILE",
                    help="one per output habit slice (same FILE may repeat)")
    a = ap.parse_args()

    specs = []
    for h in a.habits:
        sid, name, fn = h.split(":", 2)
        specs.append((int(sid), name, fn))

    ref = Dataset(specs[0][2])
    axes = {v: ref[v][:].astype(float) for v in AXIS_VARS}
    L = ref.dimensions["n_Legendre"].size
    NP = ref.dimensions["n_Phase_Elements"].size
    src_attr = getattr(ref, "Source", "")
    ref.close()

    ke = []; ka = []; kb = []; g = []; neff = []; pco = []
    ids = []; names = []; phase = []; mda = []; mdb = []
    for hid, name, fn in specs:
        d = Dataset(fn)
        for v in AXIS_VARS:
            if not np.allclose(d[v][:].astype(float), axes[v]):
                raise SystemExit("FATAL: %s axis %s mismatches reference %s" % (fn, v, specs[0][2]))
        if d.dimensions["n_Legendre"].size != L or d.dimensions["n_Phase_Elements"].size != NP:
            raise SystemExit("FATAL: %s Legendre/Phase dims mismatch" % fn)
        ke.append(d["ke"][0]); ka.append(d["ka"][0]); kb.append(d["kb"][0]); g.append(d["g"][0])
        neff.append(d["n_Legendre_Eff"][0]); pco.append(d["pcoeff"][0])
        ids.append(hid); names.append(name); phase.append(int(d["Habit_Phase"][0]))
        mda.append(float(d["mD_a"][0])); mdb.append(float(d["mD_b"][0]))
        print("  habit %2d <- %-32s %s" % (hid, name, fn))
        d.close()

    NH = len(specs)
    ke = np.stack(ke); ka = np.stack(ka); kb = np.stack(kb); g = np.stack(g)
    neff = np.stack(neff).astype("i4"); pco = np.stack(pco)

    nc = Dataset(a.output, "w", format="NETCDF4")
    for n, v in [("n_Frequency", axes["Frequency"].size), ("n_Dm", axes["Dm"].size),
                 ("n_Mu", axes["Mu"].size), ("n_Temperature", axes["Temperature"].size),
                 ("n_Habit", NH), ("n_Legendre", L), ("n_Phase_Elements", NP), ("nchar", 32)]:
        nc.createDimension(n, v)
    nc.Scheme = "CRTM-Exp"; nc.Release = np.int32(1); nc.Version = np.int32(2)
    nc.PSD = "normalized_gamma"; nc.Orientation = "random"
    nc.Source = src_attr
    nc.Note = ("multi-habit assembly via ren_assemble.py; habits: "
               + ", ".join("%d=%s" % (i, n) for i, n in zip(ids, names)))

    def cv(nm, dims, data, **k):
        x = nc.createVariable(nm, "f8", dims, **k); x[:] = data; return x

    cv("Frequency", ("n_Frequency",), axes["Frequency"]).units = "GHz"
    cv("Dm", ("n_Dm",), axes["Dm"]).units = "microns"
    cv("Mu", ("n_Mu",), axes["Mu"])
    cv("Temperature", ("n_Temperature",), axes["Temperature"]).units = "K"
    nc.createVariable("Habit_Id", "i4", ("n_Habit",))[:] = np.array(ids, "i4")
    nc.createVariable("Habit_Phase", "i4", ("n_Habit",))[:] = np.array(phase, "i4")
    cv("mD_a", ("n_Habit",), np.array(mda)); cv("mD_b", ("n_Habit",), np.array(mdb))
    nc.createVariable("Habit_Name", "S1", ("n_Habit", "nchar"))[:] = \
        stringtochar(np.array([n.ljust(32) for n in names], dtype="S32"))

    bd = ("n_Habit", "n_Dm", "n_Mu", "n_Temperature", "n_Frequency")
    for nm, arr in (("ke", ke), ("ka", ka), ("kb", kb), ("g", g)):
        cv(nm, bd, arr, zlib=True, complevel=4)
    nc.createVariable("n_Legendre_Eff", "i4", bd, zlib=True, complevel=4)[:] = neff
    cv("pcoeff", ("n_Habit", "n_Dm", "n_Mu", "n_Temperature", "n_Frequency",
                  "n_Legendre", "n_Phase_Elements"), pco, zlib=True, complevel=4)
    nc.close()
    print("wrote %s  (n_Habit=%d, ids=%s)" % (a.output, NH, ids))


if __name__ == "__main__":
    main()
