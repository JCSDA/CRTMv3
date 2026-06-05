#!/usr/bin/env bash
# Build and run the experimental Ren-snow CloudCoeff validation driver for MWR_AWS.
# Compiles test/mains/exp_aws_scatter.f90 against the in-tree libcrtm.so and runs it
# from test/ (where ./testinput/ holds the all-netCDF coeff set + the CRTM-Exp LUT).
#
# Prereqs: a built tree under build/ (cmake .. && make), and the assembled LUT linked
# into test/testinput/ (see tools/cloudcoeff_exp/ren_assemble.py).
set -euo pipefail
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
MODDIR=$(dirname "$(find "$ROOT/build/module" -name crtm_module.mod | head -1)")
INCDIR=$ROOT/test/mains/unit/Unit_Test          # Load_ECMWF84_Atm_Data.inc lives here
BIN=$ROOT/build/bin/exp_aws_scatter

echo "=== compiling exp_aws_scatter (modules: $MODDIR) ==="
# libcrtm.so is self-contained (carries its own netCDF/HDF5 RUNPATH) and the driver
# calls only CRTM symbols, so no explicit -lnetcdff/-lnetcdf is needed.
gfortran -D_REAL8_ -ffree-line-length-none -O2 -fopenmp \
  -I "$MODDIR" -I "$INCDIR" \
  "$ROOT/test/mains/exp_aws_scatter.f90" \
  -o "$BIN" \
  "$ROOT/build/lib/libcrtm.so" \
  -Wl,-rpath,"$ROOT/build/lib"
echo "compiled -> $BIN"

echo "=== running from $ROOT/test ==="
cd "$ROOT/test"
"$BIN"
