#!/usr/bin/env python3
"""Extract AWS-1 observation/geometry scenes for CRTM obs-space smoke tests.

The output CSV is intentionally model-state explicit.  AWS-1 L1B supplies
observed brightness temperatures and observation geometry; atmospheric and
surface background fields must come from a collocation pipeline.  The scalar
SST/U10/SSS arguments here are only a smoke-test convenience.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract valid AWS-1 ocean scenes to a CRTM comparison CSV."
    )
    parser.add_argument("aws_l1b", type=Path, help="AWS-1 L1B NetCDF file")
    parser.add_argument("output_csv", type=Path, help="Output scene CSV")
    parser.add_argument("--max-scenes", type=int, default=128)
    parser.add_argument("--scan-stride", type=int, default=25)
    parser.add_argument("--fov-stride", type=int, default=2)
    parser.add_argument("--geo-group", type=int, default=0)
    parser.add_argument("--ocean-only", action="store_true")
    parser.add_argument("--lat-min", type=float, default=-90.0)
    parser.add_argument("--lat-max", type=float, default=90.0)
    parser.add_argument("--sst", type=float, default=285.0, help="Smoke-test SST [K]")
    parser.add_argument("--u10", type=float, default=5.0, help="Smoke-test 10 m wind [m/s]")
    parser.add_argument("--sss", type=float, default=35.0, help="Smoke-test salinity [psu]")
    parser.add_argument(
        "--bt-flag-valid",
        type=int,
        default=0,
        help="Expected brightness-temperature flag value for this product",
    )
    parser.add_argument(
        "--surface-flag-valid",
        type=int,
        default=0,
        help="Expected surface flag value for this product",
    )
    return parser.parse_args()


def scalar(value, fill=np.nan):
    return np.ma.filled(value, fill).item()


def vector(value, fill=np.nan) -> np.ndarray:
    return np.asarray(np.ma.filled(value, fill), dtype=float)


def valid_int(value, expected: int) -> bool:
    try:
        return int(scalar(value, fill=-9999)) == expected
    except (TypeError, ValueError):
        return False


def valid_all_int(values, expected: int) -> bool:
    data = np.asarray(np.ma.filled(values, -9999), dtype=int)
    return bool(np.all(data == expected))


def main() -> None:
    args = parse_args()
    if args.max_scenes <= 0:
        raise SystemExit("--max-scenes must be positive")
    if args.scan_stride <= 0 or args.fov_stride <= 0:
        raise SystemExit("--scan-stride and --fov-stride must be positive")
    if args.lat_min > args.lat_max:
        raise SystemExit("--lat-min cannot be greater than --lat-max")

    with Dataset(args.aws_l1b) as ds:
        data = ds.groups["data"]
        nav = data.groups["navigation"]
        cal = data.groups["calibration"]
        pinfo = data.groups["processing_information"]
        qual = ds.groups["quality"]

        tb = cal.variables["aws_toa_brightness_temperature"]
        lat = nav.variables["aws_lat"]
        lon = nav.variables["aws_lon"]
        zenith = nav.variables["aws_satellite_zenith_angle"]
        azimuth = nav.variables["aws_satellite_azimuth_angle"]
        surface_type = nav.variables["aws_surface_type"]

        bt_flag = pinfo.variables["aws_brightnesstemp_flag"]
        position_flag = pinfo.variables["aws_position_flag_earthview"]
        navigation_status = pinfo.variables["aws_navigation_status"]
        surface_flag = pinfo.variables["aws_surface_flag"]
        channel_quality = pinfo.variables["aws_channel_quality_flag"]
        l1b_quality = qual.variables["L1B_quality_flag"]

        n_scans, n_fovs, n_channels = tb.shape
        if args.geo_group < 0 or args.geo_group >= lat.shape[2]:
            raise SystemExit(f"--geo-group must be in [0,{lat.shape[2] - 1}]")

        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        header = [
            "scene_id",
            "scan",
            "fov",
            "lat",
            "lon",
            "scan_angle",
            "zenith",
            "azimuth",
            "sst",
            "u10",
            "sss",
        ] + [f"obs_tb_{i}" for i in range(1, n_channels + 1)]

        scene_id = 0
        mid_fov = 0.5 * (n_fovs - 1)
        with args.output_csv.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)

            for scan in range(0, n_scans, args.scan_stride):
                if not valid_int(l1b_quality[scan], 1):
                    continue
                if not valid_int(navigation_status[scan], 0):
                    continue
                if not valid_int(surface_flag[scan, args.geo_group], args.surface_flag_valid):
                    continue
                if not valid_all_int(channel_quality[scan, :], 1):
                    continue

                for fov in range(0, n_fovs, args.fov_stride):
                    if not valid_int(position_flag[scan, fov], 1):
                        continue
                    if args.ocean_only and not valid_int(surface_type[scan, fov, args.geo_group], 0):
                        continue
                    if not valid_all_int(bt_flag[scan, fov, :], args.bt_flag_valid):
                        continue

                    obs_tb = vector(tb[scan, fov, :])
                    scene_lat = float(scalar(lat[scan, fov, args.geo_group]))
                    scene_lon = float(scalar(lon[scan, fov, args.geo_group]))
                    scene_zenith = float(scalar(zenith[scan, fov, args.geo_group]))
                    scene_azimuth = float(scalar(azimuth[scan, fov, args.geo_group]))
                    if not np.all(np.isfinite([scene_lat, scene_lon, scene_zenith, scene_azimuth])):
                        continue
                    if scene_lat < args.lat_min or scene_lat > args.lat_max:
                        continue
                    if not np.all(np.isfinite(obs_tb)):
                        continue

                    sign = -1.0 if fov < mid_fov else 1.0
                    scene_id += 1
                    writer.writerow(
                        [
                            scene_id,
                            scan,
                            fov,
                            f"{scene_lat:.6f}",
                            f"{scene_lon:.6f}",
                            f"{sign * scene_zenith:.6f}",
                            f"{scene_zenith:.6f}",
                            f"{scene_azimuth:.6f}",
                            f"{args.sst:.6f}",
                            f"{args.u10:.6f}",
                            f"{args.sss:.6f}",
                        ]
                        + [f"{value:.6f}" for value in obs_tb]
                    )

                    if scene_id >= args.max_scenes:
                        print(f"wrote {scene_id} scenes to {args.output_csv}")
                        return

        print(f"wrote {scene_id} scenes to {args.output_csv}")


if __name__ == "__main__":
    main()
