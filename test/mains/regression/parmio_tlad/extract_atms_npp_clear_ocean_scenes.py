#!/usr/bin/env python3
"""Extract clear-ocean ATMS-NPP ObsValue/GeoVaLs for standalone CRTM."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from netCDF4 import Dataset


HYDROMETEOR_NAMES = (
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "mass_content_of_cloud_ice_in_atmosphere_layer",
    "mass_content_of_rain_in_atmosphere_layer",
    "mass_content_of_snow_in_atmosphere_layer",
    "mass_content_of_graupel_in_atmosphere_layer",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract clear-ocean atms_npp ObsValue/GeoVaLs into a CRTM scene/profile CSV."
    )
    parser.add_argument("obsout_file", type=Path)
    parser.add_argument("geoval_file", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--salinity", type=float, default=35.0, help="Sea-surface salinity [psu]")
    parser.add_argument("--ocean-min", type=float, default=0.99, help="Minimum water_area_fraction")
    parser.add_argument(
        "--hydrometeor-max",
        type=float,
        default=1.0e-12,
        help="Maximum allowed absolute cloud/rain/snow/graupel mass content",
    )
    parser.add_argument(
        "--cloud-fraction-max",
        type=float,
        default=1.0e-12,
        help="Maximum allowed cloud_area_fraction_in_atmosphere_layer",
    )
    parser.add_argument(
        "--max-scenes",
        type=int,
        default=512,
        help="Maximum number of scenes to write; <=0 writes all selected scenes.",
    )
    parser.add_argument("--chunk-size", type=int, default=50000)
    return parser.parse_args()


def as_array(value, dtype=float) -> np.ndarray:
    return np.asarray(np.ma.filled(value, np.nan), dtype=dtype)


def finite_obs(obs: np.ndarray) -> np.ndarray:
    return np.isfinite(obs).all(axis=1) & (obs > 0.0).all(axis=1) & (obs < 1000.0).all(axis=1)


def wind_from_direction(eastward: float, northward: float) -> float:
    # Meteorological direction the wind is from, degrees clockwise from north.
    return float((np.degrees(np.arctan2(-eastward, -northward)) + 360.0) % 360.0)


def epoch_to_ymd(epoch_seconds: int) -> tuple[int, int, int]:
    dt = datetime.fromtimestamp(int(epoch_seconds), timezone.utc)
    return dt.year, dt.month, dt.day


def max_abs_profile(var, start: int, end: int) -> np.ndarray:
    values = np.asarray(var[start:end, :], dtype=np.float64)
    values = np.where(np.isfinite(values), np.abs(values), np.inf)
    return values.max(axis=1)


def find_candidates(args: argparse.Namespace) -> np.ndarray:
    candidates: list[np.ndarray] = []
    with Dataset(args.geoval_file) as geoval, h5py.File(args.obsout_file, "r") as obsout:
        n_locs = len(geoval.dimensions["nlocs"])
        obs_tb = obsout["/ObsValue/brightnessTemperature"]

        for start in range(0, n_locs, args.chunk_size):
            end = min(start + args.chunk_size, n_locs)
            water = as_array(geoval["water_area_fraction"][start:end, 0])
            land = as_array(geoval["land_area_fraction"][start:end, 0])
            ice = as_array(geoval["ice_area_fraction"][start:end, 0])
            snow = as_array(geoval["surface_snow_area_fraction"][start:end, 0])
            obs_ok = finite_obs(np.asarray(obs_tb[start:end, :], dtype=np.float64))

            selected = (
                (water >= args.ocean_min)
                & (land <= 1.0 - args.ocean_min)
                & (ice <= 1.0 - args.ocean_min)
                & (snow <= 1.0 - args.ocean_min)
                & obs_ok
            )

            hydro_max = np.zeros(end - start, dtype=np.float64)
            for name in HYDROMETEOR_NAMES:
                hydro_max = np.maximum(hydro_max, max_abs_profile(geoval[name], start, end))
            cloud_fraction_max = max_abs_profile(
                geoval["cloud_area_fraction_in_atmosphere_layer"], start, end
            )

            selected &= hydro_max <= args.hydrometeor_max
            selected &= cloud_fraction_max <= args.cloud_fraction_max

            idx = np.nonzero(selected)[0] + start
            if idx.size:
                candidates.append(idx.astype(np.int64))

    if not candidates:
        return np.empty(0, dtype=np.int64)
    all_candidates = np.concatenate(candidates)
    if args.max_scenes > 0 and all_candidates.size > args.max_scenes:
        pick = np.linspace(0, all_candidates.size - 1, args.max_scenes, dtype=np.int64)
        all_candidates = all_candidates[pick]
    return all_candidates


def main() -> None:
    args = parse_args()
    candidates = find_candidates(args)
    if candidates.size == 0:
        raise SystemExit("No clear-ocean ATMS scenes matched the requested thresholds")

    with Dataset(args.geoval_file) as geoval, h5py.File(args.obsout_file, "r") as obsout:
        obs_tb = obsout["/ObsValue/brightnessTemperature"]
        hofx_tb = obsout["/hofx/brightnessTemperature"]
        channels = np.asarray(obsout["/Channel"][:], dtype=int)
        n_channels = channels.size
        n_layers = len(geoval.dimensions["air_pressure_nval"])
        n_levels = len(geoval.dimensions["air_pressure_levels_nval"])

        header = [
            "loc",
            "year",
            "month",
            "day",
            "lat",
            "lon",
            "scan_position",
            "scan_angle",
            "zenith",
            "azimuth",
            "solar_zenith",
            "solar_azimuth",
            "water_fraction",
            "land_fraction",
            "ice_fraction",
            "snow_fraction",
            "sst",
            "u10",
            "wind_direction",
            "sss",
            "hydrometeor_max",
            "cloud_fraction_max",
        ]
        header += [f"obs_tb_{channel}" for channel in channels]
        header += [f"hofx_tb_{channel}" for channel in channels]
        header += [f"level_pressure_{i}" for i in range(n_levels)]
        for stem in ("pressure", "temperature", "h2o", "o3"):
            header += [f"{stem}_{i}" for i in range(1, n_layers + 1)]

        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_csv.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            for loc in candidates:
                eastward = float(geoval["eastward_wind_at_surface"][loc, 0])
                northward = float(geoval["northward_wind_at_surface"][loc, 0])
                year, month, day = epoch_to_ymd(obsout["/MetaData/dateTime"][loc])
                hydro_max = max(
                    float(np.nanmax(np.abs(geoval[name][loc, :]))) for name in HYDROMETEOR_NAMES
                )
                cloud_fraction_max = float(
                    np.nanmax(np.abs(geoval["cloud_area_fraction_in_atmosphere_layer"][loc, :]))
                )

                row = [
                    int(loc + 1),
                    year,
                    month,
                    day,
                    float(obsout["/MetaData/latitude"][loc]),
                    float(obsout["/MetaData/longitude"][loc]),
                    int(obsout["/MetaData/sensorScanPosition"][loc]),
                    float(obsout["/MetaData/sensorViewAngle"][loc]),
                    float(obsout["/MetaData/sensorZenithAngle"][loc]),
                    float(obsout["/MetaData/sensorAzimuthAngle"][loc]),
                    float(obsout["/MetaData/solarZenithAngle"][loc]),
                    float(obsout["/MetaData/solarAzimuthAngle"][loc]),
                    float(geoval["water_area_fraction"][loc, 0]),
                    float(geoval["land_area_fraction"][loc, 0]),
                    float(geoval["ice_area_fraction"][loc, 0]),
                    float(geoval["surface_snow_area_fraction"][loc, 0]),
                    float(geoval["skin_temperature_at_surface_where_sea"][loc, 0]),
                    float(np.hypot(eastward, northward)),
                    wind_from_direction(eastward, northward),
                    args.salinity,
                    hydro_max,
                    cloud_fraction_max,
                ]
                arrays = (
                    np.asarray(obs_tb[loc, :], dtype=np.float64),
                    np.asarray(hofx_tb[loc, :], dtype=np.float64),
                    np.asarray(geoval["air_pressure_levels"][loc, :], dtype=np.float64) / 100.0,
                    np.asarray(geoval["air_pressure"][loc, :], dtype=np.float64) / 100.0,
                    np.asarray(geoval["air_temperature"][loc, :], dtype=np.float64),
                    np.asarray(
                        geoval["water_vapor_mixing_ratio_wrt_dry_air"][loc, :],
                        dtype=np.float64,
                    )
                    * 1000.0,
                    np.asarray(geoval["mole_fraction_of_ozone_in_air"][loc, :], dtype=np.float64)
                    * 1.0e6,
                )
                values = row + [value for array in arrays for value in array]
                writer.writerow(
                    [f"{value:.8g}" if isinstance(value, (float, np.floating)) else value for value in values]
                )

    print(
        f"wrote {candidates.size} scenes to {args.output_csv} "
        f"(hydrometeor_max <= {args.hydrometeor_max:g}, "
        f"cloud_fraction_max <= {args.cloud_fraction_max:g})"
    )


if __name__ == "__main__":
    main()
