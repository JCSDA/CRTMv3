#!/usr/bin/env python3
"""Extract SOCA/JEDI GMI observations and GeoVaLs for CRTM comparison."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract gmi_gpm ObsValue/GeoVaLs into a CRTM scene/profile CSV."
    )
    parser.add_argument("obs_file", type=Path)
    parser.add_argument("geoval_file", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--salinity", type=float, default=35.0, help="Sea-surface salinity [psu]")
    parser.add_argument("--ocean-min", type=float, default=0.99, help="Minimum water_area_fraction")
    return parser.parse_args()


def as_array(var, dtype=float) -> np.ndarray:
    return np.asarray(np.ma.filled(var[:], np.nan), dtype=dtype)


def main() -> None:
    args = parse_args()
    with Dataset(args.obs_file) as obs, Dataset(args.geoval_file) as geoval:
        metadata = obs.groups["MetaData"]
        obs_value = obs.groups["ObsValue"]
        preqc = obs.groups["PreQC"]

        obs_tb = as_array(obs_value.variables["brightnessTemperature"])
        qc = as_array(preqc.variables["brightnessTemperature"])
        n_locs, n_channels = obs_tb.shape
        n_geoval_locs = len(geoval.dimensions["nlocs"])
        n_layers = len(geoval.dimensions["num_profile_levels"])
        if n_locs != n_geoval_locs:
            raise SystemExit(f"location mismatch: obs={n_locs} geoval={n_geoval_locs}")
        if n_channels != len(geoval.dimensions["nchans"]):
            raise SystemExit("channel-count mismatch between obs and geoval")

        obs_channels = as_array(metadata.variables["sensorChannelNumber"], dtype=int)
        geo_channels = as_array(geoval.variables["sensor_chan"], dtype=int)
        if not np.array_equal(obs_channels, geo_channels):
            raise SystemExit("sensor channel numbers differ between obs and geoval")

        lat_obs = as_array(metadata.variables["latitude"])
        lon_obs = as_array(metadata.variables["longitude"])
        lat_geo = as_array(geoval.variables["latitude"])
        lon_geo = as_array(geoval.variables["longitude"])
        if np.nanmax(np.abs(lat_obs - lat_geo)) > 1.0e-4:
            raise SystemExit("latitude arrays differ between obs and geoval")
        if np.nanmax(np.abs(lon_obs - lon_geo)) > 1.0e-4:
            raise SystemExit("longitude arrays differ between obs and geoval")

        pressure = as_array(geoval.variables["air_pressure"]) / 100.0
        level_pressure = as_array(geoval.variables["air_pressure_levels"]) / 100.0
        temperature = as_array(geoval.variables["air_temperature"])
        h2o = as_array(geoval.variables["humidity_mixing_ratio"])
        o3 = as_array(geoval.variables["mole_fraction_of_ozone_in_air"])
        co2 = as_array(geoval.variables["mole_fraction_of_carbon_dioxide_in_air"]) * 1.0e6
        clw = as_array(geoval.variables["mass_content_of_cloud_liquid_water_in_atmosphere_layer"])
        cli = as_array(geoval.variables["mass_content_of_cloud_ice_in_atmosphere_layer"])
        re_liq = as_array(geoval.variables["effective_radius_of_cloud_liquid_water_particle"])
        re_ice = as_array(geoval.variables["effective_radius_of_cloud_ice_particle"])

        water_fraction = as_array(geoval.variables["water_area_fraction"])
        land_fraction = as_array(geoval.variables["land_area_fraction"])
        ice_fraction = as_array(geoval.variables["ice_area_fraction"])
        snow_fraction = as_array(geoval.variables["surface_snow_area_fraction"])
        sst = as_array(geoval.variables["surface_temperature_where_sea"])
        u10 = as_array(geoval.variables["surface_wind_speed"])
        wind_direction = as_array(geoval.variables["surface_wind_from_direction"])

        zenith = as_array(metadata.variables["sensorZenithAngle"])
        azimuth = as_array(metadata.variables["sensorAzimuthAngle"])
        scan_position = as_array(metadata.variables["sensorScanPosition"])
        solar_zenith = as_array(metadata.variables["solarZenithAngle"])
        solar_azimuth = as_array(metadata.variables["solarAzimuthAngle"])

        header = [
            "loc",
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
        ]
        header += [f"obs_tb_{i}" for i in range(1, n_channels + 1)]
        header += [f"preqc_{i}" for i in range(1, n_channels + 1)]
        header += [f"level_pressure_{i}" for i in range(0, n_layers + 1)]
        for stem in (
            "pressure",
            "temperature",
            "h2o",
            "o3",
            "co2",
            "cloud_liquid",
            "cloud_ice",
            "re_liquid",
            "re_ice",
        ):
            header += [f"{stem}_{i}" for i in range(1, n_layers + 1)]

        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        with args.output_csv.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            for loc in range(n_locs):
                if water_fraction[loc] < args.ocean_min:
                    continue
                required = [
                    lat_geo[loc],
                    lon_geo[loc],
                    zenith[loc],
                    azimuth[loc],
                    sst[loc],
                    u10[loc],
                    wind_direction[loc],
                ]
                if not np.all(np.isfinite(required)):
                    continue
                if not np.all(np.isfinite(obs_tb[loc, :])):
                    continue
                if not np.all(np.isfinite(level_pressure[loc, :])):
                    continue
                if not np.all(np.isfinite(pressure[loc, :])):
                    continue

                row = [
                    loc + 1,
                    lat_geo[loc],
                    lon_geo[loc],
                    scan_position[loc],
                    zenith[loc],
                    zenith[loc],
                    azimuth[loc],
                    solar_zenith[loc],
                    solar_azimuth[loc],
                    water_fraction[loc],
                    land_fraction[loc],
                    ice_fraction[loc],
                    snow_fraction[loc],
                    sst[loc],
                    u10[loc],
                    wind_direction[loc],
                    args.salinity,
                ]
                arrays = (
                    obs_tb[loc, :],
                    qc[loc, :],
                    level_pressure[loc, :],
                    pressure[loc, :],
                    temperature[loc, :],
                    h2o[loc, :],
                    o3[loc, :],
                    co2[loc, :],
                    clw[loc, :],
                    cli[loc, :],
                    re_liq[loc, :],
                    re_ice[loc, :],
                )
                values = row + [value for array in arrays for value in array]
                writer.writerow([f"{value:.8g}" if isinstance(value, float) else value for value in values])
                count += 1

        print(f"wrote {count} scenes to {args.output_csv}")


if __name__ == "__main__":
    main()
