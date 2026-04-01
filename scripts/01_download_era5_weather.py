"""
Download ERA5 hourly weather variables for a given location.

Usage:
    python 01_download_era5_weather.py --city manchester --lat 53.3537 --lon -2.2749
    python 01_download_era5_weather.py --city london --lat 51.4700 --lon -0.4543
    python 01_download_era5_weather.py --city leeds --lat 53.8655 --lon -1.6609

Covers 2018-2024 (84 months). Downloads are split into instant and accumulated
variables to avoid CDS API server failures.

Requires a CDS API key configured in ~/.cdsapirc
See: https://cds.climate.copernicus.eu/api-how-to
"""

import argparse
import cdsapi
import os
import io
import glob
import time
import calendar
import zipfile
import pandas as pd
import xarray as xr
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ------------------------------------------------------------------ #
#  ERA5 variable definitions
# ------------------------------------------------------------------ #

# Instantaneous variables (hourly snapshots)
VARS_INSTANT = [
    '10m_u_component_of_wind',
    '10m_v_component_of_wind',
    '2m_temperature',
    '2m_dewpoint_temperature',
    'surface_pressure',
    'total_cloud_cover',
    'boundary_layer_height',
]

# Accumulated variables (hourly accumulations)
VARS_ACCUM = [
    'total_precipitation',
    'surface_solar_radiation_downwards',
]


# ------------------------------------------------------------------ #
#  Download a single month
# ------------------------------------------------------------------ #
def download_one_month(client, year, month, step_type, output_file, area):
    """Download one month of ERA5 data for either instant or accum variables."""
    variables = VARS_INSTANT if step_type == 'instant' else VARS_ACCUM
    _, max_day = calendar.monthrange(year, month)

    request = {
        "product_type": "reanalysis",
        "variable": variables,
        "year": str(year),
        "month": f"{month:02d}",
        "day": [f"{d:02d}" for d in range(1, max_day + 1)],
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": area,
        "data_format": "netcdf",
    }

    for attempt in range(3):
        try:
            client.retrieve("reanalysis-era5-single-levels", request, output_file)
            size_mb = os.path.getsize(output_file) / (1024 * 1024)
            if size_mb > 0.01:
                return True
            else:
                raise Exception(f"File too small ({size_mb:.3f} MB)")
        except Exception as e:
            if os.path.exists(output_file):
                os.remove(output_file)
            if attempt < 2:
                print(f"    Warning: retry {attempt + 1}/3 ({e})")
                time.sleep(10)
            else:
                print(f"    Failed: {e}")
                return False


# ------------------------------------------------------------------ #
#  Read a single .nc file (handles ZIP or plain netCDF)
# ------------------------------------------------------------------ #
def read_nc(nc_file, lat, lon):
    """Read a netCDF file, handling ZIP-wrapped files from newer CDS API."""
    try:
        with zipfile.ZipFile(nc_file, 'r') as z:
            names = z.namelist()
            nc_name = [n for n in names if n.endswith('.nc')]
            if not nc_name:
                raise ValueError(f"No .nc file found in ZIP: {names}")
            with z.open(nc_name[0]) as f:
                ds = xr.open_dataset(io.BytesIO(f.read()))
    except zipfile.BadZipFile:
        ds = xr.open_dataset(nc_file)

    df = (
        ds.sel(latitude=lat, longitude=lon, method="nearest")
        .to_dataframe()
        .reset_index()
    )
    ds.close()

    if 'valid_time' not in df.columns and 'time' in df.columns:
        df = df.rename(columns={'time': 'valid_time'})

    return df


# ------------------------------------------------------------------ #
#  Main download routine
# ------------------------------------------------------------------ #
def download_era5(city, lat, lon, data_dir, start_year=2018, end_year=2024):
    """Download ERA5 data for a given city and coordinate."""
    area = [lat + 0.3, lon - 0.3, lat - 0.3, lon + 0.3]

    print("=" * 60)
    print(f"Downloading ERA5 weather data for: {city.upper()}")
    print(f"  Coordinates: {lat}, {lon}")
    print(f"  Period: {start_year}-{end_year}")
    print("=" * 60)

    os.makedirs(data_dir, exist_ok=True)

    try:
        client = cdsapi.Client(timeout=600, verify=True)
        print("CDS API initialized successfully\n")
    except Exception as e:
        print(f"CDS API initialization failed: {e}")
        return

    total_months = (end_year - start_year + 1) * 12
    stats = {"exist_i": 0, "exist_a": 0,
             "new_i": 0, "new_a": 0,
             "fail_i": [], "fail_a": []}

    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            tag = f"{year}-{month:02d}"

            # --- Instantaneous variables ---
            f_instant = f"{data_dir}/{city}_instant_{year}_{month:02d}.nc"
            if os.path.exists(f_instant) and os.path.getsize(f_instant) > 1024 * 10:
                print(f"  {tag} instant: already exists")
                stats["exist_i"] += 1
            else:
                print(f"  {tag} instant: downloading...", end=" ", flush=True)
                ok = download_one_month(client, year, month, "instant", f_instant, area)
                if ok:
                    mb = os.path.getsize(f_instant) / 1024 / 1024
                    print(f"done ({mb:.2f} MB)")
                    stats["new_i"] += 1
                    time.sleep(1)
                else:
                    stats["fail_i"].append(tag)
                    time.sleep(3)

            # --- Accumulated variables ---
            f_accum = f"{data_dir}/{city}_accum_{year}_{month:02d}.nc"
            if os.path.exists(f_accum) and os.path.getsize(f_accum) > 1024 * 10:
                print(f"  {tag} accum:   already exists")
                stats["exist_a"] += 1
            else:
                print(f"  {tag} accum:   downloading...", end=" ", flush=True)
                ok = download_one_month(client, year, month, "accum", f_accum, area)
                if ok:
                    mb = os.path.getsize(f_accum) / 1024 / 1024
                    print(f"done ({mb:.2f} MB)")
                    stats["new_a"] += 1
                    time.sleep(1)
                else:
                    stats["fail_a"].append(tag)
                    time.sleep(3)

    # --- Summary ---
    print("\n" + "=" * 60)
    print("Download Summary:")
    done_i = stats["exist_i"] + stats["new_i"]
    done_a = stats["exist_a"] + stats["new_a"]
    print(f"  instant  completed: {done_i}/{total_months}  failed: {len(stats['fail_i'])}")
    print(f"  accum    completed: {done_a}/{total_months}  failed: {len(stats['fail_a'])}")
    if stats["fail_i"]:
        print(f"  instant failed months: {', '.join(stats['fail_i'])}")
    if stats["fail_a"]:
        print(f"  accum   failed months: {', '.join(stats['fail_a'])}")
    print("=" * 60)

    if done_i > 0 or done_a > 0:
        process_era5_data(city, lat, lon, data_dir)


# ------------------------------------------------------------------ #
#  Process and merge downloaded data
# ------------------------------------------------------------------ #
def process_era5_data(city, lat, lon, data_dir):
    """Merge monthly ERA5 files into a single weather CSV."""
    print("\n" + "=" * 60)
    print(f"Processing ERA5 data for: {city.upper()}")
    print("=" * 60)

    instant_files = sorted(glob.glob(f"{data_dir}/{city}_instant_*.nc"))
    accum_files = sorted(glob.glob(f"{data_dir}/{city}_accum_*.nc"))

    print(f"  instant files: {len(instant_files)}")
    print(f"  accum   files: {len(accum_files)}\n")

    if not instant_files and not accum_files:
        print("No files found")
        return None

    def file_key(path):
        base = os.path.basename(path)
        parts = base.replace('.nc', '').split('_')
        return f"{parts[-2]}_{parts[-1]}"

    instant_map = {file_key(f): f for f in instant_files}
    accum_map = {file_key(f): f for f in accum_files}
    all_keys = sorted(set(instant_map) | set(accum_map))

    all_data = []
    successful = 0

    for i, k in enumerate(all_keys):
        if (i + 1) % 12 == 0 or (i + 1) == len(all_keys):
            print(f"  Progress: {i + 1}/{len(all_keys)}")

        try:
            if k not in instant_map:
                print(f"  Warning: {k} missing instant file, skipping")
                continue
            df_i = read_nc(instant_map[k], lat, lon)

            df_a = read_nc(accum_map[k], lat, lon) if k in accum_map else None
            if df_a is None:
                print(f"  Warning: {k} missing accum file, filling with NaN")

            dt = pd.to_datetime(df_i['valid_time'])
            p = pd.DataFrame()

            # Time variables
            p['date'] = dt
            p['date_unix'] = dt.astype('int64') // 10**9
            p['year'] = dt.dt.year
            p['month'] = dt.dt.month
            p['week'] = dt.dt.isocalendar().week.values
            p['weekday'] = dt.dt.dayofweek + 1
            p['hour'] = dt.dt.hour
            p['day_julian'] = dt.dt.dayofyear

            # Derived instant variables
            p['air_temp'] = df_i['t2m'] - 273.15                      # K -> C
            p['ws'] = np.sqrt(df_i['u10']**2 + df_i['v10']**2)        # m/s
            p['wd'] = (270 - np.rad2deg(
                np.arctan2(df_i['v10'], df_i['u10'])
            )) % 360                                                    # degrees
            tc = p['air_temp']
            tdc = df_i['d2m'] - 273.15
            p['RH'] = 100 * (
                np.exp((17.625 * tdc) / (243.04 + tdc)) /
                np.exp((17.625 * tc) / (243.04 + tc))
            )                                                           # %
            p['atmos_pres'] = df_i['sp'] / 100                         # Pa -> hPa
            p['cloud_cover'] = df_i['tcc'] * 100                       # 0-1 -> %
            p['blh_m'] = df_i['blh']                                   # m

            # Derived accumulated variables
            if df_a is not None:
                p['precip_total'] = df_a['tp'].values * 1000           # m -> mm
                p['solar_wm2'] = (df_a['ssrd'].values / 3600).clip(min=0)  # J/m2 -> W/m2
            else:
                p['precip_total'] = np.nan
                p['solar_wm2'] = np.nan

            all_data.append(p)
            successful += 1

        except Exception as e:
            print(f"  Error processing {k}: {e}")

    print(f"\nSuccessfully processed: {successful}/{len(all_keys)} months")

    if not all_data:
        print("No data to merge")
        return None

    df_result = (
        pd.concat(all_data, ignore_index=True)
        .drop_duplicates(subset=['date'])
        .sort_values('date')
        .reset_index(drop=True)
    )

    print(f"\nFinal dataset:")
    print(f"  Rows:       {len(df_result):,}")
    print(f"  Date range: {df_result['date'].min()} to {df_result['date'].max()}")

    for col, label, unit in [
        ('air_temp', 'Temperature', 'C'),
        ('ws', 'Wind speed', 'm/s'),
        ('wd', 'Wind dir', 'deg'),
        ('RH', 'Humidity', '%'),
        ('atmos_pres', 'Pressure', 'hPa'),
        ('cloud_cover', 'Cloud', '%'),
        ('precip_total', 'Precip', 'mm'),
        ('blh_m', 'BLH', 'm'),
        ('solar_wm2', 'Solar', 'W/m2'),
    ]:
        if col in df_result.columns:
            valid = df_result[col].notna().sum()
            mn = df_result[col].min()
            mx = df_result[col].max()
            print(f"  {label:12s} valid: {valid:,}  range: {mn:.2f} - {mx:.2f} {unit}")

    out_csv = f"{data_dir}/{city}_weather_final.csv"
    df_result.to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    return df_result


# ------------------------------------------------------------------ #
#  Entry point
# ------------------------------------------------------------------ #
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download and process ERA5 weather data for a given location."
    )
    parser.add_argument("--city", required=True, help="City name (used for file naming)")
    parser.add_argument("--lat", required=True, type=float, help="Latitude")
    parser.add_argument("--lon", required=True, type=float, help="Longitude")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory (default: era5_data_{city})")
    parser.add_argument("--start-year", default=2018, type=int)
    parser.add_argument("--end-year", default=2024, type=int)

    args = parser.parse_args()
    data_dir = args.output_dir or f"era5_data_{args.city}"

    download_era5(args.city, args.lat, args.lon, data_dir,
                  args.start_year, args.end_year)
