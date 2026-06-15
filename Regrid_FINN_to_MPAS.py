# -----------------------------------------------------------------------------
# 1. Read user options and configuration from YAML file
# -----------------------------------------------------------------------------
#!/usr/bin/env python
# coding: utf-8
"""
Regrid FINNv2.5 Text files to MPAS-A compatible netcdf input files
All user options are read from config_finn_to_mpas.yaml
"""

import pandas as pd
import numpy as np
import xarray as xr
import os
import glob
import yaml
import re
from pathlib import Path
from datetime import datetime
import time
import argparse
import sys


def extract_date_from_filename(file_path: str):
    """Extract a date from filename using YYYYMMDD or YYYYDDD patterns."""
    name = Path(file_path).name
    m_ymd = re.search(r'(\d{8})', name)
    if m_ymd:
        return datetime.strptime(m_ymd.group(1), "%Y%m%d")
    m_ydoy = re.search(r'(\d{7})', name)
    if m_ydoy:
        try:
            return datetime.strptime(m_ydoy.group(1), "%Y%j")
        except ValueError:
            return None
    return None


def make_xtime_array(date_obj: datetime, hourly: bool, strlen: int = 64):
    if hourly:
        ntimes = 24
        timestrs = [f"{date_obj.strftime('%Y-%m-%d')}_{str(h).zfill(2)}:00:00" for h in range(24)]
    else:
        ntimes = 1
        timestrs = [date_obj.strftime('%Y-%m-%d_00:00:00')]
    xtime_arr = np.zeros((ntimes, strlen), dtype='S1')
    for i, s in enumerate(timestrs):
        s_padded = s.ljust(strlen, '\0')
        s_bytes = np.array(list(s_padded), dtype='S1')
        xtime_arr[i, :] = s_bytes
    return xtime_arr


def find_date_column(df_columns, preferred_column=None, candidates=None):
    """Find annual FINN date column by exact or case-insensitive match."""
    if preferred_column and preferred_column in df_columns:
        return preferred_column

    lowered = {col.lower(): col for col in df_columns}
    if preferred_column and preferred_column.lower() in lowered:
        return lowered[preferred_column.lower()]

    candidates = candidates or []
    for cand in candidates:
        if cand in df_columns:
            return cand
        if cand.lower() in lowered:
            return lowered[cand.lower()]
    return None


def parse_annual_dates(date_series, date_formats):
    """Parse mixed annual FINN date column values into pandas Timestamps."""
    s = date_series.astype(str).str.strip()
    parsed = pd.Series(pd.NaT, index=s.index)
    for fmt in date_formats:
        newly = pd.to_datetime(s, format=fmt, errors='coerce')
        parsed = parsed.fillna(newly)
    parsed = parsed.fillna(pd.to_datetime(s, errors='coerce'))
    return parsed


def extract_year_from_filename(file_path: str):
    """Extract a plausible year from filename."""
    m = re.search(r'(19\d{2}|20\d{2})', Path(file_path).name)
    if m:
        return int(m.group(1))
    return None


def parse_numeric_annual_day(date_series, fallback_year, mode='auto'):
    """Parse numeric annual day values as DOY or DOM.

    auto mode infers DOY when max value is > 31.
    """
    nums = pd.to_numeric(date_series, errors='coerce')
    parsed = pd.Series(pd.NaT, index=nums.index)
    valid = nums.notna()
    if not np.any(valid):
        return parsed

    nums_valid = nums[valid].astype(int)
    use_mode = mode
    if use_mode == 'auto':
        use_mode = 'doy' if nums_valid.max() > 31 else 'dom'

    if use_mode == 'doy':
        yyyyddd = nums_valid.map(lambda d: f"{fallback_year}{int(d):03d}")
        parsed_valid = pd.to_datetime(yyyyddd, format='%Y%j', errors='coerce')
        parsed.loc[valid] = parsed_valid
    elif use_mode == 'dom':
        # Fallback DOM interpretation: use January in fallback year.
        # Most annual FINN files use DOY; this only supports simple DOM cases.
        yyyymmdd = nums_valid.map(lambda d: f"{fallback_year}01{int(d):02d}")
        parsed_valid = pd.to_datetime(yyyymmdd, format='%Y%m%d', errors='coerce')
        parsed.loc[valid] = parsed_valid

    return parsed

# Set up the argument parser
parser = argparse.ArgumentParser(description="Process WRF-Chem data and calculate zonal stats.")
parser.add_argument("config", help="Path to the user-defined YAML configuration file")

# Parse the arguments
args = parser.parse_args()

# Open the user-defined YAML file
try:
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    print(f"Error: The file '{args.config}' was not found.")
    sys.exit(1)

# Get year and month from start_date and end_date
from datetime import datetime, timedelta
start_date = config.get('start_date', None)
end_date = config.get('end_date', None)
if not start_date or not end_date:
    raise ValueError('start_date and end_date must be specified in the config file.')
start_dt = datetime.strptime(start_date, '%Y-%m-%d')
end_dt = datetime.strptime(end_date, '%Y-%m-%d')
if end_dt < start_dt:
    raise ValueError('end_date must be greater than or equal to start_date.')

file_type = str(config.get('file_type', 'daily')).strip().lower()
if file_type not in {'daily', 'annual'}:
    raise ValueError("file_type must be either 'daily' or 'annual'.")

annual_date_column = config.get('annual_date_column', None)
# Preferred simplified option: annual_date_column may be either a string or a list.
# Backward compatibility: if legacy annual_date_columns is provided, use it as fallback.
if isinstance(annual_date_column, list):
    annual_date_candidates = annual_date_column
elif isinstance(annual_date_column, str) and annual_date_column.strip():
    annual_date_candidates = [annual_date_column.strip()]
else:
    legacy_cols = config.get('annual_date_columns', None)
    if isinstance(legacy_cols, list) and legacy_cols:
        annual_date_candidates = legacy_cols
    else:
        annual_date_candidates = ['DAY', 'DATE', 'day', 'date', 'YYYYMMDD', 'yyyymmdd', 'DATE_ID', 'date_id']
annual_date_formats = config.get('annual_date_formats', ['%Y%m%d', '%Y-%m-%d', '%Y/%m/%d', '%Y%j'])
annual_day_mode = str(config.get('annual_day_mode', 'auto')).strip().lower()
if annual_day_mode not in {'auto', 'doy', 'dom'}:
    raise ValueError("annual_day_mode must be one of 'auto', 'doy', or 'dom'.")
annual_reference_year = config.get('annual_reference_year', None)
annual_chunksize = int(config.get('annual_chunksize', 250000))

# Build emission file lists depending on input file type
emis_file_date_pairs = []  # Used for daily mode
emis_file_list = []
if file_type == 'daily':
    months = []
    dt = start_dt.replace(day=1)
    while dt <= end_dt:
        months.append((dt.year, dt.month))
        # Go to next month
        if dt.month == 12:
            dt = dt.replace(year=dt.year + 1, month=1)
        else:
            dt = dt.replace(month=dt.month + 1)

    for year, month in months:
        emis_dir = config['emis_dir'].format(year=year)
        emis_file_pattern = config['emis_file_pattern'].format(year=year, month=month)
        for f in glob.glob(os.path.join(emis_dir, emis_file_pattern)):
            file_date = extract_date_from_filename(f)
            if file_date is not None and start_dt <= file_date <= end_dt:
                emis_file_date_pairs.append((file_date, f))
    emis_file_date_pairs = sorted(emis_file_date_pairs)
    emis_file_list = [f for file_date, f in emis_file_date_pairs]
else:
    years = range(start_dt.year, end_dt.year + 1)
    found_files = set()
    for year in years:
        format_kwargs = {'year': year, 'month': '*', 'day': '*'}
        emis_dir = config['emis_dir'].format(**format_kwargs)
        emis_file_pattern = config['emis_file_pattern'].format(**format_kwargs)
        for f in glob.glob(os.path.join(emis_dir, emis_file_pattern)):
            found_files.add(f)
    emis_file_list = sorted(found_files)

if not emis_file_list:
    raise FileNotFoundError(
        f"No FINN input files found for file_type='{file_type}' between {start_date} and {end_date}."
    )
dst_file_dir = config['dst_file_dir']  # Output directory for NetCDF files
mpas_grid_file = config['mpas_grid_file']  # Path to MPAS grid NetCDF file
HOURLY = config.get('HOURLY', False)  # Whether to output hourly emissions
lt_fac = np.array(config.get('lt_fac', [.43, .43, .43, .43, .43, .43, .43, .43, .43, 3., 6., 10., 14., 17., 14., 12., 9., 6., 3., .43, .43, .43, .43, .43]))  # Diurnal profile
output_format = 'NETCDF3_64BIT'  # NetCDF output format (cdf5)
compression = config.get('compression', {'zlib': True, 'complevel': 4, 'shuffle': True})  # NetCDF compression options


# Read species mapping and types from YAML
species_map = config.get('species_map', {})  # FINN species to MPAS variable mapping
species_to_map = set(species_map.keys())  # Set of species to process
species_type = config.get('species_type', {})  # Dict specifying 'aerosol' or 'gas' for each species

# Constants for unit conversion
AVOGADRO = 6.022e23  # Avogadro's number
SECONDS_PER_DAY = 86400  # Seconds in a day
AEROSOL_MW = 12.0  # Aerosol molecular weight (g/mole)



STR_LEN = 64  # MPAS string length requirement




# --- AGGREGATE ALL EMISSIONS AND TIMES ---
print(f"Opening MPAS grid file: {mpas_grid_file}")
ds = xr.open_dataset(mpas_grid_file, engine="netcdf4", decode_times=False)
areaCell_m2 = ds["areaCell"].values
areaCell_cm2 = areaCell_m2 * 1e4
latCell_deg = np.degrees(ds["latCell"].values)
lonCell_deg = np.degrees(ds["lonCell"].values)
nCells = ds.sizes["nCells"]

# Prepare output arrays
all_xtime = []
processed_dates = []
species_arrays = {species_map[col]: [] for col in species_to_map if col in species_map and species_map[col] != "co2_biob_modis"}
scalar_std_arrays = {
    f"{species_map[col]}_std": []
    for col in species_to_map
    if col in species_map and species_map[col] != "co2_biob_modis" and species_type.get(col, None) == "scalar"
}

def process_one_day(df, current_date):
    # Modified read statement to accomodate FINNv1 fortran format strings FGL
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    for col in df.columns:
        cleaned = df[col].astype(str).str.replace('D', 'E', regex=False).str.replace('d', 'E', regex=False)
        if col not in {'LATI', 'LONGI'}:
            try:
                df[col] = pd.to_numeric(cleaned)
            except (ValueError, TypeError):
                df[col] = cleaned
        else:
            df[col] = cleaned

    if 'LATI' not in df.columns or 'LONGI' not in df.columns:
        raise KeyError("Input FINN file is missing required LATI/LONGI columns.")

    # Compute local hour offset for each cell (in hours, rounded to nearest int)
    local_hour_offset = np.round(lonCell_deg / 15.).astype(int) % 24
    lati = pd.to_numeric(df["LATI"], errors='coerce').values
    longi = pd.to_numeric(df["LONGI"], errors='coerce').values

    valid_geo = np.isfinite(lati) & np.isfinite(longi)
    if not np.any(valid_geo):
        return
    df = df.loc[valid_geo].reset_index(drop=True)
    lati = lati[valid_geo]
    longi = longi[valid_geo]

    emission_cols = [col for col in df.columns if col in species_to_map]
    if not emission_cols:
        raise ValueError(
            f"No requested species columns were found in the FINN input. Requested species: {sorted(species_to_map)}"
        )
    emission_per_cell = {col: np.zeros(nCells) for col in emission_cols}
    scalar_count_per_cell = {
        col: np.zeros(nCells, dtype=int)
        for col in emission_cols
        if species_type.get(col, None) == "scalar"
    }
    scalar_sumsq_per_cell = {
        col: np.zeros(nCells, dtype=float)
        for col in emission_cols
        if species_type.get(col, None) == "scalar"
    }
    # Ensure longitude conventions match: convert all to -180 to 180
    def to_minus180_180(lon):
        lon = np.asarray(lon)
        lon = np.where(lon > 180, lon - 360, lon)
        return lon
    lonCell_deg_match = to_minus180_180(lonCell_deg)
    longi_match = to_minus180_180(longi)
    # Assign emissions to nearest cell
    for i in range(len(df)):
        dist_sq = (float(lati[i]) - latCell_deg) ** 2 + (float(longi_match[i]) - lonCell_deg_match) ** 2
        min_index = np.argmin(dist_sq)
        for col in emission_cols:
            if col not in species_map:
                continue
            tmp = pd.to_numeric(df.iloc[i][col], errors='coerce')
            if pd.notna(tmp):
                value = float(tmp)
                emission_per_cell[col][min_index] += value
                if species_type.get(col, None) == "scalar":
                    scalar_count_per_cell[col][min_index] += 1
                    scalar_sumsq_per_cell[col][min_index] += value * value
    # For each emission species, convert units and apply diurnal profile
    for col, values in emission_per_cell.items():
        if col not in species_map:
            continue
        var_name = species_map[col]
        stype = species_type.get(col, None)
        values_conv = np.zeros_like(values, dtype=float)
        std_conv = None
        if stype == "aerosol":
            values_g = values * 1000.0
            values_mol = values_g / AEROSOL_MW
            values_mol_day = values_mol
            values_mol_s = values_mol_day / SECONDS_PER_DAY
            values_molecules = values_mol_s * AVOGADRO
            values_conv = values_molecules / areaCell_cm2
        elif stype == "gas":
            values_mol_s = values / SECONDS_PER_DAY
            values_molecules = values_mol_s * AVOGADRO
            values_conv = values_molecules / areaCell_cm2
        elif stype == "scalar":
            count = scalar_count_per_cell.get(col, np.ones_like(values, dtype=int))
            sumsq = scalar_sumsq_per_cell.get(col, np.zeros_like(values, dtype=float))
            with np.errstate(divide='ignore', invalid='ignore'):
                values_conv = np.true_divide(values, count)
                values_conv[count == 0] = 0.0
                mean_sq = np.true_divide(sumsq, count)
                mean_sq[count == 0] = 0.0
                variance = mean_sq - np.square(values_conv)
                variance = np.maximum(variance, 0.0)
                std_conv = np.sqrt(variance)
                std_conv[count == 0] = 0.0
        else:
            values_conv = values

        if HOURLY:
            wrk_lt_fac = np.array(lt_fac) / np.sum(lt_fac) * 24.
            hourly_values = np.zeros((24, values_conv.shape[0]))
            hourly_std_values = np.zeros((24, values_conv.shape[0])) if std_conv is not None else None
            for cell in range(values_conv.shape[0]):
                offset = local_hour_offset[cell]
                shifted_fac = np.roll(wrk_lt_fac, -offset)
                hourly_values[:, cell] = values_conv[cell] * shifted_fac
                if hourly_std_values is not None:
                    hourly_std_values[:, cell] = std_conv[cell] * shifted_fac
            species_arrays[var_name].append(hourly_values)
            if hourly_std_values is not None:
                scalar_std_arrays[f"{var_name}_std"].append(hourly_std_values)
        else:
            species_arrays[var_name].append(values_conv)
            if std_conv is not None:
                scalar_std_arrays[f"{var_name}_std"].append(std_conv)

    # Add xtime for this day
    xtime_arr = make_xtime_array(current_date, HOURLY, STR_LEN)
    all_xtime.append(xtime_arr)
    processed_dates.append(current_date)


if file_type == 'daily':
    for file_date, csv_file in emis_file_date_pairs:
        print(f"Processing FINN file: {csv_file}")
        df = pd.read_csv(csv_file, header=0, index_col=False)
        process_one_day(df, file_date)
else:
    for csv_file in emis_file_list:
        print(f"Processing annual FINN file: {csv_file}")
        header_df = pd.read_csv(csv_file, header=0, index_col=False, nrows=0)
        header_df = header_df.loc[:, ~header_df.columns.str.contains('^Unnamed')]
        date_col = find_date_column(header_df.columns, candidates=annual_date_candidates)
        if date_col is None:
            raise KeyError(
                f"Could not find annual date column in {csv_file}. "
                f"Set annual_date_column in YAML. Available columns: {list(header_df.columns)}"
            )

        file_year = extract_year_from_filename(csv_file)
        fallback_year = int(annual_reference_year) if annual_reference_year is not None else (file_year or start_dt.year)

        # Keep only rows for requested dates while reading in chunks.
        day_chunks = {}
        for chunk in pd.read_csv(csv_file, header=0, index_col=False, chunksize=annual_chunksize):
            chunk = chunk.loc[:, ~chunk.columns.str.contains('^Unnamed')]
            parsed_dates = parse_annual_dates(chunk[date_col], annual_date_formats)
            if parsed_dates.notna().sum() == 0:
                parsed_dates = parse_numeric_annual_day(chunk[date_col], fallback_year=fallback_year, mode=annual_day_mode)

            in_range = parsed_dates.notna() & (parsed_dates >= start_dt) & (parsed_dates <= end_dt)
            if not np.any(in_range):
                continue

            chunk = chunk.loc[in_range].copy()
            parsed_dates = parsed_dates.loc[in_range]
            chunk['__parsed_date__'] = parsed_dates.dt.normalize()

            for day_ts, day_df in chunk.groupby('__parsed_date__', sort=True):
                day_key = pd.to_datetime(day_ts).to_pydatetime()
                day_df = day_df.drop(columns=['__parsed_date__'])
                day_chunks.setdefault(day_key, []).append(day_df)

        for current_date in sorted(day_chunks.keys()):
            day_df = pd.concat(day_chunks[current_date], ignore_index=True)
            process_one_day(day_df, current_date)

if not processed_dates:
    raise RuntimeError(
        f"No FINN records found in the requested date range {start_date} to {end_date}."
    )

# --- CONCATENATE ALL DAYS/HOURS ---
print("Concatenating all emissions and times...")
if HOURLY:
    total_times = sum(arr.shape[0] for arr in all_xtime)
    # Stack all hourly arrays for each species
    for var in species_arrays:
        species_arrays[var] = np.concatenate(species_arrays[var], axis=0)  # shape: (Time, nCells)
    for var in scalar_std_arrays:
        scalar_std_arrays[var] = np.concatenate(scalar_std_arrays[var], axis=0)  # shape: (Time, nCells)
    xtime_all = np.concatenate(all_xtime, axis=0)  # shape: (Time, StrLen)
else:
    total_times = len(all_xtime)
    for var in species_arrays:
        species_arrays[var] = np.stack(species_arrays[var], axis=0)  # shape: (Time, nCells)
    for var in scalar_std_arrays:
        scalar_std_arrays[var] = np.stack(scalar_std_arrays[var], axis=0)  # shape: (Time, nCells)
    xtime_all = np.stack(all_xtime, axis=0)  # shape: (Time, StrLen)

print("Writing single output NetCDF file with unlimited Time dimension...")
from netCDF4 import Dataset
global_attrs = {
    "authors": "Forrest Lacey and Rajesh Kumar",
    "time_created": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
}
# Determine output filename from YAML pattern
output_file_pattern = config.get('output_file_pattern', 'FINNv2.5.1_modvrs_nrt_MOZART_{year}{month}_{mpas_grid_name}.grid_{freq}.nc')
first_file_date = min(processed_dates)
year = first_file_date.year
month = f"{first_file_date.month:02d}"
mpas_grid_name = os.path.splitext(os.path.basename(mpas_grid_file))[0]
freq = 'hourly' if HOURLY else 'daily'
out_file = os.path.join(dst_file_dir, output_file_pattern.format(year=year, month=month, mpas_grid_name=mpas_grid_name, freq=freq))
with Dataset(out_file, 'w', format='NETCDF3_64BIT') as dst:
    dst.createDimension('Time', None)  # Unlimited
    dst.createDimension('nCells', nCells)
    dst.createDimension('StrLen', STR_LEN)
    # Write emission variables
    for var, arr in species_arrays.items():
        v = dst.createVariable(var, arr.dtype.name, ('Time', 'nCells'), fill_value=np.nan)
        v[:, :] = arr
    for var, arr in scalar_std_arrays.items():
        v = dst.createVariable(var, arr.dtype.name, ('Time', 'nCells'), fill_value=np.nan)
        v[:, :] = arr
    # Add nCells variable (cell numbers 0 to nCells-1)
    nCells_var = dst.createVariable('nCells', 'i4', ('nCells',))
    nCells_var[:] = np.arange(nCells)
    # Write xtime variable
    xtime_var = dst.createVariable('xtime', 'S1', ('Time', 'StrLen'))
    xtime_var[:, :] = xtime_all
    xtime_var.long_name = "model times"
    xtime_var.calendar = "gregorian"
    xtime_var.cell_methods = "string1: mean"
    # Add global attributes
    for attr, val in global_attrs.items():
        dst.setncattr(attr, val)
print(f"Output file created: {out_file}\n")





