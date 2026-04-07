emis_file_pattern: 'GLOB_MOZ4_{year}{julianday:03d}.txt'# -----------------------------------------------------------------------------
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

# NetCDF output format for MPASv8+ compiled with SMIOL library
output_format = 'NETCDF3_64BIT_DATA' # NetCDF format (CDF-5) w/ 64-bit array
                                     # NETCDF3_64BIT=NETCDF3_64BIT_OFFSET (CDF-2) w/ 4GB limit

with open('config_finn_to_mpas.yaml', 'r') as f:
    config = yaml.safe_load(f)


# Get year and month from start_date and end_date
from datetime import datetime, timedelta
start_date = config.get('start_date', None)
end_date = config.get('end_date', None)
if not start_date or not end_date:
    raise ValueError('start_date and end_date must be specified in the config file.')
start_dt = datetime.strptime(start_date, '%Y-%m-%d')
end_dt = datetime.strptime(end_date, '%Y-%m-%d')

emis_file_date_pairs = []
from dateutil.rrule import rrule, DAILY
# Build list of all dates between start and end date (inclusive)
all_dates = list(rrule(DAILY, dtstart=start_dt, until=end_dt))

emis_file_date_pairs = []
for date in all_dates:
    year = date.year
    month = date.month
    day = date.day
    julianday = (date - datetime(year, 1, 1)).days + 1
    emis_dir = config['emis_dir'].format(year=year)
    # Support both {month:02d} and {julianday:03d} in pattern
    emis_file_pattern = config['emis_file_pattern'].format(year=year, month=f"{month:02d}", day=f"{day:02d}", julianday=julianday)
    for f in glob.glob(os.path.join(emis_dir, emis_file_pattern)):
        # Try to extract date from filename (YYYYMMDD or YYYYJJJ)
        m_ymd = re.search(r'(\d{8})', os.path.basename(f))
        m_jul = re.search(r'(\d{7})', os.path.basename(f))
        file_date = None
        if m_ymd:
            try:
                file_date = datetime.strptime(m_ymd.group(1), '%Y%m%d')
            except Exception:
                pass
        elif m_jul:
            try:
                file_date = datetime.strptime(m_jul.group(1), '%Y%j')
            except Exception:
                pass
        if file_date and start_dt <= file_date <= end_dt:
            emis_file_date_pairs.append((file_date, f))
# Sort by file_date
emis_file_list = [f for file_date, f in sorted(emis_file_date_pairs)]
dst_file_dir = config['dst_file_dir']  # Output directory for NetCDF files
mpas_grid_file = config['mpas_grid_file']  # Path to MPAS grid NetCDF file
HOURLY = config.get('HOURLY', False)  # Whether to output hourly emissions
lt_fac = np.array(config.get('lt_fac', [.43, .43, .43, .43, .43, .43, .43, .43, .43, 3., 6., 10., 14., 17., 14., 12., 9., 6., 3., .43, .43, .43, .43, .43]))  # Diurnal profile

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

def make_xtime_array(in_filename: str, hourly: bool, strlen: int = 64):
    # Try YYYYMMDD first, then YYYYJJJ (Julian day)
    m_ymd = re.search(r'(\d{8})', Path(in_filename).name)
    m_jul = re.search(r'(\d{7})', Path(in_filename).name)
    dt = None
    if m_ymd:
        date_str = m_ymd.group(1)
        dt = datetime.strptime(date_str, "%Y%m%d")
    elif m_jul:
        date_str = m_jul.group(1)
        dt = datetime.strptime(date_str, "%Y%j")
    else:
        raise ValueError("No YYYYMMDD or YYYYJJJ (Julian day) date found in input filename.")
    if hourly:
        ntimes = 24
        timestrs = [f"{dt.strftime('%Y-%m-%d')}_{str(h).zfill(2)}:00:00" for h in range(24)]
    else:
        ntimes = 1
        timestrs = [dt.strftime('%Y-%m-%d_00:00:00')]
    xtime_arr = np.zeros((ntimes, strlen), dtype='S1')
    for i, s in enumerate(timestrs):
        s_padded = s.ljust(strlen, '\0')
        s_bytes = np.array(list(s_padded), dtype='S1')
        xtime_arr[i, :] = s_bytes
    return xtime_arr




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
species_arrays = {species_map[col]: [] for col in species_to_map if col in species_map and species_map[col] != "co2_biob_modis"}

for csv_file in emis_file_list:
    print(f"Processing FINN file: {csv_file}")
    # Modified read statement to accomodate FINNv1 fortran format strings FGL
    df = pd.read_csv(csv_file, header=0, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    for col in df.columns:
        df[col] = (df[col].astype(str).str.replace('D', 'E', regex=False).str.replace('d', 'E', regex=False))
        df[col] = pd.to_numeric(df[col], errors='ignore')
    # Compute local hour offset for each cell (in hours, rounded to nearest int)
    local_hour_offset = np.round(lonCell_deg / 15.).astype(int) % 24
    lati = df["LATI"].values
    longi = df["LONGI"].values
    # Select all columns in df that are in species_to_map, regardless of position
    emission_cols = [col for col in df.columns if col in species_to_map]
    emission_per_cell = {col: np.zeros(nCells) for col in emission_cols}
    # For scalar type, also track count per cell
    scalar_count_per_cell = {col: np.zeros(nCells, dtype=int) for col in emission_cols if species_type.get(col, None) == "scalar"}
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
            tmp = df.iloc[i][col]
            stype = species_type.get(col, None)
            if stype == "scalar":
                emission_per_cell[col][min_index] += float(tmp)
                scalar_count_per_cell[col][min_index] += 1
            else:
                emission_per_cell[col][min_index] += float(tmp)
    # For each emission species, convert units and apply diurnal profile
    for col, values in emission_per_cell.items():
        if col not in species_map:
            continue
        var_name = species_map[col]
        stype = species_type.get(col, None)
        values_conv = np.zeros_like(values, dtype=float)
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
            # Average for scalar type
            count = scalar_count_per_cell.get(col, np.ones_like(values))
            with np.errstate(divide='ignore', invalid='ignore'):
                values_conv = np.true_divide(values, count)
                values_conv[count == 0] = 0.0
        else:
            values_conv = values

        if HOURLY:
            wrk_lt_fac = np.array(lt_fac) / np.sum(lt_fac) * 24.
            hourly_values = np.zeros((24, values_conv.shape[0]))
            for cell in range(values_conv.shape[0]):
                offset = local_hour_offset[cell]
                shifted_fac = np.roll(wrk_lt_fac, -offset)
                hourly_values[:, cell] = values_conv[cell] * shifted_fac
            species_arrays[var_name].append(hourly_values)
        else:
            species_arrays[var_name].append(values_conv)
    # Add xtime for this file
    xtime_arr = make_xtime_array(csv_file, HOURLY, STR_LEN)
    all_xtime.append(xtime_arr)

# --- CONCATENATE ALL DAYS/HOURS ---
print("Concatenating all emissions and times...")
if HOURLY:
    total_times = sum(arr.shape[0] for arr in all_xtime)
    # Stack all hourly arrays for each species
    for var in species_arrays:
        species_arrays[var] = np.concatenate(species_arrays[var], axis=0)  # shape: (Time, nCells)
    xtime_all = np.concatenate(all_xtime, axis=0)  # shape: (Time, StrLen)
else:
    total_times = len(all_xtime)
    for var in species_arrays:
        species_arrays[var] = np.stack(species_arrays[var], axis=0)  # shape: (Time, nCells)
    xtime_all = np.stack(all_xtime, axis=0)  # shape: (Time, StrLen)

print("Writing single output NetCDF file with unlimited Time dimension...")
from netCDF4 import Dataset
global_attrs = {
    "authors": "Forrest Lacey and Rajesh Kumar",
    "time_created": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
}
# Determine output filename from YAML pattern
output_file_pattern = config.get('output_file_pattern', 'FINNv2.5.1_modvrs_nrt_MOZART_{year}{month}_{mpas_grid_name}.grid_{freq}.nc')
first_file_date = emis_file_date_pairs[0][0]
year = first_file_date.year
month = f"{first_file_date.month:02d}"
mpas_grid_name = os.path.splitext(os.path.basename(mpas_grid_file))[0]
freq = 'hourly' if HOURLY else 'daily'
out_file = os.path.join(dst_file_dir, output_file_pattern.format(year=year, month=month, mpas_grid_name=mpas_grid_name, freq=freq))

with Dataset(out_file, 'w', format=output_format) as dst:
    dst.createDimension('Time', None)  # Unlimited
    dst.createDimension('nCells', nCells)
    dst.createDimension('StrLen', STR_LEN)
    # Write emission variables
    for var, arr in species_arrays.items():
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





