#!/usr/bin/env python

import numpy as np
import xarray as xr
import esmpy as ESMF  
import glob
import datetime
import os
import netCDF4
import yaml
import argparse
import sys

# Set up the argument parser
parser = argparse.ArgumentParser(description="Start UPTEMPO to Regrid Emissions")
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


from datetime import datetime
from Calc_Emis import Calc_Emis_T
from Regridding_ESMF import Add_bounds, Regridding
from postprocess_cams_anth_to_mpas import postprocess_netcdf

### Directory and filename setup (should only need to edit this section)
# === User input for year ===
year = config['year']
# Input file format and list for CAMS anthropogenic emissions (CAMS-GLOB-ANT)
inp_file_pattern = config['inp_file_format'].format(year=year).replace('SPC', '*')
CAMS_GLOB_ANT_file_list = glob.glob(inp_file_pattern)
# Destination filename format that regridded field will be saved
# "SPC" will be replaced to the real species name from CAMSv4.2 species
dst_file_format = config['dst_file_format'].format(year=year)
# SE-RR grid file (or any destination grid file)
SERR_scrip_file = config['SERR_scrip_file']
# species mapping: input name to output name
species_map = {}
for item in config['species']:
    if isinstance(item, dict):
        # If already a dict (future-proofing)
        species_map.update(item)
    elif isinstance(item, str) and ':' in item:
        k, v = [s.strip() for s in item.split(':', 1)]
        species_map[k] = v
    else:
        # fallback: treat as identity mapping
        species_map[item] = item
print('species_map:', species_map)
# uncomment below if you want to process only some fields/sectors of the file
#sectors = ['sum']
# uncomment below if you want to process all the fields/sectors in the file
sectors = config.get('sectors', [])
print(sectors)
# --- optional list of sectors to subtract from the "sum" sector during postprocessing
#     this list is read from the yaml file and passed through to the helper
sector_exclude = config.get('sector_exclude', [])
print(f"sector_exclude: {sector_exclude}")
# Exsiting grid information file for base emissions
CAMS_grid_file = config['CAMS_grid_file']
# Regridding weight file name and logical if a new one needs to be created
Regridding_weights_file = config['Regridding_weights_file']
need_weights = not os.path.exists(Regridding_weights_file)
# Input file format with SPC placeholder
inp_file_format = config['inp_file_format']
### Directory and filename setup section end (NO updates needed after this point)
# ## Make grid information NetCDF file (if "lon_bnds" and "lat_bnds" are not available for FV grid input)
def grid_file_has_bounds(grid_file):
    if not os.path.exists(grid_file):
        return False
    try:
        with netCDF4.Dataset(grid_file, 'r') as ds:
            return ('lat_bnds' in ds.variables) and ('lon_bnds' in ds.variables)
    except Exception as e:
        print(f"Error reading {grid_file}: {e}")
        return False

need_grid_bounds = not grid_file_has_bounds(CAMS_grid_file)
if need_grid_bounds:
    Add_bounds(CAMS_GLOB_ANT_file_list[0], CAMS_grid_file, creation_date=False)
    print('--- ran Add_bounds')
else:
    print(f'--- using existing grid file with bounds: {CAMS_grid_file}')
# ## Create weight file if you don't have one already. This is extremely useful especially when you have more than two files to be processed
if need_weights:
    ds_CAMS = xr.open_dataset(CAMS_GLOB_ANT_file_list[0], decode_times=False)
    print('--- begin regridding')
    #For MPAS:
    rr = Regridding( ds_CAMS.isel(time=slice(0,1)), src_grid_file=CAMS_grid_file, dst_grid_file=SERR_scrip_file, 
                   wgt_file=Regridding_weights_file, method='Conserve', save_wgt_file=True, save_results=False,
                  save_wgt_file_only=True, check_timings=True, creation_date=False )
    print('--- end regridding')
else: 
    print('--- will use existing weighting file ', Regridding_weights_file)
for input_name, output_name in species_map.items():
    print( f'Regridding: {input_name} → {output_name}', datetime.now() )
    print( '************************************************************************' )
    ds_emis = xr.open_dataset(inp_file_format.format(year=year).replace('SPC', input_name), decode_times=False)
    dst_file_base = dst_file_format.replace('SPC', output_name)
    # Only create the output file with the date suffix (no undated file)
    rr = Regridding( ds_emis, src_grid_file=CAMS_grid_file, dst_grid_file=SERR_scrip_file,
                     wgt_file=Regridding_weights_file, method='Conserve', fields=sectors,
                     dst_file=dst_file_base, save_wgt_file=False, save_results=True, check_results=False,
                     check_timings=True, creation_date=True, nc_file_format='NETCDF3_64BIT_DATA' )
    # --- Post-process output file: rename emission variables with spc_anth_ prefix, rename ncol dimension, add xtime ---
    import netCDF4
    import numpy as np
    # Find the dated output file
    import glob
    dst_file_pattern = dst_file_base.replace('.nc', '_c*.nc')
    matching_files = glob.glob(dst_file_pattern)
    if matching_files:
        dst_file = max(matching_files, key=os.path.getctime)
        # Post-process the output file
        tmp_file = dst_file + '.tmp'
        postprocess_netcdf(dst_file, tmp_file, year, output_name, sector_exclude)
