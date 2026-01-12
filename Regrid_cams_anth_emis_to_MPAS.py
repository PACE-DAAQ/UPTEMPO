#!/usr/bin/env python

import numpy as np
import xarray as xr
import esmpy as ESMF  
import glob
import datetime
import os
import netCDF4
import yaml

with open('config_cams_anth_regrid.yaml', 'r') as f:
    config = yaml.safe_load(f)

from datetime import datetime
from Calc_Emis import Calc_Emis_T
from Regridding_ESMF import Add_bounds, Regridding
from postprocess_cams_anth_to_mpas import postprocess_netcdf

### Directory and filename setup (should only need to edit this section)
# === User input for year ===
year = config['year']
# CAMS-GLOB-ANT_Glb directory
CAMS_GLOB_ANT_dir = config['CAMS_GLOB_ANT_dir']
CAMS_GLOB_ANT_file_list = glob.glob( CAMS_GLOB_ANT_dir + f'CAMS-GLOB-ANT*monthly_{year}.nc')
# Destination filename format that regridded field will be saved
# "SPC" will be replaced to the real species name from CAMSv4.2 species
dst_file_format = config['dst_file_format'].format(year=year)
# SE-RR grid file (or any destination grid file)
SERR_scrip_file = config['SERR_scrip_file']
# list of the species that will processed fron the input emission files
species = config['species']
print(species)
# uncomment below if you want to process only some fields/sectors of the file
#sectors = ['sum']
# uncomment below if you want to process all the fields/sectors in the file
sectors = config.get('sectors', [])
print(sectors)
# Exsiting grid information file for base emissions
CAMS_grid_file = config['CAMS_grid_file']
# Regridding weight file name and logical if a new one needs to be created
Regridding_weights_file = config['Regridding_weights_file']
need_weights = not os.path.exists(Regridding_weights_file)
# Input file format with SPC placeholder
inp_file_format = config['inp_file_format']
### Directory and filename setup section end (NO updates needed after this point)
# ## Make grid information NetCDF file (if "lon_bnds" and "lat_bnds" are not available for FV grid input)
Add_bounds( CAMS_GLOB_ANT_file_list[0], CAMS_grid_file, creation_date=False )
print('--- end Add_bounds')
# ## Create weight file if you don't have one already. This is extremely useful especially when you have more than two files to be processed
if need_weights:
    ds_CAMS = xr.open_dataset(CAMS_GLOB_ANT_file_list[0])
    print('--- begin regridding')
    #For MPAS:
    rr = Regridding( ds_CAMS.isel(time=slice(0,1)), src_grid_file=CAMS_grid_file, dst_grid_file=SERR_scrip_file, 
                   wgt_file=Regridding_weights_file, method='Conserve', save_wgt_file=True, save_results=False,
                  save_wgt_file_only=True, check_timings=True, creation_date=False )
    print('--- end regridding')
else: 
    print('--- will use existing weighting file ', Regridding_weights_file)
for sp1 in species:
    print( 'Regridding: ', sp1, datetime.now() )
    print( '************************************************************************' )
    #ds_emis = xr.open_dataset( '/glade/work/sshams/CAMS_glob_ant/monthly/CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_co_v5.3_monthly_2000.nc' )  
    ds_emis = xr.open_dataset(inp_file_format.format(year=year).replace('SPC', sp1))
    dst_file_base = dst_file_format.replace('SPC', sp1)
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
        postprocess_netcdf(dst_file, tmp_file, year, sp1)
