import netCDF4
import numpy as np
import os
import time

def postprocess_netcdf(src_path, dst_path, year, spc, sector_exclude=None):
    """Copy variables from *src_path* to *dst_path* renaming them
    and optionally subtracting the contributions of a few sectors from the
    ``sum`` variable before the file is replaced.

    Parameters
    ----------
    src_path : str
        Path to the ESMF-regridded netCDF file produced earlier.
    dst_path : str
        Temporary output path; the final file replaces the source on success.
    year : int
        Year used to construct the ``xtime`` variable.
    spc : str
        Species name (used to prefix variables).
    sector_exclude : list of str, optional
        Variables in the source file whose values should be subtracted from
        the ``sum`` variable.  If ``None`` or empty, no subtraction occurs.
    """
    if sector_exclude is None:
        sector_exclude = []

    exclude_vars = {'time', 'lon', 'lat', 'area', 'ncol'}
    # Open source file
    with netCDF4.Dataset(src_path, 'r') as src:
        # read any excluded sector arrays into memory in case they need to be
        # subtracted from ``sum`` later.  we'll create a dictionary so we can
        # access them by name.
        excluded_data = {}
        for sec in sector_exclude:
            if sec in src.variables:
                excluded_data[sec] = src.variables[sec][:]

        # Create new destination file
        with netCDF4.Dataset(dst_path, 'w', format='NETCDF3_64BIT_DATA') as dst:
            # Dimensions
            # Rename time to Time, ncol to nCells
            time_dim = 'Time' if 'time' in src.dimensions else 'Time'
            ncol_dim = 'nCells' if 'ncol' in src.dimensions else 'nCells'
            dst.createDimension(time_dim, src.dimensions.get('time', src.dimensions.get('Time', 12)).size if not src.dimensions.get('time', src.dimensions.get('Time', None)).isunlimited() else None)
            dst.createDimension(ncol_dim, src.dimensions.get('ncol', src.dimensions.get('nCells', None)).size)
            dst.createDimension('StrLen', 64)
            # xtime variable
            xtime_var = dst.createVariable('xtime', 'S1', (time_dim, 'StrLen'))
            months = [f"{year}-{str(m).zfill(2)}-01_00:00:00" for m in range(1,13)]
            xtime_arr = np.zeros((12, 64), dtype='S1')
            for i, s in enumerate(months):
                s_bytes = np.array(list(s), dtype='S1')
                xtime_arr[i, :len(s_bytes)] = s_bytes
            xtime_var[:, :] = xtime_arr
            xtime_var.long_name = "model times"
            xtime_var.calendar = "gregorian"
            xtime_var.cell_methods = "string1: mean"
            # Unified variable copying logic
            for name, var in src.variables.items():
                if name in exclude_vars:
                    continue

                # if this is the ``sum`` variable and we have exclusions, modify it
                if name == 'sum' and excluded_data:
                    # start with a copy of the original array
                    data = var[:].astype('f4').copy()
                    for sec, arr in excluded_data.items():
                        # only subtract if the shapes match
                        if arr.shape == data.shape:
                            data -= arr
                        else:
                            # shapes mismatch; warn the user but continue
                            print(f"warning: cannot subtract sector '{sec}' from sum "
                                  "(shape mismatch)")
                else:
                    data = var[:]

                if spc == 'isoprene':
                    new_name = f'iso_anth_{name}'
                elif spc == 'monoterpenes':
                    new_name = f'mnt_anth_{name}'
                else:
                    new_name = f'{spc}_anth_{name}'
                dims = [
                    time_dim if d in ['time', 'Time'] else ncol_dim if d in ['ncol', 'nCells'] else d
                    for d in var.dimensions
                ]
                new_var = dst.createVariable(new_name, 'f4', tuple(dims), fill_value=np.nan)
                new_var[:] = data
            # Global attributes
            dst.setncattr('authors', 'Forrest Lacey and Rajesh Kumar')
            dst.setncattr('time_created', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    # Move temporary file to destination file
    os.replace(dst_path, src_path)
# Example usage:
# postprocess_netcdf('input.nc', 'output.nc', 2024, 'monoterpenes')
# postprocess_netcdf('input.nc', 'output.nc', 2024, 'isoprene')
