
# FINN to MPAS Emissions Regridding Workflow

## 1. Overview
`Regrid_FINN_to_MPAS.py` automates the mapping of selected FINN fire emissions species to an MPAS grid, producing MPAS-compatible NetCDF files with correct units, diurnal profiles, and metadata.

## 2. Features
- Maps only user-specified FINN species to MPAS grid cells (nearest-cell assignment)
- Converts emissions to molecules/cm²/s (supports both aerosols and gases)
- Applies location-dependent diurnal profiles
- Outputs MPAS-compatible NetCDF files with proper variable names, units, and global attributes
- Fully configurable via a YAML file

## 3. Python Environment Setup
**NOTE:** If you already have the required Python libraries installed, you can skip this section and go to Section 4.

### Required Python Libraries
The workflow requires the following Python libraries:
- numpy
- pandas
- xarray
- netCDF4
- pyyaml
- pathlib

### Conda Environment Setup
To create a new conda environment and install all required libraries, run:

```bash
conda create -n finn2mpas python=3.8 numpy pandas xarray netCDF4 pyyaml
conda activate finn2mpas
```
If you need `pathlib` (usually included with Python 3.4+), it will be available by default.

## 4. YAML Configuration (`config_finn_to_mpas.yaml`)

Key entries:
- `emis_dir`: Directory containing FINN CSV files
- `emis_file_pattern`: Search pattern for FINN files (e.g., `*.csv`)
- `dst_file_dir`: Output directory for NetCDF files
- `mpas_grid_file`: Path to MPAS grid NetCDF file
- `lt_fac`: Diurnal profile array (length 24, sums to 100)
- `species_map`: Mapping of FINN species names to MPAS output variable names
- `species_type`: Specifies if each species is 'aerosol' or 'gas' (affects unit conversion)
- `HOURLY`: Set to `true` to produce hourly output files, or `false` for daily output files. The code will generate either daily or hourly emissions files depending on this flag.
- `output_format`, `compression`: NetCDF output options (Don't change because MPAS does not work with netCDF-4)

## Overview
`Regrid_FINN_to_MPAS.py` automates the mapping of selected FINN fire emissions species to an MPAS grid, producing MPAS-compatible NetCDF files with correct units, diurnal profiles, and metadata.

## Features
- Maps only user-specified FINN species to MPAS grid cells (nearest-cell assignment)
- Converts emissions to molecules/cm²/s (supports both aerosols and gases)
- Applies location-dependent diurnal profiles
- Outputs MPAS-compatible NetCDF files with proper variable names, units, and global attributes
- Fully configurable via a YAML file

## YAML Configuration (`config_finn_to_mpas.yaml`)
Key entries:
- `emis_dir`: Directory containing FINN CSV files
- `emis_file_pattern`: Glob pattern for FINN files (e.g., `*.csv`)
- `dst_file_dir`: Output directory for NetCDF files
- `mpas_grid_file`: Path to MPAS grid NetCDF file
- `lt_fac`: Diurnal profile array (length 24, sums to 100)
- `species_map`: Mapping of FINN species names to MPAS output variable names
- `species_type`: Specifies if each species is 'aerosol' or 'gas' (affects unit conversion)
- `output_format`, `compression`: NetCDF output options

Example:
```yaml
emis_dir: /path/to/emis_dir
emis_file_pattern: '*.csv'
dst_file_dir: /path/to/output_dir
mpas_grid_file: /path/to/mpas_grid.nc
lt_fac: [0.43, 0.43, ..., 3.0, ..., 0.43]  # 24 values, sum to 100
species_map:
  BC: bc_biob_modis
  OC: oc_biob_modis
  CO: co_biob_modis
  NH3: nh3_biob_modis
  SO2: so2_biob_modis
  ISOP: iso_biob_modis
  BIGENE: mnt_biob_modis
species_type:
  BC: aerosol
  OC: aerosol
  CO: gas
  NH3: gas
  SO2: gas
  ISOP: gas
  BIGENE: gas
output_format: NETCDF3_64BIT
compression:
  zlib: true
  complevel: 4
  shuffle: true
```


## 5. How to Run
1. Edit `config_finn_to_mpas.yaml` to match your input/output paths and desired species mapping.
2. Run the script:
```bash
python Regrid_FINN_to_MPAS.py config_finn_to_mpas.yaml
```

## 6. Output
- NetCDF files in the output directory, with variables named per `species_map`, units set to `molecules/cm2/s`, and global attributes for authors and creation time.

## 7. Support
For questions or issues, contact:
- Forrest Lacey (lacey@ucar.edu)
- Rajesh Kumar (rkumar@ucar.edu)
