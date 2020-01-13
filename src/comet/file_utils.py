

########################
## Standard Libraries ##
########################
from pathlib import Path

#######################
## External Packages ##
#######################
import xarray as xr
import dask.dataframe as dd
from toolz.curried import compose, curry




def get_data(files):
    return compose(read_data_files, check_files)(files)


@curry
def check_files(files):
    return compose( check_if_there_are_multiple_extensions,
                    check_if_files_exist )(files)


@curry
def check_if_files_exist(files):
    missing = [f for f in files if not Path(f).exists()]
    if len(missing) > 0:
        raise FileNotFoundError(f'ERROR: Could not find {",".join(missing)}')
    return files


@curry
def check_if_there_are_multiple_extensions(files):
    extensions = [Path(f).suffix for f in files]
    if len(set(extensions)) > 1:
        raise ValueError(f'ERROR: You cannot pass multiple file formats.')
    return files
    

@curry
def read_data_files(files):
    extension = Path(files[0]).suffix
    if extension == '.csv':
        nparts = len(files)
        return dd.concat( [dd.read_csv(f, sep=' ').astype(float)
                             for f in files] ).repartition(npartitions=nparts)
    
    elif extension in ['.nc', '.nc4', '.netcdf', '.netcdf4', '.cdf']:
        return xr.open_mfdataset(files, combine='by_coords')
    
    elif extension in ['.grb', '.grib', '.grib2', '.grb2']:
        return xr.open_mfdataset(files, combine='by_coords', engine='cfgrib')
    
    else:
        raise TypeError('Files are not the correct type... need to be netcdf, grib, or csv')

