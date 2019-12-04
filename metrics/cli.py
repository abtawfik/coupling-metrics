from pathlib import Path
from glob import glob

import click
import xarray as xr
import dask.dataframe as dd
from toolz.curried import curry, compose



@click.command()
@click.option('--metric', help='select a coupling metric to compute')
@click.option('--files' , help='specify which files to compute metrics for. can be a wildcard')
def comet(metric, files):
    '''Used to compute various land-atmosphere coupling metrics, just pass the necessary input files

    Usage:
    comet -c 

    Parameters:
    '''
    # Get all the files from the command-line
    # Figure out the file format - make sure it is netcdf, grib, csv, hdf, or parquet
    #print(metric)
    print(xr.open_mfdataset(files, combine='by_coords'))
    #myfiles = check_files(files)
    # Pass in a function to operate on data with
    return {'metric': metric}#,
            #'files' : dd.map_partitions(apply_to_partition, myfiles).compute() }


def check_files(files):
    return compose( read_multiple_files,
                    check_if_there_are_multiple_extensions,
                    check_if_files_exist )(files)


def apply_to_partition(df):
    return df['MIXR'].sum() * 1e-3


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
def read_multiple_files(files, args):
    extension = Path(files[0]).suffix
    if extension == '.csv':
        nparts = len(files)
        return dd.concat( [dd.read_csv(f, sep=' ').astype(float) for f in files] ).repartition(npartitions=nparts)
    elif extension in ['.nc', '.nc4', '.netcdf', '.netcdf4', '.cdf']:
        ds = xr.open_mfdataset(files, combine='by_coords')
        output = MixingDiagram().compute(ds,
                                         temp=args['temp'],
                                         pblh=args['hpbl'],
                                         psfc=args['pres'],
                                         sh_flux=args['shtfl'],
                                         lh_flux=args['lhtfl'],
                                         qhum=args['shum'])
        return 
    else:
        raise

    
if __name__ == '__main__':
    comet()
