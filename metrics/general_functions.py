'''
A set of functions that are not specific to any specific metric
Functions such as calculating saturation specific humidity or air density
'''

# 3rd party packages
import numpy as np
import pandas as pd
from toolz.curried import curry, compose

import constants as cnts
import xarray as xr

@curry
def air_density(p, t, q):
    '''Calculate the air density using the ideal gas law
    with the specific humidity of gas
    
    Parameters
    ----------
    p : float
        pressure [Pascals]
    t : float
        dry bulb temperature [Kelvin]
    q : float
        specific humidity [kg/kg]

    Return
    ------
    float
        air density [kg/m^3]
    '''
    ep = cnts.ep.value
    Rd = cnts.Rd.value
    return p / (Rd * t * ((1. + (q/ep)) / (1. + q)))

@curry
def temp_to_energy(t):
    '''Calculate the amount of heat energy in an air parcel
    using the specific heat capacity
    
    Parameters
    ----------
    t : float
        dry bulb temperature [Kelvin]

    Return
    ------
    float
        heat energy [J/kg]
    ''' 
    return cnts.cp.value * t


@curry
def shum_to_latent_energy(q):
    '''Calculate the amount of latent heat energy in an air parcel
    using the specific heat capacity

    Parameters
    ----------
    q : float
        specific humidity [kg/kg]
    
    Return
    ------
    float
        latent heat energy [J/kg]
    '''
    return cnts.Lv.value * q


@curry
def bowen_ratio(sensible_heat, latent_heat):
    '''Calculate the bowen ratio defined as:
    bowen = sensible heat flux / latent heat flux

    Parameters
    ----------
    sensible_heat : float
        sensible heat flux [W/m^2]

    latent_heat : float
        latent heat flux [W/m^2]

    Return
    ------
    float
        bowen ratio [unitless]
    '''
    return (sensible_heat / latent_heat).where(latent_heat != 0, 0.0)



@curry
def add_heat_to_layer(rho, heat_flux, depth, time_change):
    '''Calculate the amount of heat in a column of air over a defined time period

    Parameters
    ----------
    heat_flux : float
        sensible heat flux [W/m^2]

    rho : float
        air density [kg/m^3]
	
    depth : float
        depth of the air column to add the heat to [m]

    time_change : float
        time over which the heat was added [s]

    Return
    ------
    float
        total heat energy per unit mass [J/kg]
    '''
    return ((heat_flux * time_change) / (rho * depth)).where( rho*depth > 0, np.nan)


@curry
def check_variable_names(ds, variable_names):
    # Loop over the variable names and throw error if variable is not in dataset
    for name in variable_names:
        if not name in ds:
            raise KeyError(f'variable {name} is not in your dataset')


@curry
def to_xarray(template, data, input_name, output_name):
    return template.to_dataset().copy({input_name:data}).rename_vars({input_name:output_name})


@curry
def check_input_is_xarray(ds):
    if not isinstance(ds, xr.Dataset):
        raise TypeError(f'Input variable, ds, must be an xarray dataset type')

@curry
def check_time_is_a_dimension(ds):
    if 'time' not in list(ds.dims):
        raise TypeError(f'Xarray must have a "time" dimension. Current dimension names: {list(ds.dims)}')
        
@curry
def get_interval(ds):
	dt = np.unique(ds.time.diff(dim='time').values.astype('timedelta64[s]'))
	if dt.shape[0] > 1:
		raise ValueError(f'Time dimension has irregular time interval...must be constant interval')
	return dt[0].astype(float)