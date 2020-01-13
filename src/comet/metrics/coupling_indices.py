
######################
# 3rd party packages #
######################
import numpy as np
import xarray as xr
from toolz.curried import compose, curry

####################
# Package specific #
####################
from . import utils as gf



class CouplingIndex(object):
    '''
    Contains general equation for computing terrestrial and atmospheric 
    legs of the coupling index
    
    - Specifc Descriptions -
    
    Terrestrial coupling index
    --------------------------
    
    Equation :  TC = stddev(soil moisture) * slope(soil moisture, surface flux)
                                        OR
                TC = covariance(soil moisture, surface flux) / stddev(soil moisture)
    
    These two equations are equivalent but the 2nd one tends to be more stable
    
    Coupling index can likely be improved because this metric assumes a linear
    relationship between each variable. In the case a surface flux refers to sensible, latent,
    or total flux from the surface to the atmosphere
    
    
    Atmospheric coupling index
    --------------------------
    
    Equation: same as terrestrial but replace a soil moisture with a preferred surface flux
    and the surface fluxes with some characteristic property of the atmosphere such as boundary
    layer height, column integrated moisture content, or TKE depending on your use-case.
    
              AC = cov(surface flux, boundary layer height) / std(surface flux)
    
    General Description
    -------------------
    Physically this parameter quantifies to what extend variations in soil moisture
    correspond to change variations in surface fluxes.  If the coupling parameter
    is high then soil moisture is said to vary sufficiently to influence 
    variations in surface fluxes. Or in the atmospheric context, surface fluxes are 
    said to vary enough to influence variations in boundary layer height.
    
    One strong advantage of the coupling index is that the returned coupling index has the
    same units as the dependent variable. In the case of TC for soil moisture and latent heat
    flux the units would be in W/m2, same as latent heat flux
    
    Reference
    ---------
    Dirmeyer, The terrestrial segment of soil moisture–climate coupling (2011)
    '''
    def __init__(self):
        return None

    def compute(self,
                ds,
                xname : str,
                yname : str,
                averaging='season'):
        '''
        Compute the coupling index using 1-day averaged data.
        It is up to the user to know which variable correspond to the x and y variables

        In general for terrestrial coupling `x` will be something like soil moisture at a certain depth
        and `y` will be latent or sensible heat flux.

        Parameters
        ----------
        ds : xarray.Dataset
           an xarray dataset containing all the variable names

        xname : string
           name of the independent variable used in coupling context
        
        yname : string
           name of the dependent variable used in the coupling context

        averaging : string
           averaging time window to use when computing coupling index
           options - 'season' and 'month'

        Return
        ------
        xarray.Dataset
           returns the cooupling index between each variable across the time
           dimension. This is computed either for each month or season. Also
           The return variable can be 'lazy' or eagerly evaluted depending on 
           what type of xarray Dataset was passed.
           Output variable is returned with same units as `yname` variable
        '''
        #-----------------------------------
        # Check the incoming ds for correct
        # structure and type
        #-----------------------------------
        gf.check_input_is_xarray(ds)
        gf.check_variable_names(ds, [xname, yname] )
        #-----------------------------------
        # Get time axis name
        #-----------------------------------
        time_dim = gf.get_time_axis_name(ds)
        dim  = f'{time_dim}.{averaging}'
        #-----------------------------------
        # Make sure data are daily averages
        #-----------------------------------
        xday = ds[xname].resample({time_dim:'1D'}).mean().groupby(dim)
        yday = ds[yname].resample({time_dim:'1D'}).mean().groupby(dim)
        coupling = gf.covariance(xday,yday,dim) / xday.std()
        return coupling.rename(f'{xname}_{yname}_CI').to_dataset()


    


