
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
from comet.metrics.fortran.ctp import conv_trig_pot_mod as ctp_hi_low_fortran
ctp_fortran = ctp_hi_low_fortran.ctp_hi_low


class ConvTrig(object):
    '''  
    Calculates the convective triggering potential (CTP) given early
    morning sounding profiles and also includes a low-level humidity index (Hi_low).
    The CTP returns useful information regarding which profile is primed for convection
    under high surface senisble heat flux and which profiles are likely to trigger
    convection under enhanced latent heat flux. 
    Note that this does not apply to established boundary layer profiles because 
    assumptions are made regarding the presences of inversion.

    CTP     =  integral of curve between the moist adiabat and environmental lapse
               rate from 100 mb above the ground to 300mb above the ground
                _       _                _       _
               |         |              |         |
    Hi_low  =  | T - T_d |          -   | T - T_d |
               |_       _|50mb abg      |_       _|150mb abg
 
               where T_d = dew point temperature [K] ;  T = air temperature [K]
                     abg = above ground

 
    References:
    -----------
    Findell, K.L., Eltahir, E.A.B. 2003: Atmospheric Controls on 
    Soil Moisture-Boundary Layer Interactions. Part I: Framework
    Development. Journal of Hydrometeorology

    Findell, K.L., Eltahir, E.A.B. 2003: Atmospheric Controls on
    Soil Moisture-Boundary Layer Interactions. Part II: Feedbacks
    within the Continental United States. Journal of Hydrometeorology

    Craig R. Ferguson and Eric F. Wood, 2011: Observed Land–Atmosphere
    Coupling from Satellite Remote Sensing and Reanalysis. J. Hydrometeor,

    '''
    def __init__(self):
        return None

    def compute(self,
                ds,
                temp      : str,
                shum      : str,
                pres      : str,
                level_dim : str,
                t2m       = None,
                q2m       = None,
                psfc      = None,
                averaging ='raw'):
        '''
        Calculates the CTP and HiLow for a given atmospheric profile

        Parameters
        ----------
        ds : xarray.Dataset
           an xarray dataset containing all the variable names

        level_dim : string
           name of the level dimension.

        temp : string
           name of the air temperature profile variable. must have a vertical dimension [K]

        shum : string
           name of the specific humidity profile variable. must have a vertical dimension [kg/kg]

        pres : string
           name of the pressure profile variable. can be an array if constant dimension or variable as
           would be the case for sigma or hybrid sigma coordinates. In the case of sigma or hybrid coordinates
           pressure must have the same dimensions as temp and shum [Pa]

        t2m : string
           name of the 2-meter temperature [K] (optional)
        
        q2m : string
           name of the 2-meter specific humidity [kg/kg] (optional)

        psfc : string
           name of the surface pressure [Pa] (optional)

        averaging : string
           averaging time window to use when computing coupling index
           options - 'raw', 'season', and 'month'

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
        gf.check_variable_names(ds, [temp, psfc, shum, pres, t2m, q2m] )
        #-----------------------------------
        # Name of averaging window
        #-----------------------------------
        time_dim = gf.get_time_axis_name(ds)
        avg_dim  = None if averaging is None else f'{time_dim}.{averaging}'
        #-----------------------------------
        # Make sure data are daily averages
        #-----------------------------------
        xday = ds[xname].resample({time_dim:'1D'}).mean().groupby(dim)
        yday = ds[yname].resample({time_dim:'1D'}).mean().groupby(dim)
        coupling = gf.covariance(xday,yday,dim) / xday.std()
        return coupling.rename(f'{xname}_{yname}_CI').to_dataset()


    


