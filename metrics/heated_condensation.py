'''
Quantifies the how primed the atmosphere is to triggering moist convection.
This does not separate between shallow versus deep convection, however, or 
include non-local dynamic triggers.


Description
-----------
Provides quantities for 1) the potential temperature required to initiate moist convection
as defined by the buoyant condensation level, 2) the height and pressure the PBL needs to
reach to trigger moist convection, and 3) the potential temperature deficit between the
current 2-m potential temperature and potential temp threshold.
This method can be used as a condition for convective triggering (i.e. trigger when TDEF=0)
Note the current formulation does not distinguish between shallow and deep convection
but simply points to convective initiation.  This technique can also be used to identify
locally triggered convection versus transient events using a single profile, where 
locally triggered == where TDEF=0 and otherwise observed cloud cover is non-local

Variables returned can be categorized into two groups:
1) THRESHOLD VARIABLES -- The first set of variales are those associated with 
                          triggering convection and refer to TBM, BCLH, BCLP, and TDEF.
2) EVALUATION VARIABLE -- These variables return detailed information about the 
                          current convective regime.  Namely, the amount of sensible
                          and latent heat necessary for initiation (SHDEF_M and LHDEF_M),
                          the most energy efficient pathway among the two (EADV_M), and
                          the height, pressure, and potential temperature at which the
                          transition from one energy advantange to another occurs
                          (TRAN_H, TRAN_P, TRAN_T). These energy transition variables
                          may not be always present because some vertical profiles
                          may lie entirely within one regime. These variables will be
                          returned as missing if TDEF=0 because there is no longer a
                          physical meaning behind "energy advantage" and "transition".



References: --Methods--
            Tawfik and Dirmeyer 2014 GRL: A processed based framework for 
            quantifying the atmospheric background state of surface triggered
            convection

            Tawfik, A., P.A. Dirmeyer, J.A. Santanello Jr. (2015), The Heated
            Condensation Framework. Part I: Description and Southern Great Plains 
            Case Study, Journal of Hydromet. doi:10.1175/JHM-D-14-0117.1
 
            --Climatological Evaluation--
            Tawfik, A., P.A. Dirmeyer, J.A. Santanello Jr. (2015), The Heated 
            Condensation Framework. Part II:  Climatological behavior of convective
            initiation and land-atmosphere coupling over the Continental United States,
            Journal of Hydromet. doi:10.1175/JHM-D-14-0118.1
'''

from toolz.curried import curry, compose
import constants as cnts
import general_functions as gf



def integrate_over_column():
    pass
    

def hcf(ds,
        t,
        p,
        q,
        t2m=None,
        psfc=None,
        q2m=None,
        to_SI=None):
    '''Computes the heated condensation framework variable set

    Parameters
    ----------
    ds : xarray.Dataset or pandas.dataframe
    a dataframe containing the necessary state atmospheric profile
    state variables used to compute the HCF variable set

    t : string
    dataframe column name of the temperature variable [Kelvin]

    p : string
    dataframe column name of the pressure variable [Pascal]

    q : string
    dataframe column name of the specific humidity variable [kg/kg]

    t2m : string
    dataframe column name of the 2-meter temperature [Kelvin]
    This variable is optional. When included along with the 2m q and surface
    pressure they are appended to the bottom of the profile

    q2m : string
    dataframe column name of the 2-meter speific humidity [kg/kg]
    This variable is optional. When included along with the 2m t and surface
    pressure they are appended to the bottom of the profile

    psfc : string
    dataframe column name of the surface pressure [Pascal]
    This variable is optional. When included along with the 2m q and 2m t
    they are appended to the bottom of the profile

    to_SI : dict
    a dictionary containing variable names as keys and conversion
    function as a value. 
    Example: `to_SI = {'temp' : lambda x : x + 273.15}`

    '''
    # Check to make sure column has at least 3-values
    check_if_more_than_3layers_provided()
    # Convert to SI units if dictionary of conversions are provided
    convert_to_SI(ds,t,p,q,h,t2m,psfc,q2m,h2m)
    # Check units to make sure they are SI
    check_for_SI_units(ds,t,p,q,h,t2m,psfc,q2m,h2m)
    # If 2-meter variables are included then append them to the profile variables
    # otherwise ignore them
    append_2meter_variables_to_profile(t,p,q,h,t2m,psfc,q2m,h2m)
    # Take special care if someone includes ONLY the surface pressure. In this case
    # filter out the pressure levels greater than psfc
    remove_pressure_levels_greater_than_surface_pressure(t,p,q,h,psfc)
    # Compute integrated column specific humidity at each level
    qmix = column_integrate(ds,q,p)
    # compute saturation specific humidity
    qsat = gf.compute_saturation_shum(ds.t, ds.p)
    # Find the index of intersection
    iequilibrium = find_level_of_intersection()
    # Return the zero point
