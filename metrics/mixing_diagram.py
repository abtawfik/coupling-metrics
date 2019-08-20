'''
 Calculates the full suite of mixing diagram variables when
 given a time series of of surface fluxes, surface 2m state variables, and 
 boundary layer height. Mixing diagram variables returned are surface bowen
 ratio, entrainment bowen ratio, advection ratios, and the various latent and 
 sensible heat fluxes associated with those ratios. The utility is to characterize 
 the coupling between surface fluxes and top of boundary layer fluxes in tandem
 with knowledge regarding soil moisture state.  More details regarding motivation
 can be found in the references below.
 
  References: Santanello et al. 2009,  A Modeling and Observational Framework 
              for Diagnosing Local Land-Atmosphere Coupling on Diurnal Time Scales

              Santanello et al. 2011,  Diagnosing the Sensitivity of Local 
              Land-Atmosphere Coupling via the Soil Moisture-Boundary Layer Interaction

              ** Comprehensive Evaluation **
              Santanello et al. 2013, Diagnosing the Nature of L-A Coupling:
              A Case Study of Dry/Wet Extremes in the U.S. Southern Great Plains  
'''

# 3rd party packages
import numpy as np
import pandas as pd

import constants as cnts
import general_functions as gf


def surface_flux_increment(additional_heat, bowen, latent_energy, heat_energy):
    '''Add the new heat flux to a given initial sensible and latent heat energy 
    values. The change in latent heat in a column depends on the bowen ration

    Parameters
    ----------
    additional_heat : float
       sensible heat energy being added to the initial value [J/kg]

    bowen : float
       bowen ratio [unitless]

    latent_energy : float
       initial latent heat energy in a layer [J/kg]

    heat_energy : float
       initial sensible heat energy in a layer [J/kg]

    Return
    ------
    float, float
       updated latent and sensible heat energy after the 
       additional energy is included [J/kg]
    '''
    updated_latent_energy    =  latent_energy  +  (additional_heat/bowen)
    updated_sensible_energy  =  heat_energy    +  additional_heat
    return updated_latent_energy, updated_sensible_energy


def mixing_diagram(temp, pres, qhum, sh_flux, lh_flux, pblh, dt):
    '''
    Calculates the full suite of mixing diagram variables when
    given a time series of of surface fluxes, surface 2m state variables, and 
    boundary layer height. Mixing diagram variables returned are surface bowen
    ratio, entrainment bowen ratio, advection ratios, and the various latent and 
    sensible heat fluxes associated with those ratios.

    Parameters
    ----------
    temp : float
       dry bulb temperature [K]

    pres : float
       surface pressure [Pa]

    qhum : float
       specific humidity [kg/kg]

    lh_flux : float
       surface latent heat flux [W/m^2]

    sh_flux : float
       surface sensible heat flux [W/m^2]

    pblh : float
       planetary boundary layer height or height of the mixed layer [m]

    dt : float
       time step. also known as the sample interval [s]

    Return
    ------
    float, float
       returns the sensible heat entrainment flux (sh_ent) and 
       latent heat entrainment flux (lh_ent) [W/m^2]

    '''
    #---------------------------------------------------
    # Calculate the air density
    #---------------------------------------------------
    rho           = gf.air_density(pres, temp, qhum)
    heat_energy   = gf.temp_to_energy(temp)
    latent_energy = gf.shum_to_latent_energy(qhum)
    bowen         = gf.bowen_ratio(sh_flux, lh_flux)
    #---------------------------------------------------
    # Calculate the change in heat energy due to 
    # sensible heat flux
    #---------------------------------------------------
    dtheta        = gf.add_heat_to_layer(rho, sh_flux, pblh, dt)
    #---------------------------------------------------
    # Stagger the sensible and latent heat energy
    # so 
    #---------------------------------------------------
    initial_heat      = np.roll(heat_energy  , 1)
    initial_latent    = np.roll(latent_energy, 1)
    initial_heat[0]   = np.nan
    initial_latent[0] = np.nan

    #---------------------------------------------------
    # Calculate the air density
    #---------------------------------------------------
    updated_lh, updated_sh = surface_flux_increment( dtheta, bowen, initial_latent, initial_heat)
    #---------------------------------------------------
    # Separate all the energy components
    #---------------------------------------------------
    shf_tot = (heat_energy - initial_heat ) * (rho * pblh) / dt
    shf_ent = (heat_energy - updated_sh   ) * (rho * pblh) / dt
    shf_sfc = (updated_sh  - initial_heat ) * (rho * pblh) / dt

    lhf_tot = (latent_energy - initial_latent) * (rho * pblh) / dt
    lhf_ent = (latent_energy - updated_lh    ) * (rho * pblh) / dt
    lhf_sfc = (updated_lh    - initial_latent) * (rho * pblh) / dt
    return shf_sfc, shf_ent, shf_tot, lhf_sfc, lhf_ent, lhf_tot
