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
import xarray as xr

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
    #updated_latent_energy    =  latent_energy  +  (additional_heat/bowen)
    updated_sensible_energy  =  heat_energy    +  additional_heat
    updated_latent_energy    =  latent_energy  +  (additional_heat/bowen).where(bowen != 0, 0.0)
    return updated_latent_energy, updated_sensible_energy




def mixing_diagram(ds, 
                   temp='temperature', 
                   psfc='pressure', 
                   qhum='humidity', 
                   sh_flux='sh_flux', 
                   lh_flux='lh_flux', 
                   pblh='pbl_height'):
    '''
    Calculates the full suite of mixing diagram variables when
    given a time series of of surface fluxes, surface 2m state variables, and 
    boundary layer height. Mixing diagram variables returned are surface bowen
    ratio, entrainment bowen ratio, advection ratios, and the various latent and 
    sensible heat fluxes associated with those ratios.

    Parameters
    ----------
    ds : xarray.Dataset
        an xarray dataset containing all the variable names
    
    temp : string
        dry bulb temperature [K]

    pres : string
      surface pressure [Pa]

    qhum : string
      Specific humidity [kg/kg]

    lh_flux : string
      surface latent heat flux; ensure that the sensible heat flux is the
      average over the past timestep [W/m^2]

    sh_flux : string
      surface sensible heat flux ensure that the latent heat flux is the
      average over the past timestep [W/m^2]

    pblh : string
      planetary boundary layer height or height of the mixed layer [m]

    Return
    ------
    xarray Dataset
       returns an Dataset containing the following variables:
       * surface sensible heat flux (shf_sfc)
       * sensible heat entrainment flux (shf_ent) 
       * net sensible heat flux (shf_tot = shf_sfc + shf_ent)
       * surface latent heat flux (lhf_sfc)
       * latent heat entrainment flux (lhf_ent) [W/m^2]
       * net latent heat flux (lhf_tot = lhf_sfc + lhf_ent) 
       All output is in units of [W/m^2]

    '''
    #---------------------------------------------------
    # Check that all required variable names are in the 
    # dataset
    #---------------------------------------------------
    gf.check_input_is_xarray(ds)
    gf.check_variable_names(ds, [temp, psfc, qhum, sh_flux, lh_flux, pblh])
    gf.check_time_is_a_dimension(ds)
    dt = gf.get_interval(ds)
    #---------------------------------------------------
    # Calculate the air density
    #---------------------------------------------------
    rho           = gf.air_density(ds[psfc], ds[temp], ds[qhum])
    heat_energy   = gf.temp_to_energy(ds[temp])
    latent_energy = gf.shum_to_latent_energy(ds[qhum])
    bowen         = gf.bowen_ratio(ds[sh_flux], ds[lh_flux])
    #---------------------------------------------------
    # Calculate the change in heat energy due to 
    # sensible heat flux
    #---------------------------------------------------
    dtheta        = gf.add_heat_to_layer(rho, ds[sh_flux], ds[pblh], dt)
    #---------------------------------------------------
    # Stagger the sensible and latent heat energy
    # we can compare the change over time (delta heat and moisture)
    #---------------------------------------------------
    initial_heat    = heat_energy.shift(time=1)
    initial_latent  = latent_energy.shift(time=1)

    #---------------------------------------------------
    # Calculate the air density
    #---------------------------------------------------
    updated_lh, updated_sh = surface_flux_increment( dtheta, bowen, initial_latent, initial_heat)

    #---------------------------------------------------
    # Separate all the energy components
    #---------------------------------------------------
    shf_net = (heat_energy - initial_heat ) * (rho * ds[pblh]) / dt
    shf_ent = (heat_energy - updated_sh   ) * (rho * ds[pblh]) / dt
    shf_sfc = (updated_sh  - initial_heat ) * (rho * ds[pblh]) / dt

    lhf_net = (latent_energy - initial_latent) * (rho * ds[pblh]) / dt
    lhf_ent = (latent_energy - updated_lh    ) * (rho * ds[pblh]) / dt
    lhf_sfc = (updated_lh    - initial_latent) * (rho * ds[pblh]) / dt    
    #---------------------------------------
    # Rename output variables
    #---------------------------------------
    shf_sfc.name = 'shf_sfc'
    shf_ent.name = 'shf_ent'
    shf_net.name = 'shf_net'
    lhf_sfc.name = 'lhf_sfc'
    lhf_ent.name = 'lhf_ent'
    lhf_net.name = 'lhf_net'
    return xr.merge( [shf_sfc, shf_ent, shf_net, lhf_sfc, lhf_ent, lhf_net] )




def lcl_deficit(ds, 
                t2m ='t2m' , 
                psfc='psfc',
                q2m ='q2m' ,
                pblh='pblh'):
    '''
    Calculates the Lifted condensation level deficit
    The LCL deficit is the boundary layer height subtracted from the 
    height of the LCL.

    LCL deficit <= 0 means convective initiation is likely
    LCL deficit > 0 means convective initiation is less likely
    

    Parameters
    ----------
    ds : xarray.Dataset
        an xarray dataset containing all the variable names
    
    t2m : string
        2 meter dry bulb temperature [K]

    pres : string
      surface pressure [Pa]

    q2m : string
      2-meter specific humidity [kg/kg]

    pblh : string
      planetary boundary layer height or height of the mixed layer [m]

    Return
    ------
    xarray Dataset
       returns the LCL deficit in meters [m]

    '''
    #---------------------------------------------------
    # Check that all required variable names are in the 
    # dataset
    #---------------------------------------------------
    gf.check_input_is_xarray(ds)
    gf.check_variable_names(ds, [t2m, psfc, q2m, pblh])
    gf.check_time_is_a_dimension(ds)

    #---------------------------------------------------
    # Calculate the intermediate values
    #---------------------------------------------------
    theta = gf.potential_temperature(ds[t2m], ds[psfc])
    tsat  = gf.saturation_temperature(ds[t2m], ds[psfc], ds[q2m])
    plcl  = gf.lcl_pressure(ds[t2m], ds[psfc], ds[q2m])
    tvirt = gf.virtual_temperature(ds[t2m], ds[q2m])
    hlcl  = gf.lcl_height(ds[t2m], ds[psfc], ds[q2m])

    lcl_deficit  =  hlcl - ds[pblh]
    lcl_deficit.name = 'lcl_deficit'
    return lcl_deficit



