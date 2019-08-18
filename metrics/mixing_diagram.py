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




