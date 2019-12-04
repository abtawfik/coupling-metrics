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
import xarray as xr


def coupling_index(ds, x, y, dim='time'):
    '''
    Generic equation for computing covariance divided by standard deviation 
    across the time dimension
    '''
    return ds[x,y].cov(axis=dim) / ds[x].std(axis=dim)
