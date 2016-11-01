#Coupling Metrics Toolkit (CoMeT)
Contains several fortran 90 subroutines used to calculate various land-atmosphere coupling metrics.
For more information about each metric and use of the modules go to www.coupling-metrics.com

General use:
The modules available are meant to facilitate use of the most common land-atmosphere coupling metrics for research purposes. The modules are currently written in Fortran 90 and are individually packaged.  The most common use is to take a given module and create a library that is usable in the desired programming language.  For example if you are using the NCAR Common Langauge (NCL) then you could perform the following actions with the files in the */ncl_portable/ folder

> WRAPIT ctp_hilow.stub ctp_hilow.f90

This will create a file called ctp_hilow.so.  This shared object file can then be called from within your NCL script by adding lines like these:

> external CTPHiLow "ctp_hilow.so"

> CTP    =  new( 1, float )

> HILOW  =  new( 1, float )

> CTPHiLow::ctp_hi_low( nlev, t, q, p, t2m , q2m , psfc, CTP, HILOW , missing )

Note that you have to instantiate the output variables before you can call the routine in NCL.  A similar process can be performed in python where you create a shared object using the <b>f2py</b> command using files from the */ncl_portable/ directory.


Below is a brief description of each metric:

1) Convective Triggering Potential (aka CTP-HIlow or just CTP)

Evaluates morning atmospheric profiles to determine whether dry or wet soils are more likely to trigger convection

2) Mixing Diagrams

Uses the diurnal covariation of temperature and humidity to quantify heat and moisture fluxes into the planetary boundary layer

3) Terrestrial Coupling Index (aka terrestrial coupling parameter)

Quantifies the degree to which soil moisture variations control changes in surface energy fluxes.  Can be latent or sensible heat flux and is generally enough to apply to other surface variables such as the relationship between LAI and sensible heat flux.

4) Heated Condensation Framework

Assesses the atmospheric background state with respect to convective initiation and identifies local versus non-locally triggered moist convection.

5) Relative Humidity Tendency

Returns the contribution of surface energy fluxes, dry air entrainment, heat entrainment, and boundary layer growth to changes in top of boundary layer relative humidity.

6) Soil Moisture Memory (the statisitcal form)

Determines the timescale at which initial soil moisture anomalies are retained over time.  The lagged autocorrelation of soil moisture is used to make the determination of memory. 


