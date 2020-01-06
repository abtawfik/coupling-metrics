# Coupling Metrics Toolkit (CoMeT)
Calculate various state-of-the-art land-atmosphere coupling metrics.
For more information about specific metrics and use of the modules go to www.coupling-metrics.com

### When to Use CoMeT?
The metrics are meant to facilitate use of the most common land-atmosphere coupling metrics for research purposes. 

## Pre-requisites

## Installation
```bash
pip install coupling-metrics
```

## Usage
There are two ways to use CoMeT after installation. Using the python API or on the command-line. See both uses below.

### Command-line Interface Example
Say you are trying to compute the terrestrial coupling index for NARR data (North American Regional Reanalysis) between latent heat flux and soil moisture for a single year. You invoke this one the command-line:
```bash
comet coupling --xname=soilm --yname=lhtfl --averaging=season --outname=NARR_lhf_vs_soilm_2017.nc lhtfl.2017.nc soilm.2017.nc
```
This will create a file called NARR_lhf_vs_soilm_2017.nc containing the coupling index. If you wish to compute it across multiple years you can simply use the typical bash syntax
```bash
comet coupling --xname=soilm --yname=lhtfl --averaging=season --outname=NARR_lhf_vs_soilm_1979-2017.nc lhtfl.????.nc soilm.????.nc
```
** Note: be mindful of memory in this case because currently all the data need to be loaded into memory. This might change in the future.

If you want to see a list of metrics
```bash
comet --help
```

To see the arguments for a list of metrics
```bash
comet mixing --help
```


### Python API
Each metric can also be called directly from python.
```python
import comet as cm
import xarray as xr
from glob import glob

# Load data using xarray
latent_heat_files = glob('lhtfl.????.nc')
soilm_files       = glob('soilm.????.nc')
ds                = xr.open_mfdataset(latent_heat_files + soilm_files)

# Compute coupling index
terra_coupling = cm.CouplingIndex()
terra_coupling = terra_coupling.compute(ds, xname='soilm', yname='lhtfl', averaging='month')

# Output it to a file (may take a while depending on how much data is being processed)
terra_coupling.to_netcdf('My_new_NARR_Terra_coupling_lhf_vs_soilm.nc')
```

To see a list of metrics
```python
import comet as cm

print(cm.list_metrics)
```

## List of available metrics:

1. [Convective Triggering Potential]('http://journals.ametsoc.org/doi/abs/10.1175/1525-7541%282003%29004%3C0552%3AACOSML%3E2.0.CO%3B2') (aka CTP-HiLow or just CTP)
Evaluates morning atmospheric profiles to determine whether dry or wet soils are more likely to trigger convection

2. [Mixing Diagrams](http://journals.ametsoc.org/doi/abs/10.1175/2009JHM1066.1)
Uses the diurnal covariation of temperature and humidity to quantify heat and moisture fluxes into the planetary boundary layer

3. [Terrestrial Coupling Index](http://onlinelibrary.wiley.com/doi/10.1029/2011GL048268/abstract) (aka terrestrial coupling parameter)
Quantifies the degree to which soil moisture variations control changes in surface energy fluxes.  Can be latent or sensible heat flux and is generally enough to apply to other surface variables such as the relationship between LAI and sensible heat flux.

4. [Heated Condensation Framework](http://journals.ametsoc.org/doi/abs/10.1175/JHM-D-14-0117.1)
Assesses the atmospheric background state with respect to convective initiation and identifies local versus non-locally triggered moist convection.

5. [Relative Humidity Tendency]('http://journals.ametsoc.org/doi/abs/10.1175/1525-7541(2004)005%3C0086%3AIOSMOB%3E2.0.CO%3B2')
Returns the contribution of surface energy fluxes, dry air entrainment, heat entrainment, and boundary layer growth to changes in top of boundary layer relative humidity.

6. [Soil Moisture Memory]('http://journals.ametsoc.org/doi/abs/10.1175/1520-0442(1988)001%3C0523:TIOPEO%3E2.0.CO;2') (the statisitcal form)
Determines the timescale at which initial soil moisture anomalies are retained over time.  The lagged autocorrelation of soil moisture is used to make the determination of memory. 


