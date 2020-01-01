
####################
# Standard Library #
####################
from dataclasses import dataclass
from typing import Dict, Callable
from pathlib import Path
from datetime import datetime
######################
# Installed packages #
######################
import click

#####################
#   Local modules   #
#####################
from .file_utils import get_data
from .metrics.mixing_diagram import MixingDiagram
from .metrics.coupling_indices import CouplingIndex


    
#-------------------------------------------
# Primary entry point for the command-line 
#-------------------------------------------
@click.group()
def comet():
    '''Used to compute various land-atmosphere coupling metrics, 
    just pass the necessary input files
    '''
    click.echo('In the main comet')




    
#-------------------------------------------
# Mixing diagram interface
#-------------------------------------------
@comet.command()
@click.option('--temp',
              required=True,
              type=click.STRING,
              help='name of 2-meter temperature variable in the file(s)')
@click.option('--psfc',
              required=True,
              type=click.STRING,
              help='name of surface pressure variable in the file(s)')
@click.option('--qhum',
              required=True,
              type=click.STRING,
              help='name of 2-meter specific humidity variable in the file(s)')
@click.option('--sh_flux',
              required=True,
              type=click.STRING,
              help='name of surface sensible heat flux variable in the file(s)')
@click.option('--lh_flux',
              required=True,
              type=click.STRING,
              help='name of surface latent heat flux variable in the file(s)')
@click.option('--pblh',
              required=True,
              type=click.STRING,
              help='name of boundary layer height variable in the file(s)')
@click.option('--averaging',
              default='season',
              type=click.Choice(['month','season'], case_sensitive=False),
              help='time to average over')
@click.option('--outname',
              default=f'CoMet_MixD_{int(datetime.now().timestamp())}.nc',
              help='name of output file; default is `CoMeT_MixD_<current_timestamp>.nc`')
@click.argument('input_files', type=click.Path(exists=True), nargs=-1)
def mixing(temp, psfc, qhum, sh_flux, lh_flux, pblh, averaging, outname, input_files):
    '''Computes time averaged mixing diagram variables such as entrainment ratio and 
    diurnal cycle sensible and latent energy (`cp*T` and `Lv*q`).
    
    ----------
    References
    ----------
    Santanello et al. 2009,  A Modeling and Observational Framework
    for Diagnosing Local Land-Atmosphere Coupling on Diurnal Time Scales
    
    Santanello et al. 2011,  Diagnosing the Sensitivity of Local
    Land-Atmosphere Coupling via the Soil Moisture-Boundary Layer Interaction
    
    Santanello et al. 2013, Diagnosing the Nature of L-A Coupling:
    A Case Study of Dry/Wet Extremes in the U.S. Southern Great Plains
    '''
    click.echo('------Computing metrics for mixing diagrams-----')
    args = {'ds'        : get_data(input_files),
            'temp'      : temp,
            'psfc'      : psfc,
            'qhum'      : qhum,
            'sh_flux'   : sh_flux,
            'lh_flux'   : lh_flux,
            'pblh'      : pblh,
            'averaging' : averaging}
    ds = MixingDiagram().compute(**args).to_netcdf(Path(outname))
    click.echo(f'Done with mixing diagram computation! Output file --> {outname}')




    
#-------------------------------------------
# Coupling Index interface
#-------------------------------------------
@comet.command()
@click.option('--xname',
              required=True,
              type=click.STRING,
              help='name of the independent variable in the file(s)')
@click.option('--yname',
              required=True,
              type=click.STRING,
              help='name of dependent variable in the file(s)')
@click.option('--averaging',
              default='season',
              type=click.Choice(['month','season'], case_sensitive=False),
              help='time to average over')
@click.option('--outname',
              default=f'CoMet_CI_{int(datetime.now().timestamp())}.nc',
              help='name of output file; default is `CoMeT_CI_<current_timestamp>.nc`')
@click.argument('input_files', type=click.Path(exists=True), nargs=-1)
def coupling(xname, yname, averaging, outname, input_files):
    '''Computes the coupling index based on Dirmeyer (2011) approach averaged over
    monthly or seasonal.
    Basic formula is: coupling_index = covariance(x,y) / stddev(x)

    ---------
    Reference
    ---------
    Dirmeyer, The terrestrial segment of soil moisture–climate coupling (2011)
    '''
    click.echo(f'------Computing coupling index for {xname}-{yname}-----')
    args = {'ds'        : get_data(input_files),
            'xname'     : xname,
            'yname'     : yname,
            'averaging' : averaging}
    ds = CouplingIndex().compute(**args).to_netcdf(Path(outname))
    click.echo(f'Done with coupling index computation! Output file --> {outname}')


