
####################
# Standard Library #
####################
from dataclasses import dataclass
from typing import Dict, Callable

######################
# Installed packages #
######################
import click

#####################
#   Local modules   #
#####################
from .file_utils import get_data



#-------------------------------------------
# Metrics basic template/data structure
#-------------------------------------------
@dataclass
class BaseMetric:
    args : Dict
    func : Callable
    def __post_init__(self):
        return None
    
    def __call__(self, ds):
        return self.func(ds, **self.args)
    
    def check_necessary_variables(self, ds):
        raise NotImplementedError('Need to have a way to check variables are available')

        
# Every metric has a set of arguments and a functional interface
# Check 
    
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
@click.argument('input_files', type=click.Path(exists=True), nargs=-1)
def mixing(input_files):
    click.echo('------Computing metrics for mixing diagrams-----')


