'''Tests for mixing diagrams
'''

import numpy as np
import pandas as pd

import pytest
from expected_output import MixDiagOutput as expected

# module being tested
import sys
sys.path.append('../metrics')
from mixing_diagrams import MixDiag


#####################
# Read in test data #
#####################
df  =  pd.read_csv('sample_fluxes.csv')


def test_density():
    
