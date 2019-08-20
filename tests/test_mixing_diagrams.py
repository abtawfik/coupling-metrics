'''Tests for mixing diagrams
'''

import numpy as np
import pandas as pd

import pytest
from expected_output import MixDiagOutput as expected

# module being tested
import sys
sys.path.append('../metrics')
from mixing_diagram import mixing_diagram


#####################
# Read in test data #
#####################
df  =  pd.read_csv('./data/sample_fluxes.csv')


def test_mixing_diagram():
	dt = 8 * 60 * 60
	shf_sfc, shf_ent, shf_tot, lhf_sfc, lhf_ent, lhf_tot = mixing_diagram(df.temp.values,
																		  df.pres, 
																		  df.shum, 
																		  df.sh, 
																		  df.lh, 
																		  df.pblh, 
																		  dt)
	# check each output variable
	assert all(shf_sfc == expected.shf_sfc)
	assert all(shf_ent == expected.shf_ent)
	assert all(shf_tot == expected.shf_tot)
	assert all(lhf_sfc == expected.lhf_sfc)
	assert all(lhf_ent == expected.lhf_ent)
	assert all(lhf_tot == expected.lhf_tot)

