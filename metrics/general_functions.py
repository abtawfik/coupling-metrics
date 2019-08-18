'''
A set of functions that are not specific to any specific metric
Functions such as calculating saturation specific humidity or air density
'''

# 3rd party packages
import numpy as np
import pandas as pd
from toolz.curried import curry, compose

import constants as cnts


@curry
def air_density(p, t, q):
    ep = cnts.ep.value
    Rd = cnts.Rd.value
    return p / (Rd * t * ((1. + (q/ep)) / (1. + q)))

@curry
def temp_to_energy(t):
    return cnts.cp.value * t

@curry
def shum_to_latent_energy(q):
    return cnts.Lv.value * q

@curry
def bowen_ratio(sensible_heat, latent_heat):
    return np.where( latent_heat == 0, np.nan, sensible_heat / latent_heat )


