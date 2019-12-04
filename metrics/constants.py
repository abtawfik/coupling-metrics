'''
Contains all the constants used across functions and metrics
'''
from collections import namedtuple

# Create a constants data type that contains metadata for
# querying later
Constant = namedtuple('constant', ['long_name', 'unit', 'value'])
# Constants
p_ref = Constant(value=1e5   , long_name='reference pressure'                , unit='Pa')
Lv    = Constant(value=2.5e6 , long_name='latent heat of vaporization'       , unit='J/kg')
cp    = Constant(value=1005.7, long_name='specific heat at constant pressure', unit='J/kg')
grav  = Constant(value=9.81  , long_name='acceleration due to gravity'       , unit='m/s^2')
Rd    = Constant(value=287.04, long_name='gas constant for dry air'          , unit='J/(K * kg)')
ep    = Constant(value=0.622 , long_name='gas constant ratio'                , unit='unitless')
R_cp  = Constant(value=Rd.value/cp.value , long_name='poisson constant'      , unit='unitless')
onemp = Constant(value=1-ep.value, long_name='1 minus gas constant ratio'    , unit='unitless')
C_to_K= Constant(value=273.15    , long_name='convert from Celsius to Kelvin', unit='K')
es0   = Constant(value=6.11      , long_name='vapor pressure constant'       , unit='unitless')
a0    = Constant(value=17.269    , long_name='constant for saturation humidity', unit='unitless')
b0    = Constant(value=35.86     , long_name='constant for saturation humidity', unit='unitless')



