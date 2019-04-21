# -*- coding: utf-8 -*-
'''
Created Apr 20 2019

Contains common calculations and functions

@tawfik
'''

######################
# 3rd party packages #
######################
import numpy as np

######################
#     Constants      #
######################
from constants import cp, Rd, R_cp, C2K
from constants import one_minus_ep, ep, es0



#---------------------------------------------------------------------
# Calculates the dew point temperature following the Arden Buck Eq.
#---------------------------------------------------------------------
def dew_point(q, 
              p, 
              A=610.8, 
              B=237.3, 
              C=17.2693882):
    '''Calculates the dew point temperature following the Arden Buck Eq.

    Parameters
    ----------
    q : array
       specific humidity array [kg/kg]

    p : array
       pressure [Pa]

    A : float
      Empirically determined constant specified by Arden Beck eq.

    B : float
      Empirically determined constant specified by Arden Beck eq.

    C : float
      Empirically determined constant specified by Arden Beck eq.

    Returns
    -------
    array
       returns an array containing the dew point temperature the same size as
       q and p input arrays. Units = [K]
    '''
    #------------------------------
    # Calculate the vapor pressure
    #------------------------------
    e_vapor = shum_to_vapor_pressure(q, p)
    return ( (np.log(e_vapor/A)*B) / (C - np.log(e_vapor/A)) ) + C2K


#---------------------------------------------------------------------
# Calculates vapor pressure from specific humidity and pressure
#---------------------------------------------------------------------
def shum_to_vapor_pressure(q, p):
    ''' Calculates vapor pressure from specific humidity and pressure
    
    Parameters
    ----------
    q : array
       specific humidity array [kg/kg]

    p : array
       pressure [Pa]

    Returns
    -------
    array
       returns an array containing vapor pressure values the same size as
       q and p input arrays. Units = [Pa]
    '''
    return (q*(p/1e2))/ (0.622+0.378 * q) * 1e2


#---------------------------------------------------------------------
#
#---------------------------------------------------------------------
def saturation_shum(p, 
                    t,
                    A=17.269,
                    B=35.86 ):
    '''Calculates the saturation specific humidity [kg/kg] following the 
    AMS glossary definition
    
    

    Parameters
    ----------
    t : array or float
       dry bulb temperature [K]

    p : array or float
       pressure [Pa]

    A : float
       empirically determined constant 

    B : float
       empirically determined constant 


    Returns
    -------
    array
       returns an array containing saturation specific humidity the same size as
       t and p input arrays. Units = [kg/kg]
    '''
    #--------------------------------------------
    #--- Split out the numerator and denomenator
    #--------------------------------------------
    numerator   =  ep * (es0*exp((A*(t-C2K))/( t-B)))
    denomenator =  (p/1e2)-one_minus_ep*(es0*exp((a*( t-C2K))/( t-B)))
    esat        =  numerator/denomenator

    #--------------------------------------------
    #--- Vapor pressure and co
    #--------------------------------------------
    return  esat / (1 + esat)



