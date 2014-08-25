# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:22:00 2013

@author: blasco
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.optimize as optimize
import sys

def intercept(x, b):
    return (x + b)

   
def calculate_zeropoint(obs_magnitudes, std_magnitudes, err_magnitudes=None):
    """ Given two arrays with:
           obs_magnitudes: the observed magnitudes (corrected from extinction)
                           for standard spectrophotometric stars 
           std_magnitudes: theoretical magnitudes calculated from the spectrum 
                           convolved with the filter 
           mag_err: vector with uncertainties in the obs_magnitudes
                           
        find the zero point of the filter, i.e: 
           std_magnitudes = obs_magnitudes + zp   
           
        if we have more than one magnitude, we will find the average of the 
        values taking into account the mag_err. 
    """
    # In the unlikely case you send magnitudes without errors (do your 
    # homework, dude!) zero error will be assigned to all values. 
    if err_magnitudes == None:
         err_magnitudes = magnitudes - magnitudes
    
    guess = np.median(obs_magnitudes - std_magnitudes)
    pop, pcov = optimize.curve_fit(intercept, obs_magnitudes, std_magnitudes,
                           sigma=err_magnitudes, p0=guess)  
    return pop[0], np.sqrt(pcov[0][0])
