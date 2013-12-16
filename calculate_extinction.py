# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 12:13:00 2013

@author: blasco
"""

import numpy as np
import scipy.stats as stats
from scipy.optimize import fmin as simplex
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse

def single_slope_chi2(params, X, Y, YERR):
    # Follows: http://python4mpia.github.io/fitting_data/simplex-fitting.html

    # extract a and b from params
    a = np.array(params[0])
    b = np.array(params[1:])
    
    # Merit function is chi-square
    model = ((X * a).transpose() + b).transpose() 
    chi2 = np.sum( (model - Y)**2 / YERR**2 )
    print a, b, chi2, "\n"
    return chi2

def single_slope_curvefit(x, a, *b):
    a = np.array(a)
    b = np.array(b)
    return (((x * a).transpose() + b).transpose()).flatten()

    
def calculate_extinction(airmasses, magnitudes, err_magnitudes=None):
     """Given 2D arrays of airmasses and magnitudes of several stars. For N
     observations of M stars, we would have a MxN array for both airmasses and 
     magnitudes, so that we can iterate for each star. Then we calculate the 
     extinction. We will fit a single coefficient k to all the stars. 
     That is, for each individual star, we fit:
         
         magnitudes_corrected = magnitudes - k * airmass 
     
     where k is common for all the stars, airmass and magnitudes is known and
     magnitudes_corrected would be the magnitudes of the stars with no 
     atmospheric extinction. """
     
     # In the unlikely case you send magnitudes without errors (do your 
     # homework, dude!) zero error will be assigned to all values. 
     if err_magnitudes == None:
         err_magnitudes = magnitudes - magnitudes
     
     # Initial guess for extinction coefficient is the first value. 
     # The initial guess of the intercept for each star is the mean of the 
     # star's magnitude
     means = [np.mean(magnitudes[ii,]) for ii in range(magnitudes.shape[0])]
     guess = [0] + means    
     #result = simplex(single_slope_chi2, guess, args=(airmasses, magnitudes,
     #                                                 err_magnitudes))
     pop, pcov = curve_fit(single_slope_curvefit, airmasses, magnitudes.flatten(),
                           sigma=err_magnitudes.flatten(), p0=guess)    
     return pop[0], pcov[0,0]
