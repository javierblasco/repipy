#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Created on Sun Oct 20 23:55:05 2013

@author: blasco
"""

import numpy as np
import matplotlib.pyplot as plot
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from scipy import stats
import time

def main(images, over, margin=10):
    """ Calculate statistics of bias and overscan to compare both values. 
        Do a linear fit. 
        INPUT arguments: 
           images - names of the bias images 
           overscan - region to be considered as overscan (x0,y0,x1,y1). Careful, 
                      axis1 is the second axis for astropy 
                      (http://docs.astropy.org/en/latest/io/fits/index.html)
           margin - marging to be left around both the bias and overscan areas
                    in order to avoid border effects. 
        """
    time0 = time.ctime()
    print time0
    # Overscan and bias arrays will contain two rows for median and stddev 
    overscan_array = np.zeros((2, len(images)),dtype=np.float64)
    bias_array = np.zeros((2, len(images)), dtype=np.float64)
    for index, image in enumerate(images):
        im_data = fits.getdata(image).astype(np.float128)
        # Construct overscan and bias areas
        over_y0, over_x0, over_y1, over_x1 = np.array(over) + \
                              np.array([margin, margin, -1*margin, -1*margin])
        bias_y0, bias_x0, bias_y1, bias_x1 = \
                   margin, margin, over_y0 - margin, over_x1 - margin   
        overscan = im_data[over_x0:over_x1, over_y0:over_y1]
        bias = im_data[bias_x0:bias_x1, bias_y0:bias_y1]
        
        # Now statistics: median absolute deviation (much better than std) 
        # median and clipped mean. 
        over_MAD = np.median(np.abs(overscan - np.median(overscan))) # 1.5 * MAD ~ sigma
        over_median = np.median(overscan)
        good_over = np.where( (over_median - 4.5 * over_MAD < overscan) & 
                              (over_median + 4.5 * over_MAD > overscan))
        over_mean = np.mean(overscan[good_over])
        
        bias_MAD = np.median(np.abs(bias - np.median(bias)))
        bias_median = np.median(bias)
        good_bias = np.where( (bias_median - 4.5 * bias_MAD < bias) &
                              (bias_median + 4.5 * bias_MAD > bias) )
        bias_mean = np.mean(bias[good_bias])             
        
        overscan_array[:,index] = over_mean, over_MAD
        bias_array[:,index] = bias_mean, bias_MAD
    
#    diff = bias_array[0,:] - overscan_array[0,:]
#    print "Difference bias - overscan (median, stddev):", np.median(diff), \
#                                    np.std(diff)
#    
    # Linear fit
#    slope, intercept, r_value, p_value, std_err = \
#                        stats.linregress(overscan_array[0,:], bias_array[0,:])
#    print slope, intercept, r_value, p_value, std_err
    
    t1 = time.ctime()
    print t1    
    
    # Now we plot bias vs overscan
    plt.plot(bias_array[0,:] - overscan_array[0,:], 'o')
#    plt.plot(overscan_array[0,:], overscan_array[0,:] * slope + intercept)
    plt.show()
#    
    return None
