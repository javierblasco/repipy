# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:25:30 2013

@author: blasco
"""

################################################################################
#                   Campaign dependent                                         # 
################################################################################

# Work directory
directory = "/mnt/data/DEEP_OBS/20130606/"

# Directory where the lemon pipeline (https://github.com/vterron/lemon) is
lemon_dir = "/home/blasco/Desktop/librerias_python/lemon"
        
# Images not to be used. Early bias images show too high counts statistics. 
# And the focus image and the skyflat 1 are of different size to all the others. 
bad_bias = ["bias_20130606_" + str(ii).zfill(3) + ".fits" for ii in range(1,11)]
remove_images = ["focustelescope_20130606_sdssr_001.fits", \
                 "skyflat_20130606_sdssr_001.fits",   #saturated
                 "cig0812_20130606_sdssr_028.fits",    # bias
                 "cig0812_20130606_sdssr_018.fits",    # wrong
                 ] + bad_bias

nstars = 40  #number of stars used to align images
max_FWHM = 6     # largest reasonable FWHM. Sometimes LEMON seeing does not give
             # reasonable numbers (like 22 for the FWHM) and that messes up
             # the detection of stars. It does not need to be accurate, because 
             # seeing is used mainly to detect sources. But 22 is too much!

# Keywords in the header of images
filterk = "INSFLNAM"     # filter name
exptimek = "exptime"    # exposure time (seconds)
objectk = "object"       # name of object 
rak = "ra"               # degree
deck = "dec"             # degree
datek = "date-obs"
telescope = "CAHA2.2"
gaink = "CCDSENS"       # gain 
read_noisek = "CCDRON"  # read-out noise
pix_scale = 1         # pixel scale (arcsec)
FoV = 0.1             # rough radius of the FoV for astrometry calculations

# Directory with saved master flats. That will save calculating time later on.


