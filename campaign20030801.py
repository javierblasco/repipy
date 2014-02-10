# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:25:30 2013

@author: blasco
"""

################################################################################
#                   Campaign dependent                                         # 
################################################################################

# Work directory
directory = "/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche1/"

# Directory where the lemon pipeline (https://github.com/vterron/lemon) is
#lemon_dir = "/home/blasco/Desktop/librerias_python/lemon"
        
# Images not to be used. 
remove_images = ["bias_20030801_" + str(ii).zfill(3) + ".fits" for ii in range(4,13)]
remove_images = remove_images + [#"bd+28_20030801_rGunn_001.fits",   # standard saturated
                                 "bd+28_20030801_rGunn_002.fits",   # saturated
                                 #"bd+28_20030801_rGunn_003.fits",   # saturated
                                 "bd+25_20030801_rGunn_001.fits",   # saturated
                                 "bd+25_20030801_rGunn_002.fits",   # saturated
                                 "kop27_20030801_rGunn_001.fits",   # saturated
                                 "kop27_20030801_rGunn_003.fits"   # saturated 
                            
                 ]
# standards_campaign contains the names of the standards used for this campaign.
# They are in python Re format, and will be used during the calibration to 
# find the files that have been reduced.
standards_campaign = "bd+28|bd+25|kop27"


# Keywords in the header of images
filterk = "INSFLNAM"     # filter name
exptimek = "exptime"    # exposure time (seconds)
objectk = "object"       # name of object 
airmassk = "airmass"     # airmass keyword on image
rak = "ra"               # degree
deck = "dec"             # degree
datek = "date-obs"
telescope = "CAHA2.2"
gaink = "CCDSENS"       # gain 
read_noisek = "CCDRON"  # read-out noise
pix_scale = 0.529         # pixel scale (arcsec)
FoV = 0.25                # rough radius of the FoV for astrometry calculations


