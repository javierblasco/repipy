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


# standards_campaign contains the names of the standards used for this campaign.
# They are in python Re format, and will be used during the calibration to 
# find the files that have been reduced.
standards_campaign = "bd\+28|bd\+25|kop27"


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


# Characteristics of camera
pix_scale = 0.529         # pixel scale (arcsec)
FoV = 0.25                # rough radius of the FoV for astrometry calculations
max_counts = 55000        # consider saturated any pixel above this value
circular_FoV = True       # There is a circle exposed in the CCD, with the rest of the camera much darker, e.g.  CAHA 2.2

# Keywords for things the pipeline will calculate
skyk = "sky"
sky_stdk = "sky_std"
seeingk = "seeing"     

# Inevitably, there are some decisions that can only be done by an experienced eye. I wouldn't trust (for the moment)
# a program to tell me if a bias has structure, for example. So, while we do most things automatically, at the end
# there is always a need to run the pipeline several times, and decide upon specific things. Another example of things
# to decide is the images not to be used, that for the moment is done only after renaming the images, for simplicity.

type_of_bias_subtraction = "--median"   # options: "--median", "--mean", ""  . Mean/median to subtract the mean/median of an
                                   #          image when using repipy.arith, nothing to subtract the whole image
                                   #          because it has structure.

remove_images = ["bias_20030801_" + str(ii).zfill(3) + ".fits" for ii in range(4,13)]
remove_images = remove_images + [#"bd+28_20030801_rGunn_001.fits",   # standard saturated
                                 "bd+28_20030801_rGunn_002.fits",   # saturated
                                 #"bd+28_20030801_rGunn_003.fits",   # saturated
                                 "bd+25_20030801_rGunn_001.fits",   # saturated
                                 "bd+25_20030801_rGunn_002.fits",   # saturated
                                 "kop27_20030801_rGunn_001.fits",   # saturated
                                 "kop27_20030801_rGunn_003.fits"   # saturated
                                           ]
