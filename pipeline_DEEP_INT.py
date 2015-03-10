#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
@ author: javier blasco herrera
@ email: blasco@iaa.es, javier.blasco.herrera@gmail.com
@ adapted to INT by pablo ramírez: prm@iaa.es
Pipeline for the reduction of deep INT images.
"""

import os, shutil, re, sys, glob
import subprocess
import numpy as np
import datetime
import repipy.utilities as utilities
import repipy.combine as combine_images
import repipy.arith as arith
#import repipy.tidy_up2 as tidy_up
import repipy.create_masks as create_masks
import repipy.remove_cosmics as remove_cosmics
import repipy.find_sky as find_sky
import repipy.complete_headers as complete_headers
import repipy.calculate_airmass as calculate_airmass
import repipy.rename as rename
import repipy.median_filter as median_filter
import repipy.cross_match as cross_match
import repipy.estimate_seeing as estimate_seeing
import astropy.io.fits as fits
import dateutil.parser

from lemon import methods
import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot
    from iraf import apphot

if len(sys.argv) != 3:
    print sys.exit("Give me a CAMPAIGN FILE and the WORKING DIRECTORY...")

directory = sys.argv[2]
                   
execfile(sys.argv[1])

################################################################################
#                            Campaign "independent"
################################################################################
def tidy_pattern_dict():
    in_pattern = {}
    in_pattern["cig"] = "^(?P<name>cig)(?P<cig_num>\d{4})_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    in_pattern["bias"] = "^(?P<name>bias)_(?P<date>\d{8})_"+\
                         "(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    in_pattern["skyflats"] = "^(?P<name>skyflat)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    in_pattern["standards"] = "^(?P<name>feige34|he3|hz15)_(?P<date>\d{8})_" +\
                              "(?P<filt>.*)_(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    in_pattern["blanks"] = "^(?P<name>blank)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    return in_pattern     

def cosmic_removal_param(telescope = ''):
    if telescope == "OSN":
       cosmic_dict = {}
       cosmic_dict["gain"] = "2"
       cosmic_dict["readnoise"]= "5"
       cosmic_dict["sigclip"] = "5"            
    if telescope == "CAHA2.2": 
       cosmic_dict = {}
       cosmic_dict["gain"] = "1.43"  # from header of image
       cosmic_dict["readnoise"]= "7.4" # from header of image
       cosmic_dict["sigclip"] = "5"
    if telescope == "INT": 
       cosmic_dict = {}
       cosmic_dict["gain"] = "2.8"  # from header of image
       cosmic_dict["readnoise"]= "6.4" # from header of image
       cosmic_dict["sigclip"] = "5"
    return cosmic_dict                     
################################################################################
################################################################################

## Change names of the files to a more comprehensive structure of names. 
if not os.path.exists( os.path.join(directory, "cig"):
    print "> Changing names of fits files" 
    list_images = rename.main(arguments=["--copy", "--objectk", objectk,\
                                         "--filterk", filterk, "--datek", datek,\
                                         "--overwrite", "--exptime", exptimek,\
                                         directory])    
#
#
### Remove images in list remove_images.
#print "> Ignore images as selected by user, if any."
#try:
#    for im in remove_images:
#        index = [i for i,v in enumerate(list_images["filename"]) if str(im) in v]
#        for key in list_images.keys():  # remove that item from all the lists
#            list_images[key] = np.delete(list_images[key], index)
#except NameError: # variable remove_images not defined
#    pass
#
#
### Write down homogeneous filter names.
print "> Include homogeneous filter names into the filter keyword"
for im in list_images["filename"]:
    imfilt = utilities.get_from_header(im, filterk)
    imfilt2 = utilities.homogeneous_filter_name(imfilt)
    utilities.header_update_keyword(im, filterk+"_OLD", imfilt, "Original filter name")
    utilities.header_update_keyword(im, filterk,        imfilt2, "Revised filter name")


#
### Create masks for all images. Stars in blanks need to be removed (max_val small). 
print "> Create masks for images"
for ii,im in enumerate(list_images["filename"]):
    # Arguments to be passed to create_masks
    args = ["--max_val", str(max_counts), "--min_val", "0", "--mask_key", "mask", "--outside_val", "2", im]
    if circular_FoV:  # from the campaign file
        args = ["--circular"] + args
    create_masks.main(arguments=args )
    

### Combine bias images
print "> Combining bias images"                
whr = np.where(list_images["type"] == "bias")
bias_images = list(list_images["filename"][whr])
superbias = combine_images.main(arguments=["--average", "median", 
                                           "--all_together", "--output", os.path.join(directory, "superbias.fits"),
                                           "--mask_key", "mask", "--filterk", filterk] + bias_images[:])


### Subtract bias from all images.  
print "> Subtracting bias"
newname = arith.main(arguments=["--suffix", " -b", "--message", 
                                       "BIAS SUBTRACTED", "--mask_key", 
                                       "mask"] + list(list_images["filename"]) +
                                       [ "-", superbias["AllFilters"]])
list_images["filename"][:] = newname
 

## Combine skyflats using blocks to distinguish between sunset and sunrise flats.
print "> Combining sky flats"
skyflat_indices = np.where(list_images["type"] == "skyflats")    
times = list_images["time"][skyflat_indices]  # times of the skyflat images 
block_limits = utilities.group_images_in_blocks(times, limit=20)  
master_skyflats = {}
for ii in range(len(block_limits)-1): 
    block = list_images["filename"][skyflat_indices][block_limits[ii]:block_limits[ii+1]]
    time_block = utilities.mean_datetime(list_images["time"][skyflat_indices]
                                    [block_limits[ii]:block_limits[ii+1]] )
    skyflat = combine_images.main(arguments=["--average", "median", "--norm",
                                           "--scale", "median",
                                           "--output", os.path.join(directory, "masterskyflat{0}".format(ii) + ".fits"), 
                                           "--mask_key", "mask", "--filterk", filterk] + list(block)[:])    
    master_skyflats[time_block] = skyflat.values()


## Dividing the combined flat with a median-filtered version of itself the 
## large scale changes are removed, and only the pixel-to-pixel (p2p) 
## differences remain.         
print "> Creating pixel-to-pixel combined images"
for key,image in master_skyflats.items():
    # first filter with median
    filtered = median_filter.main(arguments= image + ["--mask_key", "mask",
    "--radius", "50"])
    # then divide "image" by "filtered"
    divided = arith.main(arguments=["--suffix", " -small_scale", 
                                           "--message", 
                                           "REMOVE LARGE SCALE STRUCT"] +
                                           image + ["/"] + filtered)
    master_skyflats[key] = divided    

## Combine blanks also in blocks. In this case, we will combine images from  
## every two consecutive blocks, because we only have three blanks per block
## and the dithering is not large enough, so too many residuals were present. 
print "> Combining blanks"
blank_indices = np.where(list_images["type"] == "blanks")[0]  #select blank images
times = list_images["time"][blank_indices]  # times of the blank images 
block_limits = utilities.group_images_in_blocks(times, limit=5)  
master_blanks = {} 
for ii in range(len(block_limits)-1): 
    block = list_images["filename"][blank_indices][block_limits[ii]:block_limits[ii+1]]
    time_block = utilities.mean_datetime(list_images["time"][blank_indices]
                                    [block_limits[ii]:block_limits[ii+1]] )
    blank = combine_images.main(arguments=["--average", "median", "--norm",
                                           "--scale", "median", "--mask_key",
                                           "mask", "--output", os.path.join(directory, "masterblank{0}".format(ii) + ".fits"), 
                                           "--nmin", "2", "--filterk", filterk] + list(block)[:])    
    master_blanks[time_block] = blank.values()


## Use the pixel-to-pixel differences in master_skyflats to correct the 
## master_blanks for this effect. 
print "> Correcting combined blanks for pixel-to-pixel (small scale) variations"
for time, image in master_blanks.items():
    # find closest flat
    time_diff = np.asarray(master_skyflats.keys()) - np.asarray(time)
    closest = np.argmin(abs(time_diff))
    # correct pixel-to-pixel differences (from skyflats)
    corrected = arith.main(arguments=["--suffix", " -sf", "--message",
                                             "REMOVE SMALL SCALE STRUCTURE",
                                             "--mask_key", "mask"]+
                                             image + ["/", master_skyflats.values()[closest]])
    smoothed = median_filter.main(arguments= [corrected , "--mask_key", "mask", "--radius", "150"])
    master_blanks[time] = smoothed



## Now we will correct each image with the closest sky flat field (for small
## scale variations) and the closest blank field (for large scale flatfielding)
print "> Correcting all images from both small scale and large scale flat."
for index in range(len(list_images["filename"])):
    time = list_images["time"][index]
    image = list_images["filename"][index]
    #print image
    # First pixel-to-pixel
    time_diff = np.asarray(master_skyflats.keys()) - time
    closest = np.argmin(abs(time_diff))
    corrected = arith.main(arguments=["--suffix", " -sf", "--message",
                                             "REMOVE SMALL SCALE STRUCTURE",
                                             image] + ["/", master_skyflats.values()[closest]])

    # Now the large scale using the blanks
    time_diff = np.asarray(master_blanks.keys()) - time 
    closest = np.argmin(abs(time_diff))
    corrected = arith.main(arguments=["--suffix", " -bf", "--message",
                                             "REMOVE LARGE SCALE STRUCTURE",
                                             corrected, "/"] + 
                                             master_blanks.values()[closest])
    list_images["filename"][index] = corrected


# Since the areas outside the FoV of a flatfield will be clearly masked with 
# the value 2 (as we indicated when creating the masks), we can use one of 
# them to zero the area outside the FoV for all the images.
if circular_FoV:
    print "> Zero the area outside the FoV (if set as True in 'campaign' file)"
    flat_indices = np.where( (list_images["type"] == "skyflats")  |
                         (list_images["type"] == "flats")     |
                         (list_images["type"] == "domeflats")   )[0]
    # Find a suitable mask, which contains number 2
    for index in flat_indices:
        print index
        print list_images["filename"][index]
        mask_name = fits.getval(list_images["filename"][index], "mask")
        mask = fits.getdata(mask_name)
        if np.any(mask == 2):
            break
    for im in list_images["filename"]:
        image = fits.open(im, mode="update")
        image[0].data[mask==2] = 0
        image.flush()
        image.close()


#### Add World Coordinate System to the images (except bias and skyflats)
#print "> Include WCS"
##print list_images["filename"]
#for index, image in enumerate(list_images["filename"]):
#    #print image
#    if list_images["type"][index] in ("cig", "unknown", "blanks", "standards"):
#        hdr = fits.getheader(image)
#        RA_header, DEC_header, time = hdr[rak], hdr[deck], hdr[datek]
#        RA_header, DEC_header = utilities.sex2deg(RA_header, DEC_header)
#        time = dateutil.parser.parse(time)
#        RA, DEC = utilities.precess_to_2000(RA_header, DEC_header, time)        
#        subprocess.call(['solve-field', "--no-plots", "--overwrite", 
#                         "--no-fits2fits","--scale-units", "arcsecperpix", 
#                         "--scale-low", str(0.98 * pix_scale), "--scale-high", 
#                         str(1.03 * pix_scale), "--quad-size-max", "1.", 
#                         "--quad-size-min", "0.05", "--ra", str(RA), "--dec", 
#                         str(DEC), "--radius", str(FoV), "--depth", "100,250", 
#                         "--solved", "solved.txt", np.str(image)]) 

                         
### Estimation of the sky level for all the images (except for bias)
print "> Estimate sky for images of CIG(s), standard(s)"
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig","standards","unknown"]:
        print image
        find_sky.main(arguments=[list_images["filename"][index]])


# Cosmic rays removal. Dependent on the telescope (campaign keyword: TELESCOP)
print "> Removing cosmic rays from images"
telescope = 'INT'
cosmic_dict = cosmic_removal_param(telescope)  # Read the parameters (gain,readout noise)
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] not in ["bias", "flats"]:        
            newname = remove_cosmics.main(arguments=["--suffix", " -c", "--gain",\
                   cosmic_dict["gain"], "--readnoise", cosmic_dict["readnoise"], \
                   "--sigclip", cosmic_dict["sigclip"], "--maxiter", "3", im])
            list_images["filename"][index] = newname


sys.exit("* * * * * * Reduction finished * * * * * *")
