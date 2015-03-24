#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 09:48:35 2013

@author: blasco
"""

import os, shutil, re, sys, glob
import subprocess
import scipy.ndimage.filters as filters
import dateutil.parser
from lemon import methods
import numpy as np
import datetime
import repipy.utilities as utilities
import repipy.combine as combine_images
import repipy.arith as arith
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
import warnings
warnings.filterwarnings("ignore")

from lemon import methods

import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot



if len(sys.argv) != 2:
    print sys.exit("Give me a campaign file.... See example in the routine "+\
                   "packages")
execfile(sys.argv[1])


def search_images(dir="."):
    """ Given the patterns for the CIGs, cluster and standard images that we
    need, read them and sort them into a dictionary"""
    in_pattern = {}
    in_pattern["bias"] = "^bias_(?P<date>\d{8})_(?P<exp_num>\d{3}).fits$"
    in_pattern["skyflats"] = "^skyflat_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"
    in_pattern["flats"] = "^flat_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"
    in_pattern["domeflats"] = "^domeflat_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"

    in_pattern["blanks"] = "^blank_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"
    in_pattern["cig"] = "^(?P<name>cig)(?P<cig_num>\d{4})_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"
    in_pattern["standards"] = "^(?P<name>" + standards_campaign +")_(?P<date>\d{8})_" +\
                              "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"
    in_pattern["blanks"] = "^(?P<name>blank)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3}).fits$"
    in_pattern["clusters"] = "^(?P<name>" + clusters_campaign +")_(?P<date>\d{8})_" +\
                             "(?P<filt>.*)_(?P<exp_num>\d{3}).fits$"
    list_images = utilities.locate_images2(dir, in_pattern)
    return list_images


if os.path.exists( os.path.join(directory, 'cig') ):
    list_images = search_images(directory)
else:
    print "Rename files"
    list_images = rename.main(arguments=["--copy", "--objectk", objectk,\
                                         "--filterk", filterk, "--datek", datek,\
                                         "--overwrite", "--exptime", exptimek,\
                                         directory])

print "Include homogeneous filter names into the filter keyword"
for im in list_images["filename"]:
    imfilt = utilities.get_from_header(im, filterk)
    imfilt2 = utilities.homogeneous_filter_name(imfilt)
    utilities.header_update_keyword(im, filterk+"_OLD", imfilt, "Original filter name")
    utilities.header_update_keyword(im, filterk,        imfilt2, "Revised filter name")
    

print "Ignore images as selected by user, if any."
try:
    for im in remove_images:
        index = [i for i,v in enumerate(list_images["filename"]) if str(im) in v]
        for key in list_images.keys():  # remove that item from all the lists
            list_images[key] = np.delete(list_images[key], index)
except NameError: # variable remove_images not defined
    pass


print "Create masks for images"
for ii,im in enumerate(list_images["filename"]):
    # Arguments to be passed to create_masks
    args = ["--max_val", str(max_counts), "--min_val", "0", "--mask_key", "mask", "--outside_val", "2", im]
    if circular_FoV:  # from the campaign file
        args = ["--circular"] + args
    create_masks.main(arguments=args )
                                 
print "Combine bias"
whr = np.where(list_images["type"] == "bias")
bias_images = list(list_images["filename"][whr])
print "Bias images", bias_images
output_bias = os.path.join(directory, "superbias.fits")
superbias = combine_images.main(arguments=["--average", "median", 
                                           "--all_together", 
                                           "--output", output_bias,
                                           "--mask_key", "mask",
                                           "--filterk", filterk] +\
                                           bias_images[:])
                                          
print "Subtract bias"
for ii, im in enumerate(list_images["filename"]):
    if list_images["type"][ii] not in ["bias", "unknown"]:
        args = ["--suffix", " -b", "--message", "BIAS SUBTRACTED", "--mask_key", "mask", im, "-", superbias["AllFilters"]]
        if type_of_bias_subtraction:
            args = [type_of_bias_subtraction] + args
        newname = arith.main(arguments=args)
        list_images["filename"][ii] = newname

print "Combine flats"
output_flats = os.path.join(directory, "masterskyflat.fits")
flat_indices = np.where(list_images["type"] == "skyflats")    
flats = combine_images.main(arguments=["--average", "median", "--norm",
                                         "--scale", "median", "--output",
                                         output_flats, "--mask_key", "mask",
                                         "--filterk", filterk] + 
                                         list(list_images["filename"][flat_indices]))   

print "Correct flat-field"
for ii,im in enumerate(list_images["filename"]):
    if list_images["type"][ii] not in ["bias", "unknown"]:
        # Guess the filter of the image from the name, find correct flat
        current_flat = [flats[kk] for kk in flats.keys() if im.count(kk) != 0]
        if len(current_flat) == 0:
            sys.exit("ERROR: Flat for image: " + im + " not found")
        else:
            current_flat = current_flat[0]            
        newname = arith.main(arguments=["--suffix", " -f", "--message", 
                                        "FLAT CORRECTED", "--mask_key", 
                                        "mask", im, "/", current_flat])
        list_images["filename"][ii] = newname 
        
        
print "Zero the area outside the FoV"
# Since the areas outside the FoV of a flatfield will be clearly masked with 
# the value 2 (as we indicated when creating the masks), we can use one of 
# them to zero the area outside the FoV for all the images.
if circular_FoV:
    flat_indices = np.where( (list_images["type"] == "skyflats")  |
                         (list_images["type"] == "flats")     |
                         (list_images["type"] == "domeflats")   )[0]
    # Find a suitable mask, which contains number 2
    for index in flat_indices:
        mask_name = fits.getval(list_images["filename"][index], "mask")
        mask = fits.getdata(mask_name)
        if np.any(mask == 2):
            break
    for im in list_images["filename"]:
        image = fits.open(im, mode="update")
        image[0].data[mask==2] = 0
        image.flush()
        image.close()
             
print "Estimate sky for images of CIG(s), standard(s) and cluster(s) "
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig","standards","clusters"]:
        find_sky.main(arguments=[list_images["filename"][index]])  
        
print "Removing cosmic rays from images"
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig", "standards", "clusters"]:        
            newname = remove_cosmics.main(arguments=["--suffix", " -c", 
                                          "--gain", str(fits.getval(im, gaink)), 
                                          "--readnoise", str(fits.getval(im, read_noisek)),
                                          "--sigclip", "5", "--maxiter", "3", im])
            list_images["filename"][index] = newname

#print "Calculate smoother version of images"
#for index, im in enumerate(list_images["filename"]):
#    if list_images["type"][index] in ["cig", "standards", "clusters"]:    
#        imdata = fits.getdata(im)
#        #im_median = np.median(imdata)
#        #im_MAD = np.median(np.abs(imdata - im_median))
#        #imdata[imdata < im_median + 6. * im_MAD] = im_median
#        imfiltered = filters.gaussian_filter(imdata, 2)
#        filtered_name = utilities.add_suffix_prefix(im, prefix="filtered-")
#        fits.writeto(filtered_name, imfiltered)

print "Include WCS"
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig", "standards", "clusters"]:
        # Get date, RA, DEC from images. 
        time, RA_current, DEC_current = utilities.get_from_header(im, datek,
                                                                  rak,deck)
        RA, DEC = utilities.precess_to_2000(RA_current, DEC_current, time)
        RA, DEC, radius = str(RA), str(DEC), str(FoV/2.5)   
    
        # Calculate WCS using sextractor
        subprocess.call(["solve-field", "--no-plots", "--no-fits2fits", 
                         "--ra", RA, "--dec", DEC, "--radius", radius,
                         "--depth", "1-30", "--depth", "1-50", "--depth", 
                         "1-100", "--depth", "10,20,30,40,50,60,70,80,90,100",
                         "--use-sextractor", "--code-tolerance", "0.01",
                         "--overwrite", im])                         

        solved = utilities.replace_extension(im, "solved")
        # Case sextractor couldn't do it, try with astrometry.net's own routine
        if not os.path.exists(solved):
            subprocess.call(["solve-field", "--no-plots", "--no-fits2fits", 
                             "--ra", RA, "--dec", DEC, "--radius", radius,
                             "--depth", "1-30", "--depth", "1-50", "--depth", 
                             "1-100", "--depth", "10,20,30,40,50,60,70,80,90,100",
                             "--code-tolerance", "0.002", "--overwrite",
                             im])                         
        output_name = utilities.replace_extension(im, "new")
        shutil.move(output_name, im)


        # Update old WCS system PC matrices and so on to avoid confusion
        utilities.update_WCS(im, im)
 
        # Build a catalogue with the coordinates of the stars found in 
        # X and Y, RA and DEC and the flux of the stars
        corr_file = utilities.replace_extension(im, "corr")
        table = fits.open(corr_file)[1]
        cat_radec = utilities.replace_extension(im, "radec")
        f = open(cat_radec, "w")
        for line in table.data:
            f.write(str(line[2]) + " " + str(line[3]) + "\n")            
        f.close()
                                                                
print "Estimate seeing for each image"
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig", "standards", "clusters"]:
        im_cat = utilities.replace_extension(im, "radec")
        estimate_seeing.main(arguments=["--cat", im_cat, "--wcs", "world", im])

