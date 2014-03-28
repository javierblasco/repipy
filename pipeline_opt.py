# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 09:48:35 2013

@author: blasco
"""

import os, shutil, re, sys, glob
import subprocess
import pyraf.iraf as iraf
import scipy.ndimage.filters as filters
import dateutil.parser
# Advice from Victor Terron in his "lemon setup.py" about how to run mkiraf 
# automatically: 
if os.path.isfile("login.cl") == False:
    p = subprocess.Popen(['mkiraf'], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    out, err = p.communicate(input = 'xgterm')
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


if len(sys.argv) != 2:
    print sys.exit("Give me a campaign file.... See example in the routine "+\
                   "packages")
execfile(sys.argv[1])

print "Rename files"
list_images = rename.main(arguments=["--copy", "--objectk", objectk,\
                                     "--filterk", filterk, "--datek", datek,\
                                     "--overwrite", "--exptime", exptimek,\
                                     directory])

print "List of files after rename in list_files.txt"
# Strip the path from the filenames and calculate the longest of them
file_list = [os.path.split(name)[1] for name in list_images["filename"]]
longest_name = max([len(name) for name in file_list])
f = open("list_files.txt", "w")
for nn, tt in zip(file_list, list_images["time"]):
    nn = nn + " " * (longest_name - len(nn))  # use spaces for padding
    f.write("{}     {}\n".format(nn, tt.isoformat()))
f.close()    

print "Include homogeneous filter names into the filter keyword"
for im in list_images["filename"]:
    imfilt, = utilities.get_from_header(im, filterk)
    imfilt2 = utilities.homogeneous_filter_name(imfilt)
    utilities.header_update_keyword(im, filterk+"_OLD", imfilt2, "Original filter name")
    utilities.header_update_keyword(im, filterk,        imfilt, "Revised filter name")
    

print "Ignore images as selected by user, if any."
for im in remove_images:
    index = [i for i,v in enumerate(list_images["filename"]) if str(im) in v]
    for key in list_images.keys():  # remove that item from all the lists
        list_images[key] = np.delete(list_images[key], index)

print "Create masks for images"
for ii,im in enumerate(list_images["filename"]):
    create_masks.main(arguments=["--max_val", "50000", "--min_val", 
                                 "0", "--mask_key", "mask", "--circular",
                                 "--outside_val", "2", im])   
                                 
print "Combine bias"
whr = np.where(list_images["type"] == "bias")
bias_images = list(list_images["filename"][whr])
print "Bias images", bias_images
superbias = combine_images.main(arguments=["--average", "median", 
                                           "--all_together", 
                                           "-o", "superbias",
                                           "--nlow", "0", 
                                           "--nhigh","0",
                                           "--mask_key", "mask",
                                           "--filterk", filterk] +\
                                           bias_images[:])

print "Superbias:", superbias["AllFilters"]
                                           
print "Subtract bias"
type_of_subtraction = "--median"   # bias has tiny structure, worth subtracting single number
for ii, im in enumerate(list_images["filename"]):
    newname = arith.main(arguments=["--suffix", " -b", "--message", "BIAS SUBTRACTED",
                                "--mask_key", "mask", type_of_subtraction, 
                                im, "-", superbias["AllFilters"]])
    list_images["filename"][ii] = newname[0]  

print "Combine flats"
flat_indices = np.where(list_images["type"] == "skyflats")    
flats = combine_images.main(arguments=["--average", "median", "--norm",
                                         "--scale", "median", "-o", 
                                         "masterskyflat", "--nhigh", "1", 
                                         "--nlow", "0", "--mask_key", "mask", 
                                         "--filterk", filterk] + 
                                         list(list_images["filename"][flat_indices]))   

print "Correct flat-field"
for ii,im in enumerate(list_images["filename"]):
    if list_images["type"][ii] != "bias":
        # Guess the filter of the image from the name, find correct flat
        current_flat = [flats[kk] for kk in flats.keys() if im.count(kk) != 0]
        if len(current_flat) == 0:
            sys.exit("ERROR: Flat for image: " + im + " not found")
        else:
            current_flat = current_flat[0]            
        newname = arith.main(arguments=["--suffix", " -f", "--message", 
                                        "FLAT CORRECTED", "--mask_key", 
                                        "mask", im, "/", current_flat])
        list_images["filename"][ii] = newname[0]  
        
        
print "Zero the area outside the FoV"
# Since the areas outside the FoV of a flatfield will be clearly masked with 
# the value 2 (as we indicated when creating the masks), we can use one of 
# them to zero the area outside the FoV. 
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
             
print "Estimate sky for images of CIG(s), standard(s) and cluster(s) "
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig","standards","clusters"]:
        find_sky.main(arguments=[list_images["filename"][index]])  
        
#print "Removing cosmic rays from images"
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
                         "--use-sextractor", "--code-tolerance", "0.002",
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
        
#        # Move all filtered images related to im to their own directory
#        directory, filename = os.path.split(im) 
#        wcs_dir = os.path.join(directory, "wcs_images")
#        utilities.if_dir_not_exists_create(wcs_dir)
#        
#        filename = utilities.replace_extension(filename, "*")
#        filt_files = glob.glob( os.path.join(directory, filename))
#        utilities.move_list(filt_files, wcs_dir)        


                                                          
print "Estimate seeing for each image"
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig", "standards", "clusters"]:
        print "\n \n"
        print im
        im_cat = utilities.replace_extension(im, "radec")
        estimate_seeing.main(arguments=["--cat", im_cat, "--wcs", "world", im])
        
print "Do photometry for each image"
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["standards", "cig", "clusters"]:
        im_cat = utilities.replace_extension(im, "radec")
        seeing = fits.getval(im, "seeing")
        sigma_sky = fits.getval(im, "sky_std")
        coords_type = "world"   #world coordinate system
        photom_file = im + ".mag.1"
        utilities.if_exists_remove(photom_file)
        iraf.noao()
        iraf.digiphot()
        iraf.apphot()
        iraf.module.phot(im, coords=im_cat, fwhmpsf=seeing, sigma=sigma_sky,
                         datamin=-100, datamax=50000, ccdread=read_noisek,
                         gain=gaink, exposure=exptimek, airmass=airmassk, 
                         filter=filterk, obstime=datek, maxshift=2, 
                         annulus=str(8*seeing), dannulus=str(2*seeing), 
                         apertures=str(4*seeing), zmag=0, radplot="no", 
                         wcsin=coords_type, verify="no", display="no", 
                         interactive="no", icommands="", output=photom_file)

        