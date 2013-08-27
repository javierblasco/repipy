#! /usr/bin/env python
"""
@ author: javier blasco herrera
@ email: blasco@iaa.es, javier.blasco.herrera@gmail.com

Pipeline of the reduction of 2013-06-06. 
The data come from the server: /mnt/DATA/GENERAL_CIG/CAHA-DEEP/IM-r/raw/20130606/*
They do not have comprehensible names, but are in sequence and contain raw data. 

"""

import os, shutil, re, sys, glob
import subprocess
# Advice from Victor Terron in his "lemon setup.py" about how to run mkiraf 
# automatically:
if os.path.isfile("login.cl") == False:
    p = subprocess.Popen(['mkiraf'], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    out, err = p.communicate(input = 'xgterm')
import numpy as np
import datetime
import repipy.utilities as utilities
import repipy.combine as combine_images
import repipy.arith as arith_images
import repipy.tidy_up2 as tidy_up
import repipy.create_masks as create_masks
import repipy.remove_cosmics as remove_cosmics
import repipy.find_sky as find_sky
import repipy.complete_headers as complete_headers
import repipy.calculate_airmass as calculate_airmass
import repipy.rename as rename
import repipy.median_filter as median_filter
import astropy.io.fits as fits

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
                 "skyflat_20130606_sdssr_001.fits"] + bad_bias

# Keywords in the header of images
filterk = "INSFLNAM"     # filter name
exptimek = "exptime"    # exposure time (seconds)
objectk = "object"       # name of object 
rak = "ra"               # degree
deck = "dec"             # degree
datek = "date-obs"
telescope = "CAHA2.2"
#timek = "date"      

# Directory with saved master flats. That will save calculating time later on.
saved_dir = "/mnt/data/DEEP_OBS/20130606/saved/"
if os.path.isfile(saved_dir + 'skyflats/masterskyflat1_sdssr-small_scale.fits'):
    print "Reading masterskyflats from:", saved_dir + "skyflats/"
    master_skyflats = {datetime.datetime(2013, 6, 7, 4, 19, 10, 297297): 
                      [saved_dir + 'skyflats/masterskyflat1_sdssr-small_scale.fits'], 
                      datetime.datetime(2013, 6, 6, 20, 2, 38, 200000): 
                      [saved_dir + 'skyflats/masterskyflat0_sdssr-small_scale.fits']}
    list_mastersky = glob.glob(saved_dir + "skyflats/*")
    if not os.path.isdir(directory + "skyflats"):
        os.mkdir(directory + "skyflats")
    for element in list_mastersky:
        shutil.copy(element, directory + "skyflats/")        
          
if os.path.isfile(saved_dir + 'blanks/masterblank0_sdssr-sf-mf.fits'):
    print "Reading masterblanks from:", saved_dir + "blanks/" 
    master_blanks = {datetime.datetime(2013, 6, 7, 2, 24, 44, 333333): 
                     [saved_dir + 'blanks/masterblank2_sdssr-sf-mf.fits'],
                     datetime.datetime(2013, 6, 6, 22, 21, 16, 333333): 
                     [saved_dir + 'blanks/masterblank0_sdssr-sf-mf.fits'], 
                     datetime.datetime(2013, 6, 7, 3, 24, 21, 333333): 
                     [saved_dir + '/blanks/masterblank3_sdssr-sf-mf.fits'], 
                     datetime.datetime(2013, 6, 6, 23, 35, 39, 333333): 
                     [saved_dir + 'blanks/masterblank1_sdssr-sf-mf.fits']}     
    list_masterblank = glob.glob(saved_dir + "blanks/*")
    if not os.path.isdir(directory + "blanks/"):
        os.mkdir(directory + "blanks/")
    for element in list_masterblank:
        shutil.copy(element, directory + "blanks/")


################################################################################
################################################################################


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
    return cosmic_dict                     
################################################################################
################################################################################


# Change names of the files to a more comprehensive structure of names. 
print "Changing names of fits files"
list_images = rename.main(arguments=["--copy", "--objectk", objectk,\
                                     "--filterk", filterk, "--datek", datek,\
                                     "--overwrite", "--exptime", exptimek,\
                                     directory])    
    
# Remove images in list remove_images.
print "Removing images as selected by user, if any."
for im in remove_images:
    index = [i for i,v in enumerate(list_images["filename"]) if str(im) in v]
    for key in list_images.keys():  # remove that item from all the lists
        list_images[key] = np.delete(list_images[key], index)

# Create masks for all images. Stars in blanks need to be removed (max_val small). 
print "Creating masks"
whr = np.where(list_images["type"] != "blanks")
create_masks.main(arguments=["--max_val", "50000", "--circular"] + 
                  list(list_images["filename"][whr]))        

whr = np.where(list_images["type"] == "blanks")
create_masks.main(arguments=["--max_val", "30000", "--min_val", "1000", "--stars",
                             "--circular"] + list(list_images["filename"][whr])) 


# In the case of  

# Combine bias images
print "Combining bias images"                
whr = np.where(list_images["type"] == "bias")
bias_images = list(list_images["filename"][whr])
superbias = combine_images.main(arguments=["--average", "median", "--all_together",
                                           "-o", "superbias", "--nlow", "0", 
                                           "--mask_key", "mask",
                                           "--filterk", filterk] + bias_images[:])    

# Subtract bias from all images.  
print "Subtracting bias"
newname = arith_images.main(arguments=["--suffix", " -b", "--message", 
                                       "BIAS SUBTRACTED", "--mask_key", 
                                       "mask"] + list(list_images["filename"]) +
                                       [ "-", superbias["AllFilters"]])
list_images["filename"] = np.asarray(newname)

 
# Combine skyflats using blocks to distinguish between sunset and sunrise flats.
try:
    dummy = len(master_skyflats)
except NameError:
    print "Combining sky flats"
    skyflat_indices = np.where(list_images["type"] == "skyflats")    
    times = list_images["time"][skyflat_indices]  # times of the skyflat images 
    block_limits = utilities.group_images_in_blocks(times, limit=20)  
    master_skyflats = {}
    for ii in range(len(block_limits)-1): 
        block = list_images["filename"][skyflat_indices][block_limits[ii]:block_limits[ii+1]]
        time_block = utilities.mean_datetime(list_images["time"][skyflat_indices]
                                        [block_limits[ii]:block_limits[ii+1]] )
        skyflat = combine_images.main(arguments=["--average", "median", "--norm",
                                               "--scale", "median", "--notest",
                                               "-o", "masterskyflat{0}".format(ii), 
                                               "--nhigh", "1", "--nlow", "0",
                                               "--mask_key", "mask", "--notest",
                                               "--filterk", filterk] + list(block)[:])    
        master_skyflats[time_block] = skyflat.values()
    
    # Dividing the combined flat with a median-filtered version of itself the 
    # large scale changes are removed, and only the pixel-to-pixel (p2p) 
    # differences remain.         
    print "Creating pixel-to-pixel combined images"
    for key,image in master_skyflats.items():
        print key
        # first filter with median
        filtered = median_filter.main(arguments= image + ["--mask_key", "mask",
        "--side", "150", "--fill_val", "0"])
        # then divide "image" by "filtered"
        divided = arith_images.main(arguments=["--suffix", " -small_scale", 
                                               "--message", 
                                               "REMOVE LARGE SCALE STRUCT"] +
                                               image + ["/"] + filtered)
        master_skyflats[key] = divided    
    print "master_skyflats ={", master_skyflats, "}"
# Combine blanks also in blocks. In this case, we will combine images from  
# every two consecutive blocks, because we only have three blanks per block
# and the dithering is not large enough, so too many residuals were present. 
try:
    dummy = len(master_blanks)
except NameError:
    print "Combining blanks"
    blank_indices = np.where(list_images["type"] == "blanks")[0]  #select blank images
    times = list_images["time"][blank_indices]  # times of the blank images 
    block_limits = utilities.group_images_in_blocks(times, limit=5)  
    master_blanks = {} 
    for ii in range(len(block_limits)-1): 
        block = list_images["filename"][blank_indices][block_limits[ii]:block_limits[ii+1]]
        time_block = utilities.mean_datetime(list_images["time"][blank_indices]
                                        [block_limits[ii]:block_limits[ii+1]] )
        blank = combine_images.main(arguments=["--average", "median", "--norm",\
                                               "--scale", "median", "--mask_key",\
                                               "mask", "-o", 
                                               "masterblank{0}".format(ii), 
                                               "--nhigh", "0", "--nlow", "0",
                                               "--nmin", "2",
                                               "--filterk", filterk] + list(block)[:])    
        master_blanks[time_block] = blank.values()
    
    # Use the pixel-to-pixel differences in master_skyflats to correct the 
    # master_blanks for this effect. 
    print "Correcting combined blanks for pixel-to-pixel (small scale) variations"
    for time, image in master_blanks.items():
        # find closest flat
        time_diff = np.asarray(master_skyflats.keys()) - np.asarray(time)
        closest = np.argmin(abs(time_diff))
        # correct pixel-to-pixel differences (from skyflats)
        print image, master_skyflats.values()[closest]
        corrected = arith_images.main(arguments=["--suffix", " -sf", "--message",
                                                 "REMOVE SMALL SCALE STRUCTURE",
                                                 "--mask_key", "mask"]+
                                                 image + ["/"] + 
                                                 master_skyflats.values()[closest])
        smoothed = median_filter.main(arguments= corrected + [ "--mask_key", "mask",
            "--side", "150", "--fill_val", "0"])
        master_blanks[time] = smoothed
    print "master_blanks = {", master_blanks, "}"    



# Now we will correct each image with the closest sky flat field (for small
# scale variations) and the closest blank field (for large scale flatfielding)
print "Correcting all images from both small scale and large scale flat."
for index in range(len(list_images["filename"])):
    time = list_images["time"][index]
    image = list_images["filename"][index] 
    # First pixel-to-pixel
    time_diff = np.asarray(master_skyflats.keys()) - time
    closest = np.argmin(abs(time_diff))
    corrected = arith_images.main(arguments=["--suffix", " -sf", "--message",
                                             "REMOVE SMALL SCALE STRUCTURE",
                                             image] + ["/"] + 
                                             master_skyflats.values()[closest])
    # Now the large scale using the blanks
    time_diff = np.asarray(master_blanks.keys()) - time 
    closest = np.argmin(abs(time_diff))
    corrected = arith_images.main(arguments=["--suffix", " -bf", "--message",
                                             "REMOVE LARGE SCALE STRUCTURE"] +
                                             corrected + ["/"] + 
                                             master_blanks.values()[closest])
    list_images["filename"][index] = corrected[0]
#print "Removing cosmic rays from images"
#cosmic_dict = cosmic_removal_param(telescope)  # Read the parameters (gain,readout noise)
#for index, im in enumerate(list_images["filename"]):
#    if list_images["type"][index] not in ["bias", "flats"]:        
#            newname = remove_cosmics.main(arguments=["--suffix", " -c", "--gain",\
#                   cosmic_dict["gain"], "--readnoise", cosmic_dict["readnoise"], \
#                   "--sigclip", cosmic_dict["sigclip"], "--maxiter", "3", im])
#            list_images["filename"][index] = newname

print "Estimate seeing from images"
print set(list_images["type"])
sys.exit()
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ("cig", "standards", "clusters"):
        # Victor Terron has promissed changing dirs will soon be unnecessary 
        curdir = os.path.abspath(os.curdir)
        os.chdir(lemon_dir)
        import lemon.seeing as seeing
        seeing.main(arguments=["--margin", "0", "--filename", '', "--suffix",
                               "-s", image, os.path.split(image)[0] ])
        newname = utilities.add_suffix_prefix(image, suffix = "-s")
        print "\n newname", newname, "\n"
        os.chdir(curdir)
        list_images["filename"][index] = newname
            
        # While running lemon.seeing a sextractor catalogue is produced. 
        catalog = fits.getheader(newname)["SEX CATALOG"]
        catalog_newname = utilities.replace_extension(newname, ".cat")
        catalog_newname = os.path.split(catalog_newname)[1]
        shutil.copy(catalog, catalog_newname)
        utilities.header_update_keyword(newname, "SEX CATALOG", catalog_newname)

print "Aligning images of the CIG(s), standards and cluster(s)"
# List of objects to be aligned. There might be several cigs, several clusters
# and several standard fields.
types_need_aligning = ["cig", "standards","cluster"]
objects_need_aligning = ()
for current_type in types_need_aligning:
    whr = np.where(list_images["type"] == current_type)
    objects_need_aligning = objects_need_aligning + tuple(set(list_images["object"][whr]))

# For each object, read x_image, y_image, mag_auto from the sextractor catalog, 
# select the top 20 brightest stars and find the translation between images
for current_object in object_need_aligning:    
    whr = np.where(list_images["object"] == current_object)[0]
    for ii in whr:
        print current_object, list_images["object"][ii]
    

sys.exit()












































## Observing date (begining of the night). This will be added to the file names
## in order to distinguish observations of the same galaxy in different nights.
#date = "20130606"
## Year/telescope/observatory/location of observations. This will be used to get 
## some values such as the names of the keywords in the header, and calculate the 
## airmasses at the moment of observations...
#year = "2013"
#telescope = "CAHA2.2"
#observatory = "CAHA"
#location = "Europe/Madrid"  # Place at the same time zone as the observatory 
#                            # and important enough to be in pytz (python time zone)
#
## Patterns of the images as they are in the raw directory
#def original_pattern_dict():
#    in_pattern = {}
#    in_pattern["cig"] = "^(?P<name>cig)(?P<cig_num>\d{3,4})-(?P<exp_num>\d{3})" +\
#                        "(?P<filt>.{3})(?P<rest>\.fit)"
#    in_pattern["bias"] = "^(?P<name>bias)-(?P<exp_num>\d{3})(?P<rest>\.fit)"
#    in_pattern["flats"] = "^(?P<name>flat)-(?P<exp_num>\d{3})" + \
#                          "(?P<filt>.{1,3})(?P<rest>\.fit)"
#    in_pattern["standards"] = "^(?P<name>feige34|he3|hz15)-(?P<exp_num>\d{3})" +\
#                               "(?P<filt>.{3})(?P<rest>\.fit)"    
#    return in_pattern
#
## CORRECT TYPO IN THE NAMES AT THE VERY BEGINING
#if os.path.isfile(directory+"feige-005H52.fit") == True:
#    shutil.move(directory+"feige-005H52.fit", "feige34-005H52.fit")
#if os.path.isfile(directory+"feige-005rGu.fit") == True:
#    shutil.move(directory+"feige-005rGu.fit", "feige34-005rGu.fit")
#
#
#
#
#
#
#################################################################################
##                   Campaign "independent"                                     #
#################################################################################
#
#def cosmic_removal_param(telescope = ''):
#    if telescope == "OSN":
#        cosmic_dict = {}
#        cosmic_dict["gain"] = "2"
#        cosmic_dict["readnoise"]= "5"
#        cosmic_dict["sigclip"] = "5"        
#    return cosmic_dict
#
# 
## Read patterns and tidy up directory. The resulting variable List_images is a 
## dictionary of numpy arrays that contains the file names, the type (cig, flats, 
## standards or bias) and the object name (e.g. cig0210, flat, feige34, ...)
#pattern = original_pattern_dict()
#list_images = tidy_up.tidy_up(directory, pattern, date)
#
## If the array "filename" in the dictionary list_images is empty it is because 
## the directory is already tidy. Just read the regular expressions for tidy 
## directories (basically "name_date_filter_expnum.fits) and read the images
## that agree with any of those patterns.  
#newpattern = tidy_pattern_dict()
#if len(list_images["filename"]) == 0:
#    list_images = utilities.locate_images2(directory, newpattern)
#
#
## In case we need to exclude images, now is the moment to do it. Bias that are 
## not stable, flats with too few or too many counts, images too saturated...
## We compare all the images with  
#print "Removing images as selected by user, if any."
#for im in remove_images: 
#        matching = [s for s in list_images["filename"] if im in s]
#        if len(matching) > 0: 
#            index = np.where(list_images["filename"] == matching[0])[0] 
#            for array in list_images:
#                list_images[array] = np.delete(list_images[array], index)
#
#
## Load a dictionary for the names and formats of the keywords in your header. 
## Make sure your telescope and year are in the routine hdr_keywords of the 
## reduction pipeline. Thus, when we say, for example, hdr_keywords["OBJECT"] 
## the dictionary will tell us the name of the keyword in the header that contains
## the object name, for example "OBJ" or "OBJNAME" or whichever that particular
## telescope and instrument included in the header. 
#hdr_keywords = header_keywords.main(telescope, year)
#
#
## Include missing details in headers. For each image this includes information like 
## the observatory, the coordinates of the object in the image, the local sidereal
## time... This info will be used to calculate the airmass and the astrometry       
#print "Including details in the header"
#for index, im in enumerate(list_images["filename"]):
#    if list_images["type"][index] not in ["bias", "flats"]:  # exclude bias/flats
#        obj_name = list_images["objname"][index]
#        complete_headers.main(arguments=["--observatory", observatory, "--object",\
#                         obj_name, "--RA_keyword", hdr_keywords["ra"][0], \
#                         "--DEC_keyword", hdr_keywords["dec"][0], im])
#
## Calculate airmass using RA, DEC, time, location... 
#print "Including/recalculating airmass in the headers"
#for index, im in enumerate(list_images["filename"]):
#    if list_images["type"][index] not in ["bias", "flats"]:  # exclude bias/flats
#        calculate_airmass.main(arguments=["--RA", "RA_hours", "--DEC", "DEC_deg",\
#                                       "--observatory", observatory, "--st", "ST",\
#                                       "--date", hdr_keywords["date"][0], "--ut",\
#                                       hdr_keywords["UT_time"][0], "--equinox",\
#                                       hdr_keywords["equinox"][0], "--location",\
#                                       location, im])
#
## Create masks for all the images of cigs or standards. Mask out, among other 
## things, saturated pixels or the outside part of the round FoV of CAHA. 
#print "Creating masks"
##whr = np.where((list_images["type"] != "bias") & (list_images["type"] != "flats"))
##create_masks.main(arguments=list_images["filename"][whr])       
#   
## Combine bias. First remove any wrong bias and then combine using your 
## favourite parameters.  
#    
## Bias has little structure, but a horizontal bar with a difference of ~1 count
## is evident. We better subtract it all instead of just making an average. Also
## update the list of images to the ones with the bias subtracted.
#print "Subtracting bias"
#for index, im in enumerate(list_images["filename"]):
#        newname = arith_images.main(arguments=["--suffix", " -b", "--message", \
#                                   "BIAS SUBTRACTED", im, "-", superbias["AllFilters"]])
#        list_images["filename"][index] = newname[0]
#            
## Now we check and combine the flats. There are flats at the beginning and at the 
## end of the night. One of them is too full of stars and another has some ~200 
## counts per pixel, so we remove them. 
#print "Combining flat images"
#whr = np.where(list_images["type"] == "flats")
#flat_images = list(list_images["filename"][whr])
#superflats = combine_images.main(arguments=["--average", "median", "--scale", \
#                                "median", "-o", "superflat_" + date , "--notest",\
#                                "--nhigh", "1", "--nlow", "0", "--filterk", "filter",\
#                                "--norm"] + flat_images[:])
#
## Perform the flat-field correction to standards and cig images. For that, loop
## through all the list of images, find (using the patterns) the appropriate 
## filter name, from it the correct superflat and then do the arith to divide 
## by it.  
#print "Flatfield correction"
#for index, im in enumerate(list_images["filename"]):
#    current_type = list_images["type"][index]
#    if current_type not in ["bias"]: 
#        data_im = re.match("^(?P<name>.*)_(?P<date>\d{8})_(?P<filt>.*)_" +\
#                         "(?P<exp_num>)(?P<rest>.*)$", os.path.split(im)[1])
#        current_filter = data_im.groupdict()["filt"]
#        current_flat = superflats[current_filter]
#        newname = arith_images.main(arguments=["--suffix", "f", "--message", \
#                                    "FLAT-FIELD CORRECTION:", im, "/", \
#                                    current_flat])
#        list_images["filename"][index] = newname[0]  
#         
## Removal of cosmic rays before doing any photometry to avoid a hit close to a 
## star distorting the photometry
#print "Cosmic ray removal"
#cosmic_dict = cosmic_removal_param(telescope)  # Read the parameters (gain,readout noise)
#for index, im in enumerate(list_images["filename"]):
#    if list_images["type"][index] not in ["bias", "flats"]:        
#            newname = remove_cosmics.main(arguments=["--suffix", "c", "--gain",\
#                   cosmic_dict["gain"], "--readnoise", cosmic_dict["readnoise"], \
#                   "--sigclip", cosmic_dict["sigclip"], "--maxiter", "3", im])
#            list_images["filename"][index] = newname
#
#
## Evaluate sky, put it in the header. 
##print "Estimation of sky value"
##for index, im in enumerate(list_images["filename"]):
##    if list_images["type"][index] not in ["bias", "flats"]:
##        sky = find_sky.main(arguments=[im])
#
#
## Find PSF for all images using the lemon pipeline (Victor Terron, vterron@iaa.es)
## to estimate the FWHM of each image individually, and copy the result (called 
## )
#print "Finding FWH"
#for index, im in enumerate(list_images["filename"]):
#    if list_images["type"][index] not in ["bias", "flats"]:
#            curdir = os.path.abspath(os.curdir)
#            os.chdir(lemon_dir)
#            import lemon.seeing as seeing
#            seeing.main(arguments=["--margin", "0", im, os.path.split(im)[0] ])
#            os.chdir(curdir)
#            # Output from lemon seeing is, by default, best_seeing.fits
#            output = os.path.join(os.path.split(im)[0], "best_seeing.fits")
#            shutil.move(output, im)
#
##print list_images["filename"]
##sys.exit()
#
#
#
## For each of the cigs/standards, calculate the offsets that will align the stars
## of the images and do photometry.  
#print " Aligning images and doing photometry"
#objects = set(list_images["objname"])
#for obj in objects:  # For each of the objects: cig0123, cig1234, he3, hz15, ...
#    if obj not in ["bias", "flats"]:
#        curdir = os.path.abspath(os.curdir)
#        os.chdir(lemon_dir) # We need to change to lemon directory
#        import lemon.photometry as photometry
#        import lemon.offsets as offsets
#        whr = np.where(list_images["objname"] == obj)
#        current_list = list(list_images["filename"][whr])
#        name_xml = os.path.join(curdir, obj + "_offsets.xml")
#        offsets.main(arguments=["--filterk", hdr_keywords["filter"][0], "--airmk", 
#                                "airmass", "--datek", "utmiddle", "--overwrite",
#                                "--output", name_xml,
#                                current_list[0]] + current_list[:])
#        name_DB = os.path.join(curdir,obj + "LEMONdbB")                      
#        photometry.main(arguments=["--output", name_DB, "--margin", "0",
#                                    "--aperture", "5", "--annulus", "5.5",
#                                    "--overwrite", "--gain", 
#                                    cosmic_dict["gain"], "--uik", "org_name", 
#                                    name_xml])
#        os.chdir(curdir) # get back to working directory.
#
## From the photometry of the standards, calculate a value for the extinction 
## coefficient, by fitting the magnitude as a function of the airmass for 
## brightest stars in the image 
#
#
#
#sys.exit()
