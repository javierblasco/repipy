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
import pyraf.iraf as iraf
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
#import repipy.tidy_up2 as tidy_up
import repipy.create_masks as create_masks
import repipy.remove_cosmics as remove_cosmics
import repipy.find_sky as find_sky
import repipy.complete_headers as complete_headers
import repipy.calculate_airmass as calculate_airmass
import repipy.rename as rename
import repipy.median_filter as median_filter
import repipy.cross_match as cross_match
import astropy.io.fits as fits

if len(sys.argv) != 2:
    print sys.exit("Give me a campaign file....")

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

# Combine bias images
print "Combining bias images"                
whr = np.where(list_images["type"] == "bias")
bias_images = list(list_images["filename"][whr])
superbias = combine_images.main(arguments=["--average", "median", 
                                           "--all_together", "-o", "superbias",
                                           "--nlow", "0", "--mask_key", "mask",
                                           "--filterk", filterk] + bias_images[:])    


# Subtract bias from all images.  
print "Subtracting bias"
newname = arith.main(arguments=["--suffix", " -b", "--message", 
                                       "BIAS SUBTRACTED", "--mask_key", 
                                       "mask"] + list(list_images["filename"]) +
                                       [ "-", superbias["AllFilters"]])
list_images["filename"][:] = newname
 
# Combine skyflats using blocks to distinguish between sunset and sunrise flats.
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
                                           "--scale", "median",
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
    "--radius", "50"])
    # then divide "image" by "filtered"
    divided = arith.main(arguments=["--suffix", " -small_scale", 
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
                                               "--filterk", filterk] + 
                                               list(block)[:])    
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
        corrected = arith.main(arguments=["--suffix", " -sf", "--message",
                                                 "REMOVE SMALL SCALE STRUCTURE",
                                                 "--mask_key", "mask"]+
                                                 image + ["/"] + 
                                                 master_skyflats.values()[closest])
        smoothed = median_filter.main(arguments= corrected + [ "--mask_key", "mask",
            "--radius", "150"])
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
    corrected = arith.main(arguments=["--suffix", " -sf", "--message",
                                             "REMOVE SMALL SCALE STRUCTURE",
                                             image] + ["/"] + 
                                             master_skyflats.values()[closest])

    # Now the large scale using the blanks
    time_diff = np.asarray(master_blanks.keys()) - time 
    closest = np.argmin(abs(time_diff))
    corrected = arith.main(arguments=["--suffix", " -bf", "--message",
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

print "Include world coordinate system"
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ("cig", "standards", "clusters"):
        hdr = fits.getheader(image)
        RA, DEC = hdr[rak], hdr[deck]
        subprocess.call(['solve-field', "--no-plots", "--overwrite", 
                         "--no-fits2fits","--scale-units", "arcsecperpix", 
                         "--scale-low", str(0.95 * pix_scale), "--scale-high", 
                         str(1.05 * pix_scale), "--quad-size-max", "0.8", 
                         "--quad-size-min", "0.1", "--ra", str(RA), "--dec", 
                         str(DEC), "--radius", str(FoV), "--depth", "20,50,100",  
                         "--solved", "solved.txt", np.str(image)]) 


sys.exit()

print "Estimate seeing from images"
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ("cig", "standards", "clusters"):
        # Victor Terron has promissed changing dirs will soon be unnecessary 
        curdir = os.path.abspath(os.curdir)
        os.chdir(lemon_dir)
        import lemon.seeing as seeing
        seeing.main(arguments=["--margin", "0", "--filename", '', "--suffix",
                               "-s", image, os.path.split(image)[0] ])
        newname = utilities.add_suffix_prefix(image, suffix = "-s")

        os.chdir(curdir)
        list_images["filename"][index] = newname
            
        # While running lemon.seeing a sextractor catalogue is produced. 
        catalog = fits.getheader(newname)["SEX CATALOG"]
        catalog_newname = utilities.replace_extension(newname, ".cat")
        shutil.copy(catalog, catalog_newname)
        # Damn FITS format and its constraints!
        short_name = os.path.split(catalog_newname)[1]
        utilities.header_update_keyword(newname, "SEX CATALOG", short_name)

print "Estimate sky for images of CIG(s), standard(s) and cluster(s) "
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig","standards","clusters"]:
        find_sky.main(arguments=[list_images["filename"][index]])

print "Detecting objects for images of CIG(s), standard(s) and cluster(s)"
# This part follows the example in the webpage:
#http://www.lancesimms.com/programs/Python/pyraf/daofind.py
for index, image in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig","standards","clusters"]:
        outfile = utilities.replace_extension(image, ".cat")
        if os.path.isfile(outfile):
            os.remove(outfile)
        # Prepare and run daofind
        hdr = fits.getheader(image)
        iraf.noao(_doprint=0)     # Load noao
        iraf.digiphot(_doprint=0) # Load digiphot
        iraf.apphot(_doprint=0)   # Load apphot        
        iraf.daofind.setParam('image',image)
        FWHM = min(max_FWHM,float(hdr["LEMON FWHM"]))
        iraf.daofind.setParam('fwhmpsf', FWHM)       
        iraf.daofind.setParam('output', outfile)
        iraf.daofind.setParam('sigma', float(hdr["SKY_STD"]))
        iraf.daofind.setParam('gain', gaink)
        iraf.daofind.setParam('readnoise', float(hdr[read_noisek]))
        iraf.daofind.setParam('roundlo', -0.3)  # Minimal roundness  (0 is round)
        iraf.daofind.setParam('roundhi', 0.3 )
        iraf.daofind.setParam('airmass', "airmass")
        iraf.daofind.setParam('filter', filterk)
        iraf.daofind.setParam('exposure', exptimek)
        iraf.daofind.setParam('obstime', "date-obs")
        iraf.daofind.setParam('verify', "no")
        iraf.daofind.setParam('datamin', 500)
        iraf.daofind.setParam('datamax', 50000)
        iraf.daofind.saveParList(filename='daofind.par')
        iraf.daofind(ParList='daofind.par')

print "Reading catalogs, calculating shifts"
# List of objects to be aligned. There might be several cigs, several clusters
# and several standard fields.
types_need_aligning = ["cig", "standards","clusters"]
objects_need_aligning = ()
for current_type in types_need_aligning:
    whr = np.where(list_images["type"] == current_type)
    objects_need_aligning = objects_need_aligning +\
                                      tuple(set(list_images["objname"][whr]))
# For each object, read x, y, mag from the sextractor catalog, 
# select the top 15 brightest stars and find the translation between images
for current_object in objects_need_aligning:    
    print "Current object", current_object
    # Open files to store data that imalign will need later on
    obj_list = open(current_object + ".lis", "w")
    shifts_list = open(current_object + ".shifts", "w")    
    coords_list = open(current_object + ".coords","w")
    output_list = open(current_object + ".out", "w")

    # All images of this object. Read fiirst as reference
    whr = np.where(list_images["objname"] == current_object)[0]
    ref_im = list_images["filename"][whr[0]]
    ref_catalog = utilities.replace_extension(ref_im, ".cat")
    
    # Txdump writes to a temporary file, then we read it
    iraf.ptools(_doprint=0)
    if os.path.isfile("temp.txt"):
        os.remove("temp.txt")
    iraf.txdump(ref_catalog, "xcenter, ycenter, mag", "yes", Stdout="temp.txt")
    x_ref=y_ref=mag_ref=np.array([],dtype=np.float16)
    with open("temp.txt") as f:
        for line in f:
            x_ref = np.append(x_ref,float(line.split()[0]))
            y_ref = np.append(y_ref,float(line.split()[1]))
            mag_ref = np.append(mag_ref, float(line.split()[2]))  
    brightest_stars = np.argsort(mag_ref)[:nstars] 
    x_ref = x_ref[brightest_stars] 
    y_ref = y_ref[brightest_stars]
    
    #Write to a file, which imalign will need later
    for ii,jj in zip(x_ref,y_ref):
        coords_list.write( str(ii) + " " + str(jj) + "\n")
    coords_list.close()
            
    # Finally, one by one, calculate the shifts  
    for index in whr:  # for all images of the current_object
        new_im = list_images["filename"][index]
        output = utilities.add_suffix_prefix(new_im, suffix="-a")
        
        # Write input and output in files for imalign 
        obj_list.write(new_im + "\n") 
        output_list.write(output + "\n")        
        
        # Catalog for the new image                
        new_catalog = utilities.replace_extension(new_im, ".cat")
        if os.path.exists("temp.txt"):
            os.remove("temp.txt")
        iraf.txdump(new_catalog, "xcenter, ycenter, mag", "yes", Stdout="temp.txt")
        x_new=y_new=mag_new=np.array([],dtype=np.float16)
        with open("temp.txt") as f:
            for line in f:
                x_new = np.append(x_new,float(line.split()[0]))
                y_new = np.append(y_new,float(line.split()[1]))
                mag_new = np.append(mag_new, float(line.split()[2]))
        brightest_stars = np.argsort(mag_new)[:nstars]
        x_new = x_new[brightest_stars]
        y_new = y_new[brightest_stars]
        os.remove("temp.txt")
        result = cross_match.main(xref=x_ref, yref=y_ref, xobj=x_new, 
                                  yobj=y_new, error=0.01, scale=1, angle=0, 
                                  flip=False, test=False)
        shifts_list.write(str(result[3][0][0]) + " " + str(result[3][0][1]) + "\n")
    shifts_list.close()
    obj_list.close()
    output_list.close()


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
#        newname = arith.main(arguments=["--suffix", " -b", "--message", \
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
#        newname = arith.main(arguments=["--suffix", "f", "--message", \
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
