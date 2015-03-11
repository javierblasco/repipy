# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:00:51 2013

@author: pablo
"""

print "############################################################################"
print "                     Observation date: 7 July 2013"
print "                      Target: CIG0756 and CIG0841"
print "############################################################################"

# Work directory:
directory = "/mnt/DATA/pablo/observations/optical/CAHA/20130707/"

# Directory where the lemon pipeline (https://github.com/vterron/lemon) is
lemon_dir = "/home/pablo/Desktop/work/routines/lemon/"

# Images to be removed:
# # # # DEPEND ON THE CAMPAIGN! # # #
# 
# BIAS: same issue with the bias, the average numbers vary during the night. I
# remove "bias_20130707_001.fits" and 002 since they have the most distant means.
# 
# FOCUS: we remove a "test_20130707_sdssr_001.fits" and "none_20130707_sdssr_001.fits"
# images that likely are test images.
#
# SKYFLAT: "skyflat_20130707_sdssr_002.fits", 003 and 069 show a very high counts
# number. We remove them.
# 
# DOMEFLATS: there are none.
#
# BLANKS: 014 presents a too high counts number. 013 shows a satellite trace, we
# remove them.
# 
# TARGET CIG0756: they all seem ok.
#
# TARGET CIG0841: from "cig0841_20130707_sdssr_021.fits" to 025 the sky turns 
# progressively brighter: from 7k-8k counts in most of the images, to ~12k and
# up to ~14k. However, the galaxy keeps having an ascending count number towards 
# the center so I do not remove them though.
#
# STANDARDS: "sao110_20130707_sdssr_006.fits" and 007 present a higher sky count 
# number but the stars are still brighter so I don't remove them.
#
# bad_skyflats = ["skyflat_20130706_sdssr_" + str(ii).zfill(3) + ".fits" for ii in range(1,8)+range(11,14)+range(17,23)+range(25,28)]

remove_images = ["bias_20130707_001.fits", "bias_20130707_002.fits", \
                "blank_20130707_sdssr_14.fits", "blank_20130707_sdssr_13.fits", "skyflat_20130707_sdssr_002.fits", \
                "skyflat_20130707_sdssr_003.fits", "skyflat_20130707_sdssr_069.fits", \
                "test_20130707_sdssr_001.fits", "none_20130707_sdssr_001.fits" ]
                
# Header keywords:
filterk = "INSFLNAM"     # filter name
exptimek = "EXPTIME"     # exposure time (seconds)
objectk = "OBJECT"       # name of object 
rak = "RA"               # degree
deck = "DEC"             # degree
datek = "DATE-OBS"       # observation date
telescope = "CAHA2.2"


# Directory with saved master flats. That will save calculating time later on.
saved_dir = "/home/pablo/Desktop/work/observations/optical/CAHA/20130707/saved/"
if os.path.isfile(saved_dir + 'skyflats/masterskyflat1_sdssr-small_scale.fits'):
    master_skyflats = {datetime.datetime(2013, 6, 7, 4, 19, 10, 297297): 
                      [saved_dir + 'skyflats/masterskyflat1_sdssr-small_scale.fits'], 
                      datetime.datetime(2013, 6, 6, 20, 2, 38, 200000): 
                      [saved_dir + 'skyflats/masterskyflat0_sdssr-small_scale.fits']}
    list_mastersky = glob.glob(saved_dir + "skyflats/*")
    if not os.path.isdir(directory + "skyflats"):
        os.mkdir(directory + "skyflats")
    for element in list_mastersky:
        shutil.copy(element, directory + "skyflats/")        
          
if os.path.isfile(saved_dir + 'blanks/masterblank2_sdssr-sf-mf.fits'):
    master_blanks = {datetime.datetime(2013, 6, 7, 2, 24, 44, 333333): 
                    [saved_dir + 'blanks/masterblank2_sdssr-sf-mf.fits'], 
                     datetime.datetime(2013, 6, 6, 22, 21, 16, 333333): 
                    [saved_dir + 'blanks/masterblank0_sdssr-sf-mf.fits'], 
                    datetime.datetime(2013, 6, 6, 23, 35, 39, 333333): 
                    [saved_dir + 'blanks/masterblank1_sdssr-sf-mf.fits']} 
    list_masterblank = glob.glob(saved_dir + "blanks/*")
    if not os.path.isdir(directory + "blanks/"):
        os.mkdir(directory + "blanks/")
    for element in list_masterblank:
        shutil.copy(element, directory + "blanks/")