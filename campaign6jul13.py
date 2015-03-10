# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:53:55 2013

@author: pablo
"""


print "############################################################################"
print "                     Observation date: 6 July 2013"
print "                            Target: CIG0756"
print "############################################################################"

# Work directory:
directory = "/mnt/DATA/pablo/observations/optical/CAHA/20130706/"

# Directory where the lemon pipeline (https://github.com/vterron/lemon) is
lemon_dir = "/home/pablo/Desktop/work/routines/lemon/"

# Images to be removed:
# # # # DEPEND ON THE CAMPAIGN! # # #
# 
# BIAS: very low difference between the 8 bias. None discarded.
# 
# FOCUS: there is none.
#
# SKYFLAT: surprisingly, there are 19 completely saturated flat fields... the
# rest have a high number of counts but can be used. We remove: 
# "skyflat_20130706_sdssr_001.fits"- 007, 011-013, 017-022, 025-027.
# 
# DOMEFLATS: there are none.
#
# BLANKS: 010 has a very high count number. 011 and 012 are saturated.
#
# TARGET CIG0756: they all seem ok.
#
# STANDARDS: we remove two "sao110_20130706_sdssr_005.fits" and 006 because they 
# present an awkward number of counts.
#
bad_skyflats = ["skyflat_20130706_sdssr_" + str(ii).zfill(3) + ".fits" for ii in range(1,8)+range(11,14)+range(17,23)+range(25,28)]
remove_images = ["blank_20130706_sdssr_010.fits", "blank_20130706_sdssr_011.fits", \
                "blank_20130706_sdssr_012.fits", "sao110_20130706_sdssr_005.fits", \
                "sao110_20130706_sdssr_006.fits"] + bad_skyflats

# Header keywords:
filterk = "INSFLNAM"     # filter name
exptimek = "EXPTIME"     # exposure time (seconds)
objectk = "OBJECT"       # name of object 
rak = "RA"               # degree
deck = "DEC"             # degree
datek = "DATE-OBS"       # observation date
telescope = "CAHA2.2"

# Directory with saved master flats. That will save calculating time later on.
#saved_dir = "/home/pablo/Desktop/work/observations/optical/CAHA/20130706/saved/"
#if os.path.isfile(saved_dir + 'skyflats/masterskyflat1_sdssr-small_scale.fits'):
#    master_skyflats = {datetime.datetime(2013, 6, 7, 4, 19, 10, 297297): 
#                      [saved_dir + 'skyflats/masterskyflat1_sdssr-small_scale.fits'], 
#                      datetime.datetime(2013, 6, 6, 20, 2, 38, 200000): 
#                      [saved_dir + 'skyflats/masterskyflat0_sdssr-small_scale.fits']}
#    list_mastersky = glob.glob(saved_dir + "skyflats/*")
#    if not os.path.isdir(directory + "skyflats"):
#        os.mkdir(directory + "skyflats")
#    for element in list_mastersky:
#        shutil.copy(element, directory + "skyflats/")        
#          
#if os.path.isfile(saved_dir + 'blanks/masterblank2_sdssr-sf-mf.fits'):
#    master_blanks = {datetime.datetime(2013, 7, 7, 2, 37, 14): 
#                    [saved_dir + 'blanks/masterblank2_sdssr-sf-mf.fits'], 
#                     datetime.datetime(2013, 7, 6, 23, 10, 19, 333333): 
#                    [saved_dir + 'blanks/masterblank0_sdssr-sf-mf.fits'], 
#                    datetime.datetime(2013, 7, 7, 0, 27, 27, 333333): 
#                    [saved_dir + 'blanks/masterblank1_sdssr-sf-mf.fits']} 
#    list_masterblank = glob.glob(saved_dir + "blanks/*")
#    if not os.path.isdir(directory + "blanks/"):
#        os.mkdir(directory + "blanks/")
#    for element in list_masterblank:
#        shutil.copy(element, directory + "blanks/")