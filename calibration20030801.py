#! /usr/bin/env python
# -*- coding: utf-8 -*-

from repipy.calculate_extinction import calculate_extinction
from repipy.calculate_zeropoint import calculate_zeropoint
import repipy.arith as arith_images
import repipy.utilities as utilities
import repipy.extract_mag_airmass_common as extract
import repipy.match_psfs as match_psfs
import repipy.utilities as utils
import repipy.combine as combine_images
import repipy.astrometry as astrometry
from repipy.objects import astronomical_object

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.optimize as optimize
import sys
import astropy.io.fits as fits
import os
import shutil
import pickle
import pyraf.iraf as iraf

import lemon.photometry as photometry

if len(sys.argv) != 2:
    print sys.exit("Give me a campaign file.... See example in the routine "+\
                   "packages")
execfile(sys.argv[1])

def search_images(dir="."):
    """ Given the patterns for the CIGs, cluster and standard images that we
    need, read them and sort them into a dictionary"""
    in_pattern = {}
    in_pattern["cig"] = "^(?P<name>cig)(?P<cig_num>\d{4})_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3})(?P<rest>.*-c\.fits)$"
    in_pattern["bias"] = "^(?P<name>bias)_(?P<date>\d{8})_"+\
                         "(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    in_pattern["skyflats"] = "^(?P<name>skyflat)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3})(?P<rest>.*-b\.fits)$"
    in_pattern["standards"] = "^(?P<name>" + standards_campaign +")_(?P<date>\d{8})_" +\
                              "(?P<filt>.*)_(?P<exp_num>\d{3})(?P<rest>.*-c\.fits)$"
    in_pattern["blanks"] = "^(?P<name>blank)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3})(?P<rest>.*-c\.fits)$"
    list_images = utilities.locate_images2(dir, in_pattern)
    return list_images



list_images = search_images(directory)

print "For each object and filter, do photometry on all images"
AllObj_set = set(x for x,t in zip(list_images["objname"], list_images["type"]) if t in ["standards", "cig", "clusters"])
for obj in AllObj_set: # for each of the standards
    obj_images = list(list_images["filename"][np.where( (list_images["objname"] == obj))])
    output_db = os.path.join(directory, obj+".db")
    utilities.if_exists_remove(output_db)
    coord_file = utilities.replace_extension(obj_images[0], "radec")
    photometry.main(arguments=["--maximum", str(max_counts), "--uik", "", "--margin", "20", "--gaink", gaink,
                              "--aperture", "4.", "--annulus", "6", "--dannulus", "2", "--individual-fwhm", "--objectk", objectk,
                              "--filterk", filterk, "--datek", datek, "--expk", exptimek, "--fwhmk", seeingk, "--airmk", airmassk,
                              "--coordinates", coord_file, obj_images[0]] + obj_images + [output_db])

print "Find extinction coefficient from photometry of standards"
standards_set = set(x for x,t in zip(list_images["objname"], list_images["type"]) if t == "standards")
for obj in standards_set:
    airmasses, magnitudes, filters = extract.main(output_db)
    for filt in set(filters.flatten()):
        # find columns with the current filter
        columns = [ii for ii in range(len(filters[0, :])) if filters[0,ii] == filt]
        # if there are more than two images
        if len(airmasses[0,columns]) > 0 and len(airmasses[0,columns]) > 2:
            ext_coeff, sigma_ext_coeff = calculate_extinction(airmasses[:, columns], magnitudes[:, columns])
            print "Extinction coefficient for ", obj, " for filter ", filt, ":", ext_coeff, "+/-", sigma_ext_coeff
        else:
            print "Not used ",  obj, " with filter " + filt + ". Too few elements."

########################################################################################################################
# After studying the values given by the different filters and fields of standards, we decide to use the value:
ext_coeff = 0.235
sigma_ext_coeff = 0.015

# From the Halpha images of Kopff 27:
airmass_kp27_Halpha = np.array([1.178, 1.735, 2.969])
magnitudes_kp27_Halpha = np.array([-10.133, -9.987, -9.696]) -\
                                              ext_coeff * airmass_kp27_Halpha
# Assuming random errors, the error of the magnitude estimation and the one 
# from the extinction correction should add quadratically.                                               
err_mag_kp27_Halpha = np.array([0.003, 0.003, 0.004])
err_extinction_kp27 = np.array([sigma_ext_coeff] * len(magnitudes_kp27_Halpha))*\
                  airmass_kp27_Halpha
err_kp27_Halpha = np.sqrt(err_mag_kp27_Halpha**2 + err_extinction_kp27**2)                  
#print "Magnitudes of Kopff 27 in Halpha:", magnitudes_kp27_Halpha
#print "Error of the magnitudes:", err_kp27_Halpha

# Now the fluxes of the standard star Kopff 27:
flux_kp27_Halpha = 1.05451e-11
m0_kp27_Halpha = - 2.5 * np.log10(flux_kp27_Halpha)

# From the Halpha images of bd+28
airmass_bd28_Halpha = np.array([1.955,1.039,1.020])
magnitudes_bd28_Halpha = np.array([-9.343,-9.547,-9.531]) -\
                                              ext_coeff * airmass_bd28_Halpha
err_mag_bd28_Halpha = np.array([0.002, 0.003, 0.002])                                              
err_extinction_bd28 = np.array([sigma_ext_coeff] * len(magnitudes_bd28_Halpha)) *\
                  airmass_bd28_Halpha 
err_bd28_Halpha = np.sqrt(err_mag_bd28_Halpha**2 + err_extinction_bd28**2)                  
#print "Magnitudes of bd28:", magnitudes_bd28_Halpha
#print "Error of the magnitudes:", err_bd28_Halpha

# Fluxes
flux_bd28_Halpha = 6.30523e-12
m0_bd28_Halpha = -2.5 * np.log10(flux_bd28_Halpha)

# Concatenate all together
obs_magnitudes = np.array(list(magnitudes_kp27_Halpha) + list(magnitudes_bd28_Halpha))
err_mag = np.array(list(err_kp27_Halpha) + list(err_bd28_Halpha))
std_magnitudes = np.array([m0_kp27_Halpha] * 3 + [m0_bd28_Halpha] * 3)
zp, sigma_zp = calculate_zeropoint(obs_magnitudes, std_magnitudes, err_mag)
#plt.plot(obs_magnitudes[0:3], std_magnitudes[0:3], 'o' )
#plt.plot(obs_magnitudes[3:], std_magnitudes[3:], 'o' )
#plt.show()
print "Zero point for Halpha:", zp, "+/-", sigma_zp
###################################################################################################################


# Equate PSFs for scientific objects
SciObj_set = set(x for x,t in zip(list_images["objname"], list_images["type"]) if t in ["cig", "clusters"])
for obj in SciObj_set:
    obj_images = list(list_images["filename"][np.where(list_images["objname"] == obj)])
    stars_images = [utils.replace_extension(ii, "radec") for ii in obj_images]
    output_images = match_psfs.main(["--input_stars"] + stars_images + ["--suffix", " -p",
                                                                        "--sigma_key", sky_stdk,
                                                                        "--gain_key", gaink,
                                                                        "--ron_key", read_noisek,
                                                                        "--expt_key", exptimek,
                                                                        "--airm_key", airmassk,
                                                                        "--FWHM_k", seeingk] + obj_images)
    list_images["filename"][np.where(list_images["objname"] == obj)] = output_images


print "Normalize using exposure time. "
for ii, im_name in enumerate(list_images["filename"]):
    if list_images["type"][ii] in ["cig", "clusters"]:
        tt = float(utils.get_from_header(im_name, exptimek))
        newname = utils.add_suffix_prefix(im_name, suffix="-t")
        mssg = "Normalized to exptime (divided by " + str(tt) + ")"
        arith_images.main(arguments=["--output", newname, "--message", mssg, "--mask_key", "MASK", im_name, "/", str(tt)])
        mssg = "Before normalizing: " + str(tt)
        # Update values for the exptime, the sky, the sky std...
        utils.header_update_keyword(newname, exptimek, 1, mssg)
        ss = float(utils.get_from_header(im_name, skyk))
        utils.header_update_keyword(newname, skyk, ss/tt)
        ss_std = float(utils.get_from_header(im_name, sky_stdk))
        utils.header_update_keyword(newname, skyk, ss_std/tt)
        list_images["filename"][ii] = newname


print "Correct for atmospheric extinction."
for ii, im_name in enumerate(list_images["filename"]):
    if list_images["type"][ii] in ["cig", "clusters"]:
        airmass = float(utils.get_from_header(im_name, airmassk))
        correcting_factor = 10**(ext_coeff * airmass / 2.5)
        newname = utils.add_suffix_prefix(im_name, suffix="-e")
        mssg = "Corrected from atmospheric extinction. Coefficient: " + str(ext_coeff)
        arith_images.main(arguments=["--output", newname,
                                     "--message", mssg,
                                     "--mask_key", "MASK",
                                     im_name, "*", str(correcting_factor)])
        mssg = "Before correcting for atmosphere: " + str(airmass)
        utils.header_update_keyword(im_name, airmassk, 0, comment=mssg)
        utils.header_update_keyword(im_name, "k_coeff", ext_coeff, comment="Extinction coefficient")
        utils.header_update_keyword(im_name, "k_err", sigma_ext_coeff, comment="Extinction coefficient 1-sigma")

# Add zero point to the header of all the Halpha images
for im_name in list_images["filename"]:
    filter = utils.get_from_header(im_name, filterk)
    if filter == "H6678":
        utils.header_update_keyword(im_name, "ZP", zp, "AB magnitude zero point." )
        utils.header_update_keyword(im_name, "ZP_err", sigma_zp, "Zero point 1-sigma. ")

print "Combine images of same object and filter"
iraf.images(_doprint=0)
iraf.immatch(_doprint=0)

keywords = dict(filterk=filterk, objectk=objectk, gaink=gaink, datek=datek, exptimek=exptimek,
                fwhmk="seeing", airmassk=airmassk)
objects_list = dict()
for obj in SciObj_set:
    obj_images = list(list_images["filename"][np.where(list_images["objname"] == obj)])
    obj_filters = [utils.get_from_header(im_name, filterk) for im_name in obj_images]
    for filt in set(obj_filters):
        input_images = [im_name for im_name, im_filt in zip(obj_images, obj_filters) if im_filt == filt]
        output_name = os.path.join(directory, obj + "_combined_" + filt + ".fits")
        input_names = ",".join(input_images)
        utilities.if_exists_remove(output_name)
        iraf.imcombine(input_names, output=output_name, combine="median", offsets="wcs")

        # Save all this in an astronomical object (see objects.py)
        if obj in objects_list.keys():
            current_object = objects_list[obj]
        else:
            current_object = astronomical_object(obj_name=obj, obj_type="galaxy", keywords=keywords)

        if filt[0:2].lower() == "ha":
            current_object.narrow_final = output_name
        elif filt[0].lower() == "r":
            current_object.cont_final = output_name
        objects_list[obj] = current_object

# 





sys.exit()






