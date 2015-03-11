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
import repipy.astroim as astroim
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.optimize as optimize
import sys
import astropy.io.fits as fits
import os
import shutil
import pickle
import subprocess
from lemon import methods
import repipy

# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import images


from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

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
coefficient_list = []
for obj in standards_set:
    input_db = os.path.join(directory, obj+".db")
    airmasses, magnitudes, filters = extract.main(input_db)
    for filt in set(filters.flatten()):
        # find columns with the current filter
        columns = [ii for ii in range(len(filters[0, :])) if filters[0,ii] == filt]
        # if there are more than two images
        if len(airmasses[0,columns]) > 0 and len(airmasses[0,columns]) > 2:
            ext_coeff, serror_ext_coeff = calculate_extinction(airmasses[:, columns], magnitudes[:, columns])
            coefficient_list.append(ext_coeff)
            #plt.plot(airmasses[:, columns], magnitudes[:, columns], 'o')
            #plt.show()
            print "Extinction coefficient for ", obj, " for filter ", filt, ":", ext_coeff, "+/-", serror_ext_coeff
        else:
            print "Not used ",  obj, " with filter " + filt + ". Too few elements."

ext_coeff, serror_ext_coeff =  np.mean(coefficient_list), np.std(coefficient_list) / np.sqrt(len(coefficient_list) - 1)
print "Mean extinction = ", ext_coeff, "+/-", serror_ext_coeff

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
        utils.header_update_keyword(newname, sky_stdk, ss_std/tt)
        list_images["filename"][ii] = newname


print "Correct for atmospheric extinction."
for ii, im_name in enumerate(list_images["filename"]):
    if list_images["type"][ii] in ["cig", "clusters", "standards"]:
        airmass = float(utils.get_from_header(im_name, airmassk))
        correcting_factor = 10**(ext_coeff * airmass / 2.5)
        newname = utils.add_suffix_prefix(im_name, suffix="-e")
        mssg = "Corrected from atmospheric extinction. Coefficient: " + str(ext_coeff)
        arith_images.main(arguments=["--output", newname,
                                     "--message", mssg,
                                     "--mask_key", "MASK",
                                     im_name, "*", str(correcting_factor)])
        mssg = "Before correcting for atmosphere: " + str(airmass)
        utils.header_update_keyword(newname, airmassk, 0, comment=mssg)
        utils.header_update_keyword(newname, "k_coeff", ext_coeff, comment="Atmospheric extinction coefficient")
        utils.header_update_keyword(newname, "k_err", serror_ext_coeff, comment="Standard error extinction coefficient")
        list_images["filename"][ii] = newname
        sky, sky_std = utils.get_from_header(im_name, "sky", "sky_std")
        utils.header_update_keyword(newname, "sky", sky * correcting_factor)
        utils.header_update_keyword(newname, "sky_std", sky_std * correcting_factor)




print "Calculate zero point from the standards"
zp = defaultdict(list)
for ii, im_name in enumerate(list_images["filename"]):
    if list_images["type"][ii] in ['standards']:
        image = astroim.Astroim(im_name)
        zp[image.filter.filter_ID].append(image.zero_point)

for kk, vv in zp.iteritems():
    zp[kk] = np.median(vv), np.median( np.abs(np.array(vv)-np.median(vv)))

print "Add zero point to headers"
for im_name in list_images["filename"]:
    filter = astroim.Astroim(im_name).filter.filter_ID
    if zp.has_key(filter):
        utils.header_update_keyword(im_name, "ZP", zp[filter][0], "AB magnitude zero point." )
        utils.header_update_keyword(im_name, "ZP_err", zp[filter][1], "Zero point 1-sigma. ")


print "Combine images of same object and filter"

# Dictionary that contains the images separated by object and filter
objects_list = defaultdict(list)

for target in set(list_images['objname']):
    obj_images = list(list_images["filename"][list_images["objname"] == target])
    obj_filters = [utils.get_from_header(im_name, filterk) for im_name in obj_images]
    for filt in set(obj_filters):
        input_images = [im_name for im_name, im_filt in zip(obj_images, obj_filters) if im_filt == filt]
        output_name = os.path.join(directory, target + "_combined_" + filt + ".fits")
        input_names = ",".join(input_images)
        utilities.if_exists_remove(output_name)
        iraf.imcombine(input_names, output=output_name, combine="median", offsets="wcs")
        objects_list[target].append(output_name)


print "Do photometry in combined images"
for target in SciObj_set:
    target_images = objects_list[target]
    print target_images
    output_db = os.path.join(directory, target+".db")
    utilities.if_exists_remove(output_db)
    photometry.main(arguments=["--maximum", str(max_counts), "--uik", "", "--margin", "20", "--gaink", gaink,
                              "--aperture", "4.", "--annulus", "6", "--dannulus", "2", "--individual-fwhm", "--objectk", objectk,
                              "--filterk", filterk, "--datek", datek, "--expk", exptimek, "--fwhmk", seeingk, "--airmk", airmassk,
                               target_images[0]] + target_images + [output_db])


print "Calculate scaling factors"

scale_factor = {}
scale_MAD = {}
for target in SciObj_set:
    scaling_factors = []
    if len(objects_list[target]) == 2: # both continuum and Halpha present
        input_db = os.path.join(directory, target+".db")
        airmasses, magnitudes, filters = extract.main(input_db)
        # Separate Halpha and rGunn magnitudes.
        Ha_index, rGunn_index = 'H-alpha' in filters[0,1], 'Gunn' in filters[0,1]
        magnitudes_Halpha, magnitudes_rGunn = magnitudes[:,Ha_index], magnitudes[:,rGunn_index]
        sf =  10 ** ( (np.array(magnitudes_Halpha) - np.array(magnitudes_rGunn)) / (-2.5))
        scaling_factors = np.append(scaling_factors, sf)

        scale_factor[target] = np.median(scaling_factors)
        scale_MAD[target] = np.median( np.abs( scaling_factors - np.median(scaling_factors)) )
        print "Scale factor for object", target, scale_factor[target]

print "Scale continuum images using stars"
for target in SciObj_set:
    if target in scale_factor.keys():
        im_cont = [im for im in objects_list[target] if 'rgunn' in im.lower()][0]
        newname = utils.add_suffix_prefix(im_cont, suffix="-scaled")
        mssg = 'Scaled to Halpha image'
        arith_images.main(arguments=["--output", newname, "--message", mssg, im_cont, "*", str(scale_factor[target])])
        sky, std_sky = utils.get_from_header(newname, "sky", "sky_std")
        utils.header_update_keyword(newname, "sky", float(sky) * scale_factor[target] )
        utils.header_update_keyword(newname, "sky_std", float(std_sky) * scale_factor[target])









