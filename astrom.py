#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Created on Tue Jul  9 11:03:50 2013

@author: javier blasco herrera
"""

import repipy.utilities as utilities
import repipy.astroim as astroim
import repipy.remove_cosmics as remove_cosmics
import subprocess
import os
import sys
import shutil
import tempfile
import argparse
import astropy.io.fits as fits
from repipy import __path__ as repipy_path
repipy_path = repipy_path[0]


def get_basename(im_name):
    """  Given a file name, return the file without the path AND without the extension
    :param im_name: input filename
    :return: the same file, but without the path and the extension, for /home/blasco/file,fits --> file
    """
    basename = os.path.splitext(os.path.basename(im_name))[0]
    return basename


def perform_cosmic_removal(im_name, output=None):
    """ Remove cosmics from fits image
    :param im_name: fits image to be cleaned from cosmic rays
    :return: name of the cleaned image
    """
    if output is None:
        output = im_name

    arguments = [im_name, "--output", output, "--maxiter", "2"]

    im = astroim.Astroim(im_name)

    needed_keys = [("--gain", im.header.gaink), ("--readnoise", im.header.ccdronk)]
    for arg, kk in needed_keys:
        if kk is None:
            continue
        arguments += [arg, kk]

    new_name = remove_cosmics.main(arguments=arguments)

    return new_name


def build_default_sex(args):
    """ Build a default.sex file in the path of repipy
    :return:
    """
    with open(os.path.join(repipy_path, "default.sex"), 'w') as fd:
        fd.write("PARAMETERS_NAME   {0}  \n".format(os.path.join(repipy_path, "default.param")))
        fd.write("FILTER_NAME       {0}  \n".format(os.path.join(repipy_path, "default.conv")))
        fd.write("STARNNW_NAME      {0}  \n".format(os.path.join(repipy_path, "default.nnw")))
        fd.write("DETECT_MINAREA         10  \nTHRESH_TYPE         RELATIVE \n"
                 "DETECT_THRESH        {0}  \nANALYSIS_THRESH     {0}  \n"
                 "CATALOG_TYPE        FITS_1.0".format(args.threshold))
    return None




def include_wcs(args):
    # Copy input names into output names
    output_names = args.input[:]
    for ii, im_name in enumerate(args.input):
        # If it is a domeflat, skyflat or bias image, calculating the astrometry makes little sense
        im = astroim.Astroim(im_name)
        if im.target.objtype in ["bias", "domeflat", "skyflat", "flat"]:
            continue

        # Prepare my output file  for the resulting image
        if args.suffix:
            new_file = utilities.add_suffix_prefix(im_name, suffix=args.suffix)
        elif args.overwrite:
            new_file = im_name

        # Copy input image into a temporary file, so that we can modify it freely, for example, remove cosmics
        # or filter it (to be done).
        basename = get_basename(im_name)
        _, input_image = tempfile.mkstemp(prefix=basename, suffix=".fits")
        shutil.copy2(im_name, input_image)

        # Remove cosmics
        if args.cosmics:
            perform_cosmic_removal(input_image)

        # The output of the WCS process of astrometry.net will go to a .new file, the coordinates to a .coor
        output_wcs = utilities.replace_extension(input_image, ".new")
        solved_file = utilities.replace_extension(input_image, ".solved")
        corrfile = utilities.replace_extension(input_image, ".cor")

        # Try first with the defaults of astrometry
        arguments_def = ["solve-field", "--no-plots", "--no-fits2fits", "--dir", "/tmp", "--overwrite",
                      "--new-fits", output_wcs, "--corr", corrfile, "--cpulimit", "1", input_image]
        try:  # Try to add the RA, DEC, Radius options to constrain the search
            ra, dec = im.header.get(im.header.RAk, im.header.DECk)
            ra, dec = utilities.sex2deg(ra, dec)
            arguments = arguments_def + ["--ra", str(ra), "--dec", str(dec), "--radius", str(args.radius),
                                      "--cpulimit", "20"]
        except:
            arguments = arguments_def
        print "Trying to find WCS with astrometry standard keywords. "
        subprocess.call(arguments)

        # Now we will try using sextractor
        build_default_sex(args)
        # To avoid having too much residual crap in the folder, the output of astrometry will go to tmp (--dir /tmp).
        arguments0 = ["solve-field", "--no-plots", "--no-fits2fits", "--use-sextractor", "--dir", "/tmp",
                      "--x-column", "X_IMAGE", "--y-column", "Y_IMAGE", "--sort-column", "MAG_AUTO",
                      "--sort-ascending", "--sextractor-config", os.path.join(repipy_path, "default.sex"),
                      "--overwrite", "--new-fits", output_wcs, "--corr", corrfile, input_image]
        arguments0 += args.extras

        try:  # Try to add the RA, DEC, Radius options to constrain the search
            ra, dec = im.header.get(im.header.RAk, im.header.DECk)
            ra, dec = utilities.sex2deg(ra, dec)
            arguments = arguments0 + ["--ra", str(ra), "--dec", str(dec), "--radius", str(args.radius),
                                      "--cpulimit", "20"]
        except:
            arguments = arguments0

        # Run astrometry, in case of not solving it on the first attempt, try fitting freely (no RA, DEC used)
        if not os.path.exists(solved_file):
            subprocess.call(arguments)

        if not os.path.exists(solved_file):
            subprocess.call(arguments0)


        # Only if we have a solution
        if os.path.exists(solved_file):
            # copy the input file into the new file if they are not the same
            if os.path.abspath(im_name) != os.path.abspath(new_file):
                shutil.copy2(im_name, new_file)

            # get WCS from the resulting file, only the part added by astrometry.net
            with fits.open(output_wcs) as file_wcs:
                hdr_wcs = file_wcs[0].header
                ind = hdr_wcs.values().index('Original key: "END"')  # new part added by astrometry
                hdr_wcs = hdr_wcs[ind:]

            # Add that header to the original one of the image
            with fits.open(new_file, 'update') as new_im:
                new_im[0].header += hdr_wcs
                new_im.flush()

            # Build a catalogue of RA, DEC for the stars found
            with fits.open(corrfile) as table:
                cat_radec = utilities.replace_extension(new_file, "radec")
                with open(cat_radec, "w") as f:
                    for ra, dec in zip(table[1].data["field_ra"], table[1].data["field_dec"]):
                        f.write(str(ra) + " " + str(dec) + "\n")
            output_names[ii] = new_file

            # If --remove_original is True
            if args.remove_orig and not args.overwrite:
                os.unlink(im_name)

    return output_names


# Create parser
parser = argparse.ArgumentParser(description=""" Find the WCS for a set of images, not necessarily of the same object.
                                          """)


# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of '
                    'images from which to compute the WCS.', nargs="+",
                    type=str)
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',
                    default='', type=str, help='suffix to be added at the end ' + \
                                               'of the image input list to generate the outputs. There ' + \
                                               'is a peculiarity with argparse: if you pass, e.g., "-c" to ' + \
                                               '--suffix, the program will understand that you want to ' + \
                                               'call the code with the flag -c, which does not exist. This ' + \
                                               'does not raise an error, but just stops execution, which is ' + \
                                               'quite annoying. One way around it is " -c" (notice the ' + \
                                               'space, since within the string is stripped. If the --overwrite is used, ' + \
                                               'the suffix will be ignored. Use wisely!')
parser.add_argument("--remove_original", dest='remove_orig', action="store_true", \
                    help=' Once the job is done, remove the original image. This will have no effect if ' + \
                         ' --overwrite is used.')
parser.add_argument("--cosmics", dest='cosmics', action="store_true", \
                    help=' Prevent some cosmic rays from being detected by astrometry.net. What happens under the hood' + \
                         'is that we copy the input image into a new one, subtract cosmic rays, solve with '+\
                         'astrometry.net and add the WCS solution to the original image, which remains otherwise '+\
                         'unaltered. ')
parser.add_argument("--overwrite", dest='overwrite', action="store_true",
                    help='Overwrite the original images. This keyword will override both --remove_original and --suffix')
parser.add_argument("--radius", dest='radius', action="store", type=float, default="1",
                    help='Search radius. If the RA and DEC are found in the header, astrometry will look for ' + \
                         'solutions within this radius of those coordinates. Default=1.0')
parser.add_argument("--threshold", dest='threshold', action="store", type=float, default=7,
                    help='Detection threshold for sextractor Search radius. If the RA and DEC are found in the header, astrometry will look for ' + \
                         'solutions within this radius of those coordinates. Default=1.0')
parser.add_argument('-o', action = 'append', dest = 'extras',
                  default = [], help = "additional options to pass to Astrometry.net's "
                  "solve-field. For example: -o ' --downsample 2' -o ' --tweak-order 2'. "
                  "This option may be used multiple times.")

def main(arguments=None):
    # Pass arguments to variable args
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)


    # Either --suffix is not "" or --overwrite is active
    if not args.suffix and not args.overwrite:
        sys.exit("Declare a --suffix or --overwrite, one of the two must be active!")


    # If overwrite is used, --suffix and --remove_original have no effect
    if args.overwrite:
        args.suffix = ""
        args.remove_original = False

    # If a suffix was passed, check that you remove any blank spaces
    args.suffix = args.suffix.strip()

    # If -o is used, make sure you split the command
    args.extras = sum( [ii.strip().split() for ii in args.extras] ,[])

    output_list = include_wcs(args)
    return output_list

if __name__ == "__main__":
    main()

