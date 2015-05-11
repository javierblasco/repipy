#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Created on Tue Jul  9 11:03:50 2013

@author: javier blasco herrera
"""

import repipy.utilities as utilities
import repipy.astroim as astroim
import subprocess
import os
import sys
import shutil
import tempfile
import argparse
import astropy.io.fits as fits

def include_WCS(args):

    # Copy input names into output names
    output_names = args.input[:]
    for ii, im_name in enumerate(args.input):
        # If it is a domeflat, skyflat or bias image, calculating the astrometry makes little sense
        im = astroim.Astroim(im_name)
        if im.target.objtype in ["bias", "domeflat", "skyflat"]:
            continue

        # Prepare my output file  for the resulting image
        new_file = im_name
        if args.suffix:
            new_file = utilities.add_suffix_prefix(new_file, suffix=args.suffix)



        # And now, to avoid having too much residual crap in the folder, the output of astrometry will go to tmp
        basename = os.path.basename(im_name)
        fd1, resultfile = tempfile.mkstemp(prefix=basename, suffix=".newfile")
        fd2, wcsfile = tempfile.mkstemp(prefix=basename, suffix=".wcs")
        fd3, xylsfile = tempfile.mkstemp(prefix=basename, suffix=".xyls")
        fd4, rdlsfile = tempfile.mkstemp(prefix=basename, suffix=".rdls")
        fd5, axyfile = tempfile.mkstemp(prefix=basename, suffix=".axy")
        fd6, solutionfile = tempfile.mkstemp(prefix=basename, suffix=".solved")
        fd7, matchfile = tempfile.mkstemp(prefix=basename, suffix=".match")
        fd8, corrfile = tempfile.mkstemp(prefix=basename, suffix=".corr")


        # Prepare the arguments of solve-field
        arguments = ["solve-field", "--no-plots", "--no-fits2fits", "--use-sextractor",
                     "--overwrite", "--new-fits", resultfile, "--wcs", wcsfile, "--index-xyls", xylsfile,
                     "--rdls", rdlsfile, "--axy", axyfile, "--solved", solutionfile, "--match", matchfile,
                     "--cpulimit", "120", "--cpulimit", "120", "--corr", corrfile, im_name]

        # Try to add the RA, DEC, Radius options to constrain the search
        with utilities.tmp_mute():
            try:
                ra, dec = im.header.get( im.header.RAk, im.header.DECk)
                ra, dec = utilities.sex2deg(ra, dec)
                radius = 1
                arguments = arguments + ["--ra", str(ra), "--dec", str(dec), "--radius", str(args.radius)]
            except:
                pass

        # Run astrometry, in case of not solving it on the first attempt, try fitting freely (no RA, DEC used)
        subprocess.call(arguments)
        if not os.path.exists(solutionfile):
            subprocess.call(arguments[:-6])

        # Only if we have a solution, move resulting files around
        if os.path.exists(solutionfile):
            # Move important files around
            shutil.move(resultfile, new_file)
            with fits.open(rdlsfile) as table:
                cat_radec = utilities.replace_extension(im_name, "radec")
                with open(cat_radec, "w") as f:
                    for line in table[1].data:
                        f.write(str(line[0]) + " " + str(line[1]) + "\n")
            output_names[ii] = new_file

        # If --remove_original is True
        if args.remove_orig:
            os.unlink(im_name)

    return output_names


# Create parser
parser = argparse.ArgumentParser(description=""" Find the WCS for a set of images, not necessarily of the same object.
                                          """)


# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'images from which to compute the WCS.', nargs="+", type=str)
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',\
                    default='', type=str, help='suffix to be added at the end '+\
                    'of the image input list to generate the outputs. There '+\
                    'is a peculiarity with argparse: if you pass, e.g., "-c" to '+\
                    '--suffix, the program will understand that you want to '+\
                    'call the code with the flag -c, which does not exist. This '+\
                    'does not raise an error, but just stops execution, which is '+\
                    'quite annoying. One way around it is " -c" (notice the '+\
                    'space, since within the string is stripped. If the --overwrite is used, '+\
                    'the suffix will be ignored. Use wisely!')
parser.add_argument("--remove_original",  dest='remove_orig', action="store_true",\
                     help=' Once the job is done, remove the original image. This will have no effect if '+\
                          ' --overwrite is used.')
parser.add_argument("--overwrite", dest='overwrite', action="store_true",
                     help='Overwrite the original images. This keyword will override both --remove_original and --suffix')
parser.add_argument("--radius", dest='radius', action="store", type=float,
                     help='Search radius. If the RA and DEC are found in the header, astrometry will look for '+\
                     'solutions within this radius of those coordinates. Default=1.0')


def main(arguments = None):
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

  output_list = include_WCS(args)
  return output_list

if __name__ == "__main__":
    main()

