#! /usr/bin/env python
# -*- coding: utf-8 -*-

import pyraf.iraf as iraf
import tempfile
import argparse
import sys
import pyraf.iraf as iraf
import repipy.astroim as astroim
import repipy.utilities as utils
import os

#fd, coords_file = tempfile.mkstemp(prefix=basename, suffix=".coords")

def apply_phot(args):
    iraf.noao(_doprint=0)
    iraf.digiphot(_doprint=0)
    basename = args.image.im_name
    hdr = args.image.header  # for short
    fd, phot_file = tempfile.mkstemp(prefix=basename, suffix=".mag.1")
    seeing, sigma = utils.get_from_header(args.image.im_name, hdr.seeingk, hdr.sigmak)
    iraf.phot(args.image.im_name, output=phot_file, coords=args.model_stars,
              wcsin=args.coords, fwhmpsf=seeing, sigma=sigma, ccdread=hdr.ccdronk,
              gain=hdr.gaink, exposure=hdr.exptimek,
              airmass=hdr.airmassk, annulus=6*seeing, dannulus=3*seeing,
              apert=3*seeing, verbose="no", verify="no", interac="no")
    print "Output phot file: ", phot_file
    return phot_file


def subtract_stars(args):
    """ Routine to subtract a set of stars defined by the user. We will model the PSF of some suitable stars in the
        images, then use that """

    # First substitute the names of the image in args.image by the corresponding astroim object, with more information
    args.image = astroim.Astroim(args.image)

    # Do photometry on the stars to be used to model the PSF
    phot_file = apply_phot(args)

# Create parser
parser = argparse.ArgumentParser(description='Extract stars from an image. ')

# Add necessary arguments to parser
parser.add_argument("image", metavar='image', action='store', help='Name ' +\
                    'of fits image to calculate extract stars from.')
parser.add_argument("--model_stars", metavar='model_stars', action='store',
                    dest="model_stars",
                    help='File containing a list of stars in two columns. ' +\
                    'Parameter "coords" will determine if the coordinates are pixels  '+\
                    'or RA DEC. ', nargs="+")
parser.add_argument("--coords", metavar='coords', action='store', default="world",
                    help='Type of coordinates that "stars" give. Allowed: '+\
                    ' "logical", tv", "physical", and "world". Default: world ')
parser.add_argument("--subt_stars", metavar='subt_stars', action='store', dest='subt_stars',
                    help='File containing the list of stars to be subtracted. ', nargs="+")
parser.add_argument("--output", metavar='output', action='store', dest='output',
                    help='Name of output file. If no output is indicated, the suffix "-e" will be added ' +\
                         'to the image name. ', nargs="+")


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)
  print "Model stars file", args.model_stars , type(args.model_stars)
  args.model_stars = os.path.abspath(args.model_stars)
  args.subt_stars = os.path.abspath(args.subt_stars)


  subtract_stars(args)
  return None

if __name__ == "__main__":
    main()


