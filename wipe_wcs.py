#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import repipy.utilities as utilities
import astropy.io.fits as fits
import sys


def wipe(args):
    for im_name in args.input:
        with fits.open(im_name, 'update') as im:
            for hdu in im:
                hdu.header = utilities.remove_WCS(hdu.header)
                im.flush()



# Create parser
parser = argparse.ArgumentParser(description='(Almost) Completely wipe the World Coordinate Sytem'+\
                                             ' from the header of a list of images')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'images from which to remove the WCS.', nargs="+",
                    type=str)


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)
  wipe(args)

if __name__ == "__main__":
    main()

