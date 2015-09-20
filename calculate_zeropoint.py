#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:22:00 2013

@author: blasco
"""
import repipy.astroim as astroim
import repipy
import numpy as np
import argparse
import sys

def calculate_zeropoint(args):
    """ Calculate the zero point for a series of images

    The images must be reduced, and contain a standard star included in the standards.csv file within your path to
    repipy.
    """
    zero_points, filters = [], []
    for im_name in args.input:
        im = astroim.Astroim(im_name)
        zero_points.append(im.zero_point())
        filters.append(im.filter.filter_name)

    if len(set(filters)) != 1:
        msg = "\n \n ERROR: The images include several filters, calculating the zero point of a mix of filters " \
              "is not allowed. Select images with the same filter! \n \n"
        sys.exit(msg)

    zp, std_zp = np.mean(zero_points), np.std(zero_points)
    print "\n Zero point for filter {0}: \n zp = {1}, std_zp = {2} \n".format(filters[0], zp, std_zp )


    return zp, std_zp

# Create parser
parser = argparse.ArgumentParser(description='\n Program to calculate the zero point of a set of images. The images '
                                             'must be reduced, and contain a standard star from the '
                                             'file: {0}/standards.csv \n'.format(repipy.__path__[0]))

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'input images for which to estimate the FWHM.',
                    nargs="+", type=str)

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)

  zp, zp_std = calculate_zeropoint(args)
  return zp, zp_std

if __name__ == "__main__":
    main()


