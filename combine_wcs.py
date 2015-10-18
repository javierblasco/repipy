#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""  Combine a list of images into a single output image.

The program takes in a list of images and returns an image with the type of combination indicated by --average,
produciendo an output file given by --output output_name. Scaling the images previous to the combination is optional.
If a mask is present in the header, the masks will be combined as well, and if more than a certain percentage of the
images contains a pixel that is masked out, the resulting pixel will be also masked out.
"""
import argparse
import sys
from lemon import methods
import repipy
from repipy import utilities
from repipy import astroim
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import images


def combine(args):
    """ Combine images using the wcs in their headers to match the pixels.
    :param args:
    :return:
    """
    utilities.if_exists_remove(args.output)
    input_names = ",".join(args.input)

    # Unfortunately the mean is considered the only average in IRAF... (sic!)
    if args.average == "mean":
        args.average = "average"

    iraf.imcombine(input_names, output=args.output, scale=args.scale.lower(), combine=args.average.lower(), offsets="wcs")










############################################################################

# Create parser
parser = argparse.ArgumentParser(description='Combine images')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', \
                 help='input images to combine.', nargs='+', type=str)
parser.add_argument("--output", metavar="output", dest='output', required=True, \
                   action='store', help='Name for output file. If not present, a name will be created by '
                                        'using the pattern: ')
parser.add_argument("--average", metavar='average', type=str, default='median', \
                   help='type of average (median, mean) to combine ' +\
                   'the images. The average will be made with the central third of the image to avoid border effects. '
                   'Default: median')
parser.add_argument("--scale", metavar='scale', type=str, default='none', \
                   help='scaling function (median, mean, mode, none) to apply' +\
                   'to	the images before combining them. As with the average, the central third of the image '
                   'will be used to avoid border effects. Default: none' )
parser.add_argument("--mask_key", metavar="mask_key", dest='mask_key', \
                    action='store', default="MASK", help=' Keyword in the header ' +\
                    'of the image that contains the name of the mask. The mask '+\
                    'will contain ones (1) in those pixels to be MASKED OUT.')
parser.add_argument("--output_mask", metavar="output_mask", dest='out_mask', \
                    action='store', default="", help=' Name of the output mask. '
                    'If none is provided, the program will add .msk '+\
                    'to the output name ( output.fits is the name of the output image. ')
parser.add_argument("--fill_val", metavar="fill_val", dest="fill_val", \
                    action='store', default='', help=' If present, this '+\
                    'keyword contains a value that substitutes the result '+\
                    'in those pixels that are masked in all images. For example,  '+\
                    'you might want to zero all values that had no valid values. ')

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  args = parser.parse_args(arguments)

	
  # Call combine, keep name of the file created
  newfile = combine(args)
  return newfile  

if __name__ == "__main__":
    main()
