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
import astropy.io.fits as fits
import tempfile
import os
import shutil
import numpy

def combine(input_images, output, scale, combine_method, offsets="wcs"):
    """ Combine images with IRAF's imcombine, using the wcs in their headers to match the pixels. """

    utilities.if_exists_remove(output)
    input_names = ",".join(input_images)
    # Unfortunately the mean is considered the only average in IRAF... (sic!)
    if combine_method == "mean":
        combine_method = "average"

    iraf.imcombine(input_names, output=output, scale=scale.lower(), combine=combine_method.lower(),
                   offsets=offsets)
    return None

def include_wcs_in_masks(input_images):
    """  Put the World Coordinate System of the images inside the masks
    """
    img_list  = [astroim.Astroim(im_name, memmap=True) for im_name in input_images]
    mask_names = [im.primary_header.get("MASK") for im in img_list]
    output = []
    for im_object, mask_name in zip(img_list, mask_names):
        with fits.open(mask_name, 'readonly') as mask:
            mask_header = im_object.chips[0].header.hdr
            mask_data = mask[0].data.copy()
            mask_data[mask_data>0] = 1
            _, path = tempfile.mkstemp(suffix=".fits")
            fits.writeto(path, mask_data * 1., mask_header, clobber=True)
            output.append(path)
    return output



def final_mask(path, output_mask, percentage=0.5):
    """ Move temporary mask to final output mask file, mask any pixel with value > percentage, unmask the rest .

    The temporary mask is an average of masks, where a masked pixel has value 1. So a value of 0.3 means that 30% of
    the input pixels where masked. The level at which a pixel is masked, because it is the combination of many masked
    pixels is established by percentage, and a value of 0.5 is currently in place.
    """
    with fits.open(path, "readonly") as temp_mask:
        mask_data = temp_mask[0].data
        mask_header = temp_mask[0].header
        mask_data[mask_data >= percentage] = 1
        mask_data[mask_data < percentage] = 0
        fits.writeto(output_mask, mask_data, mask_header, clobber=True)


def combine_images(args):
    """ Combine images by matching their World Coordinate Systems,

    Combine the set of images using imcombine. Add the masks, if present, using the same WCS as the images, mask out any
    pixel in the combined image for which more than a certain % of images where masked.
    """

    masks_with_wcs = include_wcs_in_masks(args.input)
    output_mask = args.output + ".msk"
    combine(args.input, args.output, args.scale, args.average, offsets="wcs")
    _, path = tempfile.mkstemp(suffix=".fits")
    utilities.if_exists_remove(path)
    combine(masks_with_wcs, path, "none", "average", offsets="wcs")
    final_mask(path, output_mask)
    utilities.header_update_keyword(args.output, "MASK", os.path.abspath(output_mask))
    return args.output













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
                   help='type of average (median, mean, sum) to combine ' +\
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
  newfile = combine_images(args)
  return newfile  

if __name__ == "__main__":
    main()



# TODO: Modify code to accept multi-extension fits files.