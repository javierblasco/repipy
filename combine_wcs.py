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
import scipy.stats.mstats
from reproject import reproject_interp


# def combine(input_images, output, scale, combine_method, offsets="wcs"):
#     """ Combine images with IRAF's imcombine, using the wcs in their headers to match the pixels. """
#
#     utilities.if_exists_remove(output)
#     input_names = ",".join(input_images)
#     # Unfortunately the mean is considered the only average in IRAF... (sic!)
#     if combine_method == "mean":
#         combine_method = "average"
#
#     iraf.imcombine(input_names, output=output, scale=scale.lower(), combine=combine_method.lower(),
#                    offsets=offsets)
#     return None

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
    """ Combine images using their World Coordinate Systems to match them,
    """

    # Read all images into a cube (TODO: think about the RAM)
    with fits.open(args.input[0]) as im0:
        lx, ly = im0[0].data.shape
        ref_hdr = im0[0].header

    headers = [fits.open(im_name)[0].header for im_name in args.input]
    cube = numpy.ma.zeros((len(args.input), lx, ly))
    cube.mask = numpy.zeros_like(cube.data)
    for ii, im_name in enumerate(args.input):
        with astroim.Astroim(im_name) as im:
            cube.data[ii, :,:] = im.chips[0].data
            if im.chips[0].mask is not None:
                cube.mask[ii,:,:] = im.chips[0].mask

    # Scale images
    scale_functions = {"median": numpy.ma.median,
                       "mean": numpy.ma.mean,
                       "mode": scipy.stats.mstats.mode,
                       "none": lambda x: 1}
    for ii, im_name in enumerate(args.input):
        func = scale_functions[args.scale.lower()]
        cube[ii,:,:] /= func(cube[ii,:,:])


    # Reproject all images to the ref_hdr
    for ii, _ in enumerate(args.input):
        if ii == 0:
            continue
        cube.data[ii,:,:], footprint = reproject_interp((cube.data[ii,:,:], headers[ii]), ref_hdr)
        cube.mask[ii,:,:], footprint = reproject_interp((cube.mask[ii,:,:], headers[ii]), ref_hdr)
        #whr = numpy.isnan(cube.data[ii,:,:])
        #cube.mask[ii,:,:][whr] = True

    # Do average
    average_functions = {"median": numpy.ma.median, "mean": numpy.ma.mean, "sum": numpy.ma.sum}
    func = average_functions[args.average.lower()]
    final_image = func(cube, axis=0)
    ref_hdr["NCOMBINE"] = len(args.input)

    mask_name = utilities.replace_extension(args.output, ".fits.msk")
    mask_name_header = utilities.replace_extension(os.path.basename(args.output), ".fits.msk")
    ref_hdr["MASK"] = mask_name_header
    fits.writeto(args.output, final_image.data, ref_hdr, clobber=True )
    fits.writeto(mask_name, numpy.array(final_image.mask, dtype=int), clobber=True)

    return args.output













############################################################################

# Create parser
parser = argparse.ArgumentParser(description='Combine images')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', \
                 help='input images to combine.', nargs='+', type=str)
parser.add_argument("--reference", metavar='reference', action='store', type=str, help='Reference image. The other  '
                    'images will be oriented like the reference one before combining them. ')
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
parser.add_argument("--memmap", action="store_true", dest="memmap", help='Use when many images, or when the images are '
                    'very large, to store the images in the disk instead of in RAM memmory. This will obviously slow '
                    "down the processing of the images, but it's one of the few options for big data." )

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