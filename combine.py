#!/usr/bin/env python
# -*- coding: UTF-8 -*-

###############################################################################
import os
import argparse
import numpy
from scipy.stats import mode as mode
import sys
import astropy.io.fits as fits
import collections
import repipy.utilities as utils
import repipy.arith as arith_images
import repipy.find_keywords as find_keywords

def home_made_median(map_cube, cube):
    """ map_cube indicates the number of elements of cube that are not 
        masked. For the sort that we did before, masked elements are sorted 
        as if they had the very highest values. So, if map[ii,jj] = 5
        we know that the first 5 elements of cube[:,ii,jj] are not masked
        and that the rest are masked. So, the median will be cube[2,ii,jj]
        for the pixel [ii,jj] of the images. Of course, there will be points
        of map that will be even, so the median will be the mean of elements
        N//2 and (N-1)//2, where N is the value of map for the given pixel.
        That is what we call indicesz_above and below. We will average
        the two resulting images. """

    lx, ly = map_cube.shape 
    indicesz_above = (map_cube // 2).flatten() 
    indicesz_below = ((map_cube-1) // 2).flatten()
    indicesx, indicesy = numpy.where(map_cube == map_cube)
    image_above = cube[indicesz_above, indicesx, indicesy].reshape((lx,ly))
    image_below = cube[indicesz_below, indicesx, indicesy].reshape((lx,ly))
    image = (image_above + image_below) / 2.
    return image

################################################################################
def cube_images(input_images, mask_key, scales, limits=0):
    """ From a set of images (input_images), and once divided by their scales, 
        form a cube with all the images. """
    # Set the limits to build the piece of the cube
    header = fits.getheader(input_images[0])
    if limits == 0:
        min_x, min_y = (0, 0)
        max_x, max_y = header["NAXIS2"], header["NAXIS1"]
    else:
        min_x, min_y, max_x, max_y = limits

    # Initialize the cube and the mask.
    lx, ly = max_x - min_x, max_y - min_y
    lz = len(input_images)
    mask = numpy.zeros([lz,lx,ly])  # initialized to zero, all pixels are valid
    cube = numpy.ma.array(mask, mask=mask) 

    # Now read all the images (and masks, if present), into the cube
    for index, image in enumerate(input_images):
        im = utils.read_image_with_mask(image, mask_keyword = mask_key, limits=[min_x, min_y, max_x, max_y])
        cube.data[index,:,:] = im.data/scales[index]
        cube.mask[index,:,:] = im.mask
    return cube    
    
################################################################################
def compute_scales(input_images, scale_type, mask_key):
    """ From the list of images, use the central third of the image to
        calculate the scaling factor requested by user. """
    lx, ly = fits.getdata(input_images[0]).shape
    centre = [lx/3, ly/3, lx*2/3, ly*2/3]  # min_x, min_y, max_x, max_y
    scales = []
    for image in input_images:
        im = utils.read_image_with_mask(image, mask_keyword=mask_key, limits=centre)
        if scale_type == "median":
            new_item = numpy.ma.median(im)
        elif scale_type == "mean":
            new_item = numpy.ma.mean(im)
        elif scale_type == "none":
            new_item = 1
        if new_item is numpy.ma.masked:  # Result if all elements where masked
            new_item = 1
            print "All elements of the central part of image %s seem to be masked. "+\
                  "Replacing scale calculation by 1." % image
        scales.append(new_item)
    return scales  # cases where of all elements were masked

###################################################################################
def combine(args):
    # Create the folders that do not already exist for the output file
    outdir, outfile = os.path.split(os.path.abspath(args.output))
    if outdir == "":
        outdir = "."
    utils.if_dir_not_exists_create(outdir)

    # Build a list of the filter of each image
    images_filters = utils.collect_from_images(args.input, args.filterk)

    # If user wants all images to be combined together, regardless of filter:
    if args.all_together:
        images_filters = ["AllFilters"] * len(images_filters)

    # Create a default dictionary for the resulting images
    result = collections.defaultdict(str)    
    
    # For each of the filters present combine the images (and write)
    for filt in set(images_filters):
        # list of objects with current filter (exception: allfilters is true)
        list1 = [args.input[p] for p,f in enumerate(images_filters) if f == filt ]

        # Check that all images have same dimension. If not, exit program
        if not utils.check_dimensions(list1):
            sys.exit("Dimensions of images to combine are different!")

        # Calculate scale of images
        scales = compute_scales(list1, args.scale, args.mask_key)

        # Get the sizes of the images
        lx, ly = utils.get_from_header(list1[0], "NAXIS2", "NAXIS1")

        # Now, in order to avoid loading many images in memory, we need to slice the images in pieces and combine a slice
        # of all the images at a time
        n_slices = 32          # divide the slow axis in this many pieces

        # Define the whole image and set all elements of mask to False
        whole_image = numpy.ma.zeros([lx,ly])
        whole_image.mask = numpy.zeros_like(whole_image.data)

        for xmin in range(0, lx, lx/n_slices):
            xmax = min(xmin + lx/n_slices, lx)

            # Now we can build and sort a section of the cube with all the images
            cube = cube_images(list1, args.mask_key, scales, limits=[xmin, 0, xmax, ly])
            cube.sort(axis=0)

            # Finally, average! Remember that the cube is sorted so that
            # cube[0,ii,jj] < cube[1,ii,jj] and that the highest values of all
            # are the masked elements. We will take advantage of it if the median
            # is selected, because nowadays the masked median is absurdly slow:
            # https://github.com/numpy/numpy/issues/1811
            map_cube = numpy.ma.count(cube, axis=0) # number non-masked values per pixel
            if args.average == "mean":
                image = numpy.ma.mean(cube, axis=0)
                non_masked_equivalent = numpy.mean(cube.data, axis=0)
            elif args.average == "median":
                image = home_made_median(map_cube, cube)
                non_masked_equivalent = numpy.median(cube.data, axis=0)

            # Image is a masked array, we need to fill in masked values with the
            # args.fill_val if user provided it. Also, values with less than
            # args.nmin valid values should be masked out. If user did not provide
            # a fill_val argument, we will substitute masked values with the
            # unmasked equivalent operation.
            image.mask[map_cube < args.nmin] = 1
            mask = image.mask
            if args.fill_val != '':
                image = image.filled(args.fill_val)
            else:
                image.data[mask == True] = non_masked_equivalent[mask == True]
                image = image.data

            whole_image.data[xmin:xmax, 0:ly] = image[:,:]
            whole_image.mask[xmin:xmax, 0:ly] = mask[:,:]

        # And save images. If all_together is activated, use the file name given by user. If not, we need
        # to separate by filter, so compose a new name with the one given by the user adding the filter
        if args.all_together:
            newfile = args.output
        else:
            newfile = os.path.join(outdir, utils.add_suffix_prefix(outfile, suffix="_" + filt) )

        if args.out_mask != "":
            name_mask = args.out_mask
        else:
            name_mask = newfile + ".msk"
        if os.path.isfile(newfile):
            os.remove(newfile)
        if os.path.isfile(name_mask):
            os.remove(name_mask)
        fits.writeto(newfile, whole_image.data)
        fits.writeto(name_mask, whole_image.mask.astype(numpy.int))
        result[filt] = newfile

        # Add comments to the headers
        string1 = " - Image built from the combination of the images: "+\
                 ", ".join(list1)
        string2 = " combine = " + args.average + ", scale = " + args.scale
        utils.add_history_line(newfile, string1 + string2 )
        utils.add_history_line(name_mask, " - Mask of image: " + newfile)
        if args.mask_key != "":
            utils.header_update_keyword(newfile, args.mask_key, name_mask,
                                        "Mask for this image")

        # To normalize calculate median and call arith_images to divide by it.
        if args.norm == True:
            median = compute_scales([newfile], args.scale, args.mask_key)[0]
            msg =  "- NORMALIZED USING MEDIAN VALUE:"
            arith_images.main(arguments=["--message", msg, "--output", newfile,
                                         "--mask_key", args.mask_key, "--fill_val",
                                         args.fill_val, newfile, "/", str(median)])
    return result

############################################################################

# Create parser
parser = argparse.ArgumentParser(description='Combine images')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', \
                 help='input images to combine.', nargs='+', type=str)
parser.add_argument("--output", metavar="output", dest='output', \
                   action='store', help='Name for output file.')
parser.add_argument("--average", metavar='average', type=str, default='median', \
                   help='type of average (median, mean) to combine ' +\
                   'the images. Default: median')
parser.add_argument("--scale", metavar='scale', type=str, default='none', \
                   help='scaling function (median, mean, mode, none) to apply' +\
                   'to	the images before combining them. Default: none' )
parser.add_argument("--all_together", action="store_true", dest="all_together", \
                   default=False, help=' Force all the files to be combined'+ \
                   'together, i.e. do not separate by filter (e.g. for bias)') 
parser.add_argument("--norm", action="store_true", dest="norm", default=False, \
                    help="Normalize resulting image? Default: No")
parser.add_argument("--nmin", metavar="nmin", type=int, dest="nmin", action='store', \
                    default=1, help="Minimum number of images with valid "+\
                    "(i.e. non masked, non rejected) pixels. If the valid pixels "+\
                    "are less than nmin, the fill_val value is used as result, "+\
                    "and the pixel is masked. Default: 1.")
parser.add_argument("--filterk", action="store", dest="filterk", default='filter', \
                    help="Keyword in the header that contains the name of "+\
                    "the filter. Alternatively you can provide the filter "+\
                    "name itself with the argument '--filter'")
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

  # Check that either the filterk or a config file is present:
  if args.filterk == "" and args.config == "":
      sys.exit("ERROR! Either --filterk or a --config_file are "+\
               "necessary to determine the filter of the images. Type python "+\
               "combine.py -h for help.")
	
  # Call combine, keep name of the file created
  newfile = combine(args)
  return newfile  

if __name__ == "__main__":
    main()
