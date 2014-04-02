#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Created on Wed Aug 21 09:02:07 2013

@author: blasco
"""
import numpy
import sys
import argparse
import repipy.utilities as utils
import astropy.io.fits as fits
import skimage.filter
          
def filter_image(args):
    """ Routine that uses a median filter on a list of masked images. """
    output_list = []
    for image in args.input:
        # Read image, mask and header
        im = utils.read_image_with_mask(image, mask_keyword=args.mask_key)
        hdr = fits.getheader(image)
        
        # skimage uses masks where 1 means valid and 0 invalid
        mask = (im.mask + 1) % 2
        filt_im = skimage.filter.median_filter(im.data, mask=mask, radius=args.radius)
        
        # Make name of file if not given
        if args.output == "":
            output = utils.add_suffix_prefix(image, suffix='-mf')
        else :
            output = args.output
        
        # Add history line to the header and write to file
        hdr.add_history("- Image median filtered. Radius = " + str(args.radius))
        fits.writeto(output, filt_im, header=hdr)
        output_list.append(output)
    return output_list
    
############################################################################
# Create parser
parser = argparse.ArgumentParser(description='Combine images')
# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'input images to be median filtered.', \
                    nargs="+", type=str)
parser.add_argument("--mask_key", metavar="mask_key", dest='mask_key', \
                    action='store', default="", help=' Keyword in the header ' +\
                    'of the image that contains the name of the mask. The mask '+\
                    'will contain ones (1) in those pixels to be MASKED OUT.')
parser.add_argument("--radius", metavar="radius", type=int, dest="radius", 
                    action='store', required=True, help=" Radius of circle "+\
                    "to be used to filter the image. Mandatory argument.")
parser.add_argument("--output", metavar='output', dest='output', action='store',
                   default='', help='output image in which to save the result.'+\
                    ' If not provided the suffix -mf will be added to '+\
                    'the input.')          

############################################################################                   

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  args = parser.parse_args(arguments)
  
  # Call combine, keep name of the file created
  newfile = filter_image(args)
  return newfile  

if __name__ == "__main__":
    main()
                  
