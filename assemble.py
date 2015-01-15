#!/usr/bin/env python
# -*- coding: utf-8 -*-
import astropy.io.fits as fits
import tempfile
import sys
from lemon import mosaic
import repipy.utilities as utilities
import argparse

def find_primaryHDU(HDUList):
    """ Find the primary HDU in an HDUList """
    main_hdrs = [layer for layer in HDUList if isinstance(layer, fits.hdu.image.PrimaryHDU)]
    if main_hdrs:
        return main_hdrs[0].header
    return None

def find_chips(HDUList):
    """ Find the chips within a HDUList """
    CompImageHDU = fits.hdu.compressed.CompImageHDU
    imageHDU = fits.hdu.image.ImageHDU
    chips = [layer for layer in HDUList if isinstance(layer, (CompImageHDU, imageHDU))]
    return chips


def create_images(chip_list):
    """ Create a set of images in the /tmp directory from a list of imageHDU objects. 
        Each object in the list is a chip from a multi-layered fits file. The program 
        returns the list of file names in the /tmp directory """
    output = []
    for layer in chip_list:  
        fd, path = tempfile.mkstemp(prefix="temp_im_", suffix=".fits")
        fits.writeto(path, layer.data, layer.header, output_verify="Ignore")
        output.append(path)
    return output

def build(im=None, output_file=None):
    # Read the image, separate the main header (if present) and the N chips
    im_HDUList = fits.open(im)
    main_hdr = find_primaryHDU(im_HDUList)
    chip_list = find_chips(im_HDUList)    
    
    # Create images in the /tmp directory from the N chips
    temp_images = create_images(chip_list)

    # Do the mosaic
    arguments = temp_images + [output_file, "--cores", "1", "--no-reprojection", "--overwrite"]
    mosaic.main(arguments)

    # Remove WCS from main header and plug it into the image
    main_hdr = utilities.remove_WCS(main_hdr)
    output_image = fits.open(output_file, 'update')
    output_image[0].header += main_hdr
    output_image.flush()
    output_image.close()

############################################################################

# Create parser
parser = argparse.ArgumentParser(description='Compose image')

# Add necessary arguments to parser
parser.add_argument("inputs", metavar='input', action='store', \
                 help='input images to combine.', nargs='+', type=str)
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',\
                    default=' -m', type=str, help='suffix to be added at the end '+\
                    'of the image input list to generate the outputs. There '+\
                    'is a peculiarity with argparse: if you pass, e.g., "-c" to '+\
                    '--suffix, the program will understand that you want to '+\
                    'call the code with the flag -c, which does not exist. This '+\
                    'does not raise an error, but just stops execution, which is '+\
                    'quite annoying. One way around it is " -c" (notice the space, '+\
                    'since within the string is stripped within the code.  '+\
                    ' You need to chose either a suffix or activate the --overwrite', nargs=1)
parser.add_argument("--overwrite", action="store_true", dest="overwrite", \
                    default=False, help="Allows you to overwrite the original image.")
    
def main(arguments=None):
    if arguments is None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    
    # In order to allow "  -c" to be able to use the hyphen as suffix.
    if args.suffix != "":
        args.suffix = args.suffix[0].strip()
    
    # Detecting errors
    if args.suffix == '' and args.overwrite == False:  
        sys.exit("Error! Introduce a suffix or activate the --overwrite option. \
		  For help: python arith.py -h ") 
   
    for im in args.inputs:
        if args.suffix != "":
            output_file = utilities.add_suffix_prefix(im, suffix=args.suffix)
        elif args.overwrite:
            output_file = im
        build(im, output_file)

if __name__ == "__main__":
    sys.exit(main())



