#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import tempfile
import shutil
import argparse
import sys
import os
from repipy import utilities
from lemon import methods
import repipy
import astropy.io.fits as fits

# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    import pyraf.iraf as iraf
    from iraf import imcopy


def trim(args):
    # Define region to be trimmed
    y0, y1, x0, x1 = args.region

    # args.output should be a list of output names. If they do not exist, the outputs should be the same as the
    # inputs with whatever suffix, if any, the user gave
    if not args.output:
        args.output = [utilities.add_suffix_prefix(im_name, suffix=args.suffix) for im_name in args.input]

    # Do the actual trimming. We will do it first into a temporary file, then copy it into args.output. This is just
    # in case the output and input filenames are the same, or if the output exists. IRAF will not overwrite!
    for im_name, new_name in zip(args.input, args.output):
        basename = os.path.splitext(os.path.basename(im_name))[0]
        _, temp_output = tempfile.mkstemp(prefix=basename, suffix=".fits")
        os.unlink(temp_output )
        with utilities.tmp_mute():
            imcopy(im_name + "[{0}:{1}, {2}:{3}]".format(x0, x1, y0, y1), temp_output)
            shutil.move(temp_output, new_name)

            # If a mask exists, trim exactly equally
            if args.mask_key:
                maskname = fits.getheader(im_name)[args.mask_key]
                with fits.open(maskname, 'readonly') as mask_im:
                    mask_im[0].data = mask_im[0].data[y0:y1+1, x0:x1+1]
                    mask_output = utilities.replace_extension(new_name, ".fits.msk")
                    fits.writeto(mask_output, mask_im[0].data, mask_im[0].header, clobber=True)
                utilities.header_update_keyword(new_name, "MASK", os.path.abspath(mask_output),
                                                comment="Name of mask image. ")

    return args.output

# Create parser
parser = argparse.ArgumentParser(description='Trim a given section of an image, keeping WCS information true')
parser.add_argument("input", metavar='input', action='store', nargs="+", \
                    type=str, help='Images to be trimmed.')
parser.add_argument("--output", metavar='output', action='store', dest='output', \
                    default='', help='Name of output file.')
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store', \
                    default='', help='suffix to be added at the end of the image' + \
                                     'input list to generate the output.')
parser.add_argument("--region", metavar=('x0', 'x1', 'y0', 'y1'), action='store', nargs=4, type=int,
                    required=True, dest="region", \
                    help='Region of the image to be trimmed. x0,x1 is the range in the horizontal '
                         'axis as shown by ds9, while y0, y1 mean the same in the vertical axis.'
                         'Please, note that what ds9 shows as X (i.e. the horizontal axis) is the second index in '
                         'a numpy array, the Y axis as shown by ds9 being the first. '
                         'Usually x0,x1 correspond to Naxis1 and y0,y1 to Naxis2 in the header of fits images. ')
parser.add_argument("--overwrite", dest='overwrite', action='store_true',
                    help='Overwrite the input file if no suffix or output file is given. ')
parser.add_argument("--mask_key", metavar='mask_key', action='store',
                    dest='mask_key', default="",
                   help='Key where the name of the mask image is stored in the header. The mask image should be trimmed '
                        'to be the same size as the output image.')


def main(arguments=None):
    # Pass arguments to variable
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    if args.output == '' and (not args.overwrite) and args.suffix == '':
        sys.exit("Error! Introduce a suffix, an output filename or use the '--overwrite' option. " + \
                 "For help: python trim_images.py -h ")
    if args.suffix != "":
        args.suffix = args.suffix.strip()

    newfile = trim(args)
    return newfile

if __name__ == "__main__":
    main()
