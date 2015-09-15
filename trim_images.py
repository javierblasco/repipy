#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import tempfile
import shutil
from repipy import utilities
from lemon import methods
import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import imcopy

def trim(args):
    # Define region to be trimmed
    y0, y1, x0, x1 = args.region

    # args.output should be a list of output names. If they do not exist, the outputs should be the same as the
    # inputs with whatever suffix the user gave
    if not args.output:
        args.output = [utilities.add_suffix_prefix(im_name, args.suffix) for im_name in args.input]

    # Do the actual trimming
    for im_name, new_name in zip(args.input):
        # Do the operation first into a temporary file, then copy it into args.output. This will prevent the problem
        # that iraf's imcopy does not overwrite output filenames.
        basename = os.path.splitext(os.path.basename(im_name))[0]
        _, temp_name = tempfile.mkstemp(prefix=basename, suffix=".fits")
        imcopy(args.input[x0:x1,y0:y1], temp_name)
        shutil.move(temp_name, new_name)

    return args.output

# Create parser
parser = argparse.ArgumentParser(description='Trim a given section of an image, keeping WCS information true')
parser.add_argument("input",metavar='input', action='store', nargs="+",  \
                    type=str, help='Images to be trimmed.')
parser.add_argument("--output", metavar='output', action='store', dest='output', \
                    default='', help='Name of output file.')
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store', \
                    default='',help='suffix to be added at the end of the image' +\
                    'input list to generate the output.')
parser.add_argument("--region",metavar=('x0', 'x1', 'y0', 'y1'), action='store', nargs=4, type=int,
                       required=True, dest="region", \
                    help='Region of the image to be trimmed. x0,x1 is the range in the horizontal '
                          'axis as shown by ds9, while y0, y1 mean the same in the vertical axis.'
                          'Please, note that what ds9 shows as X (i.e. the horizontal axis) is the second index in '
                          'a numpy array, the Y axis as shown by ds9 being the first. '
                          'Usually x0,x1 correspond to Naxis1 and y0,y1 to Naxis2 in the header of fits images. ')
parser.add_argument("--overwrite", metavar='overwrite', dest='overwrite', action='store_true',
                    help='Overwrite the input file if no suffix or output file is given. ')

def main(arguments=None):
    # Pass arguments to variable
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    if args.output == '' and (not args.overwrite) and args.suffix == '':
        sys.exit("Error! Introduce a suffix, an output filename or use the '--overwrite' option. "+\
                 "For help: python trim_images.py -h ")
    if args.suffix != "":
      args.suffix = args.suffix.strip()

    newfile = trim(args)
    return newfile
if __name__ == "__main__":
    main()
