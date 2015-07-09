#! /usr/bin/env python
# -*- coding: UTF-8 -*-
import shutil
import sys
import argparse
import numpy
import astropy.io.fits as fits
import repipy.utilities as utilities
from scipy.ndimage.filters import median_filter


def fit_pol(y, deg):
    """ Fit a polynomial to a vector. Return the model values
    :param y:
    :return:
    """
    x = numpy.arange(len(y))
    coeff = numpy.polyfit(x, y, deg)
    fit = numpy.poly1d(coeff)
    return fit(x)


def subtract(args):
    for im_name in args.input1:
        if args.overwrite:
            new_name = im_name
        elif args.suffix:
            new_name = utilities.add_suffix_prefix(im_name, suffix=args.suffix)

        # Read image, separate data and header
        im = fits.open(im_name)
        data = im[0].data
        hdr = im[0].header

        # Extract the overscan region. Notice that what we call x,y are the second and first axes
        y0, y1, x0, x1 = args.region
        overscan = data.copy()[x0:x1, y0:y1]


        # Average over the short axis
        if overscan.shape[0] < overscan.shape[1]:
            average = numpy.nanmedian(overscan, axis=0)
            # Fit a polynomial and return the fitted values
            fitted_overscan = fit_pol(average, args.deg)
            data[:, y0:y1] -= fitted_overscan
        else:
            average = numpy.nanmedian(overscan, axis=1)
            # Fit a polynomial and return the fitted values
            fitted_overscan = fit_pol(average, args.deg)
            data[x0:x1, :] = (data[x0:x1, :].T - fitted_overscan).T


        # Write to the output file
        hdr.add_comment("Overscan region subtracted. Region: [{0}:{1},{2}:{3}]".format(x0, x1, y0, y1))
        fits.writeto(new_name, data, hdr, clobber=True)

    return None












# Create parser
parser = argparse.ArgumentParser(description="""Subtract an average of the overscan region from the whole image.
                                                 This routine will load the image, read the region that you have
                                                 given through the parameter --region x0 x1 y0 y1 , average the short
                                                 axis and subtract a fit to the long one from the image.
                                                 For example, a region of --region 0 1000 980 1000 will subtract a
                                                 region of shape 1000x20, so it will use a median on the short axis,
                                                 perform a 3rd degree polynomial fit to the long one and subtract the
                                                 resulting fit from the corresponding region of the image.

                                             """)

# Add necessary arguments to parser
parser.add_argument("input1", metavar='input1', action='store', help='list of ' +\
                    'input images from which to subtract the overscan region', nargs="+", type=str)

mandatory = parser.add_argument_group('Mandatory argument')
mandatory.add_argument("--region",metavar=('x0', 'x1', 'y0', 'y1'), action='store', nargs=4, type=int,
                       required=True, dest="region", \
                    help='Region of the image where the overscan is situated. x0,x1 is the range in the horizontal '
                          'axis as shown by ds9, while y0, y1 mean the same in the vertical axis.'
                          'Please, note that what ds9 shows as X (i.e. the horizontal axis) is the second index in '
                          'a numpy array, the Y axis as shown by ds9 being the first. '
                          'Usually x0,x1 correspond to Naxis1 and y0,y1 to Naxis2 in the header of fits images. ')
parser.add_argument("--overwrite", action="store_true", dest="overwrite", \
                    default=False, help="Allows you to overwrite the original image.")
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',\
                    default='', type=str, help='suffix to be added at the end '+\
                    'of the image input list to generate the outputs. There '+\
                    'is a peculiarity with argparse: if you pass, e.g., "-c" to '+\
                    '--suffix, the program will understand that you want to '+\
                    'call the code with the flag -c, which does not exist. This '+\
                    'does not raise an error, but just stops execution, which is '+\
                    'quite annoying. One way around it is " -c" (notice the '+\
                    'space, since within the string is stripped.')
parser.add_argument("--deg", metavar="deg", dest='deg', action='store',\
                    default=0, type=int, help='Degree of the polynomial to fit the '+\
                    'overscan region. Default:0. ')


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)

   # Detecting errors
  if args.suffix == '' and args.overwrite == False:
      sys.exit("Error! Introduce a suffix or the --overwrite option. For help: python arith.py -h ")

  # If --suffix " -c"
  args.suffix = args.suffix.strip()

  newname = subtract(args)

  return newname

if __name__ == "__main__":
    main()
