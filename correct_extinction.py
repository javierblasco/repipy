#! /usr/bin/env python
# -*- coding: utf-8 -*-


import repipy.utilities as utilities
import repipy.astroim as astroim
import numpy
import astropy.io.fits as fits
import argparse
import sys
from scipy import spatial
import astropy.wcs.wcs as wcs

from lemon import methods
import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import obsutil


def test1(coeff):
    """ Check that either a single coefficient was passed or the names of the filters were passed

    The user passes the coefficients, and they are packed into a list of lists. This means that a valid coeff is
    coeff = [['rGunn', 0.21, 0.03]] for the name of the filter, coefficient and uncertainty, respectively.
    coeff = [[0.21]] is also valid, and will be used for all images, irrespective of their filter.

    But if two coefficients are passed, we need to check that all of them have names for the filters, otherwise they
    would be impossible to distinguish which one to apply to each image.

    :param coeff: coefficients passed by the user
    :return: True if the number of elements is one, or all the coefficients come with the name of the filter
    """

    num_elements = len(coeff)
    if num_elements == 1:
        return True

    for filt_coeff in coeff:
            # If the first thing you find in each of the lists is a number, the filter names are not there
            try:
                float(filt_coeff[0])
                message = "\n\nIf you pass several filters you MUST name them all. \n " \
                          "See correct_extinction.py -h  for help and examples. \n " \
                          "Exiting code!\n\n"
                sys.exit("\n\n {0} \n\n".format(message))
            except ValueError:
                pass
    return True


def test2(coeff):
    """ Check if all the elements of coeff have less than three elements each
    :param coeff:
    :return:
    """
    for elem in coeff:
        if len(elem) > 3: # more than three values? Nonsense
            message = "\n\nWrong number of inputs: --coeff {0}. You can pass up to three things: \n " \
                      "-filter name (optional) \n " \
                      "-extinction coefficient (mandatory) \n " \
                      "-uncertainty in the coefficient (optional) \n" \
                      "See correct_extinction.py -h  for help and examples. \n\n""".format(" ".join(elem))
            sys.exit("\n\n{0}\n\n".format(message))
    return True

def pack_coefficients(coeff):
    """
    Pack all the coefficients passed by the user into a dictionary
    :param coeff: coefficients passed by the user to be used in the program
    :return: True (passed all tests) or False (did not pass at least one of them)
    """

    # Plenty of places where the user can do it wrong, let's help him a little ;)
    assert test1(coeff) * test2(coeff) == 1

    coeff_dict = {}
    for element in coeff:
        try:
            # is it all numbers? Then it has to be a single coefficient for all filters
            coeff_dict['All'] = numpy.array(element, dtype=numpy.float)
        except ValueError:  # The first thing must be the filter name
            coeff_dict[element[0]] = numpy.array(element[1:], dtype=numpy.float)
    return coeff_dict

def build_outputs(args):
    """ From the arguments passed by the user decide on the output names for every input.
    :param args: argparse object with the parameters passed by the user
    :return: list with the names of the output files
    """
    output_files = []
    for im_name in args.input:
        new_name = im_name
        output_files.append(utilities.add_suffix_prefix(new_name, suffix=args.suffix))
    return output_files



def correct(args):
    """ Program to correct images from atmospheric extinction. The inputs
    are the images and the extinction coefficient.
    """

    args.outputs = build_outputs(args)
    for im_name, new_name in zip(args.input, args.outputs):
        im = astroim.Astroim(im_name)

        # The user might have put the string that appears in the header literal (im.filter.filter_name), the name
        # of the filter given when you print im.filter, which is im.filter.__str__ or none at all because only one
        # unnamed value was passed (hence the "All").
        possible_keys = (im.filter.filter_name, im.filter.__str__(), "All")

        # Try the possible keys until you find one
        for kk in possible_keys:
            if kk in args.coeff.keys():
                result = args.coeff[kk]

        try:
            length_of_result = len(result)
        except UnboundLocalError:
            print "\n\n{0} has filter {1}, but no coefficient was passed with that filter name. " \
                  "No action was performed for this image. \n\n".format(im_name, im.filter.filter_name)
            continue


        if len(result) == 2:
            coeff, uncertainty = result
        elif len(result) == 1:
            coeff, uncertainty = result[0], ""

        # Get the airmass from the header of the image
        airm_keyword = im.primary_header.airmassk
        airmass = im.primary_header.airmass
        if airmass == 0:
            print "{0} has zero airmass. Nothing to be done. ".format(im_name)
            continue

        # Correct the airmass
        im = fits.open(im_name, 'readonly')
        im[0].data = im[0].data * 10**(coeff * airmass / 2.5 )
        im[0].header[airm_keyword] = (0, "Airmass before correcting: {0}".format(airmass))
        im[0].header[args.extinctk] = (coeff, "Extinction coefficient used to correct")
        if uncertainty:
            im[0].header[args.extinctk + "_std"] = (uncertainty, "Uncertainty in the extinction coefficient")

        fits.writeto(new_name, im[0].data, im[0].header, clobber=True)





    return


############################################################################


# Create parser
parser = argparse.ArgumentParser(description='Program to correct from extinction using one or several extinction '
                                        'coefficients (one per filter). You can pass a list of images with several '
                                        'filters, as long as you pass an extinction coefficient for each of them. '
                                        'the extinction coefficient can be passed including an error. Example:'
                                        ''
                                        'correct_extinction.py --coeff "r Gunn" 0.21 0.03 --coeff Ha6645 0.27 *.fits '
                                        ''
                                        'correct_extinction.py --coef 0.21 --extinctk COEFF image.fits'
                                        ''
                                        'In the first example, we pass two coefficients, one for the filter rGunn, '
                                        'and one for the filter Ha6645. Those must coincide with the content of the '
                                        'keyword holding the filter in the header (many accepted keywords, see '
                                        'header.py). The coefficient of rGunn has an associated error, while the one '
                                        'of Halpha does not have one. In the second example, one image is passed, and '
                                        'the coefficient will be saved in the header under the COEFF keyword. In any '
                                        'case, as the atmospheric coefficient is being corrected, the resulting image '
                                        'must have an AIRMASS of zero. The old airmass value will be kept in the '
                                        'comments of the keyword holding the airmass. ')

# Add necessary arguments to parser
parser.add_argument('--input', metavar='input', action='store', dest="input", help='list of input images to correct '
                    'for extinction.', nargs="+", type=str)
parser.add_argument("--coeff", metavar='coeff', action='append', dest="coeff",
                    help='Extinction coefficient(s) to correct images with. If you pass a single coefficient, you '
                         'may pass it alone, but if you pass several, one per filter, you must add the name of the '
                         'filter as it appears in the header of the images (and the names, if you used rename.py '
                         'before). You can choose to send in also the error for the coefficient if you have it. '
                         ' Some examples: '
                         '"--coeff 0.23", "--coeff rGunn 0.23", "--coeff 0.21 0.03" or '
                         '"--coeff rGunn 0.21 --Ha6645 0.23 0.04 "    ', nargs='+', type=str)
parser.add_argument("--extinctk", metavar='extinctk', action='store', dest="extinctk", default="ext_coef",
                    help='Keyword for the extinction coefficient in the header of the images. Default: "Ext_coef". ',
                    nargs=1, type=str)
parser.add_argument("--overwrite", action="store_true", dest="overwrite", \
                    default=False, help="Allows you to overwrite the original image.")
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',\
                    default='', type=str, help='suffix to be added at the end '+\
                    'of the image name(s) to generate the output(s). There '+\
                    'is a peculiarity with argparse: if you pass, e.g., "-c" to '+\
                    '--suffix, the program will understand that you want to '+\
                    'call the code with the flag -c, which does not exist. This '+\
                    'does not raise an error, but just stops execution, which is '+\
                    'quite annoying. One way around it is " -c" (notice the '+\
                    'space, since within the string is stripped.')

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)


  # The --coeff keyword is very versatile, we need to make sure they are tested to check a valid option was passed
  args.coeff = pack_coefficients(args.coeff)

  # remove left spaces from the suffix
  args.suffix = args.suffix.lstrip()

  # Detecting errors
  if args.suffix == '' and args.overwrite == False:
      sys.exit("Error! Introduce a suffix or the --overwrite option. For help: python arith.py -h ")

  correct(args)
  return None

if __name__ == "__main__":
    main()
