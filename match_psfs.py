#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 10:53:55 2014

@author: blasco
"""

############################################################################
import repipy.calculate_psf_image as psf
import repipy.utilities as utils
import argparse
import sys
import shutil
import numpy

from lemon import methods
import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot

def sort_by_seeing(args):
    """ Given the list of images of the input, put as the first of the list the one with the largest seeing
     (it will be the reference image) """
    seeings = numpy.array(utils.collect_from_images(args.input, args.FWHM_key))
    arg_max = seeings.argmax()
    args.input[0], args.input[arg_max] = args.input[arg_max], args.input[0]
    args.input_stars[0], args.input_stars[arg_max] = args.input_stars[arg_max], args.input_stars[0]



def match(args):
    """ Match the PSF of a group of images with the PSF of a reference image.
    """
    sort_by_seeing(args)

    # Get seeing from header and calculate psf of reference image, the first one since we sorted the input
    ref_seeing = utils.get_from_header(args.input[0], args.FWHM_key)
    ref_psf = psf.main(arguments=[args.input[0], "--stars", args.input_stars[0],
                                  "--sigma_key", args.sigma,
                                  "--gain_key",  args.gain_key,  
                                  "--ron_key",   args.ron_key,  
                                  "--expt_key",  args.expt_key, 
                                  "--airm_key",  args.airm_key, 
                                  "--FWHM_key", args.FWHM_key])
 
    output_list = []
    for image, stars in zip(args.input, args.input_stars):
        output = utils.add_suffix_prefix(image, prefix=args.prefix, suffix=args.suffix)

        # Too small differences of seeing are not worh doing any matching
        current_seeing = utils.get_from_header(image, args.FWHM_key)
        if abs(ref_seeing - current_seeing) < 0.2:  # Not worth equating PSFs for small differences
            shutil.copy(image, output)  # copy old file into new one
        else:
            # Calculate psf of the other images
            psf_object =  psf.main(arguments=[image,
                                              "--stars",     stars,
                                              "--sigma_key", args.sigma,
                                              "--gain_key",  args.gain_key,
                                              "--ron_key",   args.ron_key,
                                              "--expt_key",  args.expt_key,
                                              "--airm_key",  args.airm_key,
                                              "--FWHM_key", args.FWHM_key])
            utils.if_exists_remove("kernel.fits", output)
            iraf.psfmatch(image, ref_psf, psf_object, "kernel.fits",
                                 convolution="psf", 
                                 filter="cosbell")                   
            iraf.psfmatch(image, ref_psf, psf_object,
                                 convolution="kernel", 
                                 kernel="kernel.fits", 
                                 output=output,
                                 verbose="no")
            utils.if_exists_remove("kernel.fits")
        mssg = "Before equating PSFs: " + str(current_seeing)
        utils.header_update_keyword(output, args.FWHM_key, ref_seeing, comment=mssg)
        output_list.append(output)
    return output_list





# Create parser
parser = argparse.ArgumentParser(description='Calculate PSF of an image')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'names of input images to match PSFs.', nargs="+", 
                    type=str) 
parser.add_argument("--input_stars", metavar='input_stars', action='store', 
                    dest="input_stars",
                    help='Files containing a list of stars in two columns. If a single ' +\
                    'is passed, it will be assumed that it contains the stars to be used '+\
                    'for all images. Otherwise, you can list a file for each image '+\
                    'in "input"', nargs="+", type=str)
parser.add_argument("--coords", metavar='coords', action='store', default="world",
                    help='Type of coordinates that input_stars give. Allowed: '+\
                    ' "logical", tv", "physical", and "world". Default: world ')                             
parser.add_argument("--prefix", metavar="prefix", dest='prefix', action='store', \
                    default='', type=str, help='prefix to be added at the '+\
                    'beginning of the image input list to generate the outputs.',
                    nargs=1)
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',\
                    default=' -p', type=str, help='suffix to be added at the end '+\
                    'of the image input list to generate the outputs. There '+\
                    'is a peculiarity with argparse: if you pass, e.g., "-c" to '+\
                    '--suffix, the program will understand that you want to '+\
                    'call the code with the flag -c, which does not exist. This '+\
                    'does not raise an error, but just stops execution, which is '+\
                    'quite annoying. One way around it is " -c" (notice the '+\
                    'space, since within the code the string is stripped.', nargs=1)
parser.add_argument("--message", metavar="hdr_message", dest='hdr_message', \
                    action='store', default="", help='Message to be added to ' +\
                    'the header via HISTORY. For example: bias subtracted.')
parser.add_argument("--sigma_key", metavar="sigma_key", dest='sigma', \
                    action='store', default="", help='Keyword in the header ' +\
                    'containing an estimate of the noise of the sky.')
parser.add_argument("--min_val", metavar="minval", dest='minval', default=0, \
                   type=float, action='store', help='Minimum allowed value. '+\
                   'Below this value, mask out. Default: 0.')
parser.add_argument("--max_val", metavar="maxval", dest='maxval', action='store',\
                    default=50000, type=float, help='Maximum allowed value. '+\
                    'Above this value, mask out. Default: 50000.')     
parser.add_argument("--gain_key", metavar="gain_key", dest='gain_key', \
                    action='store', default="", help=' Keyword in the header ' +\
                    'of the image that contains the gain of the '+\
                    'camera. If not present, it will not be used.')
parser.add_argument("--ron_key", metavar="ron_key", dest='ron_key', \
                    action='store', default="", help=' Keyword in the header ' +\
                    'of the image that contains the read-out-noise of the '+\
                    'camera. If not present, it will not be used.')
parser.add_argument("--expt_key", metavar='expt_key', action='store', \
                     default = '', help='Name of the keyword in the headers that'+\
                     ' contain the exposure time. This can also be provided '+\
                     'with the --config_file option')  
parser.add_argument("--airm_key", metavar='airm_key', action='store', \
                     default = '', help='Name of the keyword in the headers that'+\
                     ' contain the airmass. This can also be provided '+\
                     'with the --config_file option')  
parser.add_argument("--FWHM_key", metavar='seeing', dest='FWHM_key', action='store', 
                    default='SEEING', help='Keyword in which the seeing is stored ' +\
                   'in the headers. Default: "SEEING" ' )

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  
  args = parser.parse_args(arguments)
  
    # In order to allow "  -c" to be able to use the hyphen as suffix.
  if args.suffix != "":
      args.suffix = args.suffix[0].strip()
  if args.prefix != "":
      args.prefix = (args.prefix[0]).strip()  
  
  # Number of reference star files should be one or as many as input images
  if len(args.input_stars) != 1 and (len(args.input) != len(args.input_stars)):
      sys.exit("Error! You should put one input_stars file or as many as " +\
               "input images.")
  elif len(args.input_stars) == 1:
      args.input_stars =  args.input_stars * len(args.input)  

  newnames = match(args)  
  return newnames    
     
if __name__ == "__main__":
    main()
