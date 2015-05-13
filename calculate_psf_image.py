#! /usr/bin/env python
# -*- coding: UTF-8 -*-

###############################################################################
import os
import argparse
import numpy
from scipy.stats import mode as mode
import sys
import astropy.io.fits as fits
import collections
import repipy
import repipy.utilities as utils
import repipy.arith as arith_images
from lemon import methods

# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot




###################################################################################
def psf(args):
    """ Calculate the PSF of an image.
    """
    # Read the seeing and sigma of the sky from the header
    seeing, sigma = utils.get_from_header(args.input, args.FWHM_key, args.sigma)
    
    # Do photometry on the image
    #print "photometry: \n"
    photfile_name = args.input + ".mag.1"
    utils.if_exists_remove(photfile_name)
    iraf.phot(args.input, output=photfile_name, coords=args.stars,
                     wcsin=args.coords, fwhm=seeing, 
                     sigma=sigma, datamax=args.maxval, datamin=args.minval,
                     ccdread=args.ron_key, gain=args.gain_key, exposure=args.expt_key,
                     airmass=args.airm_key, annulus=36, dannulus=18,
                     apert=18, verbose="no", verify="no", interac="no")


    # Select stars on the image                 
    #print "pstselect: \n"
    pstfile_name = args.input + ".pst.1"
    utils.if_exists_remove(pstfile_name)
    iraf.pstselect(args.input, photfile=photfile_name, pstfile=pstfile_name,
                          maxnpsf=20, fwhm=seeing, sigma=sigma,
                       datamax=args.maxval, ccdread=args.ron_key, gain=args.gain_key,
                       exposure=args.expt_key,  function="auto", nclean=1, 
                       psfrad=36, fitrad=18, maxnstar=20, verbose="no", interac="no",
                       verify="no")

    # Build psf of the stars
    #print "psf: \n"
    psffile_table = args.input + ".psf.1.fits"  # iraf keeps adding the .fits :(
    psgfile_name = args.input + ".psg.1"
    pstfile_name2 = args.input + ".pst.2"   
    utils.if_exists_remove(psffile_table,psgfile_name, pstfile_name2)
    iraf.psf( args.input, photfile=photfile_name, pstfile=pstfile_name,
                     groupfile=psgfile_name, opstfile=pstfile_name2,    
                     psfimage=psffile_table,fwhm=seeing, sigma=sigma, datamax=args.maxval, 
                     datamin=args.minval, ccdread=args.ron_key, gain=args.gain_key, 
                     exposure=args.expt_key, function="moffat25", nclean=1, 
                     psfrad=36, fitrad=18, maxnstar=20, interactive="no",
                     varorder=args.varorder, verbose="no",verify="no")

    # Use seepsf to build the image of the psf
    psffile_name = args.input + ".psf.fits" 
    utils.if_exists_remove(psffile_name)
    iraf.seepsf(psffile_table, psffile_name)
     
    return psffile_name

############################################################################

# Create parser
parser = argparse.ArgumentParser(description='Calculate PSF of an image')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='Name ' +\
                    'of image to calculate PSF.') 
parser.add_argument("--stars", metavar='stars', action='store', 
                    dest="stars",
                    help='File containing a list of stars in two columns. ' +\
                    'If in wcs or in pixels is determine by the parameter '+\
                    '"coords"')                       
parser.add_argument("--coords", metavar='coords', action='store', default="world",
                    help='Type of coordinates that "stars" give. Allowed: '+\
                    ' "logical", tv", "physical", and "world". Default: world ')                             
parser.add_argument("--sigma_key", metavar="sigma_key", dest='sigma', \
                    action='store', default="", help=' Keyword in the header ' +\
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
                    default='3.', help='Keyword/value of the seeing, depending ' +\
                   'on if the provided value is a string or a float. Default '+\
                   'value is four pixels. ' )
parser.add_argument("--varorder", metavar='varorder', dest='varorder', action='store',
                    default=-1, help=' From IRAF daopars: "The order of variability ' +
                    'of the PSF model computed by the DAOPHOT PSF task." ' +\
                   'Possible values: \n '
                   '       -1: produce a single analytic model of the PSF constant over image '+\
                   '        0: produce an analytic model and a lookup table, constant over image '+\
                   '        1: analytic model and 3 lookup tables, model changes linearly over image '+\
                   '        2: analycit model and 5 lookups, model changes quadratically over image ' +\
                   ' Default: -1.')


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  
  args = parser.parse_args(arguments)
  psf_name = psf(args)  
  return psf_name    
     
if __name__ == "__main__":
    main()



"""




"""