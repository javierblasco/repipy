#! /usr/bin/env python
# -*- coding: utf-8 -*-

import tempfile
import argparse
import sys
import os


import lemon.methods as methods
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(os.path.dirname(os.path.abspath(__file__))):
    import pyraf.iraf as iraf
    from iraf import digiphot
    from iraf import daophot


import repipy.astroim as astroim
import repipy.utilities as utils
import warnings
warnings.filterwarnings("ignore")



#fd, coords_file = tempfile.mkstemp(prefix=basename, suffix=".coords")

def apply_phot(args):
    # Temporary file
    with tempfile.NamedTemporaryFile(suffix=".mag.1") as fd:
        phot_modstars = fd.name

    with tempfile.NamedTemporaryFile(suffix=".mag.1") as fd2:
        phot_subtstars = fd2.name


    hdr = args.image.header  # for short
    seeing, sigma = utils.get_from_header(args.image.im_name, hdr.seeingk, hdr.sigmak)
    iraf.phot(args.image.im_name, output=phot_modstars, coords=args.model_stars,
              wcsin=args.coords, fwhmpsf=seeing, sigma=sigma, ccdread=hdr.ccdronk,
              gain=hdr.gaink, exposure=hdr.exptimek,
              airmass=hdr.airmassk, annulus=6*seeing, dannulus=3*seeing,
              apert=3*seeing, verbose="no", verify="no", interac="no")

    iraf.phot(args.image.im_name, output=phot_subtstars, coords=args.subt_stars,
              wcsin=args.coords, fwhmpsf=seeing, sigma=sigma, ccdread=hdr.ccdronk,
              gain=hdr.gaink, exposure=hdr.exptimek,
              airmass=hdr.airmassk, annulus=6*seeing, dannulus=3*seeing,
              apert=3*seeing, verbose="no", verify="no", interac="no")
    return phot_modstars, phot_subtstars

def apply_pst(args, photfile):
    # Temporary file
    with tempfile.NamedTemporaryFile(suffix=".pst.1") as fd:
        pst_file = fd.name

    hdr = args.image.header  # for short
    seeing, sigma = utils.get_from_header(args.image.im_name, hdr.seeingk, hdr.sigmak)
    iraf.pstselect(args.image.im_name, photfile=photfile, pstfile=pst_file,
                       fwhm=seeing, sigma=sigma, ccdread=hdr.ccdronk, gain=hdr.gaink,
                       exposure=hdr.exptimek,  function="auto", nclean=1,
                       psfrad=6*seeing, fitrad=3*seeing, verbose="no",
                       verify="no")
    print "\n Pst file: ", pst_file
    return pst_file

def apply_psf(args, phot_file, pst_file):
    with tempfile.NamedTemporaryFile(suffix=".psf.1.fits") as fd0:
        psffile_table = fd0.name

    with tempfile.NamedTemporaryFile(suffix=".psf") as fd1:
        psf_name = fd1.name

    with tempfile.NamedTemporaryFile(suffix=".psg.1") as fd2:
        psg_name = fd2.name

    with tempfile.NamedTemporaryFile(suffix=".pst.2") as fd3:
        pst2_name = fd3.name

    hdr = args.image.header  # for short
    seeing, sigma = utils.get_from_header(args.image.im_name, hdr.seeingk, hdr.sigmak)
    iraf.psf( args.image.im_name, photfile=phot_file, pstfile=pst_file,
                     groupfile=psg_name, opstfile=pst2_name,
                     psfimage=psffile_table,fwhm=seeing, sigma=sigma,
                     ccdread=hdr.ccdronk, gain=hdr.gaink,
                     exposure=hdr.exptimek, function="auto", nclean=1,
                     psfrad=6*seeing, fitrad=5*seeing, interactive="no",
                     varorder=-1, verbose="no",verify="no")

    # Use seepsf to make an image of the PSF
    iraf.seepsf(psffile_table, psf_name)
    print "\n PSF file: ", psf_name
    return psf_name


def subtract_stars(args):
    """ Routine to subtract a set of stars defined by the user. We will model the PSF of some suitable stars in the
        images, then use that """

    # First substitute the names of the image in args.image by the corresponding astroim object, with more information
    args.image = astroim.Astroim(args.image)

    # Do photometry on the stars to be used to model the PSF
    phot_modstars, phot_subtstars = apply_phot(args)
    pst_file = apply_pst(args, phot_modstars)
    psf_file = apply_psf(args, phot_modstars, pst_file)

    hdr = args.image.header  # for short
    seeing, sigma = utils.get_from_header(args.image.im_name, hdr.seeingk, hdr.sigmak)
    utils.if_exists_remove(args.output)
    iraf.allstar(args.image.im_name, photfile=phot_subtstars, psfimage=psf_file, subimage=args.output,
                 psfrad=6*seeing, fitrad=3*seeing, wcsin='physical',
                 fitsky="yes", sannulus=4*seeing,wsannulus=3*seeing,
                 verify="no")



# Create parser
parser = argparse.ArgumentParser(description='Extract stars from an image. ')

# Add necessary arguments to parser
parser.add_argument("image", metavar='image', action='store', help='Name ' +\
                    'of fits image to calculate extract stars from.')
parser.add_argument("--model_stars", metavar='model_stars', action='store',
                    dest="model_stars",
                    help='File containing a list of stars in two columns. ' +\
                    'Parameter "coords" will determine if the coordinates are pixels  '+\
                    'or RA DEC.')
parser.add_argument("--coords", metavar='coords', action='store', default="world",
                    help='Type of coordinates that "stars" give. Allowed: '+\
                    ' "logical", tv", "physical", and "world". Default: world ')
parser.add_argument("--subt_stars", metavar='subt_stars', action='store', dest='subt_stars',
                    help='File containing the list of stars to be subtracted. ')
parser.add_argument("--output", metavar='output', action='store', dest='output',
                    help='Name of output file. If no output is indicated, the suffix "-e" will be added ' +\
                         'to the image name. ')


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)
  print "Model stars file", args.model_stars , type(args.model_stars)
  args.model_stars = os.path.abspath(args.model_stars)
  args.subt_stars = os.path.abspath(args.subt_stars)


  subtract_stars(args)
  return None

if __name__ == "__main__":
    main()


