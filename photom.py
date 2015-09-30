#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import lemon.photometry as photometry
import argparse
import repipy.utilities as utilities
import numpy
import sys
from astropy import log
import astropy
import astropy.io.fits as fits
import astropy.wcs as wcs
import repipy.astroim as astroim
import sep as sextractor
import warnings


def detect_sources(image, cat_name=None):
    """ A generator of (ra, dec) tuples (PROOF OF CONCEPT) """

    log.debug("Reading FITS file")
    with astropy.io.fits.open(image) as hdulist:
        data = hdulist[0].data
        header = hdulist[0].header

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        log.info("Loading WCS data")
        wcs = astropy.wcs.WCS(header)

    try:
        log.info("Estimating background")
        background = sextractor.Background(data)
    except ValueError:
        # Fix for error "Input array with dtype '>f4' has non-native
        # byte order. Only native byte order arrays are supported".
        log.debug("Converting data to native byte order")
        data = data.byteswap(True).newbyteorder()
        background = sextractor.Background(data)

    log.info("Subtracting background")
    background.subfrom(data) # in-place


    # Global "average" RMS of background
    log.info("Detecting sources")
    rms_ntimes = 1.5
    while True:
        try:
            threshold = rms_ntimes * background.globalrms
            objects = sextractor.extract(data, threshold)

        except Exception as e:

            # Fix for error "internal pixel buffer full: The limit of 300000
            # active object pixels over the detection threshold was reached.
            # Check that the image is background subtracted and the detection
            # threshold is not too low. detection threshold". If a different
            # exception is raised, just re-raise it.

            if "threshold is not too low" in str(e):
                rms_ntimes *= 1.25
                log.debug("Internal pixel buffer full")
                log.debug("Retrying with threshold {0:.2} times RMS of background".format(rms_ntimes))
            else:
                raise
        else:
            break

    pixel_coords = numpy.empty([0,2])
    print "len(objects) = ", len(objects)
    print "Image = ", image
    sys.exit()

    # https://github.com/kbarbary/sep/blob/master/sep.pyx#L555
    log.info("Transforming pixel to celestial coordinates")
    for index in range(len(objects)):
        x, y = objects[index]['x'], objects[index]['y']
        pixel_coords = numpy.row_stack(pixel_coords, (x, y))

    print pixel_coords.shape

    # If cat_name is present, save the coordinates to a file
    if cat_name:
        fd = open(cat_name, 'w')
        for ra, dec in list(wcs.all_pix2world(pixel_coords[:,0], pixel_coords[:,1], 0)):
            fd.write( "{0}  {1} \n ".format(ra, dec) )
        fd.close()

    return list(wcs.all_pix2world(pixel_coords[:,0], pixel_coords[:,1], 0))



def do_photometry(args):
    #for ii, im_name in enumerate(args.images):
        # If the catalogue exists already, use it, if it doesn't, create it using SEP
        #else:
        #    cat_name = "mierda.cat"
        #    detect_sources(im_name, cat_name)

    hdr = astroim.Astroim(args.images[0]).header
    gaink, objectk, filterk, datek = hdr.gaink, hdr.objectk, hdr.filterk, hdr.datek
    exptimek, airmassk, timek   =  hdr.exptimek, hdr.airmassk, hdr.timek



    arguments_common = ["--uik", "", "--margin", "20", "--gaink", gaink, "--cbox", args.cbox, "--individual-fwhm",
                   "--objectk", objectk, "--filterk", filterk, "--datek", datek, "--expk", exptimek, "--fwhmk", "seeing",
                   "--airmk", airmassk, "--timek", timek, "--overwrite", "--coordinates",
                   args.coordinates[0], args.images[0]] + args.images + args.output)]
    if args.unit == "FWHM":
        photometry.main(arguments = ["--aperture", args.aperture, "--annulus", "6", "--dannulus", "2"] + arguments_common)
    elif args.unit == "pix":
        photometry.main(arguments= ["--aperture-pix", args.aperture, "--annulus-pix", "6", "--dannulus-pix", "2"] +\
                                    arguments_common)

    return args.output


# Create parser
parser = argparse.ArgumentParser(description='Arithmetic operations on images')

# Add necessary arguments to parser
parser.add_argument("--images", metavar='image', action='store', dest="images", nargs="+", type=str,
                    help='list of input images from which to subtract another image or value')
parser.add_argument("--coordinates", metavar='coordinates', action='store', dest="coordinates",  nargs=1,
                    default="",
                    help='File with coordinates of stars. You can pass one catalogue for all images '+\
                         'or one per image. If no catalogue is passed, SEP will be used to find the stars in the image.')
parser.add_argument("--output", metavar='output', action='store', dest="output",  nargs=1, type=str, required=True,
                    help='Name of the output database. ')
parser.add_argument("--cbox", metavar='cbox', action='store', dest="cbox", type=float, default=10,
                    help='Size of the box within which the centre of the star will be searched. ')
parser.add_argument("--aperture", metavar='aperture', action='store', dest="aperture", type=float, default=3,
                    help='Aperture, by default in FWHM (see --unit below), to perform the photometry. '
                         'The star will be considered to extend out to this radius. ')
parser.add_argument("--annulus", metavar='annnulus', action='store', dest="annulus", type=float, default=6,
                    help='Inner circle of the ring to measure the sky. By default it is in units of FWHM '
                         '(see --unit below). Default value: 6')
parser.add_argument("--dannulus", metavar='dannnulus', action='store', dest="dannulus", type=float, default=2,
                    help='Width of the ring to measure the sky. By default it is in units of FWHM '
                         '(see --unit below). Default value: 2')
parser.add_argument("--unit", metavar='unit', action='store', dest="unit", type=float, default='FWHM',
                    help="Scale of the units in aperture, annulus and dannulus. 'Pix' and 'FWHM' are acceptable,  "
                         'depending on if the aperture, annulus and dannulus should be taken as pixels or as number '
                         'of FWHM. ')




def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)

  output_catalogue = do_photometry(args)
  return output_catalogue

if __name__ == "__main__":
    main()
