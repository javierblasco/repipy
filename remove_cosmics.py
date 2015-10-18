#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import sys
import os
import repipy.utilities as utils
import cosmics_04.cosmics as cosmics
import astropy.io.fits as fits
import repipy.astroim as astroim
""" This routine uses cosmic.py (Malte Tewes, 2010), the python version of LACOS 
    (Van Dokkum, PASP 2001) to remove cosmic rays from an astronomical image. 
     It requires for the cosmic.py module to be in the path, obviously. Some of 
     the options of cosmic rays (e.g. sigfrac, objlim) are not input by the 
     user in this wrapper. """


def remove_cosmics(args):

    for im_name in args.input:
        if args.output != '':
            newfile = args.output
        else:
            newfile = utils.add_suffix_prefix(im_name, prefix = args.prefix, \
                                            suffix = args.suffix )

        # Read the FITS :
        array, header = cosmics.fromfits(im_name)
        # Build the object :
        im = astroim.Astroim(im_name)

        # If object is a flat, bias or flat, do not remove_cosmics
        if im.target.objtype in ["bias", "domeflat", "skyflat", "flat"]:
            continue

        gain = im.primary_header.get(im.primary_header.gaink)
        readnoise = im.primary_header.get(im.primary_header.ccdronk)
        c = cosmics.cosmicsimage(array, gain = float(gain), sigfrac = 0.3, \
                                 readnoise = float(readnoise), objlim = 7.0, \
                                 sigclip = float(args.sigclip))
        # Run the full artillery :
        c.run(maxiter = int(args.maxiter))

        # Write the cleaned image into a new FITS file, conserving the header:
        cosmics.tofits(newfile, c.cleanarray, header)

        # If you want the mask, here it is :
        if args.mask == True:
            mask = utils.add_suffix_prefix(newfile, suffix="-cosmic_mask")
            cosmics.tofits(mask, c.mask, header)

        # And write info to the header:
        im = fits.open(newfile, mode='update')
        hdr = im[0].header
        hdr.add_history("COSMIC RAYS REMOVED:")
        oldname = os.path.split(im_name)[1]
        newname = os.path.split(newfile)[1]
        hdr.add_history(oldname + " --> " + newname)
        hdr.add_history("Parameters used by cosmics.py. Gain=" + str(gain) + \
                        ", sigfrac=0.3, objlim=7.0, sigclip=" + args.sigclip + \
                        ", readnoise=" + str(readnoise))
        im.flush()
        im.close()
    return newfile
    
    
    

# Create parser
parser = argparse.ArgumentParser(description='Remove cosmic rays in images')
parser.add_argument("input",metavar='input', action='store', nargs="+",  \
                    type=str, help='Images to be corrected for cosmic ray hits.')
parser.add_argument("--output", metavar='output', action='store', dest='output', \
                    default='', help='Name of output file.')                    
parser.add_argument("--prefix", metavar="prefix", dest='prefix', action='store', \
                    default='',help='prefix to be added at the beginning of the '+\
                    'image input list to generate the output.')
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store', \
                    default='',help='suffix to be added at the end of the image' +\
                    'input list to generate the output.') 
parser.add_argument("--gain", metavar='gain', action='store', dest='gain',\
                    default='2', help="Gain (either value or keyword from the header) "+\
                    "of the telescope in e-/ADU. Default: 2" )                    
parser.add_argument("--readnoise", metavar="readout noise", action='store', \
                    dest='readnoise', default="5", help="Readout noise of the"+\
                    "telescope (value or header keyword) in e-. Default: 5")
parser.add_argument("--sigclip", metavar="sigmaclip", action='store', dest='sigclip',\
                    default="5.", help=" Sigma clipping factor in order to "+\
                    "distinguish noise from cosmic rays. Higher value are "+\
                    "recommended for high S/N images. Default: 5")
parser.add_argument("--create_mask", action="store_true", dest="mask", \
                    help="Add this keyword if you want to create a mask. The "+\
                    "mask will have -cosmic_mask as suffix... ", default=False)
parser.add_argument("--maxiter", metavar='maxiter', action='store', dest='maxiter',\
                    default='3', help="Maximum number of iterations searching "+\
                    "for cosmic rays. See the documentation of cosmics or LACOS. "+\
                    " Default: 3")



def main(arguments=None):
    # Pass arguments to variable 
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    if args.output == '' and args.prefix == '' and args.suffix == '':  
        sys.exit("Error! Introduce a prefix, a suffix or an output filename. "+\
                 "For help: python remove_cosmics.py -h ")    
    if args.suffix != "":
      args.suffix = args.suffix.strip()
    newfile = remove_cosmics(args)
    return newfile             
if __name__ == "__main__":
    main()
