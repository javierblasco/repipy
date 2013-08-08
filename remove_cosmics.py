#! /usr/bin/env python
import argparse
import sys
import os
import repipy.utilities as utils
import cosmics_04.cosmics as cosmics
import pyfits
""" This routine uses cosmic.py (Malte Tewes, 2010), the python version of LACOS 
    (Van Dokkum, PASP 2001) to remove cosmic rays from an astronomical image. 
     It requires for the cosmic.py module to be in the path, obviously. Some of 
     the options of cosmic rays (e.g. sigfrac, objlim) are not input by the 
     user in this wrapper. """


def remove_cosmics(args):
    if args.output != '':
        newfile = args.output
    else:
        newfile = utils.add_suffix_prefix(args.input[0], prefix = args.prefix, \
                                        suffix = args.suffix )
    
    # Read the FITS :
    array, header = cosmics.fromfits(args.input[0])
    # Build the object :
    c = cosmics.cosmicsimage(array, gain = float(args.gain), sigfrac = 0.3, \
                             readnoise = float(args.readnoise), objlim = 5.0, \
                             sigclip = float(args.sigclip))
    # Run the full artillery :
    c.run(maxiter = int(args.maxiter))
    
    # Write the cleaned image into a new FITS file, conserving the header:
    cosmics.tofits(newfile, c.cleanarray, header)
    
    # If you want the mask, here it is :
    if args.mask == True:
        mask = utils.add_suffix_prefix(newfile, prefix="cosmic_mask")    
        cosmics.tofits(mask, c.mask, header)
                          
    # And write info to the header:
    im = pyfits.open(newfile, mode='update')
    hdr = im[0].header
    hdr.add_history("COSMIC RAYS REMOVED:")
    oldname = os.path.split(args.input[0])[1]
    newname = os.path.split(newfile)[1]
    hdr.add_history(oldname + " --> " + newname)
    hdr.add_history("Parameters used by cosmics.py. Gain=" + args.gain + \
                    ", sigfrac=0.3, objlim=5.0, sigclip=" + args.sigclip + \
                    ", readnoise=" + args.readnoise)   
    im.flush()
    im.close()                
    return newfile
    
    
    

# Create parser
parser = argparse.ArgumentParser(description='Remove cosmic rays in images')
parser.add_argument("input",metavar='input', action='store', nargs=1,  \
                    help='Image to be corrected for cosmic ray hits.')
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
                    "recommended for high S/N images.")                  
parser.add_argument("--create_mask", action="store_true", dest="mask", \
                    help="Add this keyword if you want to create a mask. The "+\
                    "mask will be called cosmic_mask-... ", default=False)
parser.add_argument("--maxiter", metavar='maxiter', action='store', dest='maxiter',\
                    default='3', help="Maximum number of iterations searching "+\
                    "for cosmic rays. See the documentation of cosmics or LACOS" )                    



def main(arguments=None):
    # Pass arguments to variable 
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    if args.output == '' and args.prefix == '' and args.suffix == '':  
        sys.exit("Error! Introduce a prefix, a suffix or an output filename. "+\
                 "For help: python remove_cosmics.py -h ")     
    newfile = remove_cosmics(args)
    return newfile             
if __name__ == "__main__":
    main()
