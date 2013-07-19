#! /usr/bin/env python
import os
import pyfits
from pyraf.iraf import imarith as imarith
import sys
import argparse
import numpy
""" Wrapper for imarith using pyraf. The program will perform an operation 
    between images"""


def arith(args):        
    output = []
    for image in args.input1:
        # Input2 or some mean or median of input2
        if args.median == True:
            im2 = pyfits.getdata(args.input2[0])
            args.input2[0] = str(numpy.median(im2))
        elif args.mean == True and (args.input2[0]).isalpha() == True:
            im2 = pyfits.getdata(args.input2[0])
            args.input2[0] = str(numpy.mean(im2))
                    
        # If output exists use it, otherwise use the input.
        if args.output != '': 
            outpt = args.output
        else:
            outpt = image
        (outdir, outfile) = os.path.split(outpt)   
        
        # Separate the .fits termination.
        outfile_root = (outfile.split(".fits"))[0]
        
        # Now construct the output file including the prefix/suffix if present
        if args.prefix != '' and args.suffix != '':
            outpt = os.path.join(outdir, args.prefix + "-" + outfile_root + "-" + \
                                 args.suffix+".fits")                             
        elif args.prefix != '':
            outpt = os.path.join(outdir, args.prefix + "-" + outfile_root + ".fits")
        elif args.suffix != '':
            outpt = os.path.join(outdir, outfile_root + "-" + args.suffix + ".fits")
        elif args.output == "": 
            outpt = os.path.join(outdir, outfile_root + "output.fits")
            
        # Iraf does not overwrite things, if input not equal to output and output 
        # exists, remove it. 
        if os.path.isfile(outpt) and outpt != image:    
            os.remove(outpt)
           
        # Actual operation
        imarith(image, args.operation[0], args.input2[0], outpt)
        output.append(outpt)        
        
        # Read the output file and update its history
        im = pyfits.open(outpt, mode='update')
        hdr = im[0].header
        name1 = os.path.split(image)[1]
        name2 = os.path.split(args.input2[0])[1]
        if args.hdr_message != None:
            hdr.add_history(args.hdr_message)
        hdr.add_history(" - Operation performed: "+ name1 + " " + args.operation[0] \
                        + " " + name2)
        im.flush()
        im.close()
    return output

########################################################################################################################


# Create parser
parser = argparse.ArgumentParser(description='Arithmetic operations on images')

# Add necessary arguments to parser
parser.add_argument("input1", metavar='input1', action='store', help='list of ' +\
                    'input images from which to subtract another image or value', \
                    nargs="+", type=str)
parser.add_argument("operation", metavar='operation', action='store', type=str, 
		   help='type of operation (+,-,*,/) to be done', nargs=1)
parser.add_argument("input2",metavar='input2', action='store', nargs=1,  \
                    help='image (or value) with which to perform the operation')
parser.add_argument("--output", metavar='output', dest='output', action='store', \
                   default='', help='output image in which to save the result.' +\
                   'If not stated, then the --prefix or --suffix must be present.')
parser.add_argument("--prefix", metavar="prefix", dest='prefix', action='store', \
                    default='',help='prefix to be added at the beginning of the '+\
                    'image input list to generate the outputs.')
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store', \
                    default='',help='suffix to be added at the end of the image' +\
                    'input list to generate the outputs.')
parser.add_argument("--message", metavar="hdr_message", dest='hdr_message', \
                    action='store', default="", help=' Message to be added to ' +\
                    'the header via HISTORY. For example: bias subtracted.')
parser.add_argument("--overwrite", action="store_true", dest="overwrite", \
                    default=False, help="Allows you to overwrite the original image.")
parser.add_argument("--mean", action="store_true", dest="mean", default=False, \
                    help='Input2 is not used as an image, but the mean is ' +\
                    'calculated and used instead.')
parser.add_argument("--median", action="store_true", dest="median", default=False, \
                    help='Input2 is not used as an image, but the median is '+\
                    'calculated and used instead.')
	

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  args = parser.parse_args(arguments)
  
  # Detecting errors
  if args.output == '' and args.prefix == '' and args.suffix == '' and args.overwrite == False:  
      sys.exit("Error! Introduce a prefif, a suffix, the --overwrite option or the --output option. \
			  For help: python arith.py -h ") 
  newname = arith(args)  
  return newname    
     
if __name__ == "__main__":
    main()