#! /usr/bin/env python
import os
import re
import pyfits
from pyraf.iraf import imarith as imarith
import sys
import argparse
import numpy
import astropy.io.fits as fits
import repipy.utilities as utils
""" Wrapper for imarith using pyraf. The program will perform an operation 
    between images"""
   
    
def arith(args):
    output_list = []        
    for image in args.input1:        
        # Read inputs as masked arrays or numbers (input2 might be a scalar)  
        im1 = utils.read_image_with_mask(image, args.mask_key)        
        try: # if it is a number        
            im2 = float(args.input2[0])
        except ValueError:
            im2 = utils.read_image_with_mask(args.input2[0], args.mask_key)        
            
        # Case of mean or median for Input2 
        if args.median == True:
            im2 = numpy.ma.median(im2)
        elif args.mean == True:
            im2 = numpy.ma.mean(im2)

        # Do the actual operation
        import operator
        operations = {"+":operator.add, 
                      "-":operator.sub, 
                      "*":operator.mul,
                      "/":operator.div,
                      "**":operator.pow}
        oper = operations[args.operation[0]]
        result = oper(im1, im2)  # Actual operation
                    
        # If output exists use it, otherwise use the input. Prefix/suffix might
        # modify things in next step.
        if args.output != '': 
            outpt = os.path.abspath(args.output)
        else:
            outpt = os.path.abspath(image)
            
        # Separate (path, file) and (file root, extensions). Then build new name. 
        outdir, outfile = os.path.split(outpt)   
        outfile_root, outfile_ext = re.match(r'(^.*?)(\..*)', outfile).groups()   
        outpt = os.path.join(outdir, args.prefix + outfile_root + args.suffix + 
                             outfile_ext)                         

        # Now name for the mask, if the name exists, use it, otherwise build it.
        if args.mask_name != "":
            mask_name = args.mask_name
        else:
            mask_name = outpt + ".msk"

        # Prepare a header starting from input1 header 
        hdr_im = fits.getheader(image)
        name2 = os.path.split(args.input2[0])[1]
        if args.median:  # if single number because median used
            name2 = "median(" + name2 + ") (" + str(im2) + ")"
        elif args.mean:  # if single number because mean used
            name2 = "mean(" + name2 + ") (" + str(im2) + ")"
            
        # Write a HISTORY element to the header    
        hdr_im.add_history(" - Operation performed: " + outfile + " " + 
                           args.operation[0] + " " + name2)
        hdr_mask = fits.PrimaryHDU(result.mask.astype(int)).header
        hdr_mask.add_history("- Mask corresponding to image: " + outpt)
        
        # Now save the resulting image and mask
        if os.path.isfile(outpt):
            os.remove(outpt)
        if os.path.isfile(mask_name):
            os.remove(mask_name)
        fits.writeto(outpt, result.filled(args.fill), header=hdr_im) 
        fits.writeto(mask_name, result.mask.astype(numpy.int0), header=hdr_mask)
        output_list.append(outpt)
    return output_list

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
                    default='', type=str, help='prefix to be added at the '+\
                    'beginning of the image input list to generate the outputs.',
                    nargs=1)
parser.add_argument("--suffix", metavar="suffix", dest='suffix', action='store',\
                    default='', type=str, help='suffix to be added at the end '+\
                    'of the image input list to generate the outputs. There '+\
                    'is a peculiarity with argparse: if you pass, e.g., "-c" to '+\
                    '--suffix, the program will understand that you want to '+\
                    'call the code with the flag -c, which does not exist. This '+\
                    'does not raise an error, but just stops execution, which is '+\
                    'quite annoying. One way around it is " -c" (notice the '+\
                    'space, since within the string is stripped.', nargs=1)
parser.add_argument("--message", metavar="hdr_message", dest='hdr_message', \
                    action='store', default="", help=' Message to be added to ' +\
                    'the header via HISTORY. For example: bias subtracted.')
parser.add_argument("--mask_key", metavar="mask_key", dest='mask_key', \
                    action='store', default="", help=' Keyword in the header ' +\
                    'of the image that contains the name of the mask. The mask '+\
                    'will contain ones (1) in those pixels to be MASKED OUT.')
parser.add_argument("--mask_name", metavar="mask_name", dest='mask_name', \
                    action='store', default="", help=' Name for the resulting ' +\
                    'mask. If the pixel of any of the images is masked out, the '+\
                    'resulting mask will obviously be masked out as well.')         
parser.add_argument("--fill", metavar="fill", dest="fill", action="store",\
                    default=0, type=int, help=' Value to fill masked pixels when '+\
                    'creating the final image. Default: 0')
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

  # In order to allow "  -c" to be able to use the hyphen as suffix.
  if args.suffix != "":
      args.suffix = args.suffix[0].strip()
  if args.prefix != "":
      args.prefix = (args.prefix[0]).strip()  
  
  # Detecting errors
  if args.output == '' and args.prefix == '' and args.suffix == '' and args.overwrite == False:  
      sys.exit("Error! Introduce a prefif, a suffix, the --overwrite option or the --output option. \
			  For help: python arith.py -h ") 
  newname = arith(args)  
  return newname    
     
if __name__ == "__main__":
    main()