#! /usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import re
import astropy.io.fits as fits
import sys
import argparse
import numpy
import repipy.utilities as utils
import operator

""" Wrapper for imarith using pyraf. The program will perform an operation 
    between images"""
   
    
def arith(args):
    output_list = []        
    for operand1, operand2 in zip(args.input1, args.input2):
        # Read inputs as masked arrays or numbers (input2 might be a scalar)  
        im1 = utils.read_image_with_mask(operand1, args.mask_key)
        try: # if second operand is a number we still need it in a masked array
            value2 = numpy.ma.array([float(operand2)],mask=[0])
        except (ValueError,TypeError):
            value2 = utils.read_image_with_mask(operand2, args.mask_key)

        # If --keyword was passed, read from the header of the image the value and use it, instead of the full image
        if args.keyword:
            value2 = numpy.ma.array( fits.getval(operand2, args.keyword) )

        # Do the actual operation. Result is a masked array which masks 
        # any pixel if it is masked either in im1 or in im2.
        operations = {"+":operator.add, 
                      "-":operator.sub, 
                      "*":operator.mul,
                      "/":operator.div,
                      "**":operator.pow}

        # Case of mean or median for Input2 
        if args.median == True:
            value2 = numpy.ma.median(value2)
        elif args.mean == True:
            value2 = numpy.ma.mean(value2)

        oper = operations[args.operation[0]]
        result = numpy.zeros_like(im1)

        result.data[:] = oper(im1.data, value2.data)  # Actual operation of images
        result.mask[:] = im1.mask | value2.mask       # If any is masked, result is masked as well


        # If args.fill_val is present, use it
        if args.fill_val:
            result.data[:] = result.filled(float(args.fill_val))
                    
        # If output exists use it, otherwise use the input. Prefix/suffix might
        # modify things in next step.
        if args.output != '': 
            outpt = os.path.abspath(args.output)
        else:
            outpt = os.path.abspath(operand1)
            
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
        hdr_im = fits.getheader(operand1)
        name2 = os.path.split(operand2)[1]
        if args.median:  # if single number because median used
            name2 = "median(" + name2 + ") (" + str(value2) + ")"
        elif args.mean:  # if single number because mean used
            name2 = "mean(" + name2 + ") (" + str(value) + ")"
            
        # Write a HISTORY element to the header    
        #hdr_im.add_history(" - Operation performed: " + image1 + " " +
        #                   args.operation[0] + " " + name2)
        try:
            hdr_mask = fits.PrimaryHDU(result.mask.astype(int)).header
        except:
            pass
           
        hdr_mask.add_history("- Mask corresponding to image: " + outpt)
        
        # Now save the resulting image and mask
        if os.path.isfile(outpt):
            os.remove(outpt)
        if os.path.isfile(mask_name):
            os.remove(mask_name)
        fits.writeto(outpt, result.data, header=hdr_im) 
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
parser.add_argument("input2",metavar='input2', action='store', nargs="+", type=str, \
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
parser.add_argument("--fill_val", metavar="fill_val", dest="fill_val", 
                    action="store", default='', help=' If present, ' +\
                    'this value will be used to fill masked pixels in the ' +\
                    'final image. By default no filling is used, and the '+\
                    'original value is used')
parser.add_argument("--overwrite", action="store_true", dest="overwrite", \
                    default=False, help="Allows you to overwrite the original image.")
parser.add_argument("--mean", action="store_true", dest="mean", default=False, \
                    help='Input2 is not used as an image, but the mean is ' +\
                    'calculated and used instead.')
parser.add_argument("--median", action="store_true", dest="median", default=False, \
                    help='Input2 is not used as an image, but the median is '+\
                    'calculated and used instead.')
parser.add_argument("--keyword", action="store", dest="keyword", default='', \
                    help=' If operand2 is a fits image and this argument is passed with the name of a keyword, '+\
                         'do the operation using the value in the header for the given keyword, instead of operating '+\
                         'using  the images. For example, if you pass "arith.py --keyword exptime im1.fits / im1.fits '+\
                         'you are telling arith.py to divide the image im1.fits by the value of the keyword exptime '+\
                         'in its header (thus normalizing to 1sec exposure time).')


	

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  
  args = parser.parse_args(arguments)

  all_together = args.input1 + args.operation + args.input2
  valid_ops = ["+", "-", "*", "/"]
  for ii, elem in enumerate(all_together):
      if elem in valid_ops:
          op_index = ii
  args.input1 = all_together[:op_index]
  args.operation = all_together[op_index]
  args.input2 = all_together[op_index+1:]

  # input2 can be a single value/image or as many as input1
  if len(args.input1) == len(args.input2):
      pass
  elif len(args.input2) == 1:
      args.input2 *= len(args.input1)
  else:
      sys.exit("ERROR! Input2 must be either 1 element or as many as input1")


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
  
  if len(newname) == 1:   # If just one element, send back an element, not a list
      newname = newname[0]
  return newname    
     
if __name__ == "__main__":
    main()
