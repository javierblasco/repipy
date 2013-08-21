# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 09:02:07 2013

@author: blasco
"""
import numpy
import sys
import argparse
import repipy.utilities as utils
import astropy.io.fits as fits

def apply_median_filter(im, side):
    """ Filter a masked np.ma array with a median filter. """
    lx, ly = im.shape
    mask = numpy.ones((lx,ly),dtype=numpy.int0)
    data = numpy.asarray(mask, dtype=numpy.float)
    filt_im = numpy.ma.array(data, mask=mask)
    for ii in range(lx):
      for jj in range(ly):
        if not im.mask[ii,jj]:  # If pixel is not masked
          minx, maxx = max([ii-side/2,0]), min([ii+side/2+1,lx])
          miny, maxy = max([jj-side/2,0]), min([jj+side/2+1,ly])
          whr = numpy.where(im.mask[minx:maxx,miny:maxy] == 0)
          filt_im.data[ii,jj] = numpy.median(im.data[minx:maxx,miny:maxy][whr])
          filt_im.mask[ii,jj] = 0
    return filt_im            
            
def filter_image(args):
    """ Routine that uses a median filter on a list of masked images. """
    output_list = []
    for image in args.input:
        im = utils.read_image_with_mask(image, mask_keyword=args.mask_key)
        filt_im = apply_median_filter(im,args.side)        
        filt_im = filt_im.filled(args.fill_val)
        if args.output == "":
            output = utils.add_suffix_prefix(image, suffix='_filtered')
        else :
            output = args.output
        fits.writeto(output, filt_im)
        output_list.append(output)
    return output_list
    
############################################################################
# Create parser
parser = argparse.ArgumentParser(description='Combine images')
# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'input images to be median filtered.', \
                    nargs="+", type=str)
parser.add_argument("--mask_key", metavar="mask_key", dest='mask_key', \
                    action='store', default="", help=' Keyword in the header ' +\
                    'of the image that contains the name of the mask. The mask '+\
                    'will contain ones (1) in those pixels to be MASKED OUT.')
parser.add_argument("--side", metavar="side", type=int, dest="side", 
                    action='store', required=True, help=" Size of the box to "+\
                    "be used to filter the image. Mandatory argument.")
parser.add_argument("--output", metavar='output', dest='output', action='store',
                   default='', help='output image in which to save the result.'+\
                    ' If not provided the suffix _filtered will be added to '+\
                    'the input.')
parser.add_argument("--fill_val", metavar="fill_val", dest="fill_val", \
                    action='store', default=0, type=int, help=' Keyword with '+\
                    'which to fill masked pixels. Default: 0')                   

############################################################################                   

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  args = parser.parse_args(arguments)

  print args.input, type(args.input)

  # Call combine, keep name of the file created
  newfile = filter_image(args)
  return newfile  

if __name__ == "__main__":
    main()
                  