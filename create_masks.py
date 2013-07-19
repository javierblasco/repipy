import pyfits
import argparse
import numpy
#import copy
from numpy import ma as ma
import os
import sys

""" Program to mask out areas of the image according to certain criteria, to prevent it from participating in 
	aspects like making statistics of the image. 
    We mask out:
	- Values below minval and above maxval. Default: (-100,55000)
     - Values outside the round FoV of the CAFOs 2.2 instrument, when the data are from there 
""" 

def mask_CAHA2(shp, rmax):
    """
    Pixels further from the centre than radius rmax are masked out (True), others are not (False). 
    """	
    # Centre of the image
    xcent = shp[0]/2
    ycent = shp[1]/2
    x = numpy.arange(shp[0])
    y = numpy.arange(shp[1])    
    
    # Distances of every pixel to the centre
    r = numpy.sqrt( (x[:,None]-xcent)**2 + (y[None,:]-ycent)**2 )
    
    # Build a mask combining the new and the old with a logical "or"
    mask = ma.getmask(ma.masked_greater(r, rmax))
    return mask 
    
def mask_image(args):
    """ 
    This program reads the images indicated by the user and creates a mask with the mask being False:
        - for pixels outside the boundaries (minval,maxval). Default: (-100,55000).
        - for pixels outside the round FOV within the square image (data from CAHA 2.2)	  
    """

    for name1 in args.input:
        # Convert to absolute path, in case the user does, i.e. "./../file"
        filename = os.path.abspath(name1)
        
        # Read the image and the header
        #print "Image: "+filename
        image = pyfits.open(filename, mode='update')
        im = image[0].data
        hdr = image[0].header

        # Minmax masking. If there is no values in the mask (no True where found) 
        # build a mask of Falses.
        im2 = ma.masked_outside(im, int(args.minval), int(args.maxval))
        mask1=ma.getmask(im2)	
        if im2.mask.shape == (): im2 = ma.array(im, mask=False)
        mask1 = ma.getmask(im2)
        #print "    Masked out for values <"+args.minval+" or > "+args.maxval		

        # CAHA masking. If telescope is CAHA 2.2 then the image is circular, 
        # mask out anything close to the border (R~800)
        mask2 = None
        if hdr.has_key("Instrume") == True and hdr["Instrume"] == 'CAFOS 2.2':
                mask2 = mask_CAHA2(im.shape, 760)   # Actually 800
                print "    Masked out values outside circular FoV (CAFOS 2.2)"
        if mask2 == None: mask2 = numpy.ma.getmask(ma.array(im, mask=False))

        # Combine both with logical "or" operator
 
        im2.mask = ma.mask_or(mask1,mask2)
        if os.path.isfile(filename+".msk") == True:
            os.remove(filename+".msk")
        im2.mask.dump(filename+".msk")
        #print "Created mask: "+filename+".msk"		

        # Modify and update header to keep history of things
        name = os.path.split(filename)[1]       
        if str(hdr.get_history()).count("-Mask created:") == 0:
            hdr.add_history('-Mask created: (see key "mask") ')
        hdr.update("mask",name + ".msk", "mask image" )
        image.flush()  
        image.close()

########################################################################################################################
# Create parser
parser = argparse.ArgumentParser(description='Create a mask from an image.')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of" +\
                    " input images for which to create a mask', nargs='+')
parser.add_argument("--maxval", metavar='maxval', action='store', default='55000', \
                    help='Maximum value allowed, mask 	above it. Default: 55000.', \
                    nargs='?')
parser.add_argument("--minval", metavar='minval', action='store', default='-100',\
                    help='Minimum value allowed, mask below it. Default: -100 ', \
                    nargs='?')

def main(arguments=None):
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    mask_image(args)

if __name__ == "__main__":
    main()
	
