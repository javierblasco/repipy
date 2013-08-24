#! /usr/bin/env python
###############################################################################
import os
import argparse
import numpy
from scipy.stats import mode as mode
import sys
import astropy.io.fits as fits
import collections
import repipy.utilities as utils
import repipy.arith as arith_images
import repipy.find_keywords as find_keywords

def home_made_median(map_cube, cube):
    """ map_cube indicates the number of elements of cube that are not 
        masked. For the sort that we did before, masked elements are sorted 
        as if they had the very highest values. So, if map[ii,jj] = 5
        we know that the first 5 elements of cube[:,ii,jj] are not masked
        and that the rest are masked. So, the median will be cube[2,ii,jj]
        for the pixel [ii,jj] of the images. Of course, there will be points
        of map that will be even, so the median will be the mean of elements
        N//2 and (N-1)//2, where N is the value of map for the given pixel.
        That is what we call indicesz_above and below. We will average
        the two resulting images. """
    lx, ly = map_cube.shape 
    indicesz_above = (map_cube // 2).flatten() 
    indicesz_below = ((map_cube-1) // 2).flatten()
    indicesx, indicesy = numpy.where(map_cube == map_cube)
    image_above = cube[indicesz_above, indicesx, indicesy].reshape((lx,ly))
    image_below = cube[indicesz_below, indicesx, indicesy].reshape((lx,ly))
    image = (image_above + image_below) / 2.
    return image


def minmax_reject(cube, nlow, nhigh):
    """ Cube is a numpy.ma.array with shape (lz, lx, ly) where lz is the number 
        of combined images, each with sizes (lx,ly). The cube was build from 
        masked images, and its own mask reflects which elements are to be used 
        (valid=0, invalid=1). We need to reject the lowest nlow and the
        highest nhigh values for each pixel. For that, we can use masks and 
        'simple' logic operations.   """
    # First, sort the cube, for each pixel cube[:,ii,jj] is in increasing order
    # with the masked values as the highest values.
    cube.sort(axis=0)    
    lz, lx, ly = cube.shape

    # Using the logic "and" operator, we can mask out the highest value of each
    # pixel quickly by looping backwards. It is not straightforward, but works!
    for ii in range(nhigh): # every time just one more value is masked!
        a = numpy.ones((lx,ly))
        # The pixels of the array a will become zeros as they mask one more 
        # element. Then they will have no effect on the rest of cube. Clever, huh?
        for ii in range(lz-1,-1,-1):
            cube.mask[ii,:,:], a = a, numpy.logical_and(a, cube.mask[ii,:,:])
            if numpy.all(a == 0): #all zeros? no need to continuethe inner loop
               break                 
    
    # Fortunately it is easier with the lowest values. Since the cube is sorted
    # with the non-valid values as the highest ones, we just need to mask the 
    # slice that we do not want from the "bottom" of the array. Actually, 
    # we will cut it away, since it is of no use and would get in the middle 
    # when calculating the median.
    cube = cube[nlow:,:,:]        
    return cube

################################################################################
def cube_images(input_images, mask_key, scales):
    """ From a set of images (input_images), and once divided by their scales, 
        form a cube with all the images. """
    # Initialize the cube and the mask. 
    lx, ly = fits.getdata(input_images[0]).shape
    lz = len(input_images)
    mask = numpy.zeros([lz,lx,ly])  # initialized to zero, all pixels are valid
    cube = numpy.ma.array(mask, mask=mask) 

    # Now read all the images (and masks, if present), into the cube
    for index, image in enumerate(input_images):
        im = utils.read_image_with_mask(image, mask_keyword = mask_key)
        cube.data[index,:,:] = im.data/scales[index]
        cube.mask[index,:,:] = im.mask
    return cube    
    
################################################################################
def compute_scales(input_images, scale, mask_key):
    """ From the list of input images use the central third of the image to 
        calculate the scale necessary to get all the images to the same level
        of flux, using the scale (median, mode, mean, none) given by user. """
    lx, ly = fits.getdata(input_images[0]).shape
    scales = []
    for image in input_images:
        im = utils.read_image_with_mask(image, mask=mask_key)
        centre = im[lx/3:lx*2/3, ly/3:ly*2/3]
        if scale == "median":
            scales.append(numpy.ma.median(centre))
        elif scale == "mean":
            scales.append(numpy.ma.mean(centre))
        elif scale == "none":
            scales.append(1) 
    return scales                                            
                                                                
################################################################################
def run_test(args, list1, filt):
    testdir = os.path.join(args.out_dir, "tests"+str(args.out_pattern))
    if os.path.isdir(testdir) == False:
        os.makedirs(testdir)
    print "PRINTING TEST FILES WITH IMAGES DIVIDED AMONG THEMSELVES: \n "
    for jj, im1 in enumerate(list1):  
        print str(jj+1) + " " + im1  # Identify to keep track of them:
        for kk, im2 in enumerate(list1[jj+1:]): # avoid duplicity
            newname =  str(args.out_pattern) + "_test_" + filt + \
                           "_"+str(jj+1)+"_"+str(kk+jj+2)+".fits"
            filename = os.path.join(testdir, newname)
            if os.path.isfile(filename) == True: 
                os.remove(filename)  # iraf does not overwrite
            arith_images.main(arguments=["--output", filename, "--mask_key",\
                                         args.mask_key, im1, "/", im2]) 		

################################################################################
def build_filter_list(args):
    """ From a set of images read the headers and build a list of the 
        filters present and another one with which filter each image has. """    
    images_filters = []
    for obj in args.in_pattern: # For every image in the input.
        hdr_obj = fits.getheader(obj)
        if args.filterk != "":
            filter_obj = hdr_obj[args.filterk]
        else:               
            keywords = find_keywords.get_keywords(hdr_obj, ["filter"], args)
            filter_obj = hdr_obj[keywords["filter"]]
              
        # Homogeneous name for the filter names!
        filter_obj = utils.homogeneous_filter_name(filter_obj)

        # Add the filter of this image to the list
        images_filters.append(filter_obj)
    filter_list = list(set(images_filters))

    # If alltogether is set, then ignore the different filters and use all the 
    # images as if they had just one filter (e.g., to combine bias images)
    if args.alltogether == True:
        filter_list = ["AllFilters"]
    return filter_list, images_filters

###################################################################################
def combine(args):
    # If output directory does not exist, create it. 
    if os.path.isdir(args.out_dir) == False:
        os.makedirs(args.out_dir)

    # Build a list of filters present and another one with the filter of each image
    filter_list, images_filters = build_filter_list(args)
    	
    # Create a default dictionary for the resulting images
    result = collections.defaultdict(str)    
    
    # For each of the filters present do two things: divide all images by each 
    # other (and write to a test folder) and combine the images (and write)
    for filt in filter_list:	   
        # If user provided a specific filter to be combined, use just that one
        if filt != args.filter:
            continue

        # list of objects with current filter (exception: allfilters is true) 
        list1 = [args.in_pattern[p] for p,f in enumerate(images_filters) if 
                     (f == filt or filt == "AllFilters") ]
 
        # Check that all images have same dimension. If not, exit program
        if not utils.check_dimensions(list1):
            sys.exit("Dimensions of images to combine are different!")

        # Divide each object by all the others and generate files called 
        # test+out_pattern+filter+num1+num2
        if args.notest == False:
            run_test(args, list1, filt)

        # In order to combine, we first need to calculate the scales, so that 
        # they are all at the same flux level
        scales = compute_scales(list1, args.scale, args.mask_key)
        
        # Now we can build a cube with all the images, masks included. 
        cube = cube_images(list1, args.mask_key, scales)
        
        # If nlow != 0 or nhigh!=0 we need to remove the necessary pixels.
        # If user was odd enough to give args.median and args.nlow = args.nhigh
        # skip the minmax rejection.
        if (args.nlow or args.nhigh):
            if not (args.average == "median" and args.nlow == args.nhigh):
               cube = minmax_reject(cube, args.nlow, args.nhigh)
        else:
            cube.sort(axis=0)  # In any case, we need cube to be sorted

        # Finally, average! Remember that the cube is sorted so that
        # cube[0,ii,jj] < cube[1,ii,jj] and that the highest values of all 
        # are the masked elements. We will take advantage of it if the median 
        # is selected, because nowadays the masked median is absurdly slow: 
        # https://github.com/numpy/numpy/issues/1811
        map_cube = numpy.ma.count(cube, axis=0) # number non-masked values per pixel                                
        if args.average == "mean":
            image = numpy.ma.mean(cube, axis=0)
        elif args.average == "median":
            image = home_made_median(map_cube, cube)


        # Image is a masked array, we need to fill in masked values with the 
        # args.fill_val. Also, values with less than args.nmin valid values 
        # should be masked out. 
        image.mask[map_cube < args.nmin] = 1 
        mask = image.mask.astype(numpy.int0) # converto to binary 
        image = image.filled(args.fill_val)
             
        # And save image
        newfile = os.path.join(args.out_dir, 
                               str(args.out_pattern) + "_" + filt + '.fits')
        if args.out_mask != "":
            name_mask = args.out_mask
        else:
            name_mask = newfile + ".msk"
        if os.path.isfile(newfile): 
            os.remove(newfile)  
        if os.path.isfile(name_mask):
            os.remove(name_mask)
        fits.writeto(newfile,image)
        fits.writeto(name_mask, mask)
        result[filt] = newfile  
        if os.path.isfile(newfile) == False:
            sys.exit(newfile + " does not exist")
 
        # Add comments to the headers
        string1 = " - Image built from the combination of the images: "+\
                 " ,".join(list1) + "\n"
        string2 = " combine = " + args.average + ", scale = " + args.scale +\
                 ", reject = minmax, nhigh = " + str(args.nhigh) + ", nlow = "+\
                 str(args.nlow) + "\n"
        utils.add_history_line(newfile, string1 + string2 )  
        utils.add_history_line(name_mask, " - Mask of image: " + newfile)
        if args.mask_key != "":
            utils.header_update_keyword(newfile, args.mask_key, name_mask, 
                                        "Mask for this image")        
       
        # To normalize calculate median and call arith_images to divide by it.
        if args.norm == True:
            im = fits.getdata(newfile)
            lx, ly = im.shape
            median = numpy.median(im[lx/3:lx*2/3,ly/3:ly*2/3])                                 
            msg =  "- NORMALIZED USING MEDIAN VALUE:"                      
            arith_images.main(arguments=["--message", msg, "--output", newfile,
                                         "--mask_key", args.mask_key, "--fill", 
                                         args.fill_val, newfile, "/", str(median)])
    return result

############################################################################

# Create parser
parser = argparse.ArgumentParser(description='Combine images')
# Add necessary arguments to parser
parser.add_argument("in_pattern", metavar='input', action='store', \
                 help='input pattern that the program will look for at the' +\
                 'beginning of the files to combine', nargs='+')
parser.add_argument("-o", metavar="output", dest='out_pattern', \
                   action='store', help='How the output file will be called')
parser.add_argument("--out_dir", metavar="out_dir", dest='out_dir', default='',\
                   action='store', help=' Directory to print output file' +\
                   'Default: same folder as input images.')
parser.add_argument("--average", metavar='average', type=str, default='median', \
                   help='type of average (median, mean) to combine ' +\
                   'the images. Default: median')
parser.add_argument("--scale", metavar='scale', type=str, default='none', \
                   help='scaling function (median, mean, mode, none) to apply' +\
                   'to	the images before combining them. Default: none' )
parser.add_argument("--all_together", action="store_true", dest="alltogether", \
                   default=False, help=' Force all the files to be combined'+ \
                   'together, i.e. do not separate by filter (e.g. for bias)') 
parser.add_argument("--norm", action="store_true", dest="norm", default=False, \
                    help="Normalize resulting image? Default: No")
parser.add_argument("--nlow", metavar="nlow", type=int, dest='nlow', action='store', \
                   default='0', help='Number of low pixels to be rejected by'+ \
                   ' the combining procedure. Default: 0')
parser.add_argument("--nhigh", metavar="nhigh", type=int, dest='nhigh', action='store', \
                   default='1', help='Number of highpixels to be rejected by' +\
                   'the combining procedure. Default: 1')
parser.add_argument("--nmin", metavar="nmin", type=int, dest="nmin", action='store', \
                    default=1, help="Minimum number of images with valid "+\
                    "(i.e. non masked, non rejected) pixels. If the valid pixels "+\
                    "are less than nmin, the fill_val value is used as result, "+\
                    "and the pixel is masked. Default: 1.")
parser.add_argument("--notest", action="store_true", dest="notest", default=False, \
                    help="Do not print the tests. Default: False")
parser.add_argument("--config_file", action="store", dest="config", default="", \
                    help="Config file to introduce the names of the keywords "+\
                    "in the headers of the images. The file should have "+\
                    "a line per keyword, where we give the common name, "+ \
                    "a '=' and the name of the keyword. E.g. for CAHA 2.2m: \n"\
                    "OBJECT = OBJECT \n FILTER = INSFLNAM \n DATE = DATE-OBS \n" +\
                    " At least the keyword FILTER should be included for this " +\
                    "routine to work. " )
parser.add_argument("--filterk", action="store", dest="filterk", default='', \
                    help="Keyword in the header that contains the name of "+\
                    "the filter. Alternatively you can provide the filter "+\
                    "name itself with the argument '--filter'")
parser.add_argument("--filter", action="store", dest="filter", default='', \
                    help="Name of the filter we want to combine. If images with "+\
                         "several filters are present, only this one will be used.")
parser.add_argument("--mask_key", metavar="mask_key", dest='mask_key', \
                    action='store', default="", help=' Keyword in the header ' +\
                    'of the image that contains the name of the mask. The mask '+\
                    'will contain ones (1) in those pixels to be MASKED OUT.')
parser.add_argument("--output_mask", metavar="output_mask", dest='out_mask', \
                    action='store', default="", help=' Name of the output mask. '
                    'If none is provided, the program will create output.fits.msk '+\
                    'where output.fits is the name of the output image. ')
parser.add_argument("--fill_val", metavar="fill_val", dest="fill_val", \
                    action='store', default=0, type=int, help=' Keyword with which to '+\
                    'fill a pixel if that pixel is masked in all the images. '+\
                    'Default: 0')

def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  args = parser.parse_args(arguments)

  # Check that either the filterk or a config file is present:
  if args.filterk == "" and args.config == "":
      sys.exit("ERROR! Either --filterk or a --config_file are "+\
               "necessary to determine the filter of the images. Type python "+\
               "combine.py -h for help.")

  # Manipulate arguments
  if args.out_dir == '':
      args.out_dir = os.path.dirname(args.in_pattern[0])
  args.out_dir = os.path.abspath(args.out_dir)
	
  # Call combine, keep name of the file created
  newfile = combine(args)
  return newfile  

if __name__ == "__main__":
    main()
