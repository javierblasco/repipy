#! /usr/bin/env python
###############################################################################
import pyfits 
import os
import argparse
import numpy
import sys
import tempfile
import StringIO
import ConfigParser
import jbh_utilities as jbh
import reduction_pipeline.arith as arith_images
import reduction_pipeline.find_keywords as find_keywords
import collections
from pyraf.iraf import imcombine as imcombine

###################################################################################
def combine(args):
    # If output directory does not exist, create it. 
    if os.path.isdir(args.out_dir) == False:
        os.makedirs(args.out_dir)

    # Keep track of the filters each image has, build a set with them.
    images_filters = []
    for obj in args.in_pattern: # For every image in the input.
        hdr_obj = pyfits.getheader(obj)
        if args.filterk != "":
            filter_obj = hdr_obj[args.filterk]
        else:               
            keywords = find_keywords.get_keywords(hdr_obj, ["filter"], args)
            filter_obj = hdr_obj[keywords["filter"]]
              
        # Homogeneous name for the filter names!
        filter_obj = jbh.homogeneous_filter_name(filter_obj)

        # Add the filter of this image to the list
        images_filters.append(filter_obj)
    filter_list = list(set(images_filters))
    
    # If alltogether is set, then ignore the different filters and use all the 
    # images as if they had just one filter (e.g., to combine bias images)
    if args.alltogether == True:
        filter_list = ["AllFilters"]
	
    # If user provided a specific filter to be combined, use just that one
    if args.filter != "":
        filter_list = args.filter
 
    # Create a default dictionary for the resulting images
    result = collections.defaultdict(str)    
    
    # For each of the filters present do two things: divide all images by each 
    # other (and write to a test folder) and combine the images (and write)
    for filt in filter_list:	   
        # list of objects with current filter (exception: allfilters is true) 
        list1 = [args.in_pattern[p] for p,f in enumerate(images_filters) if 
                     (f == filt or filt == "AllFilters") ]

        # Divide each object by all the others and generate files called 
        # test+out_pattern+filter+num1+num2
        if args.notest == False:
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
                    arith_images.main(arguments=["--output", filename, im1, "/",
                                                 im2]) 		

        # Read where the centre is to make statistics in a box around it
        hdr = pyfits.getheader(list1[0])
        centrex = int(hdr["NAXIS1"]/2)
        centrey = int(hdr["NAXIS2"]/2)
        lstatx = int(centrex/3 )
        lstaty = int(centrey/3 )
        statsection = "["+str(centrex-lstatx)+":"+str(centrex+lstatx)+","+\
			       str(centrey-lstaty)+":"+str(centrey+lstaty)+"]"

        # Create a temporary file with the names of the files to be combined
        temp_file = open("temporary.lis", "w")
        for name in list1:
            temp_file.write(name + "\n")
        temp_file.close()    
                        
        # Let's combine!
        output_name = str(args.out_pattern) + "_" + filt + '.fits'
        newfile = os.path.join(args.out_dir, output_name)
        if os.path.isfile(newfile) == True: 
            os.remove(newfile)  # iraf does not overwrite files 
        imcombine(input="@temporary.lis", output=newfile, combine=args.average, \
                           scale=args.scale, reject="minmax", mclip="yes", \
                           statsec=statsection, nhigh=args.nhigh, nlow=args.nlow)
        os.remove("temporary.lis")        
        result[filt] = newfile
  
        if os.path.isfile(newfile) == False:
            sys.exit(newfile + " no existe")
 
       # And read the image again to add comments to the header
        newimage = pyfits.open(newfile, mode="update")
        hdr = newimage[0].header
        string = " ,".join(list1)
        hdr.add_history(" - Image built from the combination of the images: "+string)
        hdr.add_history(" combine = " + args.average + ", scale = " + args.scale + \
                        ", reject = 'minmax', mclip='yes', " +\
                        "statsec = " + statsection + ", nhigh=" + args.nhigh +\
                        ", nlow = " + args.nlow)
        newimage.flush()
        newimage.close()

        # To normalize calculate median and call arith_images to divide by it.
        if args.norm == True:
            im = pyfits.getdata(newfile)
            median = numpy.median(im[centrex-lstatx:centrex+lstatx,\
                                  centrey-lstaty:centrey+lstaty])                                 
            msg =  "- NORMALIZED USING MEDIAN VALUE:"                      
            arith_images.main(arguments=["--message", msg, "--output", newfile,
                                         newfile, "/", str(median)])
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
                   help='type of average (median, mean, mode, none) to combine' +\
                   'the images. Default: median')
parser.add_argument("--scale", metavar='scale', type=str, default='none', \
                   help='scaling function (median, mean, mode, none) to apply' +\
                   'to	the images before combining them. Default: none' )
parser.add_argument("--all_together", action="store_true", dest="alltogether", \
                   default=False, help=' Force all the files to be combined'+ \
                   'together, i.e. do not separate by filter (e.g. for bias') 
parser.add_argument("--norm", action="store_true", dest="norm", default=False, \
                    help="Normalize resulting image? Default: No")
parser.add_argument("--nlow", metavar="nlow", dest='nlow', action='store', \
                   default='0', help='Number of low pixels to be rejected by'+ \
                   ' the combining procedure. Default: 0')
parser.add_argument("--nhigh", metavar="nhigh", dest='nhigh', action='store', \
                   default='1', help='Number of highpixels to be rejected by' +\
                   'the combining procedure. Default: 1')
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