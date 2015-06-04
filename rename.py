#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
@author: javier blasco herrera
@e-mail: blasco@iaa.es  javier.blasco.herrera@gmail.com 

Routine to rename files with a single homogeneous structure that actually gives
information about its content. No more 13125.fits in your data!

It also organizes your data. For example, flats in a folder, bias in another one, 
standards (if recoznized) in another... It will also return a dictionary containing
three numpy arrays. The array "filename" will obviously contain the names of the 
files that have been found. The array "type" will determine if the images are 
bias, skyflats, standards, domeflats, galaxies, ... The array object will contain
the name of the object itself, for example CIG1019 in the case of galaxies or
sao110 for standards.  

"""

#############################################################################################################################
#    TO BE FIXED
#
#   - The date is being read assuming all the observations come from the same 
#     night. Could lead to mistakes? Other options?
##############################################################################################################################

import os, glob, shutil
import collections 
import astropy.io.fits as fits
import sys
import argparse
import ConfigParser
import dateutil.parser
import datetime
import StringIO
import re
import numpy
import repipy.utilities as utils
import repipy.astroim as astroim

#############################################################################
def is_a_standar(object_name):
    """ Determines if an object is one of the standard stars included in the list.
        If so, the object is from now on identified as a standard star image, and 
        treated like it. 
        
        Input: 
            item: string containing the name of an object. 
        Output: 
            boolean return. True if item is found to be a standard, False 
                            otherwise.
    """
    # Common names of the standard fields 
    standards_list = {"bd+25":"bd+25",
                      "feige110":"feige110",
                      "feige110gunnr":"feige110",
                      "feige110ha":"feige110",
                      "f110":"feige110",
                      "grw+73+8031":"grw73",
                      "grw+73+8031gunnr":"grw73",
                      "grw+73+8031ha":"grw73",
                      "pg1708+602":"pg1708",
                      "pg1708+602gunnr":"pg1708",
                      "pg1708+602ha":"pg1708",
                      "bd+28":"bd+28", 
                      "bd+35":"bd+35", 
                      "kop27":"kop27", 
                      "grw73":"grw73", 
                      "pg1708":"pg1708", 
                      "f110":"f110", 
                      "pg170":"pg170", 
                      "wolf1346":"wolf1346", 
                      "sa113":"sa113", 
                      "hz15":"hz15", 
                      "feige34":"feige34",
                      "pg1633":"pg1633", 
                      "sa110":"sa110", 
                      "pg1323":"pg1323", 
                      "pg16":"pg16", 
                      "sao110":"sao110", 
                      "hd84937":"hd84937",
                      "f66":"f66", 
                      "h84":"hd84937", "hd84":"hd84937", "hd84937":"hd84937"}
    return standards_list.get(object_name.lower())

##############################################################################
def is_a_cluster(object_name):
    """ Sometimes we take images of clusters to use as extra flat-fielding sources
        Determines if an object is one of the standard stars included in the list.
        If so, the object is from now on identified as a standard star image, and 
        treated like it. 
        
        Input: 
            item: string containing the name of an object. 
        Output: 
            boolean return. True if item is found to be a standard, False 
                            otherwise.
    """
    # Common names of the standard fields 
    cluster_list = ["m29"]
    answer = False
    for cluster in cluster_list:
        if ((object_name).lower()).count(cluster) != 0 :
            answer = True
    return answer

###########################################################################
def distinguish_type(object_name):
    """ 
      Routine to find out the type of object from the name in the header. It 
      distinguishes between bias, skyflats, domeflats, blanks, indeterminateded 
      flats, cig galaxies, standard stars and clusters. If the galaxy is a CIG
      galaxy, the name will use four digits for the number, not less. 
      
      Standard fields and clusters (sometimes used for flatfielding purposes) are
      recognized if they are in the local lists above. It will return object name 
      and object type (the groups mentioned above). Examples of returns are: 
          "standard", "pg104" 
          "galaxy", "cig0139". 
          "bias", "bias"
          "skyflat", "skyflat"
          "unknown", "unknown"
    """
    if object_name.lower().count("bias") != 0: 
        object_name = "bias"
        object_type = "bias"
    elif object_name.count("flat") != 0 and object_name.count("sky") != 0 : 
        object_name = "skyflat"     
        object_type = "skyflats"
    elif object_name.count("flat") != 0 and object_name.count("dome") != 0 : 
        object_name = "domeflat"
        object_type = "domeflats"
    elif object_name.lower().count("blank") != 0:
        object_name = "blank"
        object_type = "blanks"
    elif object_name.lower().count("flat") != 0: # flat but no idea which type
        object_name = "flat"
        object_type = "flats"
    elif is_a_standar(object_name.lower()) is not None:
        object_name = is_a_standar(object_name.lower())
        object_type = "standards"
    elif is_a_cluster(object_name.lower()) == True:
        object_name = object_name.lower()
        object_type = "clusters"
    elif re.match(r'(cig|c)\s?(\d{1,4}).*', object_name) is not None:
        # cig usually folloed by a number
        cig, number = re.match(r'(cig|c)\s?(\d{1,4}).*', object_name).groups()[0:2]
        object_name = "cig" + "{0:04d}".format(int(number))
        object_type = "cig"
    # If "cig" is missing, but there is a number, add "cig":        
    elif object_name.isdigit() == True :
        object_name = "{0:04d}".format(int(object_name))
        object_name = "cig" + object_name
        object_type = "cig"
    else:
        object_type = "unknown"

    return object_name, object_type
    

def backup(args, fits_list):
    """ Save a list of files into a

    :param args:
    :return:
    """
    backup_dir = os.path.join(args.in_dir, "original_data")
    utils.if_dir_not_exists_create(backup_dir)
    for im in fits_list:
        shutil.copy(im, backup_dir)


def sort_by_date(fits_list):
    # Sort all images chronologically
    im_datetimes = []
    for im_name in fits_list:
        im = astroim.Astroim(im_name)
        im_date = im.header.get(im.header.datek)
        im_time = im.header.get(im.header.timek)
        if im_time:  # if not None
            im_date = "{0}T{1}".format(im_date, im_time)
        im_datetimes.append(dateutil.parser.parse(im_date))
    indices = numpy.argsort(im_datetimes)
    sorted_list = [fits_list[ii] for ii in indices]
    return sorted_list




###########################################################################
def rename(args):
    # List of fit and fits images in the directory
    fits_list1 = glob.glob(os.path.join(args.in_dir, args.in_pattern + "*.fits"))
    fits_list2 = glob.glob(os.path.join(args.in_dir, args.in_pattern + "*.fit"))
    fits_list = fits_list1+fits_list2

    # Sort all images chronologically
    fits_list = sort_by_date(fits_list)

    # If --copy was selected, copy all those files into a directory called original_data
    if args.copy == True:
        backup(args, fits_list)

    # The output of the whole code will be this dictionary, in which the images
    # are sorted in groups (bias, skyflats, domeflats, cigXXXX, ...)
    empty_array = numpy.asarray([], dtype=object)
    final_dict = {"filename":numpy.asarray([], dtype="S150"),
                  "type":empty_array,
                  "objname":empty_array,
                  "time":empty_array}
    
    # Run through all images
    for im_name in fits_list:
        # Read image and header, extract name of object and filter.
        im = astroim.Astroim(im_name)
        object_name = re.sub('[\s\-_\(\)]', "", im.target.objname.lower())
        object_type = im.target.objtype
        object_filter = im.filter.__str__()
        object_date = dateutil.parser.parse(im.header.get(im.header.datek)).date()
        object_date = re.sub('[\s\-\_\:]', "", object_date.__str__())

        # If the subfolder out_dir/object_type does not exist, create it,
        # because we will create/move the new file there.
        newdir = os.path.join(args.out_dir, object_type)
        utils.if_dir_not_exists_create(newdir)

        # Create the new name of the file
        new_name = os.path.join(newdir, "{0}_{1}".format(object_name, object_date))

        # New name for the file will be determined by the object type (for the
        # subfolder), object, date and filter. For bias frames, the filter would 
        # have no meaning, so we don't put it.
        if object_type != "bias":
            new_name = "{0}_{1}_".format(new_name, object_filter)


        # Now we need to find out which sequential number the image should have
        ans = True
        jj=1
        while ans == True:         #  until file does not exist
            newfile = os.path.join(new_name + str(jj).zfill(3)+'.fits')
            ans = os.path.isfile(newfile)
            jj += 1
        oldname_base = os.path.basename(im_name)
        newname_base = os.path.basename(newfile)


        # Add history comment into the header. If image is to be overwritten,
        # just update the image with the changes in the header and move it to 
        # its new name. Otherwise, save it to the new file immediately.
        if args.overwrite == True:
            im.flush()
            im.close()
            os.rename(im_name, newfile)
        else:
            im.writeto(newfile)

        im = fits.open(newfile, 'update')
        hdr = im[0].header
        hdr.add_history("- Image "+oldname_base+" renamed "+newname_base)
        im.flush()
        im.close()

        # Add image to the dictionary for the output.       
        final_dict["filename"] = numpy.append(final_dict["filename"], newfile)
        final_dict["objname"] = numpy.append(final_dict["objname"], object_name)
        final_dict["type"] = numpy.append(final_dict["type"], object_type)

                
    # Return the dictionary with the images sorted in groups.
    return final_dict
########################################################################################################################
# Create parser
parser = argparse.ArgumentParser(description='Rename files in a folder ' +\
                                 ' using the headers')    
# Add arguments to the parser
parser.add_argument("in_dir", metavar='input_dir', action='store', \
                      help='Directory to be sorted. ' +\
                     ' Default:     "./" ', nargs=1)
parser.add_argument("out_dir", metavar='output_dir', action='store', \
                     default = '', help='Directory to save the renamed ' +\
                     'files. Default: same as input_dir', nargs="?")  
parser.add_argument("--objectk", metavar='objectk', action='store', \
                     default = '', help='Name of the keyword in the headers that'+\
                     ' contain the name of the object. This can also be provided '+\
                     'with the --config_file option')                       
parser.add_argument("--filterk", metavar='filterk', action='store', \
                     default = '', help='Name of the keyword in the headers that'+\
                     ' contain the name of the filter. This can also be provided '+\
                     'with the --config_file option')  
parser.add_argument("--exptimek", metavar='exptimek', action='store', \
                     default = '', help='Name of the keyword in the headers that'+\
                     ' contain the exposure time. This can also be provided '+\
                     'with the --config_file option')  
parser.add_argument("--datek", metavar='datek', action='store', \
                     default = '', help='Name of the keyword in the headers that'+\
                     ' contain the date. This can also be provided '+\
                     'with the --config_file option')  
parser.add_argument("--timek", metavar='timek', action='store', \
                     default = '', help='Name of the keyword that contains the '+\
                     'observing time. If the time is already present in the '+\
                     '--datek keyword just omit this one.')                     
parser.add_argument("--in_pattern", metavar='in_pattern', action='store', \
                     default='*', help='Only sort files starting with this '+\
                     'pattern. Default: "*" ')
parser.add_argument("--out_pattern", metavar='out_pattern', action='store', \
                     default='', help='Pattern to be added at the beginning '+\
                     'of the files. Default: common part of the files or' +\
                     ' the word output. ')
parser.add_argument("--copy", action="store_true", dest="copy", default=False,\
                     help="When activated a copy of the files will be " +\
                     "created for safety in a directory called 'raw' within" +\
                     " the directory of the data") 
parser.add_argument("--overwrite", action="store_true", dest="overwrite", \
                     default=False, help="When activated the names of the " +\
                     "files will be overwritten, instead of generating new " +\
                     "ones.")
parser.add_argument("--config_file", action="store", dest="config", default="", \
                    help="Config file to introduce the names of the keywords "+\
                    "in the headers of the images. The file should have "+\
                    "a line per keyword, where we give the common name, "+ \
                    "a '=' and the name of the keyword. E.g. for CAHA 2.2m: \n"\
                    "OBJECT = OBJECT \n FILTER = INSFLNAM \n DATE = DATE-OBS \n" +\
                    " At least the keywords for OBJECT, FILTER, EXPTIME AND " +\
                    " DATE should be included for this routine to work. " )
    
def main(arguments=None):
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    args.in_dir = (args.in_dir[0])
    
    # Check arguments
    if args.in_dir[-1] != "/":    # If path passed without "/" at the end, add it
        args.in_dir = args.in_dir+"/"
    if args.out_dir == '':      # If no output dir passed, use input dir 
        args.out_dir = args.in_dir
    if args.out_dir[-1] != "/":   # If path passed without "/" at the end, add it
        args.out_dir = args.out_dir+"/"
    if args.timek == '':    # If no timek given, assume datek contains time
        args.timek = args.datek

    # Call function with the args
    fits_dict = rename(args)
    return fits_dict        
if __name__ == "__main__":
    main()
