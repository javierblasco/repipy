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
import repipy.find_keywords as find_keywords
import numpy as np
import repipy.utilities as utils

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
    standards_list = ["bd+25", "bd+28", "bd+35", "kop27", "grw73", "pg1708", \
                    "f110", "pg170", "wolf1346", "sa113", "hz15", "feige34", \
                    "pg1633", "sa110", "pg1323", "pg16", "sao110"]
    for standard in standards_list:
        if ((object_name).lower()).count(standard) != 0 :
            return True
    return False

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
    elif is_a_standar(object_name.lower()) == True:
        object_name = object_name.lower()
        object_type = "standards"
    elif is_a_cluster(object_name.lower()) == True:
        object_name = object_name.lower()
        object_type = "clusters"
    elif object_name.count("cig") != 0:
        # cig usually folloed by a number
        cig, number = re.match(r'(cig)\s?(\d{1,4}).*', object_name).groups()[0:2]
        object_name = cig + "{0:04d}".format(int(number))
	object_type = "cig"
        #if number.isdigit() == True and len(number) <= 4:
        #    number = "{0:04d}".format(int(number))
        #    object_name = "cig" + number
        #    object_type = "cig"
    # If "cig" is missing, but there is a number, add "cig":        
    elif object_name.isdigit() == True :
        object_name = "{0:04d}".format(int(object_name))
        object_name = "cig" + object_name
        object_type = "cig" + object_name
    else:
        object_type = "unknown"

    return object_name, object_type
    

###########################################################################
def rename(args):
    # List of fit and fits images in the directory
    fits_list1 = glob.glob(os.path.join(args.in_dir, args.in_pattern + "*.fits"))
    fits_list2 = glob.glob(os.path.join(args.in_dir, args.in_pattern + "*.fit"))
    fits_list = fits_list1+fits_list2

    # If the needed keywords were passed by the user, build a dictionary with 
    # them, otherwise, read from config file (if present) the names of the 
    # different keywords. 
    needed = ["object", "filter", "date", "time", "exptime"]
    if args.objectk != "" and args.filterk != "" and args.datek != "" and\
       args.exptimek != "":
           keywords = {"object":args.objectk, "filter":args.filterk,\
                       "date":args.datek, "exptime":args.exptimek,\
                       "time":args.timek}
    else:
        hdr = fits.getheader(fits_list[0])
        keywords = find_keywords.get_keywords(hdr, needed, args) 

    # The output of the whole code will be this dictionary, in which the images 
    # are sorted in groups (bias, skyflats, domeflats, cigXXXX, ...)
    empty_array = np.asarray([], dtype=object)
    final_dict = {"filename":np.asarray([], dtype="S150"), # defaults S70 too small 
                  "type":empty_array, 
                  "objname":empty_array,
                  "time":empty_array}


    # If --copy was selected, copy all those files into a directory called backup
    if args.copy == True:
        raw_dir = os.path.join(args.in_dir, "backup")
        if os.path.isdir(raw_dir) == False:
            os.makedirs(raw_dir)
        for im in fits_list:
            shutil.copy(im, raw_dir)

    # Create log file
    print "\n \n "
    print "###########################################################\n"
    print "RELEVANT INFO IS BEING WRITTEN IN: \n \n" +\
                   os.path.join(os.path.abspath(args.out_dir),"renamed.log\n")
    print "###########################################################\n "
    ff = open(os.path.join(args.out_dir, "renamed.log"), 'w')
    ff.write("CHECK WITH OBSERVING LOG: \n \n")
    ff.write("ORIG_NAME, NEW_NAME   , OBJECT   , FILTER  , DATE  , EXPTIME \n ")    

    # Look for the date and time of all images
    list_datetimes =[]
    for names in fits_list:
        im = fits.open(names) 
        hdr = im[0].header  
        date_current = dateutil.parser.parse(hdr[keywords["date"]]).date()
        time_current = dateutil.parser.parse(hdr[keywords["time"]]).time()
        datetime_current = datetime.datetime.combine(date_current, time_current) 
        # If user didn't provide any time keyword but date_current does not 
        # actually contain the time this time_current will be 00:00:00. Problem?
        list_datetimes.append(datetime_current)
               
    # Looking for the smaller of dates (i.e. before 00:00 if present)
    sort_indices = np.argsort(list_datetimes)
    min_date = list_datetimes[sort_indices[0]].date()
    date = str(min_date).replace("-","")

    # And now sort files in fits_list using the date and time
    fits_list_sorted = [fits_list[ii] for ii in sort_indices]
    fits_list = fits_list_sorted
    list_datetimes_sorted = [list_datetimes[ii] for ii in sort_indices]
    final_dict["time"] = np.asarray(list_datetimes_sorted)
    
    # Run through all images
    for image in fits_list:
        # Read image and header, extract name of object and filter.
        im = fits.open(image, mode='update')
        hdr = im[0].header
        object_name = (hdr[keywords["object"]].lower())
        remove_characters = [" ", "/", "[", "]", "_"]
        for character in remove_characters:
            object_name = object_name.replace(character,"")  


        # Homogeneous names for filters!
        object_filter = hdr[keywords["filter"]]
        object_filter = utils.homogeneous_filter_name(object_filter)

            
        # find type of objects: bias, skyflat, domeflat, blank...
        object_name, object_type = distinguish_type(object_name)

        # If the subfolder out_dir/object_type does not exist, create it, 
        # because we will create/move the new file there. 
        newdir = os.path.join(args.out_dir, object_type)
        if os.path.isdir(newdir) == False:
            os.makedirs(os.path.join(args.out_dir, object_type))
            
        # New name for the file will be determined by the object type (for the 
        # subfolder), object, date and filter. For bias frames, the filter would 
        # have no meaning, so we don't put it.
        new_name = os.path.join(newdir, object_name+"_"+date+"_")
        if object_name != "bias":
            new_name = new_name + object_filter +"_" 

        # Now we need to find out which sequential number the image should have    
        ans = True
        jj=1
        while ans == True:         #  until file does not exist
            newfile = os.path.join(new_name+str(jj).zfill(3)+'.fits')
            ans = os.path.isfile(newfile)
            jj += 1
        oldname_nodir = (os.path.split(image))[1]
        newname_nodir = (os.path.split(newfile))[1]

        # Add history comment into the header. If image is to be overwritten, 
        # just update the image with the changes in the header and move it to 
        # its new name. Otherwise, save it to the new file immediately. 
        hdr.add_history("- Image "+oldname_nodir+" renamed "+newname_nodir)
        if args.overwrite == True:
            im.flush()
            os.rename(image, newfile)
        else:
            im.writeto(newfile)

        # Add image to the dictionary for the output.       
        final_dict["filename"] = np.append(final_dict["filename"], newfile)
        final_dict["objname"] = np.append(final_dict["objname"], object_name)
        final_dict["type"] = np.append(final_dict["type"], object_type)
            
        # And write the log    
        ff.write(oldname_nodir+"  "+ newname_nodir + " " + object_name +"  "+\
                 object_filter + "  " + str(hdr[keywords["date"]]) + "  " +\
                 str(hdr[keywords["exptime"]]) + "\n")
                
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
