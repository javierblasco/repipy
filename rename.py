#! /usr/bin/env python
"""
@author: javier blasco herrera
@e-mail: blasco@iaa.es  javier.blasco.herrera@gmail.com 

Routine to rename files with a single homogeneous structure that actually gives
information about its content. No more 13125.fits in your data!

It also organizes your data. For example, flats in a folder, standards (if 
recoznized) in another... 

"""

#############################################################################################################################
#    TO BE FIXED
#
#   - Bug: When checking if the file exists the program looks for blabla.fits//b , 
#     obviously it doesn't find it. Should look for 
#     balbla-b.fits instead. Fix it whenever you have time. 
#   - Using the "telescop" keyword to determine the names of all the others is 
#     not a good idea. 
#     Better ask the user to input the file with the names of the keywords. 
#   - The date is being read assuming all the observations come from the same 
#     night. Could lead to mistakes? Other options?
#   - All the outcomes and ifs of the read_config  have not been fully tested yet
##############################################################################################################################

import os, glob
import shutil, pyfits
import sys
import argparse
import ConfigParser
import dateutil.parser
import StringIO
import re
import reduction_pipeline.find_keywords as find_keywords

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
                    "f110", "pg170", "wolf1346", "sa113", "hz15", "feige34"]
    for standard in standards_list:
        if ((object_name).lower()).count(standard) != 0 :
            return True
    return False

##############################################################################
#def read_config(hdr, needed, args):
#    """ Reads a config file that contains which keywords are used by a particular 
#        telescope/instrument in the headers. It will search for those in the 
#        variable 'needed'. 
#    """
#    # If the user gave all the needed keywords when calling the program:
#    if args.filterk != "" and args.datek != "" and args.exptimek != "" and \
#       args.objectk != "":
#           keywords = {"date": args.datek , "filter": args.filterk , \
#                       "exptime": args.exptimek, "object":args.objectk }     
#    # Check if file exists. It does not exist?
#    elif args.config == "":
#        # Check if, by any chance, ALL the needed keywords are called exactly 
#         # the same in the header.
#        ans = True
#        for key in needed:
#            ans = ans*hdr.has_key(key.upper())  
#
#        # If it found all of them, you are in luck, even with no config file 
#        # the program can build the dict. 
#        if ans == True:
#            keywords = {key:key for key in needed}
#
#        # If any of the keywords is not in the header prompt the user to create
#        # a config file or use the arguments passed to the program to identify 
#        # the needed keywords in the header. 
#        else:
#            print " \n \n"
#            print " ERROR: Unknown keywords! I do not know the names of keywords"+\
#                  " in the headers of these particular files."
#            print " Please, use the argument --config_file to pass a file with the "+\
#                  " syntaxis: \n \n OBJECT = NAME_OF_OBJECT_KEYWORD_IN_YOUR_HEADER"
#            print " FILTER = NAME_OF_FILTER_KEYWORD_IN_YOUR_HEADER "
#            print " ... = ... \n"
#            print " Alternatively, you can provide the needed keywords with the "+\
#                  " corresponding arguments. \n For that, check the help: python "+\
#                  " rename.py -h"
#            print " Either way, you need to provide how OBJECT, FILTER, EXPTIME "+\
#                  " and DATE are called in your headers.\n\n"
#            sys.exit(" Exiting program\n")
#
#    #If file exists simply read it and check the needed keywords are present. 
#    else:
#        # Trick to read a file with ConfigPaser
#        ini_str = '[root]\n'+open(args.config,'r').read()
#        ini_fp = StringIO.StringIO(ini_str)
#        config = ConfigParser.RawConfigParser()
#        config.readfp(ini_fp)
#
#        # Convert to a dictionary
#        keywords = dict(config.items("root"))
#        
#        # If any of the keys is missing, tell the user!
#        for key in needed: 
#            if keywords.has_key(key) == False: 
#                sys.exit("Error! Keyword " + key + " is missing in " + args.config) 
#    return keywords            

###########################################################################
def rename(args):
    # List of fit and fits images in the directory
    fits_list1 = glob.glob(os.path.join(args.in_dir, args.in_pattern + "*.fits"))
    fits_list2 = glob.glob(os.path.join(args.in_dir, args.in_pattern + "*.fit"))
    fits_list = fits_list1+fits_list2
    fits_list.sort()

    # If --copy was selected, copy all those files into a directory called copy
    if args.copy == True:
        raw_dir = os.path.join(args.in_dir, "backup")
        if os.path.isdir(raw_dir) == False:
            os.makedirs(raw_dir)
        for jj in range(len(fits_list)):
            shutil.copy(fits_list[jj], raw_dir)

    # Create log file
    print "\n \n "
    print "###########################################################\n"
    print "RELEVANT INFO IS BEING WRITTEN IN: \n \n" +\
                   os.path.join(os.path.abspath(args.out_dir),"renamed.log\n")
    print "###########################################################\n "
    ff = open(os.path.join(args.out_dir, "renamed.log"), 'w')
    ff.write("CHECK WITH OBSERVING LOG: \n \n")
    ff.write("ORIG_NAME, NEW_NAME   , OBJECT   , FILTER  , DATE  , EXPTIME \n ")    

    # Look for the date (the starting date if before midnight)
    date = ""
    for names in fits_list:
        im = pyfits.open(names) 
        hdr = im[0].header

        # If the needed keywords were passed by the user, build a dictionary with 
        # them, otherwise, read from config file (if present) the names of the 
        # different keywords. 
        needed = ["object", "filter", "date", "exptime"]
        if args.objectk != "" and args.filterk != "" and args.datek != "" and\
           args.exptimek != "":
               keywords = {"object":args.objectk, "filter":args.filterk,\
                           "date":args.datek, "exptime":args.exptimek}
        else:
            keywords = find_keywords.get_keywords(hdr, needed, args)
     
        # Read date in format YYYYMMDD
        date_current = dateutil.parser.parse(hdr[keywords["date"]])
        date_current = str(date_current.date())
        date_current = date_current.replace("-","")  
        
        # Looking for the smaller of dates
        if date == "": 
            date = date_current
        if date_current < date: 
            date = date_current

    # Run through all images
    for image in fits_list:
        # Read header and extract relevant information: object and filter
        im = pyfits.open(image, mode='update')
        hdr = im[0].header
        object_name = (hdr[keywords["object"]].lower())
        # lower case, no spaces, no "/", no "["...
        object_name = object_name.replace(" ","")  
        object_name = object_name.replace("/","")
        object_name = object_name.replace("[","")
        object_name = object_name.replace("]","")          
        # find type of objects
        if object_name.count("flat") != 0 and object_name.count("sky") != 0 : 
              object_name = "skyflat"     
        if object_name.count("flat") != 0 and object_name.count("dome") != 0 : 
              object_name = "domeflat"
        if object_name.lower().count("bias") != 0: 
              object_name = "bias"
        # If object is cig+number, where number is three digits, add a fourth 
         # one, with a zero at the beginning
        if object_name.count("cig") != 0:
            number = (object_name.split("cig"))[1]
            if number.isdigit() == True and len(number) < 4:
                number = "{0:04d}".format(int(number))
                object_name = "cig"+number

        # If object is just a number, it usually means the number of the cig, 
        # add "cig" at the beginning.
        if object_name.isdigit() == True :
            object_name = "{0:04d}".format(int(object_name))
            object_name = "cig"+object_name

        # Read filter from keyword. Again, no spaces, no "/"
        filter_name = hdr[keywords["filter"]].replace(" ","")     
        filter_name = filter_name.replace("/","")    

        # New name for the file will be determined by object, date and filter.
        # Except for bias frames, in which the filter has no meaning.
        new_name = object_name+"_"+date+"_"
        if object_name != "bias":
            new_name = new_name + filter_name +"_" 

        # Now we need to find out which sequential number the image should have    
        ans = True
        jj=1
        while ans == True:         #  until file does not exist
            newfile = os.path.join(args.out_dir, new_name+str(jj)+'.fits')
            ans = os.path.isfile(newfile)
            if ans == False:        # this one does not exist yet
	            # Add history comment into the header
                hdr.add_history("- Image "+image+" renamed "+newfile)
                im.flush()
                if args.overwrite == True:
                    os.rename(image, newfile)
                else:
                    shutil.copy(image, newfile)
                oldname_nodir = (os.path.split(image))[1]
                newname_nodir = (os.path.split(newfile))[1]
                ff.write(oldname_nodir+"  "+ newname_nodir + " " + object_name +"  "+\
                           filter_name + "  " + str(hdr[keywords["date"]]) + "  " +\
                           str(hdr[keywords["exptime"]]) + "\n")
            jj=jj+1

    # Copy the flats from output directory to a subdir
    list_flats = glob.glob(os.path.join(args.out_dir,"skyflat*fits"))    
    for ii in list_flats:
        skyflat_dir = os.path.join(args.out_dir, "skyflats")
        if os.path.isdir(skyflat_dir) == False: 
              os.makedirs(skyflat_dir)
        if os.path.isfile(ii) == True: 
            shutil.copy(ii,skyflat_dir)
            os.remove(ii)

    # Copy blanks from output to a subdir
    blank_dir = os.path.join(args.out_dir, "blanks")
    list_flats = glob.glob(args.out_dir + "blank*fits")    
    for ii in list_flats :    
        if os.path.isdir(blank_dir) == False: 
              os.makedirs(blank_dir)
        if os.path.isfile(ii) == True: 
            shutil.copy(ii,blank_dir)
            os.remove(ii)
        

    # Copy the domeflats from output directory to a subdir
    dome_dir = os.path.join(args.out_dir, "domeflats")
    list_domeflats = glob.glob(os.path.join(args.out_dir,"domeflat*fits"))    
    for ii in list_domeflats :    
        if os.path.isdir(dome_dir) == False: 
              os.makedirs(dome_dir)
        if os.path.isfile(ii) == True: 
            shutil.copy(ii,dome_dir)
            os.remove(ii)

    # Copy the bias
    bias_dir = os.path.join(args.out_dir,"bias")
    list_bias = glob.glob(os.path.join(args.out_dir,"bias*"))
    for ii in list_bias :    
        if os.path.isdir(bias_dir) == False: 
              os.makedirs(bias_dir)
        if os.path.isfile(ii) == True: 
            shutil.copy(ii,bias_dir)
            os.remove(ii)    
    
    # Create a folder with any cig and move there the files
    cig_list = glob.glob(args.out_dir+"cig*.fits")
    for cig in cig_list:
        nm = os.path.split(cig)[1]    # remove the dir 
        nm = (nm.split("_"))[0]       # take just before the first _
        cig_dir = os.path.join(args.out_dir, nm)
        if os.path.isdir(cig_dir) == False: 
              os.makedirs(cig_dir)
        lst = glob.glob(cig_dir+"*.fits")    
        for ii in range(len(lst)):
            shutil.copy(lst[ii],cig_dir)
            os.remove(lst[ii])        


    # Create a folder for the standards, for the remaining files check if there
    #  is any .fit or .fits
    # Then check if they are likely standard fields and move them there
    standards_dir = os.path.join(args.out_dir, "standards")
    list = glob.glob(args.out_dir+"*.fits")
    for ii in list:
        if os.path.isdir(standards_dir) == False: 
              os.makedirs(standards_dir)
        nswr = is_a_standar(ii)
        if nswr == True: 
            shutil.copy(ii, standards_dir)
            os.remove(ii)
    
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
               
    # Call function with the args
    rename(args)        
if __name__ == "__main__":
    main()