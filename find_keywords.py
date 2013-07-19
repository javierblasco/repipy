#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 13:06:55 2013
@author: blasco

This program returns how several keywords are called in the header of an image. 
The problem is that each telescope/instrument can save the information in the 
headers with different keywords. For example, a keyword in which the DATE
is stored can be called DATE or DATE-OBS. This routine is used by the programs
rename.py, combine.py... The user can either pass the path to a configuration 
file which contains:
    
    FILTER    keyword_filter_for_my_instrument
    DATE      keyword_date_for_my_instrument
    ...

and so on for all the needed keywords or, can pass when calling the program the 
necessary keywords, for example with a --filterk keyword_filter_for_my_instrument
when calling the programs. 

There is a third option, which happens when the keywords are called exactly the 
same as the alias we are using (RA, DEC, FILTER, EXPTIME, DATE, TIME, ...).
If that happens, even if the user passes nothing, the program will detect them 
and pass those keywords. 


INPUT: 
    hdr: header (pyfits) of the image to check for the existence of the keywords
    needed: list with which keywords exactly we need to find. Combine.py will need
            a set that does not coincide with the one needed by rename.py and so on
    args: parser of commands given by the user when calling combine.py, rename.py or
          any other. It includes information like the path to the configuration 
          file in which the keywords are defined or, if the user included them in 
          the call to the program all the needed keywords 
OUTPUT:
    keywords: dictionary which contains the "translation", i.e. under which keyword
              the "needed" things were saved in the header.
    

"""
import sys
import StringIO
import ConfigParser

def get_keywords(hdr, needed, args):
    """ Reads a config file that contains which keywords are used by a particular 
        telescope/instrument in the headers. It will search for those in the 
        variable 'needed'. 
    """    
    # Check if file exists. It does not exist?
    if args.config == "":
        # Check if, by any chance, ALL the needed keywords are called exactly 
         # the same in the header.
        ans = True
        for key in needed:
            ans = ans*hdr.has_key(key.upper())  

        # If it found all of them, you are in luck, even with no config file 
        # the program can build the dict. 
        if ans == True:
            keywords = {key:key for key in needed}

        # If any of the keywords is not in the header prompt the user to create
        # a config file or use the arguments passed to the program to identify 
        # the needed keywords in the header. 
        else:
            print " \n \n"
            print " ERROR: Unknown keywords! I do not know the names of keywords"+\
                  " in the headers of these particular files."
            print " Please, use the argument --config_file to pass a file with the "+\
                  " syntaxis: \n \n OBJECT = NAME_OF_OBJECT_KEYWORD_IN_YOUR_HEADER"
            print " FILTER = NAME_OF_FILTER_KEYWORD_IN_YOUR_HEADER "
            print " ... = ... \n"
            print " Alternatively, you can provide the needed keywords with the "+\
                  " corresponding arguments. \n For that, check the help: python "+\
                  " rename.py -h"
            print " Either way, you need to provide " + " ".join(needed).upper() +\
                  " and DATE are called in your headers.\n\n"
            sys.exit(" Exiting program\n")

    #If file exists simply read it and check the needed keywords are present. 
    else:
        # Trick to read a file with ConfigPaser
        ini_str = '[root]\n'+open(args.config,'r').read()
        ini_fp = StringIO.StringIO(ini_str)
        config = ConfigParser.RawConfigParser()
        config.readfp(ini_fp)

        # Convert to a dictionary
        keywords = dict(config.items("root"))
        
        # If any of the keys is missing, tell the user!
        for key in needed: 
            if keywords.has_key(key) == False: 
                sys.exit("Error! Keyword " + key + " is missing in " + args.config) 
    return keywords            
