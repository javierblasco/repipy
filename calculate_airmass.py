#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import astropy.io.fits as fits
import sys
import repipy.utilities as utils
from pyraf import iraf
import argparse
parser = argparse.ArgumentParser(description='Check header for airmass and '+\
                                'include it if missing and possible. From the '+\
                                'local time, sidereal time and universal time '+\
                                'we need at least UT. ')
parser.add_argument("input", metavar='input', action='store', nargs="+", \
                    help='Image(s) for which airmass will be checked.') 
parser.add_argument("--ut", metavar='UT', action='store', dest='UT', default='',\
                    help=" Keyword for Universal Time (hh:mm:ss).")
parser.add_argument("--st", metavar='ST', action='store', dest='ST', default='ST', \
                    help=" OUTPUT Keyword for Local Sidereal Time (hh:mm:ss)")
parser.add_argument("--lt", metavar='LT', action='store', dest='LT', default='LT',\
                    help=" OUTPUT Keyword for Local Time (hh:mm:ss) ")
parser.add_argument("--RA", metavar='Right ascension', action='store', dest='ra',\
                    default='', help= "Keyword for Right ascension of the object "+\
                    "in hours!!! (3.45 or 23.2 or 12.56)")
parser.add_argument("--DEC", metavar='Declination', action='store', dest='dec',\
                    default='', help= "Keyword for the Declination of the object "+\
                    "in degrees (80.437, 23.61, 12.6)")
parser.add_argument("--date", metavar='date', action='store', dest='date',\
                    default='', help=" Keyword for date (yyyy-mm-dd) ")
parser.add_argument("--observatory", metavar='Observatory', action='store', \
                    dest='observatory', default='', help=" Name of Observatory "+\
                    "of origin. E.g. OSN, ORM, ... ")
parser.add_argument("--location", metavar='location', action='store', dest='location',\
                    default='', help=" Time zone location. E.g.: 'Europe/Madrid' ")
parser.add_argument("--equinox", metavar='equinox', action='store', dest='equinox',\
                    default='', help=" Keyword for equinox of the coordinates ")

def estimate_airmass(args):
    for image in args.input:
        im = fits.open(image, mode="update")
        hdr = im[0].header
        
        # Read all the variables from the header
        date = hdr[args.date]        
        UT = hdr[args.UT]        
        location = args.location 
        if hdr.has_key("AIRMASS"):
            hdr.update("AIRM_0", hdr["AIRMASS"], "Original airmass (not recalcul)")
        
        # Calculate local time from UT. If user provides with a name for LT in 
        # the header save it. Otherwise no need to save it. 
        LT = utils.universal_time_to_local_time(date, UT, location)
        if args.LT != "":
            hdr.update(args.LT, LT, " Local time ")
            print "Local time ", LT
            
        # Calculate local sidereal time from local time and location of observ.
        ST = utils.local_to_sidereal_time(date, LT, args.observatory)
        hdr.update(args.ST, ST, " Local sidereal time ")    
    
        im.flush()
        im.close()
        
        # Finally, estimate airmass, including it in the header
        iraf.module.load("astutil", doprint=0)
        iraf.module.setairmass(images=image, observa=args.observatory, ra=args.ra,\
                                dec=args.dec, st=args.ST, ut=args.UT,\
                                equinox=args.equinox)        
        
    return None       
       

def main(arguments=None):
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)   
    estimate_airmass(args)

if __name__ == "__main__":
    main()
