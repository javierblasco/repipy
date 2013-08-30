# -*- coding: utf-8 -*-

import argparse
import sys
import pyraf.iraf as iraf
import astropy.io.fits as fits
import urllib2
from astropy import coordinates as coord


""" This program includes in an image the missing keywords if the user gives some
appropriate values, such as the observatory. """

parser = argparse.ArgumentParser(description='Check header and include missing' +\
                                 'keywords if possible.')
parser.add_argument("input", metavar='input', action='store', nargs=1,  \
                    help='Image for which certain parameters will be checked.')                             
parser.add_argument("--object", metavar="object", action="store", dest='object', \
                    default='', help='Name of the object. This will be used to'+\
                    ' search for coordinates in Sesame and NED.' )
parser.add_argument("--observatory", metavar='observatory', action='store', \
                    dest='observatory', default='', help=" Observatory from "+\
                    "which the data come. E.g. 'OSN', 'ORM', ...")
parser.add_argument("--longitude", metavar="longitude", action='store', default='', \
                    dest="longitude", help=" Longitude of the observatory (in "+\
                    " degrees). ")                    
parser.add_argument("--latitude", metavar="latitude", action='store', \
                    dest="latitude", default='', help=" Latitude of observatory"+\
                    " (in degrees). " )
parser.add_argument("--altitude", metavar="altitude", action='store', \
                    dest="altitude", default="", \
                    help=" Altitude of the observatory (in meters).")                    
parser.add_argument("--timezone", metavar="timezone", action="store", default="",\
                    dest="timezone", help="Timezone of the observatory (hours). "+\
                    "Positive towards the west. E.g. Mainland spain = -1.")  
parser.add_argument("--RA_keyword", metavar='RA_keyword', action='store', \
                    dest='RA_keyword', default='', help=" Keyword to access RA "+\
                    "in the header if present.")
parser.add_argument("--DEC_keyword", metavar='DEC_keyword', action='store', \
                    dest='DEC_keyword', default='', help=" Keyword to access DEC "+\
                    "in the header if present.")

def observatory_position(observatory):
    longitude = latitude = altitude = timezone = ''
    if observatory == "OSN":
        longitude = "3.38472"   # degree
        latitude = "37.06417"   # degree
        altitude = "2896"       # meters 
        timezone = "-1"
    return longitude, latitude, altitude, timezone

def add_obs_position(hdr,args):
    longitude,latitude,altitude,timezone = observatory_position(hdr["observat"])        
    # User values override the automatic ones:
    if args.longitude != "":
        longitude = args.longitude
    if args.latitude != "":
        latitude = args.latitude
    if args.altitude != "":
        altitude = args.altitude
    if args.timezone != "":
        timezone = args.timezone
    # Add to header:
    hdr.update("longit", longitude, "Longitude of observatory (degree)")
    hdr.update("latit", latitude, "Latitude of observatory (degree)")
    hdr.update("altit", altitude, "Altitude of observatory (meters)")
    hdr.update("timezone", timezone, "Time zone of observatory (>0 if West)")        
    return hdr

def add_obj_coordinates(hdr, object_name, args):
    RA_hours = DEC = redshift = equinox = ""

    # Standard star he 3 is not found in NED
    if object_name[:] == "he3":
        object_name = "WD 0644+375"

    # Search NED for the object_name
    try:
        query = urllib2.urlopen("http://ned.ipac.caltech.edu/cgi-bin/"+\
                                "nph-objsearch?extend=no&out_equinox=J2000&of="+\
                                "ascii_tab&objname=" + object_name)
        result = query.readlines()
        last_line = result[-1].split()
        RA = last_line[3]
        RA_hours = RA * 24/360.         
        DEC = last_line[4]
        equinox = "2000"
        redshift = last_line[7]             
    except:
        pass
    
    # If still no RA_obj is found, try finding them on Sesame
    if RA_hours == "":
        try:
            c = coord.ICRSCoordinates.from_name(object_name)
            RA_hours = c.ra.hours
            DEC = c.dec.degrees
            equinox = c.equinox.jyear_str[1:]
        except:
            pass

    # Maybe there is a RA, DEC in the header?
    if RA_hours == ""  and args.RA_keyword != "" and args.DEC_keyword != "":
        if hdr.has_key(args.RA_keyword):
            RA = hdr[args.RA_keyword] 
            DEC = hdr[args.DEC_keyword]
            # Format of RA?
            if RA.count(":") > 0:          # hours:min:sec
                RA = RA.split(":")         # [hours, min, sec]
                RA_hours = float(RA[0]) + float(RA[1])/60. + float(RA[2])/3600.
            elif len(RA.split()) == 3:     # hours min sec
                RA = RA.split()            # [hours, min, sec]
                RA_hours = float(RA[0]) + float(RA[1])/60. + float(RA[2])/3600.
            elif len(RA.split()) == 0:     # single number?
                try:
                    RA = float(RA)         # Assuming it is in degrees
                    RA_hours = RA*24./360
                except:
                    pass
            # Same for DEC:
            if len(DEC.split()) == 3:      # deg arcmin arcsec
                DEC = DEC.split()
                DEC = float(DEC[0]) + float(DEC[1])/60. + float(DEC[2])/3600.
            RA_hours = str(RA_hours)
            DEC = str(DEC)
            print "Not found in NED or Sesame:", object_name, RA_hours, DEC
                    
    # Write results, whatever you have found.     
    hdr.update("object", object_name, "Name of the object")
    hdr.update("RA_hours", RA_hours, "RA of the object (hours)")
    hdr.update("DEC_deg", DEC, "DEC of the object (degrees)")
    hdr.update("equinox", equinox, "equinox of coordinates")
    if redshift != "":
        hdr.update("redshift", redshift, "Redshift of the object")  
    return hdr        
    
    
    
def complete_headers(args):
    """ For each image check for all the necessary keywords, when possible
    add those not present """
    # Open image and read header in update mode so that we can include keywords. 
    im = fits.open(args.input[0], mode="update")
    hdr = im[0].header
    manipulated_keywords = ""    
    
    # Add the name of the image for future reference
    hdr.update("org_name", args.input[0], "Name of original image")

    # Add observatory keyword
    if args.observatory != "":
        hdr.update("observat", args.observatory, "Name of the observatory")
        manipulated_keywords += "Observatory name"        
    
    # Add position (long, lat, alt, timezone) of the observatory     
    if args.observatory != "" or args.longitude != "" or args.latitude != "" or \
       args.altitude != "" or args.timezone != "":
           hdr = add_obs_position(hdr, args)
           manipulated_keywords += ", observatory information" 
           
    # Now, if args.name exist or it is in the header, use it to calculate
    # the coordinates. 
    if args.object != "":
        hdr = add_obj_coordinates(hdr,args.object, args)
        manipulated_keywords += ", object coordinates"
    elif hdr.has_key("object") and hdr["object"].replace(" ","") != "": 
        hdr = add_obj_coordinates(hdr, hdr["object"], args)       
        manipulated_keywords += ", object coordinates"     
    if hdr.has_key("RA_hours") == False:  # Coordinates not added, not found?
        print "Problem: coordinates for file " + args.input[0] + " not found!"
    
    # Add history comments to record changes
    history = str(hdr.get_history())
    if history.count(manipulated_keywords) != 0 :
        hdr.add_history("HEADER UPDATED:")
        hdr.add_history(" -Updated: " + manipulated_keywords)
    im.flush()
    im.close()

def main(arguments=None):
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)   
    complete_headers(args)

if __name__ == "__main__":
    main()
