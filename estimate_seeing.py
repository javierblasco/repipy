#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Created on Mon Nov  4 01:16:30 2013

@author: blasco
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 23:14:24 2013

@author: blasco
"""
import repipy.utilities as utilities
import numpy as np
import astropy.io.fits as fits  
import argparse
import sys
from scipy import spatial 
import astropy.wcs.wcs as wcs 
import os

from lemon import methods
import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import obsutil
    from iraf import digiphot, apphot, daophot


def calculate_seeing(args):
    """ Program to estimate the seeing from an image and a list of estimates
    for the positions of stars. After calculating the seeing, some of the 
    stars might get recalculated centers. The list will be updated with the 
    """
    for im, im_cat in zip(args.input, args.cat):
        output = "output.txt"
        ignore = "ignore.txt"
        # Victor Terron, esto es horroroso, sugerencias?
        utilities.if_exists_remove("q.txt", output, ignore)
        q = open("q.txt", "w")
        q.write("q")
        q.close()
        iraf.psfmeasure(im, coords = "mark1", size = "MFWHM",
                               sbuffer = 10, swidth=10,radius=10,
                               satura = 55000, ignore = "yes", 
                               imagecur = im_cat, display = "no",
                               graphcur = "q.txt", Stdout=ignore, 
                               wcs = args.wcs, logfile=output)
        
        # Now we read the input im_cat and the output output.txt and compare
        # the location of stars. Those stars that have moved will not be 
        # trusted. 
        xout = np.array([])
        yout = np.array([])
        FWHM = np.array([])
        for line in open(output, 'r'):
            # First line is a description of the file, second is a newline \n
            # and third the description of the columns.
            if line != "\n" and line.split()[0] in im:
                xout = np.append(xout, float(line.split()[1]))
                yout = np.append(yout, float(line.split()[2]))
                FWHM = np.append(FWHM, float(line.split()[4]))
            else:
                try:
                   xout = np.append(xout, float(line.split()[0])) 
                   yout = np.append(yout, float(line.split()[1]))
                   FWHM = np.append(FWHM, float(line.split()[3]))
                except:
                   pass 
        xin, yin = np.genfromtxt(im_cat, dtype="float", unpack=True)
        xin, yin = np.array(xin), np.array(yin)  # case it is only one value


        # If args.wcs is "world" it means the input is in (RA, DEC), while 
        # the output is in pixels (X,Y). We need to convert one to the other
        # and we choose to follow the (RA,DEC) which, after all, is meaningful
        if args.wcs == "world":
            coords_xy = np.array(zip(xout, yout))
            hdr = fits.open(im)[0].header
            remove_keys = ["PC001001", "PC001002", "PC002001", "PC002002"]
            for keys in remove_keys:
                hdr.pop(keys,None)
            w = wcs.WCS(hdr)            
            coords_RADEC = w.all_pix2world(coords_xy,1)
            xout = coords_RADEC[:,0]
            yout = coords_RADEC[:,1]
           
        # find common stars using KDtrees    
        if xin.size > 1 : 
            tree_in = spatial.KDTree(zip(xin, yin))
        elif xin.size == 1:
            tree_in = spatial.KDTree([(float(xin), float(yin))])
        if xout.size > 1:
            tree_out = spatial.KDTree(zip(xout, yout))
        elif xout.size == 1:
            tree_out = spatial.KDTree([(float(xout), float(yout))])
        # If WCS use 0.001 degree (~3.6 arcsec) as limit. If not, assume 
        # pixels and say 4 pixels
        if args.wcs == "world":
            limit = 0.001
        else:
            limit = 4
        matches = tree_out.query_ball_tree(tree_in, limit) 
           
        # Two close stars in the original .cat could resolve into one  when 
        # the FWHM is calculated. This will appear as several hits in the 
        # matches with exactly the same numbers: [[0], [1,2], [1,2], [3]]
        # One solution is to erase one of them 
        for index, value in enumerate(matches):      
            if matches.count(value) > 1:
                matches[index] = []
                       
        # Now restrict to the common objects
        remove_indices = [ii for ii,jj in enumerate(matches) if jj == []]
        xout = np.delete(xout, remove_indices)
        yout = np.delete(yout, remove_indices)
        FWHM = np.delete(FWHM, remove_indices)
 
        # Finally, calculate the median FWHM of the image and rewrite valid
        # stars to the im_cat file. 
        median_fwhm = np.median(FWHM)
        utilities.header_update_keyword(im, "seeing", median_fwhm, "FWHM of image")
        f = open(im_cat, 'w') # write the "good" stars in the catalogue
        for ii in range(len(xout)):
            f.write(str(xout[ii]) + "  " + str(yout[ii]) + "\n")
        f.close()

        # And clean after yourself!
        utilities.if_exists_remove("q.txt", output, ignore)


############################################################################


# Create parser
parser = argparse.ArgumentParser(description='Program to find stars in an image. ')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'input images for which to estimate the FWHM.',
                    nargs="+", type=str)
parser.add_argument("--cat", metavar='cat', action='store', dest="cat",
                    help='list of catalogues of the position of stars for the input images. ' +\
                        'You can decide to give 1 catalogue per image, 1 catalogue for all the images '+\
                        '(assming all images are of the same object) or None if catalogues exist for all the ' +\
                        'images and are named exactly like the image but with .radec extension.',
                    nargs=1, type=str)
parser.add_argument("--wcs", metavar="wcs_in", action="store", dest="wcs", 
                    default="logical",
                    help = "System in which the input coordinates are. Can "
                    "be 'logical', 'tv', 'physical' and 'world' accordind to "
                    "IRAF. If your coordinates are, for example, in RA and DEC "
                    "you should provide 'world'. Default: logical.")


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)

  # Option 1: the user did not provide any catalogues, because they are named like the images, but with '.radec' extension
  # Option 2: the user provided a single catalogue, because all the images contain the same object
  # Option 3: the user provided an invalid number of catalogues for the number of images.
  if args.cat is None:
      args.cat = [utilities.replace_extension(im_name, ".radec") for im_name in args.input]
      if not all([os.path.exists(cat_name) for cat_name in args.cat]):
          sys.exit("Catalogues not found! Check the --cat option in the description: type estimate_seeing.py -h")
  elif len(args.cat) == 1:
      args.cat *= len(args.input)
  elif len(args.input) != len(args.cat):
      sys.exit("\n\n number of star catalogues and input images do not coincide \n ")


  calculate_seeing(args)  
  return None    
     
if __name__ == "__main__":
    main()
