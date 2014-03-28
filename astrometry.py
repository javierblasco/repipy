#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Created on Tue Jul  9 11:03:50 2013

@author: javier blasco herrera
"""
import numpy as np
import repipy.cross_match as cross_match
import sys
from astropy.io import fits

input_usno = open("./hz15_ra_dec_mag.dat", "r")
input_image = open("/home/blasco/Desktop/mierda/20050110.reg", "r")
ra_array = dec_array = mag_array = xx = yy = np.array([], dtype=np.float64)
for line in input_usno:
     if line.split()[0] != "#":
         ra, dec, mag = line.split()[:]
         ra_array = np.append(ra_array, float(ra))
         dec_array = np.append(dec_array, float(dec))
         mag_array = np.append(mag_array, float(mag))

if len(ra_array) > 20: # if more than 20 stars, select the 20 brightest
    bright20 = np.argsort(mag_array)[0:20]
    ra_array = ra_array[bright20]
    dec_array = dec_array[bright20]
    mag_array = mag_array[bright20]

for line in input_image:
    xim, yim = line.split()
    xx = np.append(xx, float(xim))
    yy = np.append(yy, float(yim))

    
result = cross_match.main(xref=ra_array, yref=dec_array, xobj=xx, 
                           yobj=yy, error=0.008, scale=0.232/3600.)
#result = cross_match.main(test=True, error=0.01)
scale = result[0][0]
flip = result[1]
angle = result[2][0][0]
deltax = result[3][0][0]
deltay = result[3][0][1] 


print "scale =", scale, "angle =", angle, "flip =", flip, "deltax =",\
      deltax, "deltay =", deltay, "precission = ", result[4]
coords = np.asarray([xx,yy])



image = fits.open("/home/blasco/Desktop/mierda/prueba.fits",
                  mode="update")
hdr = image[0].header

hdr.update("CRPIX1",  1., "X reference pixel")
hdr.update("CRPIX2",  1., "Y reference pixel")
hdr.update("CRVAL1", deltax, "RA of reference pixel")
hdr.update("CRVAL2", deltay, "DEC of reference pixel")
hdr.update("CD1_1",  np.cos(angle) * scale, "rotation matrix coefficient")
hdr.update("CD2_2",  np.cos(angle) * scale, "rotation matrix coefficient")
hdr.update("CD1_2",  -np.sin(angle) * scale, "rotation matrix coefficient")
hdr.update("CD2_1",  np.sin(angle) * scale, "rotation matrix coefficient")
hdr.update("CUNIT1", "deg", "degrees")
hdr.update("CUNIT2", "deg", "degrees")
hdr.update("WCSNAME", "DSS", "local WCS approximation")
hdr.update("CTYPE1", "RA---TAN", "RA-Gnomic projection")
hdr.update("CTYPE2", "DEC---TAN", "DEC-Gnomic projection")
image.close()

