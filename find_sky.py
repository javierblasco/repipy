#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import scipy
import numpy
import astropy.io.fits as fits
import os
import sys
import argparse
from scipy.optimize import curve_fit
from matplotlib import pyplot

""" This program estimates the sky in a fits image, by assuming that the mode
of the image represents the sky. It also estimates de standard deviation of 
the sky pixels and introduces the relevant keywords in the header of the 
image. Finally, it adds a "history" line to the header.  """

def gauss(x, *p):
    A,mu,sigma,cont = p
    return cont + A*numpy.exp(-(x-mu)**2/(2.*sigma**2))    

def find_sky(args):

    for im_name in args.input:
        # Separate path and file
        imdir, imname = os.path.split(im_name)

        # Read images
        im = fits.open(im_name, mode='update')
        hdr = im[0].header
        data = fits.getdata(im_name).astype(numpy.float64)

        # If mask exist read it, otherwise build it with all values unmasked
        try:
            maskname = hdr[args.mask_key]
            mask = fits.getdata(maskname)
        except KeyError:
            mask = numpy.ma.make_mask_none(data.shape)

        # Make a copy of the array, but only with the unmasked pixels
        whr2 = numpy.where( (mask == 0) & (numpy.isfinite(data)))
        data2 = data[whr2]
        median = numpy.median(data2)
        MAD = numpy.median( numpy.abs( data2 - median))


        # Do a histogram, with limits of +- 7.5 MAD (~+-5 sigma for a normal distribution)
        minx, maxx = median-7.5*MAD, median+7.5*MAD
        nbins = 15
        n, bins = numpy.histogram(data2, bins=nbins, range=[minx, maxx])
        bincenters = 0.5 * (bins[1:]+bins[:-1])

        # Find max position and value
        maxpos = n.argmax()
        maxvalue = bincenters[maxpos]

        # First guess for a fit
        p0 = [n[maxpos], maxvalue, 1.5 * MAD, numpy.min(n)]
        coeff, varmatrix = curve_fit(gauss, bincenters,n,p0=p0)

        # Plot the histogram
        pyplot.plot(bincenters, n, 'o')

        # Fit the histogram to calculate centre and standard deviation
        hist_fit = gauss(bincenters, *coeff)

        # Plot the resulting fit
        if args.plot == True:
            pyplot.plot(bincenters,hist_fit)
            pyplot.draw()
            pyplot.show()

        # Including (or updating) sky values in the header of the image
        hdr.add_history("- Added sky value and std dev estimated from histogram of" +\
                        " image. See sky and sky_std keywords.")
        hdr["sky"] = (str(coeff[1]), "Sky value")
        hdr["sky_std"] = (str(coeff[2]), "Standard deviation of sky")
        im.flush()
        im.close()
            
    return coeff

def main(arguments=None):
    if arguments == None:
        arguments = sys.argv[1:]
    args = parser.parse_args(arguments)
    find_sky(args)



# Create parser
parser = argparse.ArgumentParser(description='Find sky for a selection of images.')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', nargs="+", type=str, \
                   help='list of input images for which to estimate sky.')
parser.add_argument("--plot", dest='plot', action='store_true', default=False, \
                    help='Plot histogram of image to find sky.')
parser.add_argument("--mask_key", metavar='mask_key', action='store', 
                    dest='mask_key', default="",
                   help='Key where the name of the mask image is stored in the header.')


if __name__ == "__main__" :
        main()





