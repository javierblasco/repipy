#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Created on Fri Aug  9 16:44:52 2013

@author: blasco
"""
import sys
import os
import argparse
import astropy.io.fits as fits
import numpy
import scipy
import repipy.utilities as utils
from scipy import ndimage
from scipy import optimize
from scipy.ndimage.filters import median_filter
import matplotlib.pyplot as plt

from skimage.util import img_as_ubyte
from skimage.draw import circle_perimeter
from skimage.feature import peak_local_max, canny
from skimage.transform import hough_circle
from skimage import data, color

from repipy import astroim

"""  Program to mask fits images.
+
+    For the moment, three criteria are used to mask pixels:
+
+        - minmax:  masks anything below a given minimum value or above a maximum one
+        - circular:  masks anything outside a circular area, that is detected from the image. This method is implemented
+                     because some telescopes have their circular Field of View within a square image, and the part
+                     outside the circle has not been exposed.
+        - stars:    Masks stars from the image. This is useful, for example, when you want to detect the sky and/or
+                    exclude possible stars from a flat image.
+    This program will open the image, mask according to the criteria explained above, save the mask and put the name
+    of the mask fits file into the header, under the keyword MASK. If a WCS is present, it will be copied into the mask
+    image.
+
+    Assumptions:
+        - This routine assumes that, if a multi-layer image is given, the first layer will contain a Primary HDU, with
+          the main header of the image (and possibly no data associated). It is in this layer where the MASK keyword
+          will be stored.
+
+"""

def gauss(x, *p):
    A,mu,sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))    

def apply_sobel_filter(image):
    """ Apply sobel filter on an image, return the filtered object. 
        This routine roughly follows the solution provided in:
        http://stackoverflow.com/questions/7185655/applying-the-sobel-filter-using-scipy         
    """
    dx = ndimage.sobel(image, 0)  # horizontal derivative
    dy = ndimage.sobel(image, 1)  # vertical derivative
    mag = numpy.hypot(dx, dy)     # magnitude
    return mag
    
def zero_edges(image, edge=5):
    """ Return same image with the edges reset to zero """
    image[0:edge, :] = 0.
    image[:, 0:edge] = 0.
    image[image.shape[0]-edge:, :] = 0.
    image[:, image.shape[1]-edge:] = 0.
    return image    

def detect_circular_FoV(data, args):
    ''' A Sobel edge detection algorithm is used to detect a sharp circular edge 
        within a rectangular 2D numpy array. The array is the only input, while
        the centre of the image and the radius of the circle are provided as 
        outputs.     
    '''
    mag = apply_sobel_filter(data)  # Filter image

    # Exclude edges of image. Sobel uses a filter of size 3.  
    mag = zero_edges(mag, edge=3)
    
    # Find those values in the highest 0.5%, the sharpest edge
    percentile = numpy.percentile(mag, args.contrast)
    indices = numpy.where((mag > percentile))
    x,y = indices
           
    # Now fit resulting points to a circle 
    xc, yc, radius, radius_MAD = fit_to_circle(x, y)

#    # Create an image of those points that were used to fit 
#    fitted_points = mag * 0.
#    fitted_points[indices] = 1
#    if os.path.isfile("fitted_points.fits"):
#        os.remove("fitted_points.fits")   
#    fits.writeto("fitted_points.fits", fitted_points )


    # If radius_MAD > 5% of the radius, data was not originally a circle
    if (radius_MAD / radius * 100) < 5:
        result = xc, yc, radius 
    else:
        result = None
    return result

def fit_to_circle(x, y, xc=None, yc=None):
    """ Fit to a circle using a variant from the method shown by the scipy 
        cookbook:
        http://wiki.scipy.org/Cookbook/Least_Squares_Circle """
    if not xc:
        estimate = numpy.median(x), numpy.median(y)  # first guess for centre

    ii = 0
    # until convergence
    while True:    
        # optimize centre 
        (xc, yc), ier = optimize.leastsq(f_2, estimate, args=(x,y)) # fitted xc, yc

        # calculate radii of points, median for best value, median absolute 
        # deviation (MAD) for error estimates.
        radii_fit = calc_R(x, y, xc, yc)   
        radius = numpy.median(radii_fit)   
        radius_MAD =  numpy.median(numpy.abs(radii_fit-radius)) 

        # Use only those values within 5 times the median absolute deviation. 
        # For a Gaussian distribution this would be > 3 sigma. 
        whr_good = numpy.where((radii_fit > radius - 5 * radius_MAD) &
                               (radii_fit < radius + 5 * radius_MAD))[0]
        # convergence if all values are good                       
        if len(whr_good) == len(radii_fit):
            break
        else:
            x = x[whr_good]
            y = y[whr_good]
            ii += 1
    return xc, yc, radius, radius_MAD        
            

def calc_R(x, y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return numpy.sqrt((x-xc)**2 + (y-yc)**2)
   
def f_2((xc, yc), x, y):
    """ calculate the algebraic distance between the data points and the mean 
        circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, xc, yc)
    return Ri - numpy.mean(Ri)

def mask_circle(image, xc, yc, radius, value=0):
    ''' Mask with value (defect, value=0) a circle within an image '''
    lx, ly = image.shape
    # Create an array like image with the separation from xc of each pixel.
    # Then, same for y axis. 
    x = numpy.arange(lx).reshape(lx,1) * numpy.ones(ly) - xc
    y = numpy.arange(ly) * numpy.ones(lx).reshape(lx,1) - yc
    
    # Calculate radius
    r = numpy.sqrt(x**2 + y**2)
    image[numpy.where(r > radius)] = value
    return image


def detect_Hough(data):
    image = data.copy()
    edges = canny(image, sigma=10, low_threshold=60, high_threshold=90)

    # Detect circles between 80% and 100% of image semi-diagonal
    lx, ly = data.shape
    sizex, sizey = lx/2., ly/2.
    max_r = numpy.sqrt(sizex**2 + sizey**2) 
    hough_radii = numpy.linspace(0.5*max_r, 0.9 * max_r, 20)
    hough_res = hough_circle(edges, hough_radii)


    centers = []
    accums = []
    radii = []
    for radius, h in zip(hough_radii, hough_res):
        # For each radius, extract two circles
        num_peaks = 2
        peaks = peak_local_max(h, num_peaks=num_peaks)
        centers.extend(peaks)
        accums.extend(h[peaks[:, 0], peaks[:, 1]])
        radii.extend([radius] * num_peaks)

    # Use the most prominent circle
    idx = numpy.argsort(accums)[::-1][:1]
    center_x, center_y = centers[idx]
    radius = radii[idx]
    return center_x, center_y, radius

def get_mesh(data):
    " Get a mesh grid the same shape of the data"
    lx, ly = data.shape
    x = numpy.linspace(0, ly, ly)
    y = numpy.linspace(0, lx, lx)
    xv, yv = numpy.meshgrid(x, y)
    return xv, yv

def get_center(data):
    " Find the centroid of the image"
    xv, yv = get_mesh(data)
    cy = numpy.sum( data * yv ) / numpy.sum(data) 
    cx = numpy.sum( data * xv ) / numpy.sum(data)
    return cy, cx

def get_radius(data):
    rr, avg = means(data)
    der_avg = derivative(rr, avg)
    whr = der_avg.argmin()
    #plt.plot( (rr[0:-1] + rr[1:]) / 2, der_avg, 'o')
    #plt.show()
    return rr[whr]

def distance(data):
    " For each pixel, tell me how far it is from the centre"
    xv, yv = get_mesh(data)
    cy, cx = get_center(data)
    distance = numpy.sqrt( (xv - cx)**2 + (yv - cy)**2 )
    return distance 


def means(data):
    nbins = 40
    dist = distance(data)
    radii = numpy.linspace(1, dist.max(), nbins+1)  
    rr, avg, std, suma = numpy.zeros(nbins), numpy.zeros(nbins), numpy.zeros(nbins), numpy.zeros(nbins)
    for ii in range(0, nbins):
        rr[ii] = (radii[ii] + radii[ii+1]) / 2.    
        whr = numpy.where( (radii[ii] <= dist) & (dist <= radii[ii+1]) )
        avg[ii] = numpy.mean(data[whr])    
    return rr, avg


def derivative(x,y):
    return (y[1:] - y[0:-1]) / ( x[1:] - x[0:-1])

def cutre_detect(data):
    
    yc, xc = get_center(data)
    radius = get_radius(data) 
    return yc, xc, radius


def mask(args):
    for image in args.image:
        im = fits.open(image, mode='update')
        data = im[0].data.astype(numpy.float64)
        header = im[0].header
        mask = numpy.ones(data.shape, dtype=numpy.int) * args.true_val #create mask
        
        # If circular field of view within rectangular image:
        if args.circular:
            result = cutre_detect(data)
            if result:
                xc, yc, radius = result
                radius = radius - args.margin   # avoid border effects
                mask = mask_circle(mask, xc, yc, radius, value=args.outside_val)
    
        # Maxval, minval masking
        bad_pixels = numpy.where((data < args.minval) | (data > args.maxval))
        mask[bad_pixels] = args.false_val
        
        # Star masking fitting sky
        if args.stars:  # if stars in the image
            unmasked = data[mask == 0].flatten()
	    try:         
                n, bins = numpy.histogram(unmasked, bins=range(int(min(unmasked)),
                                                               int(max(unmasked)),50))
	    except TypeError:
                print "Error in image: ", image
                raise

            bincenters = 0.5*(bins[1:]+bins[:-1])   
            max_pos = n.argmax()
            max_value = bincenters[max_pos]  #value of the sky
            p0 = [n[max_pos], max_value, 50]
            coeff, varmatrix = optimize.curve_fit(gauss, 
                                                  bincenters[max_pos-5:max_pos+5],
                                                  n[max_pos-5:max_pos+5],p0=p0)
            max_sky = coeff[1] + 3. * coeff[2] # 2 sigma above sky level
            mask[data > max_sky] = 1
        
        # Save mask image
        maskname = args.output
        if not maskname:
            maskname = image + ".msk"
        if os.path.isfile(maskname):
            os.remove(maskname)
        maskname = os.path.abspath(maskname)    
        fits.writeto(maskname, mask)

        # Add message to image header
        header.add_history("- Created mask of image, see mask keyword")
        header[args.mask_key] = (maskname, "Mask of original image")
        im.flush()
        im.close()        
            
def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]
  args = parser.parse_args(arguments)
  
  # True_val and false_val can not be the same
  if args.true_val == args.false_val:
      sys.exit("\n\n ERROR: true_val and false_val are the same value: " + \
               str(args.false_val) + " Use --true_val and --false_val \n\n")

  
  # Call combine, keep name of the file created
  masknames = mask(args)
  return masknames   



# Create parser
parser = argparse.ArgumentParser(description='Create masks for images')
# Add necessary arguments to parser
parser.add_argument("image", metavar='image', action='store',\
                    help='Image(s) from which to create masks.', nargs='+')
parser.add_argument("--max_val", metavar="maxval", dest='maxval', action='store',\
                    default=50000, type=float, help='Maximum allowed value. '+\
                    'Above this value, mask out. Default: 50000.')
parser.add_argument("--min_val", metavar="minval", dest='minval', default=0, \
                   type=float, action='store', help='Minimum allowed value. '+\
                   'Below this value, mask out. Default: 0.')
parser.add_argument("--output", metavar="output", dest='output', action='store',\
                    default='', help='Name of the output mask. ' +\
                    'Default: same as input images, but ending in .fits.msk.')
parser.add_argument("--circular", action="store_true", dest="circular", \
                   default=False, help=' Use if the field of view is circular, '+\
                   ' while the image is a rectangle. Mask is set to zero_value '+\
                   'the circle. ')
parser.add_argument("--stars", action="store_true", dest="stars", \
                   default=False, help=' Use if you want to mask stars in '+\
                   ' an image. ')                  
parser.add_argument("--true_val", metavar="true_val", dest='true_val', default=0, \
                   type=int, action='store', help='Value for the VALID points. '+\
                   'those you DO NOT want to mask out. Default: 0 ')
parser.add_argument("--false_val", metavar="false_val", dest='false_val', default=1, \
                   type=int, action='store', help='Value for the INVALID points. '+\
                   'those you DO want to mask out. Default: 1 ')
parser.add_argument("--outside_val", metavar="outside_val", dest="outside_val", 
                    type=int, action="store", default=2,
                    help="If --circular is used, this is the value to be used "+\
                    "to mask out the points outside the circular FoV. Default:2 ")
parser.add_argument("--margin", metavar="margin", dest='margin', default=10, \
                   type=int, action='store', help='Margin around the edges for '+\
                   'which the mask is set to zero. If --circular is used '+\
                   'this margin will be reduced from the calculated radius. '+\
                   'Default: 10')
parser.add_argument("--mask_key", metavar="mask_key", dest='mask_key', action='store',\
                    default='mask', help='Name of the keyword in the header ' +\
                    'that contains the name of the mask. Default: mask')
parser.add_argument("--contrast", metavar="contrast", dest="contrast", action='store',\
                    default=99.5, type=float,\
                    help="When a Sobel filter is applied, the " +\
                    "result is an image with the contrast of the pixels "+\
                    "respect the others in a box of 3x3. --contrast is set "+\
                    "as the minimum contrast to be taken into account as "+\
                    "the detection of a border. Default: 99.5, i.e. top 0.5. "+\
                    "percent. This value is good for flats, images with long "+\
                    "exposures... This method will not work for images with "+\
                    "less contrast, between the exposed area and the bias "+\
                    "level, such as short exposures with narrow filters")
                   

if __name__ == "__main__":
    main()
