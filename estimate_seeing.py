#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from astropy.modeling import models, fitting
import astropy.io.fits as fits
import astropy.wcs as wcs
import numpy as np
import matplotlib.pyplot as plt
import repipy.utilities as utilities
import argparse
import sys
import os

class Star(object):
    """ Define a star by its coordinates and modelled FWHM

        Given the coordinates of a star within a 2D array, fit a model to the star and determine its
        Full Width at Half Maximum (FWHM).The star will be modelled using astropy.modelling. Currently
        accepted models are: 'Gaussian2D', 'Moffat2D'
    """


    _GAUSSIAN2D = 'Gaussian2D'
    _MOFFAT2D = 'Moffat2D'
    _MODELS = set([_GAUSSIAN2D, _MOFFAT2D])

    def __init__(self, x0, y0, data, model_type=_GAUSSIAN2D, box=40):
        """ Instantiation method for the class Star.

        The 2D array in which the star is located (data), together with the pixel coordinates (x0,y0) must be
        passed to the instantiation method. .
        """
        self.x = x0
        self.y = y0
        self._box = box
        self._XGrid, self._YGrid = self._grid_around_star(x0, y0, data)
        self.data = data[self._XGrid, self._YGrid]
        self.model_type = model_type

    @property
    def model(self):
        """ Fit a model to the star. """
        return self._fit_model()

    @property
    def model_psf(self):
        """ Return a modelled PSF for the given model  """
        return self.model(self._XGrid, self._YGrid)

    @property
    def fwhm(self):
        """ Extract the FWHM from the model of the star.

            The FWHM needs to be calculated for each model. For the Moffat, the FWHM is a function of the gamma and
            alpha parameters (in other words, the scaling factor and the exponent of the expression), while for a
            Gaussian FWHM = 2.3548 * sigma. Unfortunately, our case is a 2D Gaussian, so a compromise between the
            two sigmas (sigma_x, sigma_y) must be reached. We will use the average of the two.
        """
        model_dict = dict(zip(self.model().param_names, self.model().parameters))
        if self.model_type == self._MOFFAT2D:
            gamma, alpha = [model_dict[ii] for ii in ("gamma_0", "alpha_0")]
            FWHM = 2. * gamma * np.sqrt(2 ** (1/alpha) -1)
        elif self.model_type == self._GAUSSIAN2D:
            sigma_x, sigma_y = [model_dict[ii] for ii in ("x_stddev_0", "y_stddev_0")]
            FWHM = 2.3548 * np.mean([sigma_x, sigma_y])
        return FWHM

    @utilities.memoize
    def _fit_model(self):
        fit_p = fitting.LevMarLSQFitter()
        model = self._initialize_model()
        p = fit_p(model, self._XGrid, self._YGrid, self.data)
        return p

    def _initialize_model(self):
        """ Initialize a model with first guesses for the parameters.

        The user can select between several astropy models, e.g., 'Gaussian2D', 'Moffat2D'. We will use the data to get
        the first estimates of the parameters of each model. Finally, a Constant2D model is added to account for the
        background or sky level around the star.
        """
        max_value = self.data.max()
        if self.model_type == "Gaussian2D":
            model = models.Gaussian2D(x_stddev=1, y_stddev=1, x_mean=self.x, y_mean=self.y, amplitude=max_value)
            # Establish reasonable bounds for the fitted parameters
            model.x_stddev.bounds = (0, self._box/4)
            model.y_stddev.bounds = (0, self._box/4)

        elif self.model_type == "Moffat2D":
            model = models.Moffat2D(x_0=self.x, y_0=self.y, gamma=2, alpha=2, amplitude=max_value)
            model.alpha.bounds = (1,6)
            model.gamma.bounds = (0, self._box/4)
        model += models.Const2D(self.fit_sky())
        model.amplitude_1.fixed=True
        print "Sky is fixed to value: ", model.amplitude_1
        return model

    def fit_sky(self):
        """ Fit the sky using a Ring2D model in which all parameters but the amplitude are fixed.
        """
        min_value = self.data.min()
        ring_model = models.Ring2D(min_value, self.x, self.y, self._box * 0.4, width=self._box * 0.4)
        ring_model.r_in.fixed = True
        ring_model.width.fixed = True
        ring_model.x_0.fixed = True
        ring_model.y_0.fixed = True
        fit_p = fitting.LevMarLSQFitter()
        return fit_p(ring_model, self._XGrid, self._YGrid, self.data).amplitude


    def _grid_around_star(self, x0, y0, data):
        """ Build a grid of side 'box' centered in coordinates (x0,y0). """
        lenx, leny = data.shape
        xmin, xmax = max(x0-self._box/2, 0), min(x0+self._box/2+1, lenx-1)
        ymin, ymax = max(y0-self._box/2, 0), min(y0+self._box/2+1, leny-1)
        return np.mgrid[int(xmin):int(xmax), int(ymin):int(ymax)]

    def plot_resulting_model(self):
        """ Make a plot showing data, model and residuals. """
        data = self.data
        model = self.model(self._XGrid, self._YGrid)
        residuals = data - model
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        plt.imshow(data, origin='lower', interpolation='nearest', vmin=data.min(), vmax=data.max())
        plt.colorbar()
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(model, origin='lower', interpolation='nearest', vmin=data.min(), vmax=data.max())
        plt.colorbar()
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(residuals, origin='lower', interpolation='nearest')
        plt.colorbar()
        plt.title("Residual")
        plt.show()


class StarField(object):
    """ Fit a model to a list of stars from an image, estimate their FWHM, write it to the image.

    To initialize the object you need the name of an astronomical image (im_name) and the name of a text file
     with the coordinates of the stars to be used to estimate the FWHM. The file must contain two columns with the
     RA and DEC of the stars, one row per star. If image pixels are given instead of RA and DEC, the wcs=False flag
     must be passed.
    """
    def __init__(self, im_name, coords_file, model_type, wcs=True):
        self.im_name = im_name
        self.im_data = fits.open(im_name)[0].data
        self.coords_file = coords_file
        self._wcs = wcs
        self.model_type = model_type

    def __iter__(self):
        """ Iterate over all the stars defined in the coords_file."""
        return iter(self._star_list)

    @property
    def star_coords(self):
        """ Read from coords_file the coordinates of the stars.

        The file must have two columns, with one star per row. If the coordinates are in (RA,DEC) we will transform
        them into image pixels.
        """
        x, y = np.genfromtxt(self.coords_file, unpack=True)
        if self._wcs:  # if x,y are not pixels but RA,DEC
            with fits.open(self.im_name, 'readonly') as im:
                w = wcs.WCS(im[0].header, im)
            y, x = w.all_world2pix(x, y, 1)
        return zip(x,y)

    @property
    def _star_list(self):
        """ Return a list of Star objects from the image data and the coordinates of the stars."""
        return [Star(x0, y0, self.im_data, self.model_type) for (x0,y0) in self.star_coords]

    @property
    def FWHM(self):
        """ Determine the FWHM (seeing) of the image. """
        return np.median([star.fwhm for star in self._star_list])

    def _write_FWHM(self):
        """ Write the FWHM to the header of the fits image."""
        with fits.open(self.im_name, 'update') as im:
            im[0].header["seeing"] = (self.FWHM, 'Seeing estimated in pixels')
        return None


def calculate_seeing(args):
    """ Program to estimate the seeing from an image and a list of estimates
    for the positions of stars.
    """
    for im_name, star_cat in zip(args.input, args.cat):
        im = StarField(im_name, star_cat, args.model, wcs=args.wcs)
        im._write_FWHM()



############################################################################


# Create parser
parser = argparse.ArgumentParser(description='Program to calculate the FWHM (seeing) of an image '
                                             'by fitting its stars. ')

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
parser.add_argument("--wcs", metavar="wcs_in", action="store", dest="wcs", type=bool,
                    default=True, help = "If coordinates are in pixels, set this flag to False, if in RA,DEC "
                                         "set to True. Default: True ")
parser.add_argument("--model", metavar="model", action="store", dest="model", default="Moffat2D",
                    help="Model to be fit to the stars: Moffat2D or Gaussian2D. Both models will contain a "
                         "background value. Default: 'Moffat2D' ")


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


