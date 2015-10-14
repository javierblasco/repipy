import numpy
from sklearn import linear_model
from sklearn.cross_validation import train_test_split


class Imstats(object):
    """ Do statistics in an astronomical image with or without an associated mask.

    We will do several statistics (mean, median, mode, fit to plane, statistical tests...) in an astronomical 2D image.
    A mask can be used: the mask should be another fits image, with 0 in the non-masked pixels, and any integer
    different from 0 in those pixels not to be used. The mask, if it exists, must be saved into a keyword in the header.
    That keyword is, by default, "MASK", but if the keyword has any other name, you must pass it when creating the
    Imstats object.

    By default, the whole image is used to do statistics. If you want a square section of the image, you need
    to provide region = ((x0, y0), (x1, y1)), where (x0,y0) and (x1, y1) the coordinates of two non-consecutive
    corners of the square, with Y corresponding to the NAXIS1 in the header and X to the NAXIS2 of the header.
    Notice that this convention is the one numpy arrays use, and it has flipped axes with respect to the way DS9 and
    matplotlib shows data, i.e., the point (100, 300) in ds9/matplotlib will be (300,100) in this convention.
    If you prefer a circular aperture, just give region = ((x0, y0), radius), again following the convention of Y
    for NAXIS1 and X for NAXIS2.
    """

    def __init__(self, astroim, mask_keyword="MASK", region=None):
        self.im_name = astroim.im_name
        self.region = region
        mask_name =  astroim.header._get_value("MASK")
        self.im_mask = fits.open(mask_name)
        self.im_data = astroim.data
        self.user_mask = self._get_user_mask()
        self.final_mask = numpy.logical_or(self.im_mask, self.user_mask)
        self.masked_data = numpy.ma.array(self.im_data, mask=self.final_mask)

    def __str__(self):
        """ Print the whole statistical set.
        """
        variables = [self.im_name, self.mean(), self.std(), self.median(), self.MAD()]
        output_message = "Statistics of image {0}: \n Mean = {1}, Std = {2}, Median = {3}, MAD = {4}".format(*variables)
        return output_message

    def _get_user_mask(self):
        """ From the coordinates passed by the user, build a mask of valid pixels within the region provided.

        If user passed ( (x0, y0), (x1, y1) ) type of region the region should be a rectangle, if
        """
        # Default: mask nothing!
        mask = numpy.zeros_like(self.im_mask)
        if self.region is None:
            return mask

        try:
            ((x0, y0), (x1, y1)) = self.region
            mask = self._rectangular_mask()
        except TypeError:
            pass

        try:
            ((x0, y0), r0) = self.region
            mask = self._circular_mask()
        except TypeError:
            pass

        return mask

    def _rectangular_mask(self):
        """ From two points representing two non-consecutive vertices of a rectangle, return a mask

            The points within the rectangle will have value zero, outside the rectangle the mask will have ones.
            :return: ndarray with the mask
        """
        ((x0, y0), (x1, y1)) = self.region
        if x0 > x1:
            x0, x1 = x1, x0
        if y0 > y1:
            y0, y1 = y1, y0

        mask = numpy.ones_like(self.im_mask)
        mask[x0:x1+1, y0:y1+1] = 0
        return mask

    def _circular_mask(self):
        """ From a point in the image and a radius, build a circular mask

            Inside the circle, the values are zero, outside the mask has ones.
            :return:  ndarray with the mask
        """
        (x0, y0), radius = self.region
        lenx, leny = self.im_mask.shape
        y, x = numpy.mgrid[-x0:lenx-x0, -y0:leny-y0]

        mask = numpy.ones_like(self.im_mask)
        mask[x * x + y * y <= radius * radius] = 0
        return mask

    def mean(self):
        """ Calculate the mean of the non-masked area of an image
        """
        return numpy.ma.mean(self.masked_data)

    def std(self):
        """ Calculate the standard deviation of the non-masked area of an image
        """
        return numpy.ma.std(self.masked_data)

    def median(self):
        """  Calculate the median of the non-masked area of an image
        """
        return numpy.ma.median(self.masked_data)

    def MAD(self):
        """  Calculate the median absolute deviation of the non-masked area of an image
        """
        return numpy.ma.median( numpy.ma.abs(self.masked_data - self.median()) )

    def fit_plane(self, model="linear"):
        """  Fit a plane to the non-masked area of an image.
        :param model: the algorithm used to fit: at the moment only "Ridge" and "Linear" are valid, for the Ridge and
                      Linear Regression, respectively.
        :return:  the model of the fitted plane.
        """
        lenx, leny = self.im_data.shape
        y, x = numpy.mgrid[0:lenx, 0:leny]

        # Select non-masked elements only
        x = x[self.final_mask == 0]
        y = y[self.final_mask == 0]
        data = self.im_data[self.final_mask == 0]

        # Do the linear regression fit
        X = numpy.array([x, y]).T
        Y = data

        # Split train/test samples with a 80% / 20% split
        X_train, X_test, Y_train, Y_test =  train_test_split(X, Y, test_size=0.2)

        if model.lower() == "ridge":
            model = linear_model.RidgeCV(alphas = [0.001, 0.01, 0.1, 1])
        elif model.lower() == "linear":
            model = linear_model.LinearRegression()

        model.fit(X_train,Y_train)

        # Check residuals with the test sample and print the results
        residuals = numpy.abs( model.predict(X_test) - Y_test)
        output_coefficients = list(model.coef_) + [model.intercept_, numpy.median(residuals), numpy.mean(residuals)]
        output_message = """
                        Fit to a plane:
                               counts = {0:3f} * X + {1:3f} * Y + {2:3f}.
                               Median residual: {3:3f}
                               Mean residual: {4:3f}
                         """.format(*output_coefficients)
        print output_message
        return model
