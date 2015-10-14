import numpy
from sklearn import linear_model
from sklearn.cross_validation import train_test_split


class Imstats(object):
    """ Point statistics, statistical fits and tests for images.

    By default, the whole image is used to do statistics, but the initializer accepts a parameter region to define
    circular or rectangular regions of the image to be used instead. The regions must be an iterable of lists/tuples
    with three elements (x0, y0, radius) for circular regions or four elements (x0, y0, x1, y1)
    """

    def __init__(self, data, mask=None):
        self.im_data = numpy.ma.array(data, mask=mask)

    def __str__(self):
        """ Print some statistical info.
        """
        variables = [self.mean(), self.std(), self.median(), self.mad()]
        output_message = "Statistics of image {0}: \n Mean = {1}, Std = {2}, Median = {3}, MAD = {4}".format(*variables)
        return output_message

    def mean(self):
        """ Calculate the mean of the non-masked area of an image
        """
        return numpy.ma.mean(self.im_data)

    def std(self):
        """ Calculate the standard deviation of the non-masked area of an image
        """
        return numpy.ma.std(self.im_data)

    def median(self):
        """  Calculate the median of the non-masked area of an image
        """
        return numpy.ma.median(self.im_data)

    def mad(self):
        """  Calculate the median absolute deviation of the non-masked area of an image
        """
        return numpy.ma.median(numpy.ma.abs(self.im_data - self.median()))

    def fit_plane(self, model="linear"):
        """  Fit a plane to the non-masked area of an image.
        :param model: the algorithm used to fit: at the moment only "Ridge" and "Linear" are valid, for the Ridge and
                      Linear Regression, respectively.
        :return:  the model of the fitted plane.
        """
        lenx, leny = self.im_data.shape
        y, x = numpy.mgrid[0:lenx, 0:leny]

        # Select non-masked elements only
        x = x[self.im_data.mask == 0]
        y = y[self.im_data.mask == 0]
        data = self.im_data[self.im_data.mask == 0]

        # Do the linear regression fit
        X = numpy.array([x, y]).T
        Y = data

        # Split train/test samples with a 80% / 20% split
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

        if model.lower() == "ridge":
            model = linear_model.RidgeCV(alphas=[0.001, 0.01, 0.1, 1])
        elif model.lower() == "linear":
            model = linear_model.LinearRegression()

        model.fit(X_train, Y_train)

        # Check residuals with the test sample and print the results
        residuals = numpy.abs(model.predict(X_test) - Y_test)
        output_coefficients = list(model.coef_) + [model.intercept_, numpy.median(residuals), numpy.mean(residuals)]
        output_message = """
                        Fit to a plane:
                               counts = {0:3f} * X + {1:3f} * Y + {2:3f}.
                               Median residual: {3:3f}
                               Mean residual: {4:3f}
                         """.format(*output_coefficients)
        print output_message
        return model
