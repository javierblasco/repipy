import repipy.header as header
import repipy.target as target
import repipy.filter as imfilter
import astropy.io.fits as fits
import astropy.wcs as wcs
import repipy.imstats as imstats

class chip(object):
    """ Each of the CCD Header/Data Units (HDUs) of an astronomical image
    """
    def __init__(self, HDU):
        self.data = HDU.data
        self.header  = header.Header(HDU.header)
        self.mask = None
        self.wcs = self._get_wcs()

    def _get_wcs(self):
        """ Get the World Coordinate System solution from the header, if present.

        The WCS module of astropy gives always a solution. If there is no WCS information in the header, it will return
        False in the attribute "is_celestial", and return a set of parameters, which are essentially equivalent to
        assuming that the image is pointed at Ra, Dec = 0,0  and with a pixel scale of 1 degree per pixel. Whenever that
         happens we will return None instead.
        """

        w = wcs.WCS(self.header.hdr)
        if not w.is_celestial:
            w = None
        return w

    def has_coords(self, ra, dec):
        """ Check if given coordinates are within the chip boundaries.
        """

        # With no valid wcs in the header, the rest is useless
        if not self.wcs:
            return None

        ly, lx = self.data.shape
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        if (0 < x < (lx-1)) and (0 < y < (ly-1)):
            return True
        else:
            return False



class Astroim(object):
    """
    """
    def __init__(self, image):
        self.im_name = image
        self._HDUList = fits.open(self.im_name, 'readonly')
        self.header = header.Header(self._HDUList[0].header)
        self.filter = imfilter.Filter(self.header)
        self.target = target.Target(self.header, self.filter)
        self.imstats = imstats.Imstats(self)

    def zero_point(self, aperture=None):
        return self.filter.zero_point(self.target, aperture=aperture)

    @property
    def data(self):
        """ Create an array (or a list of arrays, if several chips are present) with all the data in the image

        For each of the objects in the fits image, collect the associated data, if present. At the end, if the image
        contains only one chip, do not return an list of a single element, return simply the array with the data.

        :return:
        """
        data = [ii.data for ii in self._HDUList if ii.data is not None]
        if len(data) == 1:
            data = data[0]
        return data
