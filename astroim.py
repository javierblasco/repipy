import repipy.header as header
import repipy.target as target
import repipy.filter as imfilter
import astropy.io.fits as fits
import astropy.wcs as wcs
import repipy.imstats as imstats
import numpy


class Chip(object):
    """ Each of the CCD Header Data Units (HDUs) of an astronomical image which contains data (i.e., no main headers)
    """
    def __init__(self, hdu, mask=None):
        self.im_data = numpy.ma.array(hdu.data, mask=mask)
        self.header = header.Header(hdu.header)
        self.wcs = self._get_wcs()

    @property
    def imstats(self):
        """ Imstat object, as a @property to allow for changes in imstats following those in mask and data
        :return:
        """
        return imstats.Imstats(self.im_data.data, self.im_data.mask)

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

    def contains_coords(self, ra, dec):
        """ Check if given coordinates are within the chip.
        """

        # With no valid wcs in the header, the rest is useless
        if not self.wcs:
            return None

        ly, lx = self.im_data.data.shape
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
        self._HDUmask = self._getmask()
        self.primary_header = self._get_primary_header()
        self.chips = self._get_chips()
        self.filter = imfilter.Filter(self.primary_header)
        self.target = target.Target(self.primary_header, self.filter)


    def __iter__(self):
        """ Iterate over all the chips of the image.  """
        return iter(self.chips)

    def _get_primary_header(self):
        """ Return the main header of a fits file.

        A main header is the header with no data unit associated. It usually contains data about the telescope and
        observations, with the chip information being stored in their own headers. If not present in our data, the
        first header (Primary header) is passed, since it will contain all the info available about the general
        conditions.
        :return: main header, if present, otherwise the first header in the fits file.
        """
        main_header = [hdu.header for hdu in self._HDUList if not hdu.data]
        if not main_header:
            main_header = self._HDUList[0].header
        return main_header

    def _get_chips(self):
        """ Build Chip class objects for all the chips present in the image
        """
        hdu_with_data = [hdu for hdu in self._HDUList if hdu.data is not None]
        hdu_masks = [mask for mask in self._HDUmask]
        chip_objects = []
        for hdu, mask in zip(hdu_with_data, hdu_masks):
            if hdu.data is not None:
                chip_objects.append(Chip(hdu, mask.data))
        return chip_objects

    def _getmask(self):
        """ Get the mask, if any is present in the header, for the image.

        :return: the mask if there is a fits file under the keyword MASK, None otherwise
        """
        try:
            return fits.open(self.header.get("MASK"))
        except IOError:  # Image not found
            return None

    def zero_point(self, aperture=None):
        return self.filter.zero_point(self.target, aperture=aperture)
