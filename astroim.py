import repipy.header as header
import repipy.target as target
import repipy.filter as imfilter
import repipy.utilities as utilities
import astropy.io.fits as fits
import astropy.wcs as wcs
import repipy.imstats as imstats
import sys

class Chip(object):
    """ Each of the CCD Header Data Units (HDUs) of an astronomical image which contains data (i.e., no main headers)
    """
    def __init__(self, hdu, mask=None):
        self.data = hdu.data
        self.mask = mask
        self.header = header.Header(hdu.header)
        self.wcs = self._get_wcs()

    @property
    def imstats(self):
        """ Return an object with some basic statistical information about the image.
        :return: repipy's Imstat object
        As it is a @property, any change in the data or the mask will be reflected in the Imstats object.
        """
        return imstats.Imstats(self.data, self.mask)

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
        ly, lx = self.data.shape
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        if (0 < x < (lx-1)) and (0 < y < (ly-1)):
            return True
        else:
            return False

class Astroim(object):
    """
    """
    def __init__(self, image, memmap=False):
        self.memmap = memmap
        self.im_name = image
        self.HDUList = fits.open(self.im_name, memmap=self.memmap)
        self.primary_header = self._get_primary_header()
        self.mask_name = self.primary_header.get("MASK")
        self.HDUList_mask = self._read_mask()
        self.chips = self._get_chips()
        self.filter = imfilter.Filter(self.primary_header)
        self.target = target.Target(self.primary_header, self.filter, self.chips)

    def __iter__(self):
        """ Use an iterator to go through all the chips of the image.  """
        return iter(self.chips)

    def _read_mask(self):
        """ Read the image under the keyword "MASK" in the header, if present.

        :return: the astropy HDUList of the image of None
        """
        mask = None
        if self.mask_name:
            try:
                mask = fits.open(self.mask_name, memmap=self.memmap)
            except IOError:
                mssg = "\n In the header of image {1} appears image {0} as the mask. The image is not present! " \
                       "Remove the keyword or put the image where it should be!  ".format(self.mask_name, self.im_name)
                sys.exit(mssg)
        return mask

    def _copy_wcs_to_mask(self):
        """ Copy all the WCS-related keywords from one list of HDUs to another.

        :param hdu_list_NoWCS: HDUList without WCS, so destination HDUList
        :param hdi_list_WCS:  HDUList with WCS, this will be the origin of the WCS copy
        :return: HDUList with the resulting target HDUList
        """
        for hdu_origin, hdu_target in zip(self.HDUList, self.HDUList_mask):
            hdu_target.header = utilities.copy_WCS(hdu_target.header, hdu_origin.header)
        return None

    def zero_point(self, aperture=None):
        return self.filter.zero_point(self.target, aperture=aperture)

    def _get_primary_header(self):
        """ Return the main header of a fits file.

        A main header is the header with no data unit associated. It usually contains data about the telescope and
        observations, with the chip information being stored in their own headers. If not present in our data, the
        first header (Primary header) is passed, since it will contain all the info available about the general
        conditions.
        :return: main header, if present, otherwise the first header in the fits file.
        """
        for hdu in self.HDUList:
            if hdu.data is None:
                primary_header = hdu.header
                break
        else:
            primary_header = self.HDUList[0].header
        return header.Header(primary_header, im_name=self.im_name)

    def _get_chips(self):
        """ Build Chip class objects for all the chips present in the image
        """
        if self.HDUList_mask:
            return [Chip(hdu, msk.data) for hdu, msk in zip(self.HDUList, self.HDUList_mask) if hdu.data is not None ]
        else:
            return [Chip(hdu) for hdu in self.HDUList if hdu.data is not None]

    def write(self, output_name=None):
        """Write output to a fits file.
        :param output_name: name of the output fits file
        :return: None
        """
        if output_name is None:
            output_name = self.im_name

        self.HDUList.writeto(output_name, clobber=True)
        return None
