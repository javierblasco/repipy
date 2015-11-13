import repipy.header as header
import repipy.target as target
import repipy.filter as imfilter
import repipy.utilities as utilities
import astropy.io.fits as fits
import astropy.wcs as wcs
import repipy.imstats as imstats
import numpy
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
        As it is a @property, any change in the data will be reflected in the Imstats object.
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
        ly, lx = self.im_data.data.shape
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
        self.HDUList_mask = self._get_HDUList_mask()
        self.chips = self._get_chips()
        self.filter = imfilter.Filter(self.primary_header)
        self.target = target.Target(self.primary_header, self.filter, self.chips)

    def __iter__(self):
        """ Use an iterator to go through all the chips of the image.  """
        return iter(self.chips)

    def _get_HDUList_mask(self):
         """ HDUs corresponding to the mask image.

          Load the mask image from file, if it exists. If it doesn not exit, create it from scratch.
         :return: astropy HDUList
         """
         mask = self._read_mask()
         if mask is None:
             mask = self._build_mask()
         return mask

    def _read_mask(self):
        """ Read the image under the keyword "MASK" in the header, if present.

        :return: the astropy HDUList of the image of None
        """
        mask_name = self.primary_header.get("MASK")
        mask = None
        if mask_name:
            try:
                mask = fits.open(mask_name, memmap=self.memmap)
            except IOError:
                mssg = "\n In the header of image {1} appears image {0} as the mask. The image is not present! " \
                       "Remove the keyword or put the image where it should be!  ".format(mask_name, self.im_name)
                sys.exit(mssg)
        return mask

    def _build_mask(self):
        """ Using the original image as a model, build a similar image for the mask.

        :return: HDUlist with as many HDUs as the original image, with arrays of zeros wherever the original
        image had data, and with the WCS information, if the original had it.
        """
        HDUList_mask = fits.HDUList( [fits.PrimaryHDU()] + [fits.ImageHDU() for _ in self.HDUList[1:]])
        for ii, hdu in enumerate(self.HDUList):
            hdr = HDUList_mask[ii].header
            hdr = utilities.copy_WCS(hdr, hdu.header)   # copy WCS into the new header!
            if hdu.data is not None:
                HDUList_mask[ii].data = numpy.zeros_like(hdu.data)
        return HDUList_mask

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
        chip_objects = []
        with fits.open(self.im_name, memmap=self.memmap) as HDUList:
            hdu_with_data = [hdu for hdu in HDUList if hdu.data is not None]

            # Read the mask, if present in the header, create a mask of None otherwise.
            mask_name = self.primary_header.get("MASK")
            if mask_name is not None:
                hdu_masks = [mask for mask in fits.open(mask_name, memmap=self.memmap)]
            else:
                hdu_masks = [None for hdu in hdu_with_data]

            for hdu, mask in zip(hdu_with_data, hdu_masks):
                if mask:
                    chip_objects.append(Chip(hdu, mask.data))
                else:
                    chip_objects.append(Chip(hdu))
        return chip_objects
