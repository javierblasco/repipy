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
            if hdu.data is not None:
                HDUList_mask[ii].data = numpy.zeros_like(hdu.data)
        #Include WCS into the headers
        HDUList_mask = self._copy_wcs_to_mask(HDUList_mask, self.HDUList)
        return HDUList_mask

    def _copy_wcs_to_mask(self, hdu_list_NoWCS, hdu_list_WCS):
        """ Copy all the WCS-related keywords from one list of HDUs to another.

        :param hdu_list_NoWCS: HDUList without WCS, so destination HDUList
        :param hdi_list_WCS:  HDUList with WCS, this will be the origin of the WCS copy
        :return: HDUList with the resulting target HDUList
        """
        for hdu_origin, hdu_target in zip(hdu_list_WCS, hdu_list_NoWCS):
            hdr_target = utilities.copy_WCS(hdr_target, hdu_origin)
        return hdu_list_NoWCS



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

        return [Chip(hdu, msk.data) for hdu, msk in zip(self.HDUList, self.HDUList_mask) if hdu.data is not None ]

    def write(self, output_name=None, mask_name=None, no_mask=False):
        """Write output to a fits file.
        :param output_name: name of the output fits file
        :param mask_name: name of the output mask file. If None, it will be same as input with .fits.msk
        :param no_mask: if True, no mask will be written.
        :return: None
        """
        if output_name is None:
            output_name = self.im_name

        if mask_name is None:
            mask_name = os.path.abspath( utilities.replace_extension(output_name, ".fits.msk") )

        self.HDUList.writeto(output_name, clobber=True)
        if no_mask is not True:
            # update and record mask
            self.HDUList_mask = self._copy_wcs_to_mask(self.HDUList_mask, self.HDU)
            self.HDUList_mask[0].header.add_comment("Mask for image {0}".format(output_name))
            self.HDUList_mask.writeto(mask_name, clobber=True)
            self.primary_header.hdr["MASK"] = (mask_name, "Name of mask image.")
        return None
