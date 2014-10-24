import utilities as utils
import os
import numpy

class Filter(object):
    def __init__(self, header):
        self.header = header

    @property
    def filter_ID(self):
        """ Identify the filter ID.

        Some telescopes identify the filter ID in the header of the images. Unfortunately, the filters in
        CAHA don't use the same ID in the filter curves in the webpages and in the headers, so it's difficult to link
        them. """
        if self.header.telescope.lower() == "caha":
            return self.header.telescope.lower() + "_" + str(self.filter_wavelength) + "." + str(self.filter_width)
        else:
            return self.header.telescope.lower() + "_" + str(self._get_filterID())

    @property
    def filter_wavelength(self):
        """ IF the wavelength of the filter is present, try to find it.

        The wavelengths are usually in nanometers with one decimal digit.
         """
        wav = self._get_filterwav()[0]
        if wav: # If not None
            wav = int(round(float(wav)))
        return wav

    @property
    def filter_width(self):
        """ Find width of filter. """

        width = self._get_filterwav()[1]
        if width: # If not None
            width = int(round(float(width)))
        return width

    @property
    def filter_curve(self):
        """ Read the filter curve from the collection in the repipy/filters folder"""
        dir = '/home/blasco/Desktop/librerias_python/repipy/filters/'
        file = os.path.join(dir, self.filter_ID)
        wav, trans = numpy.genfromtxt(file).transpose()
        # In case the transmissivity is in % instead of normalized to 1
        if trans.max() > 1:
            trans /= 100
        return numpy.array([wav, trans]).transpose()

    @utils.memoize
    def _get_filterwav(self):
        """ Try to find the central wavelength of the filter used for the image.

        I know, right? Crazy to uniquely identify a filter in the header of the image...
        """

        wavelength = self.header._get_value(self.header._KEYWORDS_ALIASES['FILTER_WAVELENGTH'])
        width = self.header._get_value(self.header._KEYWORDS_ALIASES['FILTER_WIDTH'])
        return wavelength, width