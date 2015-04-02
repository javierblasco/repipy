import utilities as utils
import os
import re
import numpy
import repipy.target as target
from repipy import __path__ as repipy_path
repipy_path = repipy_path[0]

class Filter(object):
    def __init__(self, header):
        self.header = header

    def zero_point(self, target):
        """ Return the zero point, given this target """
        if target.objtype == 'standard':
            return  2.5 * (numpy.log10(target.counts / self.header.hdr[self.header.exptimek]) - numpy.log10(target.flux))

    @property
    def filter_ID(self):
        """ Identify the filter ID.

        Some telescopes identify the filter ID in the header of the images. Unfortunately, the filters in
        CAHA don't use the same ID in the filter curves in the webpages and in the headers, so it's difficult to link
        them. """
        if self.header.telescope.lower() == "caha":
            return self.header.telescope.lower() + "_" + str(self.filter_wavelength) + "." + str(self.filter_width)
        else:
            if self.header.telescope and self._get_filterID():
                return self.header.telescope.lower() + "_" + str(self._get_filterID())

    @property
    def filter_sys(self):
        """ Return the filter system of the filter: Gunn, SDSS, Johnson, Harris, ...
        :return:
        """
        return self._get_filterSYS()

    @property
    def filter_name(self):
        """ Give the name of the filter in the header.
        :return:
        """
        return self.header.hdr[self.header.filterk]

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
        dir = os.path.join(repipy_path, 'filters')
        file = os.path.join(dir, self.filter_ID)
        wav, trans = numpy.genfromtxt(file).transpose()
        # In case the transmissivity is in % instead of normalized to 1
        if trans.max() > 1:
            trans /= 100
        # In case the wavelength is in nanometers, not Angstroms
        if wav.max() < 1000:
            wav *= 10
        return numpy.array([wav, trans]).transpose()

    @property
    def filter_integral(self):
        """ Integral of the filter curve"""



    @utils.memoize
    def _get_filterwav(self):
        """ Try to find the central wavelength of the filter used for the image.

        I know, right? Crazy to uniquely identify a filter in the header of the image...
        """

        wavelength = self.header._get_value(self.header._KEYWORDS_ALIASES['FILTER_WAVELENGTH'])
        width = self.header._get_value(self.header._KEYWORDS_ALIASES['FILTER_WIDTH'])
        return wavelength, width

    @utils.memoize
    def _get_filterID(self):
        """ Try to find the ID of the filter used in the image.

        """
        return self.header._get_value(self.header._KEYWORDS_ALIASES['FILTER_ID'])

    @utils.memoize
    def _get_filterSYS(self):
        """ Find the filter system: SDSS, Gunn, Johnson, Harris
        :return:
        """
        filter_sys = self.header._get_value(self.header._KEYWORDS_ALIASES['FILTER_SYS'])
        filter_name = self.filter_name

        system_dict = {'Har': '.*har(ris)?.*',
                       '': '.*(H(a(lpha)?)?|H(a)?)\d{4}.*',
                       'Joh': '.*j(oh(nson)?)?.*',
                       'Gunn': '.*gun(n)?.*',
                       'sdss': '.*(sdss|sloan).*'
                       }
        # The order of the filter checking is important, if 'Har' is present, the filter is Harris, not Halpha
        # and if sloan is pressent, the filter is SDSS, no matter if the header says "sloan Gunn r" for historical
        # reasons. Therefore, we sort the search to beat those problems.
        order = ['Har', '', 'sdss', 'Gunn', 'John']
        for system in order:
            regexp = system_dict[system]
            string = filter_sys + filter_name
            if re.match(regexp, string, re.I):
                return system
        return ''