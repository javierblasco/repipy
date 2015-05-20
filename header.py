"""



 """

import astropy.io.fits as fits
import sys
import repipy.utilities as utils
import numpy as np
import os

class NotFoundKeyword(KeyError):
    mssg = ("ERROR! None of the known keywords for {} is found in the header of image {}."
            " Add your own keyword to the file keywords_aliases.json or include in the header a valid keyword: {}")

    def __init__(self, keyword, image, valid):
        self.keyword = keyword
        self.im_name = image
        self.valid = valid

    def __str__(self):
        return self.mssg.format(self.keyword, self.im_name, ", ".join(self.valid))


class Header(object):
    """ Object with the information we can gather from the header of the image about observatory, telescope, instrument,
        and the keywords for the most important parameters: airmass, object, exposure time, ...
    """
    _KEYWORDS_ALIASES = dict(
        FILTER = ('INSFLNAM', 'FILTER', 'JAGFBAND', 'ALFLTNM', 'WFFBAND'),
        EXPTIME = ('EXPTIME',),
        OBJECT = ('OBJECT',),
        DATE = ('DATE-OBS',),
        TIME = ('TIME-OBS',"UTSTART"),
        GAIN = ('CCDSENS', 'GAIN' ),
        RON = ('CCDRON', 'READNOIS', 'READNOISE'),
        AIRMASS = ('AIRMASS',),
        FILTER_WAVELENGTH = ('INSFLWL',),
        FILTER_WIDTH = ('INSFLDWL',),
        FILTER_ID = ('ALFLTID', 'FAFTLID', 'JAGFID', 'INSFLID', 'WFFID'),
        FILTER_SYS = ('JAGFSYS', 'WFFPSYS'),
        TELESCOPE = ('TELESCOP', 'INSTRUME', 'ORIGIN', 'INSTRID'),
        SEEING = ('SEEING', 'FWHM', 'LEMON FWHM'),
        SKY = ("SKY",),
        SIGMA = ("SKY_STD", "SIGMA"),
        RA = ("RA",),
        DEC = ("DEC",)
    )

    _TELESCOPES_ALIASES =  { 'NOT' : ['ALFOSC', 'NOT'],
                             'CAHA' : ['DSAZ', 'CAFOS', 'CA-2.2', 'CAHA'],
                             'OSN' : ['OSN', 'OSN 1.5m'],
                             'JKT' : ['JKT'],
                             'INT' : ['INT']
    }


    def __init__(self, image):
        self.im_name = image
        # Some fits images will have more than one extension, the header beeing split between the 0 and the others
        hdulist = fits.open(self.im_name)
        if len(hdulist) == 1:
            self.hdr = hdulist[0].header
            self.ndetectors = 1
        else:
            self.hdr = hdulist[0].header + hdulist[1].header
            self.ndetectors = utils.number_of_chips(hdulist)

    @property
    def telescope(self):
        return self._get_telescope()[1]

    @property
    def filterk(self):
        """ Determine the keyword that keeps the name of the filter in the header."""
        return self._get_filterk()

    @property
    def airmass(self):
        return self.hdr[self.airmassk]

    @utils.memoize
    def _get_telescope(self):
        """ Try to find the telescope name from which the image comes.

         Since every observatory, telescope and/or instrument has different keywords, we need to go into some effort
         to do this search. The information of the telescope could be hidden in any/all/none of the keywords under the
         TELESCOPE key of the dictionary _KEYWORDS_ALIASES above. And the value could be any of a long list of aliases
         such as putting the instrument instead of the telescope in the keyword 'TELESCOP' (sigh). Thus the selection 
         of aliases used. 

        """

        key, value = self.find_in_header(self._KEYWORDS_ALIASES['TELESCOPE'], self._TELESCOPES_ALIASES)
        return key, value


    @utils.memoize
    def _get_filterID(self):
        """ Try to find the ID of the filter used for the image.

        I know, right? Crazy to uniquely identify a filter in the header of the image...
        """
        value = self._get_value(self._KEYWORDS_ALIASES['FILTER_ID'])
        return value

    @utils.memoize
    def _get_keyword(self, keywords):
        """ Try keywords from the dictionary until you find one that exists. Return the keyword. """
        return self._key_and_value(keywords)[0]

    @utils.memoize
    def _get_value(self, keywords):
        """ Try keywords from the dictionary until you find one that exists. Return the value. """
        return self._key_and_value(keywords)[1]


    def get(self, *keywords):
        """ Return the value of a keyword in the header.
        :param keywords: keywords whose value we want to return
        :return: value of the keyword given by user
        """
        try:
            result = []
            for ii in keywords:
                result.append(self.hdr[ii])
            if len(result) == 1:
                return result[0]
            else:
                return result
        except:
            return None


    def find_in_header(self, list_keywords, dict_patterns):
        """ Find the keyword in which elements from both lists coincide.

        """
        k, p = None, None
        dict_targets = {key: self.hdr.get(key) for key in list_keywords}
        for pattern_key, pattern_value in dict_patterns.iteritems():
            for target_key, target_value in dict_targets.iteritems():
                if target_value in pattern_value:
                    k, p =  target_key, pattern_key
        return k, p

    def _key_and_value(self, keywords):
        """ Find if any of the unique keywords is in the header. Return its value

        Look fot all the keywords in the header, if any is present and is not empty, return its value
        """
        k, v = None, None
        for key in keywords:
            if self.hdr.get(key) is not None:
                k,v = key, self.hdr.get(key)
                break
        return k, v


def _add_property(name, keyword):
    """ Dynamically add a property to the 'header' class.

    Add to the 'header' class the property 'name', whose getter function is
    header._get_keyword(), called with header._KEYWORDS_ALIASES[keyword] as
    its sole argument.

    """

    def getter(self):
        return self._get_keyword(self._KEYWORDS_ALIASES[keyword])
    setattr(Header, name, property(getter))

_add_property('filterk', 'FILTER')
_add_property('filtersysk', 'FILTER_SYS')
_add_property('airmassk', 'AIRMASS')
_add_property('exptimek', 'EXPTIME')
_add_property('objectk', 'OBJECT')
_add_property('datek', 'DATE')
_add_property('timek', 'TIME')
_add_property('seeingk', 'SEEING')
_add_property('gaink', 'GAIN')
_add_property('skyk', 'SKY')
_add_property('sigmak', 'SIGMA')
_add_property('ccdronk', 'RON')
_add_property('RAk', 'RA')
_add_property('DECk', 'DEC')
_add_property('timek', 'TIME')