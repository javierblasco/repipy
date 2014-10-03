"""



 """

import astropy.io.fits as fits
import sys
import repipy.utilities as utils

class NotFoundKeyword(KeyError):
    mssg = ("ERROR! None of the known keywords for {} is found in the header of image {}."
            " Add your own keyword to the file keywords_aliases.json or include in the header a valid keyword: {}")

    def __init__(self, keyword, image, valid):
        self.keyword = keyword
        self.image = image
        self.valid = valid

    def __str__(self):
        return self.mssg.format(self.keyword, self.image, ", ".join(self.valid))


class header(object):
    """ Object with the information we can gather from the header of the image about observatory, telescope, instrument,
        and the keywords for the most important parameters: airmass, object, exposure time, ...
    """
    _KEYWORDS_ALIASES = dict(
        FILTER = ['INSFLNAM', 'FILTER', 'JAGFBAND'],
        EXPTIME = ['EXPTIME'],
        OBJECT = ['OBJECT'],
        DATE = ['DATE-OBS'],
        TIME = ['TIME-OBS'],
        AIRMASS = ['AIRMASS'],
        FILTER_WAVELENGTH = ['INSFLWL'],
        FILTER_WIDTH = ['INSFLDWL'],
        FILTER_ID = ['ALFLID', 'FAFTLID', 'JAGID', 'INSFLID'],
        TELESCOPE = ['TELESCOP', 'INSTRUME', 'ORIGIN', 'INSTRID']
    )

    _TELESCOPES_ALIASES = dict( NOT = ['ALFOSC', 'NOT'],
                               CAHA = ['DSAZ', 'CAFOS', 'CA-2.2', 'CAHA']
    )


    #_MANDATORY = ['FILTER', 'EXPTIME', 'OBJECT', 'DATE', 'AIRMASS']

    def __init__(self, image):
        self.image = image
        self.hdr = fits.getheader(self.image)


    @property
    def telescope(self):
        return self._get_telescope()[1]

    @property
    def filter_ID(self):
        return self._get_filterID()

    @property
    def filter_wavelength(self):
        """ IF the wavelength of the filter is present, try to find it. """
        pass

    @property
    def filter_width(self):
        pass

    @property
    def filter_curve(self):
        """ Find the filter used """
        pass




    @utils.memoize
    def _get_telescope(self):
        """ Try to find the telescope name from which the image comes.

         Since every observatory, telescope and/or instrument has different keywords, we need to go into some effort
         to do this search. The information of the telescope could be hidden in any/all/none of the keywords under the
         INSTRUMENT_TELESCOPE key of the dictionary _KEYWORDS_ALIASES above.

        """

        key, value = self.find_in_header(self._KEYWORDS_ALIASES['TELESCOPE'], self._TELESCOPES_ALIASES)
        return key, value

    @utils.memoize
    def _get_filterID(self):
        """ Try to find the ID of the filter used for the image.

        I know, right? Crazy to uniquely identify a filter in the header of the image...
        """
        print "Trying to get filter ID"
        value = self.value_in_header(self._KEYWORDS_ALIASES['FILTER_ID'])
        return value


    def find_in_header(self, list_keywords, dict_patterns):
        """ Find the keyword in which elements from both lists coincide.

        list_keywords is a list of keywords you might expect in normal headers, list(['keyword1', 'keyword2']).
        dict_patterns is something of the sort dict( name_a = [a1, a2, a3], name_b = [b1, b2]).
        What we want is to search for all the aliases of all the items of the dictionary within the elements in the
        header under the keywords of list_targets. In our example, we want to check  a1, a2, a3, b1, b2 agains the list
        formed by header[keyword1], header[keyword2]. If a match is found, say, between the alias b2 and the content of
        keyword keyword1, we want to return both keyword1 and the key name_b

        """

        dict_targets = {key: self.hdr.get(key) for key in list_keywords}
        for pattern_key, pattern_value in dict_patterns.iteritems():
            for target_key, target_value in dict_targets.iteritems():
                if any( ii in target_value for ii in pattern_value ):
                    return target_key, pattern_key
        return None, None

    def value_in_header(self, keywords):
        """ Find if any of the unique keywords is in the header. Return its value

        Look fot all the keywords in the header, if any is present and is not empty, return its value
        """

        for key in keywords:
            if self.hdr.get(key):
                return self.hdr.get(key)





hdr1 = header("./pruebas/NOT1.fits")
print hdr1.telescope

hdr2 = header("./CAHA1.fits")
print hdr2.telescope

#
# def guess_keywords(image):
#     """ Find those keywords present for each of the interesting keywords above.
#
#     :param image: FULL PATH to the image from which we want to learn the keywords
#     """
#     hdr = fits.getheader(image)
#     newdict = dict()
#     # For each of the keywords loop searching for the different options until finding the correct one
#     for k, v in keywords_aliases.iteritems():
#         present = [ii for ii in v if ii in hdr.keys()]
#         if len(present) == 0 and k in mandatory:  # no keywords found and it is a mandatory keyword!
#             raise NotFoundKeyword(k, image, v)
#         newdict[k] = present
#     return newdict
#
#
# def guess_filter(image):
#     """
#     Find out the wavelength and FWHM of filters for different observatories
#
#     :param image: path to the image we want to scan to find the filter
#     :return: (wavel, width), where wavel is the central wavelength of the pixel and width its width
#     """
#
#     NOT_filterIDs = dict()
#     JKT_filterIDs = dict()
#
#     # Read header and use the routine guess_keywords to find the correct keywords in the header
#     hdr = fits.getheader(image)
#     key_dict = guess_keywords(image)
#
#     # Wouldn't it be beautiful if every header specified with an ID exactly which filter was used?
#     # Modern NOT and JKT do it!
#
#
#
#
#     # Try and see if the telescope, wavelenth and FWHM of the filter were found in the keywords
#     wavek, fwhmk = key_dict['FILTER_WAVELENGTH'], key_dict['FILTER_WIDTH']
#     if wavek and fwhmk:
#         wavel = hdr.get(wavek)
#         width = hdr.get(fwhmk)
#
#     return (wavel, width)


