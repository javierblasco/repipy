import lemon.photometry as photometry
import lemon.passband as passband
import repipy.extract_mag_airmass_common as extract
import repipy.utilities as utilities
import os
import astropy.io.fits as fits
import repipy.header as header
import numpy
import re
import astropy.wcs as wcs
import scipy
from repipy import __path__ as repipy_path


regexp_dict = {'.*BIAS.*'             : 'bias',
               '(?=.*SKY)(?=.*FLAT)'  : 'skyflat',
               '(?=.*DOME)(?=.*FLAT)' : 'domeflat',
               '(.*FLAT)'           :  'flat',
               '(.* BLANK)'         : 'blank',
               '(?P<name>C(?:IG)?)(?P<number>\d{1,4})'    : 'cig'
               }

stds = numpy.genfromtxt("standards.csv", delimiter=",", dtype=None, autostrip=True, names=['std_names', 'ra', 'dec'])


class target(object):
    def __init__(self, image):
        self.header = header.header(image)
        self.image = image

    def __str__(self):
        return re.sub('[-+\s_]', "", self.objname.lower())

    @property
    def objtype(self):
        """ Find out what type of image (bias, flat, cig, standard, ...) this object is"""
        return self._get_object()[0]

    @property
    def objname(self):
        """ Using the header object, find out the name of the target in the image """
        return self._get_object()[1]

    @property
    def RA(self):
        return float(self._get_RaDec()[0])

    @property
    def DEC(self):
        return float(self._get_RaDec()[1])


    @property
    def spectra(self):
        """ Read in the spectra  of the standards from repipy/standard_spectra/"""
        dir = os.path.join(repipy_path[0], "standard_spectra")
        file = os.path.join(dir, self.__str__())

        # Determine how long is the header:
        with open(file, 'r') as ff:
            for ii, line in enumerate(ff):
                try:
                    a, b = line.split()
                    a, b = float(a), float(b)
                    nn = ii  # number of lines to skip
                    break
                except ValueError:
                    continue

        return numpy.genfromtxt(file, skip_header=nn)




    @utilities.memoize
    def _get_RaDec(self):
        index =  numpy.where(stds['std_names'] == self.objname)[0]
        return stds['ra'][index], stds['dec'][index]

    @utilities.memoize
    def _get_object(self):
        """ From the header read the object name, and try to find the type of object and name using regular expressions.
        Failing that, try to find in the image one of the list of standard stars: 'standards.csv' """
        name_in_header = self.header.hdr[self.header.objectk]

        # Search for the patterns above in the "object field" of the header. If the name corresponds to a bias, for
        # example, both the type and the name will be 'bias'. If it is a CIG galaxy, the type will be 'cig', but the
        # name is 'cig' followed by the number of the cig
        for key, value in regexp_dict.iteritems():
            match = re.match(key, name_in_header.replace(" ",""), re.I)
            if match:
                type = value
                name = value
                try:
                    name += match.group('number')   # case of CIGs
                except IndexError:
                    pass
                break
        else:  # If none of the above was found, try to match the coordinates of the image with the list of standards
            type, name = self._name_using_coordinates()
        return type, name

    def _name_using_coordinates(self):
        """ Check if any of the standard stars is within the image.

        Create a wcs object from the header and search for each of the stars within the image. If any is present,
        return the type 'standards' and the name of the standard star in the field of view.
        """
        type, name = 'Unknown', 'Unknown'
        w = wcs.wcs.WCS(self.header.hdr)
        ly, lx = fits.getdata(self.image).shape
        # ra_min, dec_min will be the coordinates of pixel (0,0)
        # ra_max, dec_max the coordinates of the upper right corner of the image.
        # But beware, the orientation could be such that ra_min > ra_max if the image is not oriented North
        ra_min, dec_min = w.all_pix2world(numpy.array(zip([0],[0])), 1)[0]
        ra_max, dec_max = w.all_pix2world(numpy.array(zip([lx], [ly])), 1)[0]
        for ii, star in enumerate(stds['std_names']):
            ra, dec = stds['ra'][ii], stds['dec'][ii]
            if min(ra_min, ra_max) < ra < max(ra_min, ra_max) and min(dec_min, dec_max) < dec < max(dec_min, dec_max):
                    type, name = 'standard', star
                    break
        return type, name























    # def scale_cont(self):
    #     dir = os.path.dirname(self.narrow_final)
    #     if not dir:
    #         dir = "."
    #     filename = self.Name + "_scaling.db"
    #     output_db = os.path.join(dir, filename)
    #     utilities.if_exists_remove(output_db)
    #     photometry.main(arguments=["--maximum", "55000", "--uik", "", "--margin", "20", "--gaink", self.keywords["gaink"],
    #                           "--aperture", "4.", "--annulus", "6", "--dannulus", "2", "--individual-fwhm", "--objectk",
    #                           self.keywords["objectk"], "--filterk", self.keywords["filterk"], "--datek", self.keywords["datek"],
    #                           "--expk", self.keywords["exptimek"], "--fwhmk", "seeing", "--airmk", self.keywords["airmassk"],
    #                           self.cont_final, self.narrow_final,  self.cont_final, output_db])
    #     airmasses, magnitudes, filters = extract.main(output_db)
    #     narrow_filter, cont_filter = utilities.collect_from_images([self.narrow_final, self.cont_final], self.keywords["filter"])
    #     narrow_filter, cont_filter = (passband.Passband(filt).__str__() for filt in (narrow_filter, cont_filter))
    #
    #     magnitudes_narrow = magnitudes[filters == narrow_filter]
    #     magnitudes_cont = magnitudes[filters == cont_filter]
    #     #scaling_factor =
    #
    #     print "\n"
    #     for a, mag, ff in zip(airmasses, magnitudes, filters):
    #         print a, mag, ff
    #
