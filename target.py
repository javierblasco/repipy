import lemon.photometry as photometry
import lemon.passband as passband
import repipy.extract_mag_airmass_common as extract
import repipy.utilities as utilities
import os
import sys
import astropy.io.fits as fits
import repipy.header as header
import numpy
import re
import astropy.wcs as wcs
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import simps
from repipy import __path__ as repipy_path
import subprocess
import tempfile
from collections import OrderedDict
from astropy.coordinates import SkyCoord
import astropy.units as u
from lemon import methods
import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import digiphot
    from iraf import daophot
    from iraf import apphot


regexp_dict = OrderedDict()
items = (('.*BIAS.*', 'bias'),
         ('(?=.*SKY)(?=.*FLAT)?', 'skyflat'),
         ('(?=.*DOME)(?=.*FLAT)', 'domeflat'),
         ('(.*FLAT.*)', 'flat'),
         ('(.* BLANK)', 'blank'),
         ('(?P<name>C(?:IG)?)(?P<number>\d{1,4})', 'cig')
         )
for key, val in items:
    regexp_dict[key] = val

standards_file = os.path.join(repipy_path[0], "standards.csv")
stds = numpy.genfromtxt(standards_file, delimiter=",", dtype=None, autostrip=True, names=['std_names', 'ra', 'dec'])


class Target(object):
    def __init__(self, hdr, filt, chips):
        self.header = hdr
        self.filter = filt
        self.chips = chips

    def __str__(self):
        return re.sub('[-+\s_]', "", self.objname.lower())

    @property
    def objtype(self):
        """ Find out the type of image, (bias, flat, cig, standard). """
        return self._type_and_name()[0]

    @property
    def objname(self):
        """ Find out the name of the target in the image """
        return self._type_and_name()[1]

    @utilities.memoize
    def _type_and_name(self):
        """ Find out what type of image (bias, flat, cig, standard, ...) this object is

        Several types of search are combined to find out the sort of object observed (bias, flat, cig, ...) and the name
        of the object. A match using regular expressions is tried for those objects with easy patterns. As a second
        option, using the coordinates of the object it is determined whether it is a standard star and which one. If
        both fail, the type and name of the object are set to 'Unknown'.

        """
        if self.regex_match():
            objtype, match = self.regex_match()
            objname = objtype
        elif self.is_a_standard():
            objtype = "standards"
            objname = self.is_a_standard()
        else:
            objtype = 'Unknown'
            name = self.header.get(self.header.objectk)  #whatever the header says
            objname = re.sub("[_:\s\t\|\\\\]","", name)

        # CIGs are usually followed by a number, CIG23, CIG815, CIG1020
        if objtype == 'cig':
            objname += match.group('number').zfill(4)
        return objtype, objname

    @utilities.memoize
    def regex_match(self):
        """ Use regular expressions to guess the name of the object

        Use the regular expressions in the regexp_dict to match the field that contains the objext name in the header.
        :return: match of the regular expression
        """
        name = self.header.get(self.header.objectk)
        if name is None:
            return None

        name = re.sub("[_:\s\t\|\\\\]","", name)
        for key, value in regexp_dict.iteritems():
            match = re.match(key, name, re.I)
            if match:
                return value, match


    @utilities.memoize
    def is_a_standard(self):
        """ Find if the coordinates of a set of spectrophotometric standards is in any of the chips of the image.
        :return: the name of the standard, if any is found. None otherwise.
        """
        for ii, star in enumerate(stds['std_names']):
            ra, dec = stds['ra'][ii], stds['dec'][ii]
            for chip in self.chips:
                try:
                    if chip.contains_coords(ra, dec):
                        return star
                except wcs.wcs.NoConvergence:
                    continue
        return None

    @property
    def RA(self):
        """ RA of the object targeted by the observation, if known. """
        return float(self._get_RaDec()[0])

    @property
    def DEC(self):
        """ DEC of the object targeted by the observation, if known."""
        return float(self._get_RaDec()[1])


    def counts(self, aperture=None):
        """ Do photometry in the object to get the counts/sec of the source. """
        return self._get_photometry(aperture=aperture)

    def flux(self):
        """ If you have both a spectra for the object and the filter curve, calculate the flux under the filter"""
        if self.spectra is not None and self.filter.filter_curve is not None:
            return self._get_flux()

    @utilities.memoize
    def _get_RaDec(self):
        """ Get the RA and DEC from the header of an image.

        If the object is in the catalogue of standars, forget about anything else, just use the coordinates of the
        standard. If the object is not a standard, try to get it from the keywords of the header.
        :return:
        """
        ra, dec = None, None
        if self.objtype == 'standards':
            index =  numpy.where(stds['std_names'] == self.objname)[0]
            ra, dec =  stds['ra'][index], stds['dec'][index]
        else:  # object not a standard
            try: # try reading it from the header.
                ra = self.header.hdr[self.header.RAk]
                dec = self.header.hdr[self.header.DECk]
                # If ra is in hours, it could be "16:34:23.345" or "16h34m23.345s" or even "16 34 23.345"... Difficult!
                regexp_hours = "\d{1,2}[:,h,\s+]\d{1,2}[:,m,\s+]\d{1,2}\.?\d?"
                regexp_degree = "\d{1,3}\.?\d"
                if re.search(regexp_hours, ra):
                    ra_units = u.hourangle
                elif re.search(regexp_degree):
                    ra_units = u.deg
                else:
                    sys.exit("Units of RA in the header not understood! Read _get_RaDec method in target.py. ")
                c = SkyCoord(ra, dec, unit=(ra_units, u.deg))
                ra, dec = c.ra.deg, c.dec.deg
            except:  # not able to read RA, DEC from the header
                raise
        return ra, dec

    @property
    def spectra(self):
        """ Read in the spectra  of the standards from repipy/standard_spectra/"""
        if self.objtype == 'standards':
            dir = os.path.join(repipy_path[0], "standard_spectra")
            file = os.path.join(dir, self.__str__())

            # Determine how long is the header of the file, how many lines to skip.
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

    def _get_photometry(self, seeing=None, aperture=None):
        """ Get the photometry for the target.

        If the target is a standard star, aperture photometry will be performed. For the moment nothing is done with
        the others, but in due time (TODO) photometry.py will be included here. """

        if self.objtype is not "standards":
            return None

        basename = "standards"
        fd, coords_file = tempfile.mkstemp(prefix=basename, suffix=".coords")
        os.write(fd, "{0} {1} \n".format(self.RA, self.DEC))
        os.close(fd)


        # If aperture was not given by the user, try and use the seeing as a reference for default values
        if not aperture:
            try:
                seeing = self.header.hdr[self.header.seeingk]
                aperture = 3 * seeing
            except ValueError:  # keyword was not correctly guessed
                pass

        # If aperture exists, do the photometry using it
        if aperture:
            annulus = 2 * aperture
            dannulus = max([aperture, 3])  # minimum of 3 pixels thickness for the sky annulus
            fd, photfile_name = tempfile.mkstemp(".mag.1")
            utilities.if_exists_remove(photfile_name)
            kwargs =  dict(output=photfile_name, coords=coords_file, salgorithm='median',
                      wcsin='world', fwhm=seeing, gain=self.header.gaink, exposure=self.header.exptimek,
                      airmass=self.header.airmassk, annulus=annulus, dannulus=dannulus,
                      apertures=aperture, verbose="no", verify="no", interac="no")
            iraf.phot(self.header.im_name, **kwargs)
            [counts] = iraf.txdump(photfile_name, 'FLUX', 'yes', Stdout=subprocess.PIPE)
        else:
            sys.exit("\n \n Sorry, no aperture was passed by you, and a seeing keyword was not "
                     "found in the header. \n\n ")

        utilities.if_exists_remove(coords_file)
        return float(counts)

    @utilities.memoize
    def _get_object(self):
        """ From the header read the object name, and try to find the type of object and name using regular expressions.
        Failing that, try to find in the image one of the list of standard stars in the file 'standards.csv' of the
        folder in repipy. """

        # Default type is Unknown, default name, whatever comes in the header
        type, name = "Unknown", "Unknown"

        try:
            name = self.header.hdr[self.header.objectk]
            if name.replace(" ","") == "":   # no name in the OBJECT field
                name = "Unknown"
        except TypeError: # no keyword OBJECT of any sort in header
            pass

        # Search for the patterns above in the "object field" of the header. If the name corresponds to a bias, for
        # example, both the type and the name will be 'bias'. If it is a CIG galaxy, the type will be 'cig', and the
        # name is 'cig' followed by the number of the cig
        for key, value in regexp_dict.iteritems():
            match = re.match(key, re.sub("[_:\s\t\|\\\\]","", name), re.I)
            if match:
                type = value
                name = value
                try:
                    name += match.group('number')   # case of CIGs
                except IndexError:
                    pass
                break
        else:  # If none of the above was found, try to match the coordinates of the image with the list of standards
            result = self._name_using_coordinates()
            if result:  # if not None
                type, name = result
        return type, name

    @utilities.memoize
    def _get_flux(self):
        """ Get the flux of the object under the filter by convolving the filter curve with the spectra of the object

        This program follows one in IDL called convol.pro from Jorge Iglesias, IAA.
        """
        if self.spectra is not None:  # For standards, where the spectra are used to calibrate in flux
            # Get the spectra of the star and the filter curve. The spectra must be in Angstrom and magnitudes AB
            wavelength_star = self.spectra.transpose()[0]
            magnitude_star = self.spectra.transpose()[1]
            wavelength_filter = self.filter.filter_curve.transpose()[0]
            transmittance_filter = self.filter.filter_curve.transpose()[1]

            # We will use the area under the filter, just right of the first point, left of the last
            wav_min = int(numpy.ceil(wavelength_filter.min()))
            wav_max = int(numpy.floor(wavelength_filter.max()))
            wavelength = numpy.arange(wav_min, wav_max)

            # Now interpolate the magnitude of the star and the transmittance of the filter
            f = InterpolatedUnivariateSpline(wavelength_filter, transmittance_filter, k=3)
            transmittance_interpolated = f(wavelength)
            g = InterpolatedUnivariateSpline(wavelength_star, magnitude_star, k=3)
            magnitude_interpolated = g(wavelength)


            # We convert from AB magnitudes to flux in erg/s/cm2/Hz.
            flux_star = 10**((-48.6-magnitude_interpolated)/2.5)
            delta_lambda = wavelength[1] - wavelength[0]
            flux_star = flux_star * 3e18 / wavelength ** 2


            # Finally, we integrate the multiplication of both curves
            y = transmittance_interpolated * flux_star
            x = wavelength
            # Assuming delta_lambda is sufficiently small that there is no large changes in y
            result = sum(y * delta_lambda)
            return result


    def _name_using_coordinates(self):
        """ Check if any of the standard stars is within the image.

        Create a wcs object from the header and search for each of the stars within the image. If any is present,
        return the type 'standards' and the name of the standard star in the field of view.
        """
        try:
            w = wcs.wcs.WCS(self.header.hdr)
        except ValueError:   # problems, for example, with the SIP polynomials
            return None

        if not w.is_celestial:
            return None

        for ii, star in enumerate(stds['std_names']):
            ra, dec = stds['ra'][ii], stds['dec'][ii]
            if (ra_min < ra < ra_max) and (dec_min < dec < dec_max):
                type, name = 'standards', star
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
