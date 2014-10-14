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
from scipy.interpolate import interp1d
from scipy.integrate import simps
from repipy import __path__ as repipy_path
import pyraf.iraf as iraf
import subprocess



regexp_dict = {'.*BIAS.*'             : 'bias',
               '(?=.*SKY)(?=.*FLAT)'  : 'skyflat',
               '(?=.*DOME)(?=.*FLAT)' : 'domeflat',
               '(.*FLAT)'           :  'flat',
               '(.* BLANK)'         : 'blank',
               '(?P<name>C(?:IG)?)(?P<number>\d{1,4})'    : 'cig'
               }

stds = numpy.genfromtxt("standards.csv", delimiter=",", dtype=None, autostrip=True, names=['std_names', 'ra', 'dec'])


class target(object):
    def __init__(self, image, hdr, filt):
        self.header = hdr
        self.filter = filt
        self.im_name = image

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
    def counts(self):
        """ Do photometry in the object to get the counts/sec of the source. """
        return self._get_photometry()

    @property
    def flux(self):
        """ If you have both a spectra for the object and the filter curve, calculate the flux under the filter"""
        if self.spectra != None and self.filter.filter_curve != None:
            return self._get_flux()


    @utilities.memoize
    def _get_RaDec(self):
        index =  numpy.where(stds['std_names'] == self.objname)[0]
        return stds['ra'][index], stds['dec'][index]


    @property
    def spectra(self):
        """ Read in the spectra  of the standards from repipy/standard_spectra/"""
        if self.objtype == 'standard':
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
    def _get_photometry(self):
        """ Get the photometry for the target.

        If the target is a standard star, aperture photometry will be performed. For the moment nothing is done with
        the others, but in due time (TODO) photometry.py will be included here. """

        coords_file = "remove_this.txt"
        with open(coords_file, 'w') as ff:
            ff.write(str(self.RA) + " " + str(self.DEC))

        if self.objtype == "standard":
            iraf.noao(_doprint=0)
            iraf.digiphot(_doprint=0)
            iraf.apphot(_doprint=0)
            seeing = self.header.hdr[self.header.seeingk]
            #seeing = 2.5
            print "seeing =", seeing
            photfile_name = self.im_name + ".mag.1"
            utilities.if_exists_remove(photfile_name)
            iraf.module.phot(self.im_name, output=photfile_name, coords=coords_file,
                      wcsin='world', fwhm=seeing, gain=self.header.gaink, exposure=self.header.exptimek,
                      airmass=self.header.airmassk, annulus=6*seeing, dannulus=3*seeing,
                      apert=2*seeing, verbose="no", verify="no", interac="no")
            [counts] = iraf.module.txdump(photfile_name, 'FLUX', 'yes', Stdout=subprocess.PIPE)
            utilities.if_exists_remove(coords_file)
            return counts


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

    @utilities.memoize
    def _get_flux(self):
        """ Get the flux of the object under the filter by convolving the filter curve with the spectra of the object

        This program follows one called convol.pro from Jorge Iglesias, IAA.
        """

        # Get the spectra of the star and the filter curve. The spectra must be in A
        wavelength_star = self.spectra.transpose()[0]
        magnitude_star = self.spectra.transpose()[1]
        wavelength_filter = self.filter.filter_curve.transpose()[0]
        transmittance_filter = self.filter.filter_curve.transpose()[1]

        # We will use only the area where filter and star overlap, tipically the spectra of the star is much larger
        wav_min = max(wavelength_filter.min(), wavelength_star.min())
        wav_max = min(wavelength_filter.max(), wavelength_star.max())


        # Now use those cuts
        magnitude_star = magnitude_star[numpy.where(   (wav_min < wavelength_star) & (wavelength_star < wav_max))[0]]
        wavelength_star = wavelength_star[numpy.where( (wav_min < wavelength_star) & (wavelength_star < wav_max) )[0]]
        transmittance_filter = transmittance_filter[numpy.where( (wav_min < wavelength_filter) &
                                                                 (wavelength_filter < wav_max) )]
        wavelength_filter = wavelength_filter[numpy.where( (wav_min < wavelength_filter) &
                                                           (wavelength_filter < wav_max) ) ]

        # We convert from AB magnitudes to flux in erg/s/cm2/Hz.
        flux_star = 10**((-48.6-magnitude_star)/2.5)

        # Then from F_nu we need to convert to F_lambda
        delta_lambda = wavelength_star[1] - wavelength_star[0]
        flux_star = flux_star * 3e18 / wavelength_star ** 2


        # Then we need to sample the transmittance of the filter at the wavelengths of the star.
        f = interp1d(wavelength_filter, transmittance_filter, kind='cubic', fill_value=0, bounds_error=False)
        transmisttance_interpolated = f(wavelength_star)

        # Finally, we integrate the multiplication of both curves
        y = transmisttance_interpolated * flux_star
        x = wavelength_star
        result = sum(y * delta_lambda)
        print "Integral = ", result
        return result






    def _name_using_coordinates(self):
        """ Check if any of the standard stars is within the image.

        Create a wcs object from the header and search for each of the stars within the image. If any is present,
        return the type 'standards' and the name of the standard star in the field of view.
        """
        type, name = 'Unknown', 'Unknown'
        w = wcs.wcs.WCS(self.header.hdr)
        ly, lx = fits.getdata(self.im_name).shape
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
