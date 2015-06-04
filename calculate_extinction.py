#! /usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Given some databases from lemon photometry (.db files), each from a different field,
calculate the extinction coefficient for each of the fields and each of the filters



Created on Mon Dec 16 12:13:00 2013

@author: blasco
"""

import numpy
import sys
import scipy.stats as stats
from scipy.optimize import fmin as simplex
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse
import repipy.extract_mag_airmass_common as extract
import lemon.snr as snr

# def single_slope_chi2(params, X, Y, YERR):
# # Follows: http://python4mpia.github.io/fitting_data/simplex-fitting.html
#
#     # Split params into a and b
#     a = numpy.array(params[0])
#     b = numpy.array(params[1:])
#
#     # Merit function is chi-square
#     model = ((X * a).transpose() + b).transpose()
#     chi2 = numpy.sum( (model - Y)**2 / YERR**2 )
#     print a, b, chi2, "\n"
#     return chi2

def single_slope_curvefit(x, a, *b):
    a = numpy.array(a)
    b = numpy.array(b)
    return (((x * a).transpose() + b).transpose()).flatten()


def extinction_fit(airmasses, magnitudes, err_magnitudes=None):
    """Given 2D arrays of airmasses and magnitudes of several stars. For N
     observations of M stars, we would have a MxN array for both airmasses and 
     magnitudes, so that we can iterate for each star. Then we calculate the 
     extinction. We will fit a single coefficient k to all the stars. 
     That is, for each individual star, we fit:
         
         magnitudes_corrected = magnitudes - k * airmass 
     
     where k is common for all the stars, airmass and magnitudes is known and
     magnitudes_corrected would be the magnitudes of the stars with no 
     atmospheric extinction. """

    # In the unlikely case you send magnitudes without errors (do your
    # homework, dude!) same error will be assigned to all values, so weights are
    # not used.
    if err_magnitudes is None:
        err_magnitudes = numpy.ones_like(magnitudes) * 0.1

    # Initial guess for extinction coefficient is the first value.
    # The initial guess of the intercept for each star is the mean of the
    # star's magnitude
    means = [numpy.mean(magnitudes[ii,]) for ii in range(magnitudes.shape[0])]
    guess = [4] + means
    pop, pcov = curve_fit(single_slope_curvefit, airmasses, magnitudes.flatten(),
                          sigma=err_magnitudes.flatten(), p0=guess)
    return pop[0], numpy.sqrt(pcov[0, 0])


def select_brightest(magnitudes, *args, **kwargs):
    """ Select the nstars brightest objects from the array magnitudes. Select the same positions for the rest of arrays.
    magnitudes:  magnitudes of the stars in the different images. Array of shape NxM with N the number
     of stars and M the number of images.
    nstars: select this number of the brightest in the field
    :return: return the same arrays in the same order, just using the brightest stars.
    """

    # Default value for nstars is None, and args shoule be a list so we can modify it
    nstars = kwargs.pop("nstars", None)  # Default value for nstars
    args = list(args)  # so that we can manipulate it

    # Sort according to the average magnitude
    average_magnitude = numpy.nanmean(magnitudes, axis=1)
    srt_args = numpy.argsort(average_magnitude)

    # If brightest Nstar stars should be returned, cut the sorted arguments accordingly
    if nstars is not None:
        srt_args = srt_args[0:nstars]

    # Cut all the arrays passed
    magnitudes = magnitudes[srt_args, :]
    for ii, arr in enumerate(args):
        args[ii] = arr[srt_args, :]

    if args:
        return [magnitudes] + args
    else:
        return magnitudes


def calculate_extinctions(args):
    """ Iterate through the databases and find the extinction coefficient for every filter
    """
    results = {}  # Dictionary for the outputs
    for DB in args.inputdb:
        print "\nStudying database {0}:".format(DB)
        airmass, magnitudes, filters, SNR = extract.main(DB)

        # SNR is the signal-to-noise-ratio, we want to use the errors in the magnitudes
        magnitude_errors_minus, magnitude_errors_plus = snr.snr_to_error(SNR)
        mag_errors = (-magnitude_errors_minus + magnitude_errors_plus) / 2.

        # Which filters are present?
        filter_set = set(filters.flatten().tolist())

        # Prepare the figure for the plots.
        if args.plot:
            figure = plt.figure()
            figure.set_label("BD+25")
            ax0 = plt.subplot2grid((2,2),(0, 0))
            ax1 = plt.subplot2grid((2,2),(0, 1))
            ax2 = plt.subplot2grid((2,2),(1, 0))
            ax = [ax0, ax1, ax2]
            plt.subplots_adjust(hspace=0.3)
            plt.subplots_adjust(wspace=0.12)
            plt.subplots_adjust(top=0.95)
            plt.subplots_adjust(bottom=0.05)
            plt.subplots_adjust(left=0.05)
            plt.subplots_adjust(right=0.95)



        for ii, filt in enumerate(filter_set):
            # Select only those data for the particular filter
            whr_filt = numpy.where(filters[0, :] == filt)[0]


            # If the airmass difference is not, at least, 0.2, skip this test because the fit will be meaningless
            airm_diff = numpy.max(airmass[:, whr_filt]) - numpy.min(airmass[:, whr_filt])
            if airm_diff < 0.2:
                print "Airmasses are too close for filter {0} ".format(filt)
                plt.close()
                continue

            # Select the brightest objects and do the fit
            mags_filt, airm_filt, errors_filt = select_brightest(magnitudes[:, whr_filt],
                                                                 airmass[:, whr_filt],
                                                                 mag_errors[:, whr_filt],
                                                                 nstars=args.nstars)
            ext, ext_err = extinction_fit(airm_filt, mags_filt, errors_filt)
            print "Extinction coefficient for filter {0}: {1} +- {2}".format(filt, ext, ext_err)
            results[(DB, filt)] = (ext, ext_err)

            # Do the plot of the selected stars
            if args.plot:
                for star_airmass, star_mag in zip(airm_filt, mags_filt):
                    ax[ii].plot(star_airmass, star_mag, 'o')
                    ax[ii].set_title(filt)
                    ax[ii].set_xlabel("Airmass")

        if args.plot:
            plt.show()

    return results


########################################################################################################################


# Create parser
parser = argparse.ArgumentParser(description='''Calculate the extinction coefficients from one (or more) lemon
                                                databases (the result of lemon photometry). The program will print
                                                on the screen the values, and will return as many dictionaries as
                                                databases where in the input. The keys of each dictionary will be
                                                the names of the filters present, the values will be the extinction
                                                coefficient and its uncertainties. ''')

# Add necessary arguments to parser
parser.add_argument("inputdb", metavar='inputdb', action='store', help='list of ' + \
                                                                       'input databases from which to calculate the extinction coefficient. ', \
                    nargs="+", type=str)
parser.add_argument("--plot", action="store_true", dest="plot", default=False,
                    help="Shows the plots of magnitude vs airmass for all db and filters. Default: False")
parser.add_argument("--silent", action="store_true", dest="silent", default=False,
                    help="Activate this argument if you don't want the program to print on the screen.")
parser.add_argument("--nstars", action="store", dest="nstars", default=None, type=int,
                    help="Use only the N brightest stars to the fit the extinction coefficient. By default all stars "
                         "are used. ")


def main(arguments=None):
    # Pass arguments to variable args
    if arguments == None:
        arguments = sys.argv[1:]

    args = parser.parse_args(arguments)

    dictionaries = calculate_extinctions(args)

    return dictionaries


if __name__ == "__main__":
    main()
