#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Created on Tue Jul  9 11:03:50 2013

@author: javier blasco herrera
"""

import repipy.utilities as utilities
import subprocess
import os

def main(im_name, ra=None, dec=None, FoV=None):


    if ra == None or dec == None:
        if FoV == None:
            FoV = 20  # Largest image I've seen, 20 degrees, just in case

        # User passed the keywords?
        try:
            ra, dec = utilities.get_from_header(im_name, ra, dec)
            argunents = ["--ra", ra, "--dec", dec, "--radius", float(FoV)/2.]
        except KeyError:
            arguments = []

        # User passed the numbers
        try:
            ra, dec = float(ra), float(dec)
            argunents = ["--ra", ra, "--dec", dec, "--radius", float(FoV)/2.]
        except ValueError:
            arguments = []


    arguments = arguments + ["solve-field", "--no-plots", "--no-fits2fits",
                             "--use-sextractor", "--code-tolerance", "0.01",
                             "--overwrite", im_name]
    subprocess.call(arguments)

    # Sort things out, clean up...
    solved = utilities.replace_extension(im_name, "solved")
    output_name = utilities.replace_extension(im_name, "new")
    assert os.path.exists(output_name)
    shutil.move(output_name, im_name)

    # Update old WCS system PC matrices and so on to avoid confusion
    utilities.update_WCS(im_name, im_name)

    # Build a catalogue with the coordinates of the stars found in
    # X and Y, RA and DEC and the flux of the stars
    corr_file = utilities.replace_extension(im_name, "corr")
    table = fits.open(corr_file)[1]
    cat_radec = utilities.replace_extension(im_name, "radec")
    f = open(cat_radec, "w")
    for line in table.data:
        f.write(str(line[2]) + " " + str(line[3]) + "\n")
    f.close()
