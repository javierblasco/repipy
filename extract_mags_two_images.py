#! /usr/bin/env python

# Author: Victor Terron (c) 2014
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

""" Extract the magnitudes of the stars common to two images """

# LEMON modules
import lemon.database as database
import lemon.fitsimage as fitsimage


def main(db_name, im1_name, im2_name):

#    path = "phot.LEMONdB"
#    db = database.LEMONdB("phot.LEMONdB")
#    img1 = fitsimage.FITSImage("ngc2264_sa/ferM_061sa.fits")
#    img2 = fitsimage.FITSImage("ngc2264_sa/ferM_135sa.fits")

    db = database.LEMONdB(db_name)
    img1 = fitsimage.FITSImage(im1_name)
    img2 = fitsimage.FITSImage(im2_name)

    utime1 = img1.date()
    utime2 = img2.date()

    # Extract (star ID, mag) for a given Unix time
    query = """
            SELECT p.star_id, p.magnitude
            FROM photometry AS p
            INNER JOIN images AS i
            ON p.image_id = i.id
            WHERE i.unix_time = ?
            """

    # Map each ID to its magnitude
    db._execute(query, (utime1,))
    img1 = dict(db._rows)
    db._execute(query, (utime2,))
    img2 = dict(db._rows)

    # IDs of the stars common to both images
    img1_ids = set(img1.iterkeys())
    img2_ids = set(img2.iterkeys())
    common_ids = img1_ids.intersection(img2_ids)


    # For each star, get magnitude in both images
    mag1 = []
    mag2 = []
    for star_id in sorted(common_ids):
        mag1.append(img1[star_id])
        mag2.append(img2[star_id])
        
    return (mag1, mag2)


