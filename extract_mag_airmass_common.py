#! /usr/bin/env python
# Routine created by Victor Terron, modified by Javier Blasco to suit his needs.
# LEMON modules
import lemon.database as database
import matplotlib.pyplot as plt
import numpy as np

def main(db_name):
    db = database.LEMONdB(db_name)

    # (1) Identify stars that have been observed in the same number of images
    # that the star that has been observed in the maximum number of images :-)
    # I could have simply used the image count of the images, but that would
    # discard all the stars if there were an image without stars.

    query = """
            SELECT s.id
            FROM stars AS s
            INNER JOIN photometry AS p
            ON s.id = p.star_id
            GROUP BY p.star_id
            HAVING COUNT(*) = (SELECT MAX(nimages)
                               FROM (SELECT COUNT(*) AS nimages
                               FROM photometry
                               GROUP BY star_id))
            """


    db._execute(query)

    # IDs of the stars observed in all the images
    star_ids = set(x[0] for x in db._rows)
    num_stars = len(star_ids)


    # (2) Identify the number of images involved so that we can return airmass and magnitude of stars in all images
    query = """
           SELECT id
           FROM images
           WHERE sources = 0
           """
    db._execute(query)
    num_airmasses = len(set(_ for _ in db._rows))

    # Create numpy arrays to store airmasses and magnitudes
    airmasses = np.zeros([num_stars, num_airmasses], dtype=np.float64)
    magnitudes = np.zeros([num_stars, num_airmasses], dtype=np.float64)
    filters = np.array([""] * num_stars * num_airmasses, dtype="S15").reshape(num_stars, num_airmasses)

    for ii, id_ in enumerate(star_ids):
        query = """
                SELECT p.star_id, f.name, p.magnitude, i.airmass
                FROM photometry AS p
                INNER JOIN images AS i
                ON p.image_id = i.id
                INNER JOIN photometric_filters AS f
                ON i.filter_id = f.id
                WHERE p.star_id = ?
                """
        t = (id_, )
        db._execute(query, t)

        for jj, row in enumerate(db._rows):
            star_id, filter_name, magnitude, airmass = row
            airmasses[ii,jj] = airmass
            magnitudes[ii,jj] = magnitude
            filters[ii, jj] = filter_name

    return airmasses, magnitudes, filters