import lemon.photometry as photometry
import repipy.extract_mag_airmass_common as extract
import repipy.utilities as utilities
import os

class astronomical_object(object):
    def __init__(self, obj_name="Unknown", obj_type="Unknown", keywords=None):
        self.Name = obj_name
        self.obj_type = obj_type
        self.keywords = keywords
        self.objtype = None
        self.narrow_images = None
        self.cont_images = None
        self.narrow_final = None
        self.cont_final = None

    def scale_cont(self):
        dir = os.path.dirname(self.narrow_final)
        if not dir:
            dir = "."
        filename = self.Name + "_scaling.db"
        output_db = os.path.join(dir, filename)
        utilities.if_exists_remove(output_db)
        photometry.main(arguments=["--maximum", "55000", "--uik", "", "--margin", "20", "--gaink", self.keywords["gaink"],
                              "--aperture", "4.", "--annulus", "6", "--dannulus", "2", "--individual-fwhm", "--objectk",
                              self.keywords["objectk"], "--filterk", self.keywords["filterk"], "--datek", self.keywords["datek"],
                              "--expk", self.keywords["exptimek"], "--fwhmk", "seeing", "--airmk", self.keywords["airmassk"],
                              self.cont_final, self.narrow_final,  self.cont_final, output_db])
        airmasses, magnitudes, filters = extract.main(output_db)
        print "\n"
        for a, mag, ff in zip(airmasses, magnitudes, filters):
            print a, mag, ff

