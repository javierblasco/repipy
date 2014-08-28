import lemon.photometry as photometry
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
        filename = self.Name + "_final.db"
        output_db = os.path.join(dir, filename)
        photometry.main(arguments=["--maximum", "55000", "--uik", "", "--margin", "20", "--gaink", keywords["gaink"],
                              "--aperture", "4.", "--annulus", "6", "--dannulus", "2", "--individual-fwhm", "--objectk",
                              keywords["objectk"], "--filterk", keywords["filterk"], "--datek", keywords["datek"],
                              "--expk", keywords["exptimek"], "--fwhmk", "seeing", "--airmk", keywords["airmassk"],
                              self.cont_final, self.narrow_final,  self.cont_final, output_db])



