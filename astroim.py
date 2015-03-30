import repipy.header as header
import repipy.target as target
import repipy.filter as filter
import repipy.utilities as utilities
import numpy

class Astroim(object):
    """
    """
    def __init__(self, image):
        self.im_name = image
        self.header = header.Header(self.im_name)
        self.filter = filter.Filter(self.header)
        self.target = target.Target(self.header, self.filter)

    def zero_point():
        return self.filter.zero_point(self.target)





