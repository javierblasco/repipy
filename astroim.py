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
        self.target = target.Target(self.im_name, self.header, self.filter)
        #self.zero_point = self._get_zero_point()

    def _get_zero_point(self):  # esto no deberia ir aqui. Victor.
        exptime = self.header.hdr[self.header.exptimek]
        return 2.5 * (numpy.log10(self.target.counts / exptime) - numpy.log10(self.target.flux))



