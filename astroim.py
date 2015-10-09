import repipy.header as header
import repipy.target as target
import repipy.filter as imfilter


class Astroim(object):
    """
    """
    def __init__(self, image):
        self.im_name = image
        self.header = header.Header(self.im_name, self._HDUList)
        self.filter = imfilter.Filter(self.header)
        self._HDUList = fits.open(self.im_name, 'readonly')
        self.target = target.Target(self.header, self.filter)
        self.imstats = imstats.Imstats(self)

    def zero_point(self, aperture=None):
        return self.filter.zero_point(self.target, aperture=aperture)





