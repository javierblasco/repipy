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

    @property
    def data(self):
        """ Create an array (or a list of arrays, if several chips are present) with all the data in the image

        For each of the objects in the fits image, collect the associated data, if present. At the end, if the image
        contains only one chip, do not return an list of a single element, return simply the array with the data.

        :return:
        """
        data = [ii.data for ii in self._HDUList if ii.data is not None]
        if len(data) == 1:
            data = data[0]
        return data
