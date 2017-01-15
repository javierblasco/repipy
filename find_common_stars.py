import scipy.spatial as spatial
import numpy

def main(coord_files, outputfile, error=0.0006):
    x1, y1 = numpy.genfromtxt(coord_files[0], unpack=True)
    for cfile in coord_files[1:]:
        tree_ref = spatial.KDTree( zip(x1, y1))
        x2, y2 = numpy.genfromtxt( cfile, unpack=True )
        tree_obj = spatial.KDTree( zip(x2, y2))
        matches = tree_ref.query_ball_tree(tree_obj, error)
        matches_list = [elem[0] for elem in matches if len(elem) ==1]
        xcommon = [x2[match] for match in matches_list]
        ycommon = [y2[match] for match in matches_list]
        x1, y1 = xcommon, ycommon

    with open(outputfile, 'w') as fd:
        for x, y in zip(xcommon, ycommon):
            fd.write( str(x) + "  " + str(y) + "\n")