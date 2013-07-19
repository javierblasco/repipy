import numpy as np
import scipy.spatial as spatial
from scipy import optimize as optimize
import sys
import matplotlib.pyplot as plt
# Papers about matching used for this routines:
# Groth (1986), AJ 91, 1244
# Valdes et al (1995), PASP 107, 1119
#
# If imrovements in the voting matrix are needed, check:
# Marszalek & Roquita 2004, "Pattern matching with differential voting and median
#                            transformation derivation" (chapter of book "Computer
#                            vision and graphics")

def calculate_distances(xx, yy):
    """ Routine to calculate the distances between all the pairs of points. It 
    recieves the coordinates xx and yy as numpy.ndarrays, runs through the points 
    in a loop and calculates the distances of all pairs of objects. Gives back a 
    2D array in which the value [i,j] represents the distance between stars i and j. """
    nn = len(xx)   # Number of points
    distances = np.zeros([nn,nn], dtype=np.float64)  # matrix of distances

    # Run through all the stars, calculating the distances of the rest of points 
    # respect to that particular one. Then add to the matrix the distances. We 
    # avoid doing (1,2) and (2,1), which have the same distance. 
    for index_star in range(nn):
        dist_to_star = np.sqrt( (xx[index_star+1:] - xx[index_star])**2 +
                                (yy[index_star+1:] - yy[index_star])**2  )                            
        # Now, the matrix will be symmetric, so this distance represents both 
        # the row "index_star" (or part of it) AND the column "index_star". Thus:
        distances[index_star, index_star + 1:] = dist_to_star
        distances[index_star + 1:, index_star] = dist_to_star
    return distances


def create_triangles(distances):
    """ This routine creates all the possible triangles given a set of points.
        The matrix distances is a numpy array that contains the distances between
        the different points. So, distances[0,1] contains the distance between 
        points 0 and 1, and obviously distances[0,1] = distances[1,0]. 
        The triangles are represented as proposed by Valdes et al (1995), PASP 107
        and Stetson 1989 http://nedwww.ipac.caltech.edu/level5/Stetson/frames.html. 

        This means representing each triangle as a pair of coordinates: 
            x = b/a an y = c/a  
        where a, b and c are the sides of the triangle sorted by decreasing size.
        When using this representation, the (x,y) coordintes of a triangle is 
        independent of rotation, translation, flipping and scaling. 
        
        In order to facilitate the identification of triangles, we have added a 
        third parameter, which is just a, the longest of the three sides. The 
        division of the third parameter between matched triangles will give the 
        relative scale between the images.
            
        There will be nn * (nn-1) * (nn-2) combinations of three points, with nn
        the number of points. But the triangle formed by points (1,2,3), (2,1,3)
        (3,2,1) ... they are all the same triangle, so at the end only 
        nn * (nn-1) * (nn-2) /6 need to be calculated. 
        
        The three coordinates (x,y,a) of each triangle will be saved in an 
        array, while the indices of the stars forming each triangle are stored 
        in the array points_indices. These indices are always sorted as [C,B,A] 
        where A indicates the vertix point opposite the largest side and so on.
        This will allow 1 to 1 identification of the points of the triangles 
        when a match between two similar triangles are found. """
    nn = distances.shape[0]  # matrix is necessarily squared nn x nn.
    triangles = np.zeros([3, nn * (nn - 1) * (nn - 2) / 6])
    points_indices = np.zeros([3, nn * (nn - 1) * (nn - 2) /6], dtype=np.int64) 
    iter_index = 0
    for index1 in range(0, nn):
        for index2 in range(index1 + 1, nn):
            for index3 in range(index2 + 1, nn):
                sides = np.asarray([distances[index1, index2], 
                         distances[index2, index3], 
                         distances[index1, index3] ])
                points = np.asarray([index3, index1, index2])# points opposite 
                # Sort both sides and points
                indices_sorted_sides = np.argsort(sides)
                sides = sides[indices_sorted_sides] # Sorted increasing
                points = points[indices_sorted_sides] # also the points
                # Now we include the information (x,y, a) and (points) in the 
                # arrays that correspond
                triangles[:, iter_index] = [sides[1]/sides[2], sides[0]/sides[2],\
                                            sides[2]]
                points_indices[:, iter_index] = points
                iter_index = iter_index + 1                                          
    return triangles, points_indices  

def sort_by_first_index(input_array1, input_array2):
    """ Routine to sort two arrays in such a way that the elements of array[0,:]
        are sorted, and the elements of array[1,:] follow the others.
        The array input_array2 will also suffer the same modifications, so the 
        relationship between the arrays is kept safely. Example, the array:  
        [[1,4,2],[3,5,1]] would be sorted as [[1,2,4],[3,1,5]] , so that 1 is 
        still paired with 3, 4 with 5 and 2 with 1. The second array would be 
        also sorted using the order derived from the first one."""
    indices_sorted_array = np.argsort(input_array1[0,:])
    input_array1 = input_array1[:, indices_sorted_array]
    input_array2 = input_array2[:, indices_sorted_array]
    return input_array1, input_array2

def scale_rotation_translation(coords, scale, ang, deltax, deltay):
    """ Calculate translation, rotation and scaling necessary to match the 
        coordinates of one coordinate system coords_obj to a reference system,
        coords_ref. Receives a numpy array which contains [x_ref, y_ref, x_obj, 
        y_obj], and the parameters to fit. This function is called by 
        scipy.optimize.fit_curve below."""
    # First separate reference and object coordinates
    coords_ref = coords[0:2,:]
    coords_obj = coords[2:,:]
    # Convert coordinates to an affin matrix by adding 1s.
    len_ref = len(coords_ref[0,:])
    coords_ref = np.asarray([coords_ref[0,:], coords_ref[1,:], np.ones([len_ref])])    
    len_obj = len(coords_obj[0,:])
    coords_obj = np.asarray([coords_obj[0,:], coords_obj[1,:], np.ones([len_obj])])
    # Scale, rotation and translation matrices.    
    scalemat = np.asarray([[scale, 0, 0], [0, scale, 0], [0, 0, 1]])
    rotmat = np.asarray([[np.cos(ang), -np.sin(ang), 0],
                         [np.sin(ang),  np.cos(ang), 0],
                         [ 0 , 0 , 1 ] ])             
    transmat = np.asarray([[1, 0, deltax], [0, 1, deltay], [0, 0, 1]])        
    # First scale transformation, then rotation and finally translation
    transformed = np.dot(np.dot(transmat, rotmat), np.dot(scalemat, coords_obj))                 
    dist_sq = (coords_ref[0,:] - transformed[0,:])**2 +\
              (coords_ref[1,:] - transformed[1,:])**2
    return np.sqrt(dist_sq)

def calculate_scale(coords):
    """ Calculate the scale factor between two images using a set of common 
        points in both images. """
    # Coords contains the coordinates in reference and objective image:
    [xx_ref, yy_ref, xx_obj, yy_obj] = [coords[ii,:] for ii in range(4)]
      
    # Calculate distances between objects and the ratio between both images
    # Ideally that ratio would be a matrix with all elements = scale
    distances_ref = calculate_distances(xx_ref, yy_ref)
    distances_obj = calculate_distances(xx_obj, yy_obj)
    
    # Minimize errors: only the top 50% distances are used  
    large_dist = np.where((distances_ref > np.median(distances_ref)) &\
                          (distances_obj > np.median(distances_obj)))
    ratio = distances_ref[large_dist] / distances_obj[large_dist]

    # First estimate, use the largest 50% of the distances
    estimate = np.median(ratio)
    std_estimate = np.std(ratio)

    # Then we select those values within 3*sigma of the median
    cut_3sigma = np.where(abs(ratio - estimate) < 3 * std_estimate)
    scale = np.median(ratio[cut_3sigma])
    std_scale = np.std(ratio[cut_3sigma])
    
    plt.hist(ratio, bins=10)
    plt.show() 
    return scale, std_scale

def calculate_rotation(coords):
    """ Routine to calculate the rotation angle between two sets of points in 
        different reference systems, assumed both have the same scale. Also, if
        there is a flip between the images, we want to know about it. 
        The best rotation angle is that which minimizes the dispersion in the 
        distances between the objects in one frame respect to the same objects 
        in the other. Of course, this means that in the same step we could 
        actually calculate the translation, but we prefer to leave that to a 
        separate step.  
    """
    # Do a fit with no flipping of the object coordinate system.
    ang, ang_cov = optimize.curve_fit(minimize_angle, coords,
                                          np.zeros(len(coords[0,:])),
                                          maxfev=200)  
                                          
    # Check if it is a good fit by applying the rotation. The distance of the 
    # objects will be, in general, different from zero, since no translation has
    # yet been calculated. But the dispersion between those distances will be 
    # minimized if only a translation is necessary for both sets to match. 
    rot_matrix = np.asarray([[np.cos(ang[0]), -np.sin(ang[0])],
                             [np.sin(ang[0]), np.cos(ang[0])]], dtype=np.float64)
    coords_ref = coords[0:2,:]
    coords_obj = coords[2:,:]
    coords_obj = np.dot(rot_matrix, coords_obj)    
    distx = coords_obj[0,:] - coords_ref[0,:]
    disty = coords_obj[1,:] - coords_ref[1,:]
    dist = np.sqrt(distx**2 + disty**2)        

    # Now flip the coordinates in the object image and repeat the fit.     
    flip = np.asarray([[-1,0],[0,1]])
    coords_flipped = np.dot(flip, coords[2:,:])
    coords = np.asarray([coords[0,:], coords[1,:], 
                         coords_flipped[0,:], coords_flipped[1,:]])
    ang2, ang2_cov = optimize.curve_fit(minimize_angle, coords,
                                          np.zeros(len(coords[0,:])),
                                          maxfev=15200)  
                                          
    # And check again how good the rotation fits the data                                          
    rot_matrix = np.asarray([[np.cos(ang2[0]), -np.sin(ang2[0])],
                             [np.sin(ang2[0]), np.cos(ang2[0])]], dtype=np.float64)
    coords_obj = np.dot(rot_matrix, coords_flipped)
    distx = coords_obj[0,:] - coords_ref[0,:]
    disty = coords_obj[0,:] - coords_ref[1,:]
    dist2 = np.sqrt(distx**2 + disty**2)
    
    # The one with smaller dispersion in the distances is the good one
    if np.std(dist) < np.std(dist2):
        angle = ang
        std_angle = np.sqrt(ang_cov)            
        flip = False
    else:
        angle = ang2
        std_angle = np.sqrt(ang2_cov)
        flip = True  # flip is always expressed as np.asarray([[-1,0],[0,1]])
        
    return angle, std_angle, flip

def minimize_angle(coords, ang): #, xflip, yflip):    
    # Recover the two sets of coordinates from coords
    coords_ref = coords[0:2,:]
    coords_obj = coords[2:,:]
    
    # Build rotation matrix
    rot_matrix = np.asarray([[np.cos(ang), -np.sin(ang)],
                             [np.sin(ang), np.cos(ang)]], dtype=np.float64)
                             
    # Rotate object coordinates to reference ones                          
    coords_obj = np.dot(rot_matrix, coords_obj)  
    
    # In order to get a good estimate of the covariance matrix, we need to 
    # actually subtract the values for the translation. As we said, we don't 
    # want to save them here, just estimate them. 
    deltax = np.median(coords_obj[0,:] - coords_ref[0,:])
    deltay = np.median(coords_obj[1,:] - coords_ref[1,:])
    coords_obj = coords_obj - np.asarray([[deltax],[deltay]])

    # Finally, we compare the distances of the objects
    distx = coords_obj[0,:] - coords_ref[0,:]
    disty = coords_obj[1,:] - coords_ref[1,:]
    dist = np.sqrt(distx**2 + disty**2)      
    return dist

def transform_coordinates(coords_obj, scale=1.e0, angle=0., deltax=0., deltay=0., 
                          flip=False):
    """ Transform a set of coordinates into another system of coordinates given
        a scale, rotation angle around point (0,0) of the object system and 
        translation. 
        
        Default values are basically no transformation, so it is 
        useful even if only some of them are necessary. Since 
        the order of the transformations completely changes the final result, 
        we perform the individual transformations in this order: scaling (with 
        flip if flip = 1 or scale <0), rotation and translation. 
        Also, the flip is meant to produce a mirror image with the 
        y-axis as central axis. The difference between that ant a mirror in the 
        x-axis can be reproduced just by rotating 180 degree extra. A flip
        in each axis is also equivalent to no flip and 180 degree rotation, so 
        we will not consider that case either.
        
        In any case, if you calculated your corrections flipping about another 
        axis, or translating first and rotating later or things like that, then 
        you are screwed ;). You could always use this routine once for each of 
        the transformations, I guess. 
    """
    # If flip is active and scale is positive, then scale = -scale. This is 
    # because of the way we write the scale matrix calculated below.
    if flip == True and scale > 0:
        scale = -scale
    
    # Convert coordinates to an affin matrix by adding 1s. By writting it like 
    # this, the translation can be expressed as a multiplication of matrices.
    length = len(coords_obj[0,:])
    coords_obj = np.asarray([coords_obj[0,:], coords_obj[1,:], np.ones([length])])    

    # Scale, rotation and translation matrices.    
    scalemat = np.asarray([[scale, 0, 0], [0, abs(scale), 0], [0, 0, 1]])
    rotmat = np.asarray([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle),  np.cos(angle), 0],
                         [ 0 , 0 , 1 ] ])             
    transmat = np.asarray([[1, 0, deltax], [0, 1, deltay], [0, 0, 1]])
    
    # First scale transformation, then rotation and finally translation
    transformed = np.dot(np.dot(transmat, rotmat), np.dot(scalemat, coords_obj))                 
    return transformed[0:2,:]

def calculate_translation(coords):
        """ Given a set of coordinates [x0, y0, x1, y1], where all are numpy 
        arrays or lists of numbers, calculate the translation that converts from 
        (x1,y1) to (x0, y0). """
        # Calculate the differences between both sets of coordinates
        distx = coords[0,:] - coords[2,:]
        disty = coords[1,:] - coords[3,:]
        
        # First estimate
        delta = (np.median(distx), np.median(disty))
        sigma= (np.std(distx), np.std(disty))
        
        # Sigma clipping with 3 sigma cut
        good = np.where((abs(distx - delta[0]) < 3 * sigma[0]) &
                        (abs(disty - delta[1]) < 3 * sigma[1]))
                        
        # Second estimate of delta and sigma
        distx = distx[good]
        disty = disty[good]
        delta = (np.median(distx), np.median(disty))
        sigma= (np.std(distx), np.std(disty))
        return delta, sigma 
        

def main(xref="", yref="", xobj="", yobj="", error=0.01, scale="", angle="", 
         flip=False, test=False):
    """ Routine to calculate the transformation from a set of coordinates, 
        (xobj, yobj) to a reference system of coordinates (xref,yref). The routine
        follows the methods described in Groth (1986), AJ 91, 1244 and Valdes et 
        al (1995), PASP 107, 1119. There has to be some overlap between the two 
        sets of points, about 25% for about 25 stars, according to Groth (1986). 
        
        The transformation include all or any of: scaling, flipping, rotating and 
        translating. It will identify the stars by matching triangles, and use 
        the most probable matches to estimate the best transformation. 
        USAGE:
        scale, flip, angle, delta, precission = cros-match.main(xref=xref, yref=yref, 
                                                                xobj=xobj, yobj=yobj,
                                                                error=error, flip=flip,
                                                                scale=scale, angle=angle)
        INPUTS: 
            xref: array with the x values of the points in the reference system
            yref: same as xref for the y values
            xobj: array with the x values of the system to be converted to the 
                  reference system
            yobj: same as xobj for the y values
            error: triangles are described by two numbers, representing the ratio
                   between the largest side and the other two. For two triangles 
                   to be identified as the same, both numbers must agree to better
                   than error. 
            scale: Set this to a value if you do not want the routine to estimate
                   the scale. 
            angle: Set this keyword to a value if you do not want the routine to 
                   fit the rotation. 
            flip: Set this keyword to False if you do not wish the routine to 
                  calculate if a flip is needed. 
            test: Boolean keyword. If test = True, xref, yref, xobj and yobj will
                  be substituted by values from the test runs of this routine. 
                  The results should be a scale of 2, no flipping, rotation of 
                  0 degree and translations of -441 and 172 pixels. 
        OUTPUTS:
            scale: Tuple with the scale factor in the first element and the 
                   error of the scale in the second. The scale was
                   calculated dividing the distances between matched points in 
                   both systems of reference. The error corresponds to the 
                   standard deviation of that distribution of scales. 
            flip: Boolean variable. It will be True if the data in the objective 
                  system of coordinates need to be flipped about the y-axis before 
                  any rotation or translation is applied. 
            angle: Tuple with the angle of rotation (counterclokwise) and its 
                   error. This was calculated to minimize the dispersion in the 
                   distances between matching stars. 
            delta: Tuple of tuples. The first tuple corresponds to deltax and 
                   deltay, which needs to be added to the objective coordinate 
                   system AFTER the scaling, flipping and rotating have been 
                   performed. The second tuple contains the errors of the two 
                   parameters. 
            precission: A float which contains the median distance between the 
                   matched points, as an estimate of to which precission the 
                   transformation provides a good match. 
                   """
    if test == True:
        """ Objects obtained from images in Jan 2005, field hz15, images rGunn 001 
            for two consecutive days: 20050110 and 20050111. One of them has a 
            rebin 2x2, and has clear offsets respect to the other. Some of the 
            objects are in common, and both lists have objects that are not present 
            in the other.""" 
        xref = np.asarray([ 377,  421,   945, 1089,  433,  813,  1157, 1605, 
                            353,  681,   661,  221, 1357,  333,  1889, 1793,
                            475,  639,  1301,  481, 1237, 1187,  1149, 1905,
                            433,  591,   975, 1374,  925,  890,  540])
        yref = np.asarray([1524, 1732,  1672,  992,   52,  364,  1656, 1904, 
                           1648, 1480,   920,  676,  632,  336,  1736, 1032,
                           1817, 1687,  1293,  673,  857,  667,   462,  829,
                           188,   156,   164,  540,  914, 1639,  1557])
        xobj = np.asarray([ 409,  399,   431,  561,  693,  551,   763,  331, 
                            439,  795,   813,  841,  211,   93,   683,  799,
                            458,  490,   871,  461,  901,  628,   386,  106,
                            437,  516,   708,  907,  113,  540,   667])
        yobj = np.asarray([848,   912,   954,  826,  924,  548,   582,  424, 
                           114,   318,   422,  516,  668,  462,   546,  914,
                           996,   867,   735,  423,  403,  270,   255,  683,
                           182,   167,   170,  357,  232,  931,   905]) 
        xref = xref[0:18]
        yref = yref[0:18]
        xobj = xobj[0:18]
        yobj = yobj[0:18]    
                          
    n_ref = len(xref)
    n_obj = len(xobj)
    
    # Calculate the distances between all pairs of points in both sets of coordinates. 
    distances_ref = calculate_distances(xref, yref)
    distances_obj = calculate_distances(xobj, yobj)
    
    # Create all possible triangles with any three points for each set. Define 
    # each triangle by two quantities x=b/a and y=c/a where (a,b,c) are the sides
    # of the triangle in decreasing order. Triangle[0,:] contains all the x while 
    # Triangle[1,:] contains y. Notice that x and y are independent of rotation, 
    # translation, flipping and scale changes. The array called indices_ contains
    # the points that formed each triangle, sorted as [C,B,A] where C is the point
    # opposite to side c and so on.
    triangle_ref, indices_ref = create_triangles(distances_ref)
    triangle_obj, indices_obj = create_triangles(distances_obj)
    
    # Sort the triangles coordinates and the array of indices according to the 
    # first of the two coordinates in triangle space, x. This will require to 
    # sort triangle[0,:], and apply the same order to triangle[1,:] and 
    # indices[0,:], indices[1,:] and indices[2,:]
    triangle_ref, indices_ref = sort_by_first_index(triangle_ref, indices_ref)
    triangle_obj, indices_obj = sort_by_first_index(triangle_obj, indices_obj)
    
    # Points very close together an image will produce that many of the triangles 
    # between those two points and a third one will have x ~ 1 and y ~ 0. 
    # Distinguishing them is very difficult, so it is best to remove them. Notice
    # that we will not remove a triangle in which the third object will also be 
    # very close. 
    valid_triangles = np.where(triangle_ref[0,:] < 0.9)[0]
    triangle_ref = triangle_ref[:, valid_triangles]
    indices_ref = indices_ref[:, valid_triangles]
   
    valid_triangles = np.where(triangle_obj[0,:] < 0.9)[0]
    triangle_obj = triangle_obj[:, valid_triangles]
    indices_obj = indices_obj[:, valid_triangles]
    
    # And now we convert into arrays of tuples of coordinates. Now the arrays
    # contain tuples with the two numbers that define each triangle. If you wonder
    # why we didn't start by doing tuples in the first place, the answer is that 
    # it was easier to sort this way. 
    coords_obj = zip(triangle_obj[0,:], triangle_obj[1,:])
    coords_ref = zip(triangle_ref[0,:], triangle_ref[1,:]) 
    
    
    # We search for matches between the triangles using kdtree (k-dimensional tree)
    # from scipy.spatial.
    tree_ref = spatial.KDTree(coords_ref)
    tree_obj = spatial.KDTree(coords_obj)
    matches = tree_ref.query_ball_tree(tree_obj, error)
    #for ref_id, obj_id in enumerate(matches):
    #    if obj_id != []:
    #        print "Matched triangles:", ref_id, obj_id
    #        print coords_ref[ref_id], coords_obj[obj_id[0]]
    
    # And now, following the method by Groth 1986 (PASP, 107, 1119), we are going 
    # to build a "voting matrix" with the stars of the reference image in one 
    # axis of the matrix and those of the object image in the other. 
    # For every matched triangle, three stars
    # are identified. Assuming star "i" in one image is identified as 
    # star "j" in the other for a particular triangle match, the matrix increases 
    # linked together, at the end that element of the matrix will have long number
    # as compared with random missmatches. 
    voting_matrix = np.zeros([n_ref, n_obj], dtype = np.int)    
    for index_ref, index_list in enumerate(matches): # for every triangle in tree_ref
        for index_obj in index_list:  # run through all the matches, ideally 1.
            for ii in range(3):     # each of the points of a triangle
                point_ref = indices_ref[:, index_ref][ii]
                point_obj = indices_obj[:, index_obj][ii]
                voting_matrix[point_ref, point_obj] += 1                

    # Since false matches are a little bit irksome, we estimate a rough
    # value of the scale of the image using the four matches with most votes.
    # Then we will repeat the process of the voting matrix not using anything 
    # that differs more than error from the estimated value.
    minimum = np.sort(voting_matrix.flatten())[-4]    
    largest_votes = np.where(voting_matrix >= minimum)#
    scales = []
    for ii in range(len(largest_votes[0])):
        xhighest_ref = xref[largest_votes[0][:]]
        yhighest_ref = yref[largest_votes[0][:]]
        xhighest_obj = xobj[largest_votes[1][:]]
        yhighest_obj = yobj[largest_votes[1][:]]
        distances_highest_ref = calculate_distances(xhighest_ref, yhighest_ref)
        distances_highest_obj = calculate_distances(xhighest_obj, yhighest_obj)
        non_zero = np.where((distances_highest_ref !=0) & (distances_highest_obj !=0))
    scale_estimate = np.median(distances_highest_ref[non_zero] /
                               distances_highest_obj[non_zero])
    scale_std_estimate = np.std(distances_highest_ref[non_zero] /
                                distances_highest_obj[non_zero])

    # And we repeat the voting matrix for those triangles with simimlar scale
    voting_matrix = np.zeros([n_ref, n_obj], dtype = np.int)    
    for index_ref, index_list in enumerate(matches): # for every triangle in tree_ref
        for index_obj in index_list:  # run through all the matches, ideally 1.
            sc = triangle_ref[2, index_ref] / triangle_obj[2, index_obj]
            if abs(sc - scale_estimate) < 3 * scale_std_estimate:          
                for ii in range(3):     # each of the points of a triangle
                    point_ref = indices_ref[:, index_ref][ii]
                    point_obj = indices_obj[:, index_obj][ii]
                    voting_matrix[point_ref, point_obj] += 1                
    print voting_matrix
       
    # Now we will accept as valid those matches that are 3 sigma above the "noise" 
    # of the voting matrix. Some might still be missmatches, but the process to 
    # calculate rotation, translation, flipping and scaling will cope well with 
    # a handful of mismatches.
    median_matrix = np.median(voting_matrix)
    stddev_matrix = np.std(voting_matrix)
    matched_objects = np.where(voting_matrix > median_matrix + 3. * stddev_matrix) 
        
        
    # Now we start a new part of the process . We want to calculate the 
    # transformation between the two coordinate systems. For that we can start 
    # with a fresh set of coordinates containing the objects that are supposedly 
    # matched:
    xmatched_ref = xref[matched_objects[0]]
    ymatched_ref = yref[matched_objects[0]]
    xmatched_obj = xobj[matched_objects[1]]
    ymatched_obj = yobj[matched_objects[1]]
    coords = np.asarray([xmatched_ref, ymatched_ref, xmatched_obj, ymatched_obj])

                 
    #***************************************************************************
    # We calculate affin transformations for these two coordinate systems, including 
    # scaling, rotation around (0,0) in the objective reference system and 
    # translation to be performed to make coincide both images.
    # result, cov_result = optimize.curve_fit(scale_rotation_translation, coords, 
    #                                         np.zeros(len(xmatched_ref)),
    #                                         p0=np.array([2,0,0,0]))
    # print result                                               
    # print np.sqrt(cov_result[0,0]), np.sqrt(cov_result[1,1]),\
    #       np.sqrt(cov_result[2,2]), np.sqrt(cov_result[3,3])
    #***************************************************************************
    
    
    # Now we calculate the scale factor if user has not provided one
    if scale == "":
        scale, scale_std = calculate_scale(coords)
    else:
        scale, scale_std = scale, None
 
    # Recalculate coords to correct for scale. Note that flipping is not yet
    # taken into account. This will be done when calculating rotation.
    coords[2:,:] = transform_coordinates(coords[2:,:], scale=scale)
                         
    # Fit for the angle if user has not provided one:
    if angle != "":
        angle, angle_std, flip = calculate_rotation(coords)
        angle = angle
    else:
        angle, angle_std, flip = calculate_rotation(coords)
    
    # Now we calculate the rotation and flip if present.
    coords[2:,:] = transform_coordinates(coords[2:,:], angle=angle, flip=flip)
    
    # And now calculate the translation of the resulting coordinate systems.
    delta, delta_std = calculate_translation(coords)
    
    # Finally, check how accurate the transformation is:
    coords[2:,:] = transform_coordinates(coords[2:,:], deltax=delta[0], 
                                         deltay=delta[1])
    dist = np.sqrt((coords[0,:] - coords[2,:])**2 +
                  (coords[1,:] - coords[3,:])**2) 
    print "dist:", distx
    precission = np.median(dist)
    return (scale, scale_std), flip, (angle, angle_std), (delta, delta_std), precission
        
if __name__ == "__main__":
  main(test=True)