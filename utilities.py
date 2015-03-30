 # -*- coding: utf-8 -*-
""" List of functions to help in general life """

import os
import re
import collections
import datetime
import pytz
import sys
import numpy as np
import datetime
import astropy.io.fits as fits
import astropy.units as u
import astropy.coordinates as coords
from astropy.time import Time
import dateutil.parser
import shutil
import functools
from lemon import methods

import repipy
# Change to the directory where repipy is installed to load pyraf
with methods.tmp_chdir(repipy.__path__[0]):
    from pyraf import iraf
    from iraf import astutil

def sex2deg(RA, DEC):
    try:
        hh, mm, ss = RA.split(":")
        RA = (float(hh) + float(mm) / 60. + float(ss) / 3600.) * 15
    except ValueError, AttributeError:
        pass
    try:
        dd, mm, ss = DEC.split(":")
        DEC = (float(dd) + float(mm) / 60. + float(ss) / 3600.)
    except ValueError, AttributeError:
        pass
    return RA, DEC


def number_of_chips(hdu_list):
    """ From an astropy fits object (i.e. an HDUList), find out how many chips form the image. INT has 4 chips per image,
        VST 32."""
    HDU_components = fits.HDUList(hdu_list)
    return sum([isinstance(hdu, (fits.ImageHDU, fits.hdu.compressed.CompImageHDU)) for hdu in HDU_components])

def memoize(f):
    """ Minimalistic memoization decorator.
    http://code.activestate.com/recipes/577219-minimalistic-memoization/ """

    cache = {}
    @functools.wraps(f)
    def memf(*x):
        if x not in cache:
            cache[x] = f(*x)
        return cache[x]
    return memf

def collect_from_images(image_list, keyword):
    """ From a list of images collect a single keyword """
    try:
        return [fits.getval(im, keyword) for im in image_list]
    except KeyError:
        sys.exit("Keyword %s does not exist in image %s" % (keyword, im))

def remove_WCS(header):
    """ Wipe any sign of a WCS in a header object """
    hdr = header.copy()  # To avoid modifying the original header 
    wcs_keywords = ["wcsaxes", "ctype1", "ctype2", "equinox", "lonpole", 
              "latpole", "crval1", "crval2", "crpix1", "crpix2", 
              "cunit1", "cunit2", "cd1_1", "cd1_2", "cd2_1", "cd2_2", 
              "PC001001", "PC001002", "PC002001", "PC002002"]
    for keyword in wcs_keywords:
        try:
            hdr.remove(keyword)
        except:
            pass
    return hdr

def update_WCS(image_without_wcs, image_with_wcs):
    """ Export the WCS information from an image to another one, by updating 
        certain keywords in the header of the target image"""
    image = fits.open(image_without_wcs, mode='update')
    hdr_im = image[0].header
    hdr_wcs = fits.getheader(image_with_wcs)
    
    # WCS parameters are:
    change = ["wcsaxes", "ctype1", "ctype2", "equinox", "lonpole", 
              "latpole", "crval1", "crval2", "crpix1", "crpix2", 
              "cunit1", "cunit2", "cd1_1", "cd1_2", "cd2_1", "cd2_2"]          
    for key in change:
        hdr_im[key] = hdr_wcs[key]

    # Remove unnecessary and confusing ones if they exist
    delete_these = ["PC001001", "PC001002", "PC002001", "PC002002"]
    for key in delete_these:
        del hdr_im[key]

    image.flush()
    image.close()
    return None

def get_from_header(image_name, *args):
    """ From the header of an image, get the values corresponding to the 
        keywords passed in args"""
    if len(args) == 1:
        return fits.getval(image_name, args[0])
    else:
        return (fits.getval(image_name, x) for x in args)

def precess_to_2000(RA, DEC, time):
    """ From the actual coordinates of an object in the sky for a certain 
        time, recalculate the J2000 coordinates. This is util to look 
        within catalogues. RA and DEC are in degree, time is a 
        datetime.datetime object (or will be converted into it at the 
        begining of the program)."""
    try:
        time = dateutil.parser.parse(time)
    except: 
        pass
    fk5 = coords.FK5(ra=RA, dec=DEC, unit=(u.degree, u.degree),
                                equinox=Time(time.year, format="jyear", 
                                             scale="utc"))
    fknew = coords.FK5(equinox='J2000')
    fk_2000 = fk5.transform_to(fknew)
    return fk_2000.ra.degree, fk_2000.dec.degree

def check_dimensions(image_list):
    """ Check that fits images in a list have all the same dimensions, so that 
        operations can be performed. If one is not a file (or can not be read), 
        check if it is a float. If it is a float, operations can be performed, 
        so ignore the exception. Otherwise, raise it."""
    # Check and save dimensions of every image
    dimensions = []
    for image in image_list:
        try:  # try to open
            dimensions.append(fits.open(image, mode="readonly")[0].shape)
        except IOError:
            try: #check if float
                float(image)
            except ValueError: # if not, reraise original error
                raise IOError ("[Errno 2] No such file or directory: " + "'" +\
                               image + "'")
    # If len(set(dimensions))  == 1 return True (all images are the same size)       
    if len(set(dimensions)) == 1:
        return True
    else: # print sizes of images
        print "\n Sizes of images are not the same! \n"
        print "Image: " + "\t" * 4 + "Size"    
        for image, size in zip(image_list, dimensions):
            print image, size        
        return False

def read_image_with_mask(image, mask_keyword=None, limits = 0, header=None):
    """ Read an image and a mask (from a keyword in the image), save it into a numpy.ma array.
    The mask should contain 1 for pixels to be masked out. Limits allows to read only a part of an image, avoiding
    the need to read it all into memory. Limit should be limits = (min_x, min_y, max_x, max_y), use limit=0 to use all
    image (DEFAULT: all image)."""

    # Set the limits
    header = fits.getheader(image)
    if limits == 0:
        min_x, min_y = (0, 0)
        max_x, max_y = header["NAXIS2"], header["NAXIS1"]
    else:
        min_x, min_y, max_x, max_y = limits


    # Open both image and mask (if present)
    with fits.open(image, memtype=True) as im:
        data = im[0].data[min_x:max_x, min_y:max_y]
        if mask_keyword:
            mask = fits.open(header[mask_keyword], memtype=True)[0].data[min_x:max_x, min_y:max_y]
        else:
            mask = np.zeros_like(data)
    return np.ma.array(data, mask=mask, dtype=np.float64)

def mean_datetime(datetimes):
    """ This function returns the average datetime from a given set of datetime 
        objects. We have to calculate the differences between each time and 
        the first value, then calculate the average of those differences and 
        add that to the first value. """
    delta_times = [(times - datetimes[0]).total_seconds() for times in datetimes]
    average_delta = sum(delta_times) / float(len(delta_times))
    average_time = datetimes[0] + datetime.timedelta(0,average_delta)
    return average_time
    
def group_images_in_blocks(times, limit=3):
    """ In a night at the telescope, we can observe blocks of images on the same 
        field of the sky. For example, five blank fields, then the object, then 
        another 4 blank fields, then the object... We might want to distinguish 
        those blocks, in order to, e.g., combine the blank fields of each block, 
        correct a block of imafges of the object with a certain blank or bias
        image... This routine gets the datetime.datetime objects that indicate
        the date and time of observations, and separates them in blocks, giving 
        back an array:
            indices = [ind0,ind1,ind2,ind3]
        so that, each slice [ind0:ind1], [ind1:ind2] and [ind2:ind3] gives a block
        of the incoming images """
    delta_times = np.asarray( [(times[ii+1] - times[ii]).seconds for ii in 
                               range(len(times)-1)] )
    # Median and median absolute deviation of delta_times as first guesses
    median_delta = np.median(delta_times)
    MAD = np.median(abs(delta_times - median_delta))

    # We have found a boundary between blocks it the time between images is larger
    # than limit. The limit will be assigned to the upper image, hence the +1
    block_limits = np.where(delta_times > median_delta + limit * MAD)[0] +1
    
    # Now we have the limits between blocks, we need to add the first image (where 
    # the first block starts) and the last image (where last block ends)
    block_limits = np.insert( np.append(block_limits, len(times)+1) , 0, 0) 
    return block_limits

def universal_time_to_local_time(date, time, location):
    """ Routine to convert universal time to local time for a given location. 
        The input consists of three strings with the date, time and location, 
        with the format "yyyy-mm-dd", "hh:mm:ss" and location (e.g. 
        Europe/Madrid) """
    st_tz = pytz.timezone(location)
    date_time = date + " " + time
    utc_dt = datetime.datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S')
    utc_dt = utc_dt.replace(tzinfo=pytz.utc)
    st_dt = st_tz.normalize(utc_dt.astimezone(st_tz))
    return st_dt.strftime("%H:%M:%S")        

def local_to_sidereal_time(date, time, observatory):
    """ From date, time and observatory name (e.g. OSN) calculate the sidereal 
        local time """
    year, month, day = date.split("-")
    hour, minute, second = time.split(":")   
    # iraf's routine asttimes needs the time in hours!      
    time_hours = float(hour) + float(minute)/60. + float(second)/3600.
    screen = sys.stdout
    sys.stdout = open("temporal.txt", "w")
    iraf.asttimes(observatory=observatory, year=year, month=month, day=day, \
                  time=time_hours)
    sys.stdout = screen
    for lines in open("temporal.txt"):
        if len(lines) > 2:
            sidereal_time = lines.split()[-1]
    os.remove("temporal.txt")
    return sidereal_time

def add_suffix_prefix(filename, prefix='', suffix=''):
    (outdir, outfile) = os.path.split(filename)   
    """ Routine to add a prefix or a suffix to the name of a file.
      E.g.: /home/user/text.doc might become /home/user/text-b.doc """    
    
    # Separate the extension(s).
    outfile_root = (outfile.split(os.extsep))[0]
    out_extensions = ""
    for extension in outfile.split(os.extsep)[1:]:
        out_extensions += "." + extension
        
    # Now construct the output file including the prefix/suffix if present
    if prefix != '' and suffix != '':
        outpt = os.path.join(outdir, prefix + outfile_root + \
                             suffix + out_extensions)                             
    elif prefix != '':
        outpt = os.path.join(outdir, prefix + outfile_root + out_extensions)
    elif suffix != '':
        outpt = os.path.join(outdir, outfile_root + suffix + out_extensions)

    return outpt
    
def replace_extension(filename, new_extension):
    """ From the name of a file (possibly full path) change the extension 
        by another one given by user"""
    if new_extension[0] != ".":
        new_extension = "." + new_extension
    path, name = os.path.split(filename)
    name_root, name_ext = os.path.splitext(name)
    return os.path.join(path, name_root + new_extension)

def homogeneous_filter_name(filt):
    """ Find a common name for all those thousands of differents ways of 
       writting the filter names.  """
    filt = filt.lower()
    remove_characters = [" ", "/", "[", "]", "_"]
    for character in remove_characters:
        filt = filt.replace(character,"")
    filt_dict = {"rgu": "rGunn", "rgunn":"rGunn", "gunnr":"rGunn",\
                 "johv":"VJoh", "vjoh":"VJoh", "johnv":"VJoh",\
                 "b":"B",\
                 "johr":"RJoh", "rjoh":"RJoh", "johnr":"RJoh",\
                 "cousr":"RCous", "rcous":"RCous", \
                 "sdssr":"sdssr", "rsdss":"sdssr", "r":"R", \
                 "sdssi":"sdssi", "isdss":"sdssi",\
                 "sdssg":"sdssg", "gsdss":"sdssg",\
                 "h6607":"H6607", "h07":"H6607", "6607":"H6607", \
                 "ha6607":"H6607",\
                 "h6652":"H6652", "h52":"H6652", "6652":"H6652", \
                 "h6633":"H6633", "6633":"H6633",\
                 "h6625":"H6625", "6625":"H6625", \
                 "h6650":"H6650", "h50":"H6650", "6650":"H6650", \
                 "h6678":"H6678", "h78":"H6678", "6678":"H6678", \
                 "ha6678":"H6678",\
                 "h6613":"H6613", "6613":"H6613",\
                 "67418":"H6741", "h6741":"H6741",\
                 "clear":"Clear", "cle":"Clear", "clear":"Clear",\
                 "free":"Clear",\
                 "none":"Clear", "No":"Clear" , "i":"I", "v":"V"}
    try:
        return filt_dict[filt]
    except KeyError:
        print "len(filt)", len(filt)
        print "\n\n Probably you have not defined the filter " + filt + \
              " in the dictionary! Include it in homogeneous_filter_name"+\
              " in repipy/utilities.py .  \n \n"
        sys.exit("Exiting program")

def locate_images(directory, pattern):
    """ Given a directory, save all the files that fit any of the patternsof a 
        dictionary. Save in a dictionary, whose keys are the same as the ones of
        patterns. """
    list_files = collections.defaultdict(list)
    for dd, ss, ff in os.walk(directory):
        for filename in ff:
            for key in pattern.keys():
                if re.search(pattern[key], filename, re.I):
                    list_files[key].append(os.path.join(dd, filename))
    return list_files
        
def locate_images2(directory, pattern):
    """ 
        Given a directory, save all the files that fit any of the patterns of a 
        dictionary. Save in a dictionary, whose keys are 3 numpy arrays that 
        contain the filenames, type (cig, standards, flats) and the object name
        according to the given patterns. 
    """
    # Empty array of dtype=object is the best way I found to be able to save 
    # variable length strings in an array. With dtype=string_ it just truncates
    # the strings if you try to make them longer than the original size.
    empty_array = np.asarray([], dtype=object)   
    final_dict = {"filename":empty_array, 
                  "type":empty_array, 
                  "objname":empty_array
                   }
    for dd, ss, ff in os.walk(directory):
        for filename in ff:
            #print "Checking if image {0} follows a pattern".format(filename)
            for key in pattern.keys():
                if re.match(pattern[key], filename, re.I):
                    #print "It does!"
                    match = re.search(pattern[key], filename, re.I)
                    match = match.groupdict()
                    # Find the name of the object 
                    if key == "cig": 
                        name = "cig" + str(format(int(match["cig_num"]), "04d"))
                    elif key == "bias":
                        name = "bias"
                    elif key == "standards":
                        name = match["name"]
                    elif key == "skyflats":
                        name = "skyflat"
                    elif key == "domeflats":
                        name = "domeflat"
                    elif key == "blanks":
                        name = "blank"
                    elif key == "flats":
                        name = "flat"
                    elif key == "clusters":
                        name = match["name"]
                    else:
                        name = key

                    final_dict["filename"] = np.append(final_dict["filename"],\
                                                       os.path.join(dd, filename))
                    final_dict["type"] = np.append(final_dict["type"], key)
                    final_dict["objname"] = np.append(final_dict["objname"], \
                                                      name)
    return final_dict

def create_lists(pattern, directory):
    """ Given a list of patterns for re.search make a list of files according 
    to them in a certain directory. """  
    lista = {}
    for key in pattern.keys():
        lista[key] = []
    for f in os.listdir(directory): 
        for key in pattern.keys():
            filename = os.path.join(directory,f)
            if re.search(pattern[key], filename, re.I) != None:
                lista[key].append(f)
    return lista

def create_dirs(dirs, directory):
    """ Create in a directory the list of subfolders contained in dirs"""
    for subfolder in dirs:
        if os.path.isdir(os.path.join(directory, subfolder)) == False:
            os.makedirs(os.path.join(directory,subfolder))

def name_pattern(pattern, lista):
    """ Extract info from the patterns in a list of filenames and saves it in 
    dictionaries. Pattern needs to be a dictionary where the result of the 
    pattern generates
    a dictionary (e.g: "cig(?P<cig_num>\d{3,4}).fits" would be a valid pattern, 
    "cig(\d{3,4}).fits" would not). lista be a list of names. """
    dictionary = {}
    for filename in lista:
        directory, name = os.path.split(filename)
        for key in pattern.keys():
            result = re.search(pattern[key], name, re.I)
            if result != None:
                dict_result = result.groupdict()
                dict_result["key"] = key
                dictionary[filename] = dict_result
    return dictionary

def list_dir(directory):
    """ Routine to list a dir with the full path """
    lista = os.listdir(directory)
    lista = [os.path.join(directory, name) for name in lista]  
    return lista
        
def add_history_line(image, text):
    """ Add a history line to the image with the text given """
    im = fits.open(image, mode="update")
    hdr = im[0].header
    hdr.add_history(text)
    im.flush()
    im.close()

def header_update_keyword(image, keyword, value, comment=""):
    """ Update a header keyword, or create it if not present """
    im = fits.open(image, mode="update")
    hdr = im[0].header
    hdr.update(keyword, value, comment)
    im.flush()
    im.close()

def read_from_sextractor_catalogue(filename, *keys):
    """ Read from a sextractor catalogue the given keys, for example 
        keys can be "X_IMAGE", "Y_IMAGE", "MAG_AUTO" """
    indices = [0] * len(keys)
    with open(filename, "r") as f:
        line = f.readline().split()
        # Read which row the keys have
        while line[0] == "#": 
            if line[2] in keys:
                indices[keys.index(line[2])] = int(line[1]) - 1
            line = f.readline().split()
        # If not all keys found
        if 0 in indices:
            which_are_zero = [keys[ii] for ii in indices if indices == 0]
            sys.exit("Impossible to read sextractor catalog, some keys are "+\
                     "not present: " + ", ".join(which_are_zero))
                 
        values = [np.asarray([])] * len(keys) 
        while line != []:  # until end of file
            for ii in range(len(keys)):
                values[ii] = np.append(values[ii], float(line[indices[ii]]))
            line = f.readline().split()
        return values
         
def if_exists_remove(*filename):
    """ Check if a file exists. If so, remove it"""
    for f in filename:
        if os.path.isfile(f):
            os.remove(f)

def if_dir_not_exists_create(folder):
    """ If a directory does not exist, create it. Recursively create any intermediate subfolders. """
    abs_path = os.path.abspath(folder)
    previous_folder = os.path.split(abs_path)[0]
    if not os.path.isdir(folder):  # Folder does not exist
        # Create the previous folder first
        if_dir_not_exists_create(previous_folder)
        os.mkdir(folder)



def move_list(file_list, target_dir):
    """ Move each of the elements of a list to a given folder """ 
    for item in file_list:             
        shutil.move(item, target_dir)
    return None
    

