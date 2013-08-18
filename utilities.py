 # -*- coding: utf-8 -*-
""" List of functions to help in general life """

import os, re
import pyfits
import collections
import datetime
import pytz
import sys
import numpy as np
from pyraf import iraf
import datetime
import astropy.io.fits as fits

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
    
def read_image_with_mask(image, mask_keyword=None):
    """ Read an image and a mask (from a keyword in the image), save it into a 
        numpy.ma array. The mask should contain 1 for pixels to be masked out. """
    data = fits.getdata(image)
    if mask_keyword:
        header = fits.getheader(image)
        mask_name = header[mask_keyword]
        mask = fits.getdata(mask_name)
    else:
        mask = np.zeros_like(data)
    return np.ma.array(data, mask=mask)         

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
        correct a block of images of the object with a certain blank or bias 
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
    iraf.module.load("astutil", doprint=0)
    screen = sys.stdout
    sys.stdout = open("temporal.txt", "w")
    iraf.module.asttimes(observatory=observatory, year=year, month=month, day=day, \
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
        outpt = os.path.join(outdir, prefix + outfile_root + "-" + \
                             suffix + out_extensions)                             
    elif prefix != '':
        outpt = os.path.join(outdir, prefix + outfile_root+ "-" + out_extensions)
    elif suffix != '':
        outpt = os.path.join(outdir, outfile_root + "-" + suffix + out_extensions)

    return outpt

def homogeneous_filter_name(filt):
    """ Find a common name for all those thousands of differents ways of 
       writting the filter names.  """
    filt = filt.lower()
    remove_characters = [" ", "/", "[", "]", "_"]
    for character in remove_characters:
        filt = filt.replace(character,"")
    filt_dict = {"rgu": "rGunn", "rgunn":"rGunn", "gunnr":"rGunn",\
                 "sdssr":"sdssr", "rsdss":"sdssr", "r":"R", \
                 "sdssi":"sdssi", "isdss":"sdssi",\
                 "h6607":"H6607", "h07":"H6607", "6607":"H6607", \
                 "h6652":"H6652", "h52":"H6652", "6652":"H6652", \
                 "h6650":"H6650", "h50":"H6650", "6650":"H6650", \
                 "h6678":"H6678", "h78":"H6678", "6678":"H6678", \
                 "clear":"Clear", "cle":"Clear", "clear":"Clear",\
                 "none":"Clear", "No":"Clear"  , "i":"I", "v":"V"}
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
                if re.search(pattern[key], filename) != None:
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
            for key in pattern.keys():
                if re.match(pattern[key], filename) != None:
                    match = re.search(pattern[key], filename)
                    match = match.groupdict()
                    # Find the name of the object 
                    if key == "cig": 
                        name = "cig" + str(format(int(match["cig_num"]), "04d"))
                    elif key == "standards":
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
    dictionaries. Pattern needs to be a dictionary where the result of the pattern generates
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
    im = pyfits.open(image, mode="update")
    hdr = im[0].header
    hdr.add_history(text)
    im.flush()
    im.close()

def check_sizes(imagelist):
    """ Uses pyraf to check if a list of images have the same size """
        