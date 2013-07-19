# -*- coding: utf-8 -*-
import jbh_utilities as jbh
import re, shutil, os
import collections
import numpy as np

def tidy_up(directory, pattern, date):
    """
        Routine to tidy up the directory we are going to analyze. It also returns
        a dictionary with three numpy arrays. 
            1. final_dict["filename"] contains the names of the files found. 
            2. final_dict["type"] contains the type of image (cig, standard, bias,
               flat...depends on the name of the patterns introduced). 
            3. final_dict["objname"] contains the specific name of the object 
               (e.g. cig0345, feige34, flat, NGC2030). For this the pattern must
               have identified the name of the object. 
        Examples of patterns: 
                
           in_pattern = {}
           in_pattern["cig"] = "^(?P<name>cig)(?P<cig_num>\d{3,4})-"+\
                               "(?P<exp_num>\d{3})(?P<filt>.{3})(?P<rest>\.fit)"
           in_pattern["bias"] = "^(?P<name>bias)-(?P<exp_num>\d{3})(?P<rest>\.fit)"    
    """
        
    # List all files in directory and extract pattern information
    lista = jbh.list_dir(directory)
    dictionary_images = jbh.name_pattern(pattern, lista)

    # Empty array of dtype=object is the best way I found to be able to save 
    # variable length strings in an array. With dtype=string_ it just truncates
    # the strings if you try to make them longer than the original size.
    empty_array = np.asarray([], dtype=object)
    final_dict = {"filename":empty_array, 
                  "type":empty_array, 
                  "objname":empty_array}
    dirs = []  # keep track of the folders created
    # Loop through the items that corresponded to a pattern
    for filename, fileinfo in dictionary_images.items():
        # For each type find the name of the object and the folder to move to 
        if fileinfo["key"] == "cig": 
            name = "cig" + str(format(int(fileinfo["cig_num"]), "04d"))
            folder = name 
        elif fileinfo["key"] == "standards":
            name = fileinfo["name"]
            folder = "standards"
        else:
            name = fileinfo["key"]
            folder = fileinfo["key"]

        # Create the folders if not present
        if folder not in dirs:
            jbh.create_dirs([folder], directory)
            dirs.append(folder)

        # Find filter
        try:  
            filt = jbh.homogeneous_filter_name(fileinfo["filt"])
        except KeyError:   # bias has no filter!
            filt = "none"
            
        # Build new name    
        newname = os.path.join(directory, folder, name) + "_" + date + "_" + \
                             filt + "_" + fileinfo["exp_num"] + "_" + fileinfo["rest"]
        newname = newname.replace("_.", ".")   # blabla_001_.fit is painful to see
        if re.search("fit(?!s)", newname):  # if .fit and not .fits
            newname = newname.replace(".fit", ".fits")
        # Add file name, type and objname to the dictionary
        final_dict["filename"] = np.append(final_dict["filename"], newname)         
        final_dict["type"] = np.append(final_dict["type"], fileinfo["key"])       
        final_dict["objname"] = np.append(final_dict["objname"], name)        
        
        # Now move and update the header adding a HISTORY comment. 
        shutil.move(filename, newname)
        text = "File renamed: " + os.path.split(filename)[1] + " --> " + os.path.split(newname)[1]  
        jbh.add_history_line(newname, text)                
    return final_dict