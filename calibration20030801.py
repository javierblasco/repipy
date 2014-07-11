# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.optimize as optimize
from repipy.calculate_extinction import calculate_extinction
from repipy.calculate_zeropoint import calculate_zeropoint
import repipy.utilities as utilities
import sys
import astropy.io.fits as fits

if len(sys.argv) != 2:
    print sys.exit("Give me a campaign file.... See example in the routine "+\
                   "packages")
execfile(sys.argv[1])

def search_images(dir="."):
    """ Given the patterns for the CIGs, cluster and standard images that we 
    :param dir: directory where data should be
    need, read them and sort them into a dictionary"""
    in_pattern = {}
    in_pattern["cig"] = "^(?P<name>cig)(?P<cig_num>\d{4})_(?P<date>\d{8})_" +\
                        "(?P<filt>.*)_(?P<exp_num>\d{3})(?P<rest>.*-c\.fits)$"
    in_pattern["bias"] = "^(?P<name>bias)_(?P<date>\d{8})_"+\
                         "(?P<exp_num>\d{3})(?P<rest>\.fits)$"
    in_pattern["skyflats"] = "^(?P<name>skyflat)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3})(?P<rest>.*-b\.fits)$"
    in_pattern["standards"] = "^(?P<name>" + standards_campaign +")_(?P<date>\d{8})_" +\
                              "(?P<filt>.*)_(?P<exp_num>\d{3})(?P<rest>.*-c\.fits)$"
    in_pattern["blanks"] = "^(?P<name>blank)_(?P<date>\d{8})_(?P<filt>.*)_" +\
                         "(?P<exp_num>\d{3})(?P<rest>.*-c\.fits)$"
    list_images = utilities.locate_images2(dir, in_pattern)
    return list_images

    
magnitudes = np.array( [[-11.619, -11.469, -11.162], # kop27, H6678, (265.9581, 5.4206)
                        [-11.452, -11.300, -10.997], # kop27, H6678, (265.9206, 5.3809)

                        [-9.514, -9.370, -9.074],    # kop27, H6678, (265.9096, 5.3395)
                        [-7.443, -7.312, -7.021],    # kop27, H6678, (265.9909, 5.3791)                        

                        [-7.459, -7.307, -7.002],    # kop27, H6678, (265.9374, 5.3193)                                                
                        [-11.561, -11.420, -11.117], # kop27, rGunn, (265.9784, 5.4226) 

 
                        [-9.244, -9.119, -8.776],    # kop27, rGunn, (265.9099, 5.4154)
                        [-9.270, -9.487, -9.469],    # bd+28, H6678, (327.7022, 28.8661)

                        [-8.264, -8.477, -8.456],    # bd+28, H6678, (327.7656, 28.8997)
                        [-6.710, -6.912, -6.938],    # bd+28, H6678, (327.8318, 28.8364)
                        
                        [-9.991, -10.242, -10.232],  # bd+28, rGunn, (327.7914, 28.9180)
                        [-10.286, -10.544, -10.534]  # bd+28, rGunn, (327.7039, 28.8450)
                        ] )  

err_magnitudes = np.array( [[0.002, 0.002, 0.002],
                            [0.002, 0.002, 0.002],

                            [0.004, 0.005, 0.006],
                            [0.016, 0.015, 0.022],
 
                            [0.016, 0.016, 0.023], 
                            [0.002, 0.002, 0.003],

 
                            [0.007, 0.008, 0.010],
                            [0.002, 0.003, 0.002],

                            [0.004, 0.005, 0.003],
                            [0.010, 0.011, 0.008],

                            [0.002, 0.003, 0.004],
                            [0.002, 0.003, 0.003]  ] )
                  
airmasses = np.array( [ [1.178, 1.735, 2.969],
                        [1.178, 1.735, 2.969],

                        [1.178, 1.735, 2.969],
                        [1.178, 1.735, 2.969],
 
                        [1.178, 1.735, 2.969], 
                        [1.177, 1.708, 2.876],

                        [1.177, 1.708, 2.876],
                        [1.955, 1.039, 1.020],

                        [1.955, 1.039, 1.020],
                        [1.955, 1.039, 1.020],

                        [2.035, 1.044, 1.018],
                        [2.035, 1.044, 1.018]   ] )

ext_coeff, sigma_ext_coeff = calculate_extinction(airmasses,magnitudes,\
                                                  err_magnitudes=err_magnitudes) 
print "Extinction coefficient for Halpha: ", ext_coeff, "+/-", sigma_ext_coeff, "\n"


# Playing around with different combinations of star fields and filters we 
# decide to use the value 0.253 +/- 0.006 for the extinction coefficient. 
ext_coeff = 0.253
sigma_ext_coeff = 0.006

    
# From the Halpha images of Kopff 27:
airmass_kp27_Halpha = np.array([1.178, 1.735, 2.969])
magnitudes_kp27_Halpha = np.array([-10.133, -9.987, -9.696]) -\
                                              ext_coeff * airmass_kp27_Halpha
                                              
# Assuming random errors, the error of the magnitude estimation and the one 
# from the extinction correction should add quadratically.                                               
err_mag_kp27_Halpha = np.array([0.003, 0.003, 0.004])
err_extinction_kp27 = np.array([sigma_ext_coeff] * len(magnitudes_kp27_Halpha))*\
                  airmass_kp27_Halpha
err_kp27_Halpha = np.sqrt(err_mag_kp27_Halpha**2 + err_extinction_kp27**2)                  
print "Magnitudes of Kopff 27 in Halpha:", magnitudes_kp27_Halpha
print "Error of the magnitudes:", err_kp27_Halpha, "\n\n"


# Now the fluxes of the standard star Kopff 27:
flux_kp27_Halpha = 1.05451e-11
m0_kp27_Halpha = - 2.5 * np.log10(flux_kp27_Halpha)


# From the Halpha images of bd+28
airmass_bd28_Halpha = np.array([1.955,1.039,1.020])
magnitudes_bd28_Halpha = np.array([-9.343,-9.547,-9.531]) -\
                                              ext_coeff * airmass_bd28_Halpha
err_mag_bd28_Halpha = np.array([0.002, 0.003, 0.002])                                              
err_extinction_bd28 = np.array([sigma_ext_coeff] * len(magnitudes_bd28_Halpha)) *\
                  airmass_bd28_Halpha 
err_bd28_Halpha = np.sqrt(err_mag_bd28_Halpha**2 + err_extinction_bd28**2)                  
print "Magnitudes of bd28:", magnitudes_bd28_Halpha
print "Error of the magnitudes:", err_bd28_Halpha, "\n\n"

# Fluxes
flux_bd28_Halpha = 6.30523e-12
m0_bd28_Halpha = -2.5 * np.log10(flux_bd28_Halpha)


# Concatenate all together
obs_magnitudes = np.array(list(magnitudes_kp27_Halpha) + list(magnitudes_bd28_Halpha))
err_mag = np.array(list(err_kp27_Halpha) + list(err_bd28_Halpha))
std_magnitudes = np.array([m0_kp27_Halpha] * 3 + [m0_bd28_Halpha] * 3)
zp, sigma_zp = calculate_zeropoint(obs_magnitudes, std_magnitudes, err_mag)
print "Zero point for Halpha:", zp, "+/-", sigma_zp

zp = 37.829
sigma_zp = 0.0003

list_images = search_images(directory)


print "Correct from atmosphere"
for index, im in enumerate(list_images["filename"]):
    if list_images["type"][index] in ["cig", "standards", "clusters"]:  
        airm, = utilities.get_from_header(im, airmassk)
        image = fits.open(im)
        image[0].data = image[0].data * 10**(ext_coeff * airm/2.5)
        newname = utilities.add_suffix_prefix(im, suffix="-e")
        utilities.if_exists_remove(newname)
        fits.writeto(newname, image[0].data, image[0].header)
        # Update header
        utilities.header_update_keyword(newname, airmassk, 0, "Corrected for atmosphere")
        list_images["filename"][index] = newname
        image.close()

print "Align images"
#print "Objects:", set(list_images["object"])

 
# Para ser más eficiente con astrometry.net 
#solve-field cig0895_20030801_H6678_003-b-f-c.fits --overwrite --no-fits2fits 
#--ra 315.2901 --dec 10.303 --radius 0.5 --no-plots --use-sextractor --resort 
#--depth 30 --code-tolerance 0.002