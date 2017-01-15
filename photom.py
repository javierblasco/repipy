import matplotlib.pyplot as plt
import pyregion
import ds9
import astropy.io.fits as fits
import numpy as np
import sys
import Tkinter
import tkMessageBox
import repipy.utilities as utils
import repipy.extract_stars as extract_stars
import astropy.wcs.wcs as wcs

class myWindow:

    def __init__(self, message):
        self.mw = Tkinter.Tk()
        self.mw.withdraw()
        self.mw.option_add("*font", ("Arial", 15, "normal"))
        self.mw.geometry("+250+200")
        self.message = message
        
    def box(self):
        box = tkMessageBox.showinfo("Action required", self.message)
        
  
        
def mask_from_ds9(image_name, message):
    # Display image
    d = ds9.ds9()
    d.set("file " + image_name)

    # Tell the user to select the regions
    answer = False
    while answer == False:
        window = myWindow(message)
        window.box()
        try:
            region = pyregion.parse(d.get("region"))
            answer = True
        except ValueError:  # No regions defined
            print "\n There is no region loaded in the image! \n"
    return region

def positions_of_stars_from_ds9(image_name, message, catalogue_name):
    # Display image
    d = ds9.ds9()
    d.set("file " + image_name)

    # Tell the user to select the regions
    answer = False
    while answer == False:
        window = myWindow(message)
        window.box()
        try:
            region = pyregion.parse(d.get("region"))
            answer = True
        except ValueError:  # No regions defined
                print "\n There is no region loaded in the image! \n"

    coords_RADEC = [cc.coord_list[0:2] for cc in region]

    #hdr = fits.getheader(image_name)
    #w = wcs.WCS(hdr)
    #coords_RADEC = w.all_pix2world(xyout,1)
    with open(catalogue_name, 'w') as fd:
        for reg in coords_RADEC:
            fd.write(" {0}  {1} \n".format(*reg))
    return region


def main(Ha_name, rgunn_name, scaling_factor, zp_Ha, T_Ha, T_rGunn):
    # Tell the user to load regions for the galaxy
    galaxy_messg = "\n Create/load region(s) in ds9 to enclose the galaxy. "+\
              " Save it for latter use. Then hit the 'OK' button!"
    galaxy_region = mask_from_ds9(rgunn_name, galaxy_messg)

    # Sky regions
    sky_messg = "\n Now do the same for sky region (or regions!) "
    sky_region =   mask_from_ds9(rgunn_name, sky_messg)

    # Stars to model PSF
    mod_stars_catalogue = utils.replace_extension(rgunn_name, ".model_stars")
    mod_stars_messg = "\n Select nice isolated non-saturated stars to model the PSF. "
    mod_stars_region = positions_of_stars_from_ds9(rgunn_name, mod_stars_messg, mod_stars_catalogue)

    # Stars to be subtracted
    subt_stars_catalogue = utils.replace_extension(rgunn_name, ".subt_stars")
    subt_stars_messg = "\n Finally, select the stars to be subtracted from the images"
    subt_stars_region = positions_of_stars_from_ds9(rgunn_name, subt_stars_messg, subt_stars_catalogue)

    # Now remove the stars from the images
    output_Ha = utils.add_suffix_prefix(Ha_name, suffix='-s')
    extract_stars.main(arguments=['--model_stars', mod_stars_catalogue,
                                  "--subt_stars", subt_stars_catalogue,
                                  '--coords', 'world',
                                  '--output', output_Ha,
                                  Ha_name])
    output_rgunn = utils.add_suffix_prefix(rgunn_name, suffix='-s')
    extract_stars.main(arguments=['--model_stars', mod_stars_catalogue,
                                  "--subt_stars", subt_stars_catalogue,
                                  '--coords', 'world',
                                  '--output', output_rgunn,
                                  rgunn_name])

    # Read images
    image_Ha = fits.open(output_Ha)
    Ha_data = np.array(image_Ha[0].data, dtype=np.float64)
    Ha_header = image_Ha[0].header
    Ha_shape = Ha_data.shape

    image_R = fits.open(output_rgunn)
    R_data = np.array(image_R[0].data, dtype=np.float64)
    R_header  = image_R[0].header
    R_shape = R_data.shape

    # Halpha image: get galaxy mask, sky mask, subtract sky from galaxy flux
    Ha_gal_mask = galaxy_region.as_imagecoord(Ha_header).get_mask(shape=Ha_shape)
    Ha_sky_mask = sky_region.as_imagecoord(Ha_header).get_mask(shape=Ha_shape)
    Ha_sky_flux = np.median(Ha_data[Ha_sky_mask == 1])
    Ha_data_nosky = Ha_data - Ha_sky_flux
    Ha_and_R_counts = np.sum( Ha_data_nosky[Ha_gal_mask == 1])
    print "Narrow Ha filter counts: ", Ha_and_R_counts

    # Same for R image: get galaxy mask, sky mask, subtract sky from galaxy flux
    R_gal_mask = galaxy_region.as_imagecoord(R_header).get_mask(shape=R_shape)
    R_sky_mask = sky_region.as_imagecoord(R_header).get_mask(shape=R_shape)
    R_sky_flux = np.median(R_data[R_sky_mask == 1])
    R_data_nosky = R_data - R_sky_flux
    R_counts = np.sum( R_data_nosky[R_gal_mask == 1])
    print "R filter counts: ", R_counts


    # Halpha is basically the subtractiong of the counts in Halpha filter minus the scaled R counts
    # but more precisely, the contrinbution of Halpha line to the filter rGunn should be corrected.
    # For:
    #   - a transmittance T(Gunn) of redshifted Halpha in the rGunn filter
    #   - T(Ha_filter) the transmittance of redshifted Halpha in the Halpha narrow filter
    #   - scaling_factor the factor of scale between counts in filter rGunn and in the Ha filter (typically from stars)
    #   - rgunn_scaled the number of counts in the Gunn filter, already scaled to match Halpha
    #   - Ha_filter is the counts in the narrow Ha filter, which obviously contain a contribution from R.
    #   - Halpha the number of counts of the Halpha line in the Halpha filter (the contribution of Halpha alone)
    #
    # we would need to subtract the Halpha contribution to rGunn before we scale it, but since we are starting
    # already from a scaled version of rGunn, we need to divide the contribution from Halpha by the scaling factor.
    # The system of equations would be:
    #
    #  R_counts = rgunn_scaled - Halpha * T(Gunn)/T(Ha_filter) / scaling_factor
    #  Halpha = Ha_filter - R_counts
    #
    # with solution:
    #  Halpha = (Ha_filter - rgunn_scaled) * scaling_factor/(scaling_factor - T(Gunn)/T(Halpha)
    #
    # which gives, for a typical scaling factor of ~11 and similar transmittances in both filters a correction
    # of the order of ~10%.

    Halpha = (Ha_and_R_counts - R_counts * scaling_factor) / (1 - scaling_factor * T_rGunn / T_Ha)
    print "Halpha counts: ", Halpha
    print Halpha  / T_Ha * 10**(-zp_Ha/2.5)



    mask = np.logical_or(Ha_gal_mask, Ha_sky_mask)
    d = ds9.ds9()
    d.set_np2arr(Ha_data * mask, dtype=np.float64)


#########################################
image_list = [("/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche1/cig0895_combined_H6678.fits",
               "/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche1/cig0895_combined_rGunn.fits"),
              ("/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche2/cig0828_combined_H6678.fits",
               "/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche2/cig0828_combined_rGunn.fits"),
              ("/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche2/cig0854_combined_H6678.fits",
               "/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche2/cig0854_combined_rGunn.fits"),
              ("/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche2/cig0633_combined_H6678.fits",
               "/mnt/data/OPTICAL_DATA/OPTICO/CAHA2.2/2003/AUG_03/Noche2/cig0633_combined_rGunn.fits")]

cigs = ["cig0895", "cig0828", "cig0854", "cig0633" ]

# From the filter transmittance curve of 667.8  for Halpha at redshift 0.0161 (for CIG0812)
T_Ha = [0.774, 0.774, 0.774, 0.744]
T_rGunn = [0.73, 0.73, 0.73, 0.73]
zp_Ha = [37.68, 37.873, 37.873, 37.873]
scaling_factor = [0.109849984506, 0.102706989603, 0.11220184543, 0.107102698809]

#########################################

def test():
    for ii, pair in enumerate(image_list):
        print "CIG ", cigs[ii]
        main(pair[0], pair[1], scaling_factor[ii], zp_Ha[ii], T_Ha[ii], T_rGunn[ii])
        print "\n" * 2

