'''
Scrit to download colored images from database
'''
# Import the necessary packages 
import splusdata 
import pandas as pd
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
from astropy.wcs import WCS
import os
import argparse
import sys
from pathlib import Path
ROOT_PATH = Path("cluster-hdbscan")

############################################################
# Definition to decompress the images ######################
############################################################
def fz2fits(image):
    """
    It converts SPLUS images
    from .fz to .fits
    """
    datos = fits.open(image)[1].data
    heada = fits.open(image)[1].header
    imageout = image[:-2] + 'fits'
    print ('Creating file: ')
    print (imageout)
    fits.writeto(imageout, datos, heada, overwrite=True)
############################################################

parser = argparse.ArgumentParser(
    description="""Get colored image and cut image in the r-band""")

parser.add_argument("source", type=str,
                    default="known-PN-jplus-idr",
                    help="Name of source, taken the prefix ")

parser.add_argument("--radi", type=float, default=None,
                    help="""Size of the images in pixel""")

parser.add_argument("--name", type=str, default=None,
                    help="""Name of the object""")

cmd_args = parser.parse_args()
file_ = cmd_args.source + ".csv"
try:
    data = pd.read_csv(ROOT_PATH / file_)[100:200]
except FileNotFoundError:
    data = pd.read_csv(file_)[100:200]

ra = data["RA"]
dec = data["DEC"]

# Radius
rad = int(cmd_args.radi)

# Nome of the object if has
Name = cmd_args.name

# Connect
conn = splusdata.connect('Luis', 'plutarco*80')

# Getting the colored imge
for key, value in data.iterrows():
    img = conn.twelve_band_img(value.RA, value.DEC, radius=rad, noise=0.15, saturation=0.15)
    # Getting the Fits image in the r-band
    hdu = conn.get_cut(value.RA, value.DEC, rad, 'R')

    # Save the image, note that the output image is compress
    hdu.writeto('{}_{}_r.fz'.format(value.ID.split("R3.")[-1], rad), overwrite=True) # write to fits


    # Decompress
    hdufits = fz2fits('{}_{}_r.fz'.format(value.ID.split("R3.")[-1], rad))

    # Read the FITS file
    hdul = fits.open('{}_{}_r.fits'.format(value.ID.split("R3.")[-1], rad))[0]
    wcs = WCS(hdul.header)

    print(wcs)                 

    f = plt.figure(figsize=(18,9))

    ax1 = aplpy.FITSFigure(hdul, figure=f, subplot=(1, 1, 1))#, north=True)
    plt.imshow(img, origin='lower', cmap='cividis', aspect='equal')
                 
    ax1.add_scalebar(20.0/3600)
    ax1.scalebar.set_label('20 arcsec')
    ax1.scalebar.set(color='yellow', linewidth=4, alpha=0.9)
    ax1.scalebar.set_font(size=35, weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')

    ax1.axis_labels.set_font(size=22, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    ax1.axis_labels.hide()
    ax1.tick_labels.hide()
    #ax1.axis_labels.hide_y()

    ax1.tick_labels.set_font(size=22, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    #ax1.list_layers()
    #ax1.show_markers(ra, dec, layer='marker', edgecolor='green', facecolor='none', marker='o', s=10, alpha=0.9, linewidths=60)#, layer='marker_set_1', edgecolor='black', facecolor='none', s=30, alpha=     0.5, linewidths=20)

    # ax1.axis_labels.hide_y()
    # ax1.tick_labels.hide_y()

    #ax2.colorbar.set_box([0.95, 0.1, 0.015, 0.8])
    ax1.set_theme('publication')
    #f.tight_layout()
    #f.savefig("-".join([image_name, "images.pdf"]))

    plt.savefig("Plots-hdbscan/image_"+str(value.ID.split("R3.")[-1]).replace(".", "-")+".pdf")
    f.clear()
    plt.close(f)
            
    os.remove('{}_{}_r.fz'.format(value.ID.split("R3.")[-1], rad))
    os.remove('{}_{}_r.fits'.format(value.ID.split("R3.")[-1], rad))
    print("FZ and FITS Files Removed!")
