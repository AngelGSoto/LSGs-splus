'''
Make photo-spectra from observed SPLUS objects. This program is an updated version of the program: photo-spectra-SPLUSDR2.py.
I madified this one to work with SPLUS SMC catalog
'''
from __future__ import print_function
import numpy as np
import glob
import json
import matplotlib.pyplot as plt
from astropy.table import Table
import seaborn as sns
import sys
import argparse
import pandas as pd
import os
from colour import Color
from pathlib import Path
ROOT_PATH = Path("cluster-hdbscan")

Number = []

wl = [3485, 3785, 3950, 4100, 4300, 4803, 5150, 6250, 6600, 7660, 8610, 9110]
color = ["#CC00FF", "#9900FF", "#6600FF", "#0000FF", "#009999", "#006600", "#DD8000", "#FF0000", "#CC0066", "#990033", "#660033", "#330034"]
marker = ["s", "o", "o", "o", "o", "s", "o", "s", "o", "s", "o", "s"] ### tienen todos los filtros


parser = argparse.ArgumentParser(
    description="""Write wave and magnitude of a spectrum""")

parser.add_argument("source", type=str,
                    default="known-PN-jplus-idr",
                    help="Name of source, taken the prefix ")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
# file_ = args.source + ".ecsv"


# try:
#     data = Table.read(ROOT_PATH / file_, format="ascii.ecsv")
# except FileNotFoundError:
# file_ = args.source + ".dat"
# data = Table.read(ROOT_PATH / file_, format="ascii"):

file_ = args.source + ".csv"
try:
    data = pd.read_csv(ROOT_PATH / file_)[:100]
except FileNotFoundError:
    data = pd.read_csv(file_)[:100]


# Number of objects
n = len(data)


# Magnitudes list
Number = []
mag_auto  = [[] for _ in range(n)]
mag_petro = [[] for _ in range(n)]
mag_aper = [[] for _ in range(n)]

# Error lists
mag_auto_err  = [[] for _ in range(n)]
mag_petro_err  = [[] for _ in range(n)]
mag_aper_err  = [[] for _ in range(n)]

print("The numbre of objects for plotting is: %d" % n)

for i in range(n):
    mag_aper[i].append(data["u_auto"][i]) 
    mag_aper[i].append(data["J0378_auto"][i])
    mag_aper[i].append(data["J0395_auto"][i])
    mag_aper[i].append(data["J0410_auto"][i])
    mag_aper[i].append(data["J0430_auto"][i])
    mag_aper[i].append(data["g_auto"][i])
    mag_aper[i].append(data["J0515_auto"][i]) 
    mag_aper[i].append(data["r_auto"][i]) 
    mag_aper[i].append(data["J0660_auto"][i])
    mag_aper[i].append(data["i_auto"][i]) 
    mag_aper[i].append(data["J0861_auto"][i]) 
    mag_aper[i].append(data["z_auto"][i])

    #ERRO Aper
    mag_aper_err[i].append(float(data["e_u_auto"][i]))
    mag_aper_err[i].append(float(data["e_J0378_auto"][i]))
    mag_aper_err[i].append(float(data["e_J0395_auto"][i]))
    mag_aper_err[i].append(float(data["e_J0410_auto"][i]))
    mag_aper_err[i].append(float(data["e_J0430_auto"][i]))
    mag_aper_err[i].append(float(data["e_g_auto"][i]))
    mag_aper_err[i].append(float(data["e_J0515_auto"][i])) 
    mag_aper_err[i].append(float(data["e_r_auto"][i])) 
    mag_aper_err[i].append(float(data["e_J0660_auto"][i])) 
    mag_aper_err[i].append(float(data["e_i_auto"][i]))
    mag_aper_err[i].append(float(data["e_J0861_auto"][i]))
    mag_aper_err[i].append(float(data["e_z_auto"][i]))
   
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
    ##########################################################################################
    # Plotting -- auto  ######################################################################
    ##########################################################################################
    if file_.endswith("csv"):
        plotfile = "Plots-lsg-confirmed/photopectrum_splus_"+str(data["ID"][i].split("R3.")[-1]).replace(".", "-")+"_{}_auto.pdf".format(file_.split('.csv')[0])
    else:
        plotfile = "Plots-lsg-confirmed/photopectrum_splus_"+str(data["ID"][i].split("R3.")[-1]).replace(".", "-")+"_{}_auto.pdf".format(file_.split('.dat')[0])
    fig = plt.figure(figsize=(15.5, 9.5))
    ax = fig.add_subplot(1,1,1)
    plt.tick_params(axis='x', labelsize=42) 
    plt.tick_params(axis='y', labelsize=42)
    ax.set_xlim(left=3000, right=9700)
    #ax.set_ylim(ymin=17.5,ymax=23)
    #ax1.set_xlabel(r'$\lambda$')
    ax.set_xlabel(r'Wavelength $[\mathrm{\AA]}$', fontsize = 44)
    ax.set_ylabel(r'Magnitude [AB]', fontsize = 44)

    mask = [mag_aper[i][m] != 99.0 for m in range(len(mag_aper[0]))]
    mask_err = [mag_aper_err[i][m] <= 0.9 for m in range(len(mag_aper_err[0]))]
    mask_t = np.array(mask) & np.array(mask_err)
    ax.plot(np.array(wl)[mask_err], np.array(mag_aper[i])[mask_err], '-k', alpha=0.2)#, label='Auto')
    for wl1, mag, mag_err, colors, marker_ in zip(np.array(wl), np.array(mag_aper[i]), np.array(mag_aper_err[i]), color, marker):
        if mag_err <= 0.9:
            ax.scatter(wl1, mag, color = colors, marker=marker_, s=600, zorder=10)
            ax.errorbar(wl1, mag, yerr=mag_err, marker='.', fmt='.', color=colors, ecolor=colors, elinewidth=5.9, markeredgewidth=5.2,  capsize=20)
    # plt.text(0.06, 0.1, "Fr 2-21",
    #          transform=ax.transAxes, fontsize=48,  fontdict=font)
    #plt.subplots_adjust(bottom=0.19)
    plt.legend(fontsize=20.0)
    plt.tight_layout()
    plt.gca().invert_yaxis()
    if args.debug:
        print("Finished plot S-spectra of:", data["ID"][i].split("R3.")[-1].replace(".", "-"))
    #save_path = '../../../Dropbox/JPAS/paper-phot/'
    #file_save = os.path.join(save_path, plotfile)
    plt.savefig(plotfile)
    #plt.clf()
    fig.clear()
    plt.close(fig)
