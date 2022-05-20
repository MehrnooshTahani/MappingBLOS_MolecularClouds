"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os

from astropy.wcs import WCS
from astropy.io import fits

import math
import pandas as pd

import matplotlib.pyplot as plt
import adjustText

import LocalLibraries.config as config
from LocalLibraries.RegionOfInterest import Region
import LocalLibraries.RefJudgeLib as rjl
import LocalLibraries.PlotTemplates as pt

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
saveFilePath_ALlPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_allPotRefPoints + config.cloud + '.txt')
saveFigureDir_RefPointMap = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots)
# -------- DEFINE FILES AND PATHS. --------

# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- LOAD ALL POTENTIAL REFERENCE POINTS --------
AllPotentialRefPoints = pd.read_csv(saveFilePath_ALlPotentialRefPoints)
print('---------------------\n')
# -------- LOAD ALL POTENTIAL REFERENCE POINTS. --------

# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS --------
n_AllRef = list(AllPotentialRefPoints['ID#'])
Ra_AllRef = list(AllPotentialRefPoints['Ra(deg)'])
Dec_AllRef = list(AllPotentialRefPoints['Dec(deg)'])
# ---- Convert Ra and Dec of reference points into pixel values of the fits file
x_AllRef, y_AllRef = rjl.RADec2xy(Ra_AllRef, Dec_AllRef, wcs)
# ---- Convert Ra and Dec of reference points into pixel values of the fits file.
# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS. --------

# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP --------
fig, ax = pt.extinctionPlot(hdu, regionOfInterest)

plt.title('All Potential Reference Points' + ' in the ' + cloudName + ' region\n', fontsize=12, y=1.08)
plt.scatter(x_AllRef, y_AllRef, marker='o', facecolor='green', linewidth=.5, edgecolors='black', s=50)

# ---- Annotate the chosen reference points
text = []
for i, number in enumerate(n_AllRef):
    # Each point is labelled in order of increasing extinction value
    # To label with ID number use: txt = ax.text(x_AllRef[i], y_AllRef[i], str(number), size=9, color='w')
    txt = ax.text(x_AllRef[i], y_AllRef[i], str(i + 1), size=9, color='w')
    text.append(txt)
adjustText.adjust_text(text)
# ---- Annotate the chosen reference points

# ---- Display or save the figure
saveFigurePath_RefPointMap = saveFigureDir_RefPointMap + os.sep + 'RefPointMap_AllPotentialRefPoints.png'
plt.savefig(saveFigurePath_RefPointMap)
plt.show()
plt.close()
# ---- Display or save the figure.
print('Saving the map of all potential reference points to '+saveFigurePath_RefPointMap)
# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP. --------