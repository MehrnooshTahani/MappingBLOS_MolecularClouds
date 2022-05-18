"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os
import pandas as pd
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import math
import matplotlib.pyplot as plt
from Classes.RegionOfInterest import Region
from Classes.FindAllPotentialRefPoints import FindAllPotentialReferencePoints
from Classes.FindOptimalRefPoints import FindOptimalRefPoints
import adjustText
from Classes.CalculateB import CalculateB
import Classes.config as config

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

# -------- CHOOSE THE THRESHOLD EXTINCTION --------
print('\n---------------------')

print('All potential reference points will be taken to be all points with a visual extinction value less than the '
      'extinction threshold.')
if abs(regionOfInterest.cloudLatitude) < config.offDiskLatitude:
    Av_threshold = config.onDiskAvThresh
    print('\t-For clouds that appear near the disk, such as {}, an appropriate threshold value is {}.'
          .format(cloudName, Av_threshold))
else:
    Av_threshold = config.offDiskAvThresh
    print('\t-For clouds that appear off the disk, such as {}, an appropriate threshold value is {}.'
          .format(cloudName, Av_threshold))

print("Given this information, the threshold extinction has been set to the suggested {}".format(Av_threshold))
# -------- CHOOSE THE THRESHOLD EXTINCTION. --------

# -------- FIND ALL POTENTIAL REFERENCE POINTS --------
AllPotenitalRefPoints = FindAllPotentialReferencePoints(cloudName, Av_threshold, saveFilePath=saveFilePath_ALlPotentialRefPoints)
print('Based on this threshold extinction, a total of {} potential reference points were found.'.format(AllPotenitalRefPoints.numAllRefPoints))

print(AllPotenitalRefPoints.AllRefPoints)
print('---------------------\n')
# -------- FIND ALL POTENTIAL REFERENCE POINTS. --------

# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS --------
n_AllRef = list(AllPotenitalRefPoints.AllRefPoints['ID#'])
Ra_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Ra(deg)'])
Dec_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Dec(deg)'])
RM_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Rotation_Measure(rad/m2)'])
Av_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Extinction_Value'])
# ---- Convert Ra and Dec of reference points into pixel values of the fits file
x_AllRef = []  # x pixel coordinate of reference
y_AllRef = []  # y pixel coordinate of reference
for i in range(len(Ra_AllRef)):
    pixelRow, pixelColumn = wcs.wcs_world2pix(Ra_AllRef[i], Dec_AllRef[i], 0)
    x_AllRef.append(pixelRow)
    y_AllRef.append(pixelColumn)
# ---- Convert Ra and Dec of reference points into pixel values of the fits file.
# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS. --------

# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP --------
fig = plt.figure(figsize=(8, 8), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection=wcs)

plt.title('All Potential Reference Points' + ' in the ' + cloudName + ' region\n', fontsize=12, y=1.08)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
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

# ---- Style the main axes and their grid
if not math.isnan(regionOfInterest.xmax) and not math.isnan(regionOfInterest.xmin):
    ax.set_xlim(regionOfInterest.xmin, regionOfInterest.xmax)
if not math.isnan(regionOfInterest.ymax) and not math.isnan(regionOfInterest.ymin):
    ax.set_ylim(regionOfInterest.ymin, regionOfInterest.ymax)

ra = ax.coords[0]
dec = ax.coords[1]
ra.set_major_formatter('d')
dec.set_major_formatter('d')
ra.set_axislabel('RA (degree)')
dec.set_axislabel('Dec (degree)')

dec.set_ticks(number=10)
ra.set_ticks(number=20)
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_minor_frequency(10)

ra.grid(color='black', alpha=0.5, linestyle='solid')
dec.grid(color='black', alpha=0.5, linestyle='solid')
# ---- Style the main axes and their grid.

# ---- Style the overlay and its grid
overlay = ax.get_coords_overlay('galactic')

overlay[0].set_axislabel('Longitude')
overlay[1].set_axislabel('Latitude')

overlay[0].set_ticks(color='grey', number=20)
overlay[1].set_ticks(color='grey', number=20)

overlay.grid(color='grey', linestyle='solid', alpha=0.7)
# ---- Style the overlay and its grid.

# ---- Style the colour bar
if regionOfInterest.fitsDataType == 'HydrogenColumnDensity':
    cb = plt.colorbar(im, ticklocation='right', fraction=0.02, pad=0.145, format='%.0e')
    cb.ax.set_title('Hydrogen Column Density', linespacing=0.5, fontsize=12)
elif regionOfInterest.fitsDataType == 'VisualExtinction':
    cb = plt.colorbar(im, ticklocation='right', fraction=0.02, pad=0.145)
    cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
# ---- Style the colour bar.

# ---- Display or save the figure
saveFigurePath_RefPointMap = saveFigureDir_RefPointMap + os.sep + 'RefPointMap_AllPotentialRefPoints.png'
plt.savefig(saveFigurePath_RefPointMap)
plt.show()
plt.close()
# ---- Display or save the figure.
print('Saving the map of all potential reference points to '+saveFigurePath_RefPointMap)
# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP. --------