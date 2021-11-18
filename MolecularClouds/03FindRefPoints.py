"""
This file is for finding reference points

    - Set the number of reference points and the file will find that many reference points in the specified region of
    interest and save their information as a text file. It will then produce a plot of the reference points scattered
    on the visual extinction file

This file was made for debugging / testing code - I do not think it is needed anymore
See instead 0XOptimalNumRefPoints
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from Classes.refPointFinder import RefPointFinder
from MolecularClouds.Classes.RegionOfInterest import Region

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = 'OrionB'
# cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- SET NUMBER OF REFERENCE POINTS --------
NRef = 2
# -------- SET NUMBER OF REFERENCE POINTS. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/MatchedRMExtinction'+cloudName+'.txt')
saveFilePath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/RefPoints/NRef'+str(NRef)+'.txt')
saveFigurePath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/Plots/ReferencePoints'+str(NRef)+'.png')
# -------- DEFINE FILES AND PATHS. --------

# -------- FIND REFERENCE POINTS --------
RefPointFinder(MatchedRMExtincPath, NRef, saveFilePath)
# -------- FIND REFERENCE POINTS. --------

# -------- READ AND UNPACK REFERENCE POINT DATA --------
referenceData = pd.read_csv(saveFilePath)
Ra = list(referenceData[referenceData.columns[2]])
Dec = list(referenceData[referenceData.columns[3]])
RM = list(referenceData[referenceData.columns[4]])
# -------- READ AND UNPACK REFERENCE POINT DATA. --------

# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- PREPARE TO PLOT reference POINTS --------
# ---- Convert Ra and Dec of reference points into pixel values of the fits file
x = []  # x pixel coordinate of reference
y = []  # y pixel coordinate of reference
for i in range(len(Ra)):
    pixelRow, pixelColumn = wcs.wcs_world2pix(Ra[i], Dec[i], 0)
    x.append(pixelRow)
    y.append(pixelColumn)
# ---- Convert Ra and Dec of reference points into pixel values of the fits file.

# -------- CREATE A FIGURE --------
fig = plt.figure(num=None, figsize=(12, 12), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection=wcs)

plt.title('Reference Rotation Measure Points \n NumPoints = {}'.format(NRef), fontsize=12, y=1.08)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
plt.scatter(x, y, marker='o', facecolor='green', linewidth=.5, edgecolors='black')

# ---- Style the main axes and their grid
# ax.set_xlim(xmin_pix, xmax_pix)
# ax.set_ylim(ymin_pix, ymax_pix)

ra = ax.coords[0]
dec = ax.coords[1]
ra.set_major_formatter('d')
dec.set_major_formatter('d')
ra.set_axislabel('Ra (degree)')
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
cb = plt.colorbar(im, orientation='vertical', ticklocation='right', fraction=0.02, pad=0.145)
cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
# ---- Style the colour bar.

# ---- Style the legend
# marker1 = plt.scatter([], [], s=10, facecolor=(1, 1, 1, 0.7), edgecolor='black')
# marker2 = plt.scatter([], [], s=50, facecolor=(1, 1, 1, 0.7), edgecolor='black')
# marker3 = plt.scatter([], [], s=100, facecolor=(1, 1, 1, 0.7), edgecolor='black')
# marker4 = plt.scatter([], [], s=200, facecolor=(1, 1, 1, 0.7), edgecolor='black')
# marker5 = plt.scatter([], [], s=100, facecolor=(1, 0, 0, 0.7), edgecolor='black')
# marker6 = plt.scatter([], [], s=100, facecolor=(0, 0, 1, 0.7), edgecolor='black')
# legend_markers = [marker1, marker2, marker3, marker4, marker5, marker6]
#
# labels = [
#     str(10) + ' rad m' + r'$^{-2}$',
#     str(50) + ' rad m' + r'$^{-2}$',
#     str(100) + ' rad m' + r'$^{-2}$',
#     str(200) + ' rad m' + r'$^{-2}$',
#     'Negative RM',
#     'Positive RM', ]
#
# plt.legend(handles=legend_markers, labels=labels, scatterpoints=1)
# ---- Style the legend.

# ---- Display or save the figure
# plt.show()
plt.savefig(saveFigurePath)
# ---- Display or save the figure.
# -------- CREATE A FIGURE. --------
