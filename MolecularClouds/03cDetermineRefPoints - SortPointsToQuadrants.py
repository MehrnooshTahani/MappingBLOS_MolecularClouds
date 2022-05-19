"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os
import pandas as pd
from astropy.wcs import WCS
from astropy.io import fits
import math
import matplotlib.pyplot as plt
from LocalLibraries.RegionOfInterest import Region
import adjustText
import LocalLibraries.config as config
import LocalLibraries.RefJudgeLib as rjl

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
saveFilePath_ALlPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_allPotRefPoints + config.cloud + '.txt')
# -------- DEFINE FILES AND PATHS. --------

# -------- LOAD ALL POTENTIAL REFERENCE POINTS --------
AllPotentialRefPoints = pd.read_csv(saveFilePath_ALlPotentialRefPoints)
print('---------------------\n')
# -------- LOAD ALL POTENTIAL REFERENCE POINTS. --------

# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- FIND REGIONS TO SPLIT THE CLOUD INTO. --------
cloudCenterX, cloudCenterY = rjl.findWeightedCenter(hdu.data, regionOfInterest.xmin, regionOfInterest.xmax, regionOfInterest.ymin, regionOfInterest.ymax)
m, b = rjl.getDividingLine(hdu.data, regionOfInterest.xmin, regionOfInterest.xmax, regionOfInterest.ymin, regionOfInterest.ymax)
mPerp, bPerp = rjl.getPerpendicularLine(cloudCenterX, cloudCenterY, m)

# ---- Get points to graph these lines
x = range(0, hdu.data.shape[0])
if not math.isnan(regionOfInterest.xmax) and not math.isnan(regionOfInterest.xmin):
    x = range(int(regionOfInterest.xmin), int(regionOfInterest.xmax))
y = m * x + b
y2 = mPerp * x + bPerp
# ---- Get points to graph these lines
# -------- FIND REGIONS TO SPLIT THE CLOUD INTO. --------

# -------- SORT REF POINTS INTO THESE REGIONS. --------
Q1 = []
Q2 = []
Q3 = []
Q4 = []
for i in range(len(AllPotentialRefPoints)):
    idNum = AllPotentialRefPoints['ID#'][i]
    px = AllPotentialRefPoints['Extinction_Index_x'][i]
    py = AllPotentialRefPoints['Extinction_Index_y'][i]

    # ---- Sort into quadrant
    aboveCloudLine = rjl.isPointAboveLine(px, py, m, b)
    aboveCloudPerpLine = rjl.isPointAboveLine(px, py, mPerp, bPerp)

    if aboveCloudLine and aboveCloudPerpLine:
        Q1.append(i+1)
    elif aboveCloudLine and not aboveCloudPerpLine:
        Q2.append(i+1)
    elif not aboveCloudLine and aboveCloudPerpLine:
        Q3.append(i+1)
    elif not aboveCloudLine and not aboveCloudPerpLine:
        Q4.append(i+1)
    # ---- Sort into quadrant
# -------- SORT REF POINTS INTO THESE REGIONS. --------

# -------- OUTPUT RESULTS. --------
print("The potential reference points, sorted by quadrant, are:")
print("Q1: {}".format(Q1))
print("Q2: {}".format(Q2))
print("Q3: {}".format(Q3))
print("Q4: {}".format(Q4))
# -------- OUTPUT RESULTS. --------

# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS --------
n_AllRef = list(AllPotentialRefPoints['ID#'])
Ra_AllRef = list(AllPotentialRefPoints['Ra(deg)'])
Dec_AllRef = list(AllPotentialRefPoints['Dec(deg)'])
RM_AllRef = list(AllPotentialRefPoints['Rotation_Measure(rad/m2)'])
Av_AllRef = list(AllPotentialRefPoints['Extinction_Value'])
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

plt.title('Quadrant division of the ' + cloudName + ' region\n', fontsize=12, y=1.08)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')

plt.plot(x, y)
plt.plot(x, y2)

ax.set_xlim(0, hdu.data.shape[1])
ax.set_ylim(0, hdu.data.shape[0])

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
plt.show()
plt.close()
# ---- Display or save the figure.
# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP. --------