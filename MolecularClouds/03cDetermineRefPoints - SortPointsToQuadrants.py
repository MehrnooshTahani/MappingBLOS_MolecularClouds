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
import LocalLibraries.PlotTemplates as pt

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
saveFilePath_ALlPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_allPotRefPoints + config.cloud + '.txt')
saveFigurePath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots, config.cloud + "QuadrantDivision.png")
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
Q1, Q2, Q3, Q4 = rjl.sortQuadrants(list(AllPotentialRefPoints.index), AllPotentialRefPoints['Extinction_Index_x'], AllPotentialRefPoints['Extinction_Index_y'], m, b, mPerp, bPerp)
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
# ---- Convert Ra and Dec of reference points into pixel values of the fits file
x_AllRef, y_AllRef = rjl.RADec2xy(Ra_AllRef, Dec_AllRef, wcs)
# ---- Convert Ra and Dec of reference points into pixel values of the fits file.
# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS. --------

# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP --------
fig, ax = pt.extinctionPlot(hdu, regionOfInterest)

plt.title('Quadrant division of the ' + cloudName + ' region\n', fontsize=12, y=1.08)

plt.plot(x, y)
plt.plot(x, y2)
ax.set_xlim(0, hdu.data.shape[1])
ax.set_ylim(0, hdu.data.shape[0])
if not math.isnan(regionOfInterest.xmax) and not math.isnan(regionOfInterest.xmin):
    ax.set_xlim(regionOfInterest.xmin, regionOfInterest.xmax)
if not math.isnan(regionOfInterest.ymax) and not math.isnan(regionOfInterest.ymin):
    ax.set_ylim(regionOfInterest.ymin, regionOfInterest.ymax)

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
plt.savefig(saveFigurePath)
plt.show()
plt.close()
# ---- Display or save the figure.
# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP. --------