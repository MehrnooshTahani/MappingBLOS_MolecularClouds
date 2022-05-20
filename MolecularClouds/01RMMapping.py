"""
This file maps the RMs on the extinction files as a first step to get an understanding of the rotation measure coverage
of the region of interest.
"""
from LocalLibraries.DataFile import DataFile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.wcs import WCS
from astropy.io import fits
import os
from LocalLibraries.RegionOfInterest import Region
import LocalLibraries.config as config
import math
import LocalLibraries.PlotTemplates as pt

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
RMCatalogPath = os.path.join(config.dir_root, config.dir_data, config.file_RMCatalogue)
saveFigurePath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots, config.prefix_rmMapping + config.cloud + '.png')
# -------- DEFINE FILES AND PATHS. --------


# -------- FUNCTION DEFINITION --------
def rm2RGB(rm):
    """
    Takes rotation measure values and assigns them a marker colour and size for use in plotting rotation measure data

    :param rm: The rotation measure, or list of rotation measures
    :return:  A tuple of (colour, size) corresponding to the rotation measure. Note "colour" is a tuple of (RBG,alpha)
    """
    c = []  # Marker colour
    s = []  # Marker size

    for item in rm:
        s.append(abs(item))

        alpha = 1  # Optional: set the transparency
        if int(np.sign(item)) == -1:
            c.append((1, 0, 0, alpha))  # Negative rotation measures assigned red
        elif int(np.sign(item)) == 1:
            c.append((0, 0, 1, alpha))  # Positive rotation measures assigned blue
        elif np.sign(item) == 0:
            c.append((0, 1, 0, alpha))  # Zero-value rotation measures assigned green

    # return the list of RGBA tuples and sizes
    return c, s
# -------- FUNCTION DEFINITION. --------


# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- READ ROTATION MEASURE FILE --------
# Get all the rm points within the region of interest
rmData = DataFile(RMCatalogPath, regionOfInterest.raHoursMax, regionOfInterest.raMinsMax, regionOfInterest.raSecMax,
                  regionOfInterest.raHoursMin, regionOfInterest.raMinsMin, regionOfInterest.raSecMin,
                  regionOfInterest.decDegMax, regionOfInterest.decDegMin)
# -------- READ ROTATION MEASURE FILE. --------

# -------- PREPARE TO PLOT ROTATION MEASURES --------
# ---- Convert Ra and Dec of RMs into pixel values of the fits file
x = []  # x pixel coordinate of RM
y = []  # y pixel coordinate of RM
for i in range(len(rmData.targetRaHourMinSecToDeg)):
    pixelRow, pixelColumn = wcs.wcs_world2pix(rmData.targetRaHourMinSecToDeg[i], rmData.targetDecDegArcMinSecs[i], 0)
    x.append(pixelRow)
    y.append(pixelColumn)
# ---- Convert Ra and Dec of RMs into pixel values of the fits file.

color, size = rm2RGB(rmData.targetRotationMeasures)
# -------- PREPARE TO PLOT ROTATION MEASURES. --------

# -------- CREATE A FIGURE --------
fig, ax = pt.extinctionPlot(hdu, regionOfInterest)

plt.title('Rotation Measure Data' + ' in the '+cloudName+' region\n', fontsize=12, y=1.08)
plt.scatter(x, y, marker='o', s=size, facecolor=color, linewidth=.5, edgecolors='black')

# ---- Style the legend
marker1 = plt.scatter([], [], s=10, facecolor=(1, 1, 1, 0.7), edgecolor='black')
marker2 = plt.scatter([], [], s=50, facecolor=(1, 1, 1, 0.7), edgecolor='black')
marker3 = plt.scatter([], [], s=100, facecolor=(1, 1, 1, 0.7), edgecolor='black')
marker4 = plt.scatter([], [], s=200, facecolor=(1, 1, 1, 0.7), edgecolor='black')
marker5 = plt.scatter([], [], s=100, facecolor=(1, 0, 0, 0.7), edgecolor='black')
marker6 = plt.scatter([], [], s=100, facecolor=(0, 0, 1, 0.7), edgecolor='black')
legend_markers = [marker1, marker2, marker3, marker4, marker5, marker6]

labels = [
    str(10) + ' rad m' + r'$^{-2}$',
    str(50) + ' rad m' + r'$^{-2}$',
    str(100) + ' rad m' + r'$^{-2}$',
    str(200) + ' rad m' + r'$^{-2}$',
    'Negative RM',
    'Positive RM', ]

legend = plt.legend(handles=legend_markers, labels=labels, scatterpoints=1)

frame = legend.get_frame()
frame.set_facecolor('1')
frame.set_alpha(0.4)
# ---- Style the legend.

plt.savefig(saveFigurePath, bbox_inches='tight')
print('Saving figure to '+saveFigurePath)
plt.show()
# -------- CREATE A FIGURE. --------
