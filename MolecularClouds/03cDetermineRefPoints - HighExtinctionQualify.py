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
print('---------------------\n')
# -------- FIND ALL POTENTIAL REFERENCE POINTS. --------

# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION --------
print('---------------------')
print('We will now check if any of the potential reference points are near a region of high extinction.')
# -------- Define the range
# The distance the point can be from a region of high extinction and still be thought to sample the background
cloudDistance = regionOfInterest.distance  # [pc]
cloudJeansLength = config.cloudJeansLength  # [pc]
minDiff = np.degrees(np.arctan(cloudJeansLength / cloudDistance))  # [deg]

minDiff_pix = minDiff / abs(hdu.header['CDELT1'])
NDelt = config.pixelCheckMultiplier * 5 * math.ceil(minDiff_pix)  # Round up
print("\t-A region around the point has been defined to the suggested {} pixels".format(NDelt))

# Choose the minimum extinction value which you want to correspond to an "on" position
highExtinctionThreshold = config.highExtinctionThreshMultiplier * Av_threshold
print("\t-A region of high extinction has been defined to the suggested suggested Av={} (y/n)".format(highExtinctionThreshold))

# -------- Define the range.

# -------- For each potential reference point
nearHighExtinctionRegion = []
for i in range(len(AllPotenitalRefPoints.AllRefPoints)):
    idNum = AllPotenitalRefPoints.AllRefPoints['ID#'][i]
    px = AllPotenitalRefPoints.AllRefPoints['Extinction_Index_x'][i]
    py = AllPotenitalRefPoints.AllRefPoints['Extinction_Index_y'][i]

    # ---- Find the extinction range for the given point
    ind_xmax = px + NDelt + 1  # add 1 to be inclusive of the upper bound
    ind_ymax = py + NDelt + 1  # add 1 to be inclusive of the upper bound
    ind_xmin = px - NDelt
    ind_ymin = py - NDelt
    # ---- Find the extinction range for the given point.

    # ---- Cycle through extinction values within the range
    # If an extinction value within this range is too high, then it cannot be considered as a reference point
    highExtinction = False
    for pxx in range(ind_xmin, ind_xmax):
        for pyy in range(ind_ymin, ind_ymax):
            if 0 <= pxx < hdu.data.shape[1] and 0 <= pyy < hdu.data.shape[0]:
                extinction_val = hdu.data[pyy, pxx]
                if extinction_val > highExtinctionThreshold:
                    highExtinction = True
    if highExtinction == False:
        nearHighExtinctionRegion.append(i + 1)  # To identify points numbered in order of increasing extinction
    # ---- Cycle through extinction values within the range.
# -------- For each potential reference point.
print('The potential reference point(s) {} are far from a region of high extinction'.format(nearHighExtinctionRegion))
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION. --------
