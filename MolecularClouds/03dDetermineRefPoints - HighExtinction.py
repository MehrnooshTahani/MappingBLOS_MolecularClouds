"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os

import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

import math
import pandas as pd

from LocalLibraries.RegionOfInterest import Region
import LocalLibraries.config as config
import LocalLibraries.RefJudgeLib as rjl

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

# -------- LOAD ALL POTENTIAL REFERENCE POINTS --------
AllPotentialRefPoints = pd.read_csv(saveFilePath_ALlPotentialRefPoints)
print('---------------------\n')
# -------- LOAD ALL POTENTIAL REFERENCE POINTS. --------

# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION --------
print('---------------------')
print('We will now check if any of the potential reference points are near a region of high extinction.')
# -------- Define the range
# The distance the point can be from a region of high extinction and still be thought to sample the background
cloudDistance = regionOfInterest.distance  # [pc]
cloudJeansLength = config.cloudJeansLength  # [pc]
minDiff = np.degrees(np.arctan(cloudJeansLength / cloudDistance))  # [deg]

minDiff_pix = minDiff / abs(hdu.header['CDELT1'])
NDelt = config.pixelCheckMultiplier * math.ceil(minDiff_pix)  # Round up
print("\t-A region around the point has been defined to the suggested {} pixels".format(NDelt))

# Choose the minimum extinction value which you want to correspond to an "on" position
highExtinctionThreshold = config.highExtinctionThreshMultiplier * Av_threshold
print("\t-A region of high extinction has been defined to the suggested suggested Av={}".format(highExtinctionThreshold))

# -------- Define the range.

# -------- For each potential reference point
nearHighExtinctionRegion = []
farHighExtinctionRegion = []
for i in range(len(AllPotentialRefPoints)):
    idNum = AllPotentialRefPoints['ID#'][i]
    px = AllPotentialRefPoints['Extinction_Index_x'][i]
    py = AllPotentialRefPoints['Extinction_Index_y'][i]

    # ---- Find the extinction range for the given point
    if rjl.nearHighExtinction(px, py, hdu.data, NDelt, highExtinctionThreshold):
        nearHighExtinctionRegion.append(i + 1)
    if not rjl.nearHighExtinction(px, py, hdu.data, NDelt * 2, highExtinctionThreshold):
        farHighExtinctionRegion.append(i + 1)
    # ---- Find the extinction range for the given point.
# -------- For each potential reference point.
print('The potential reference point(s) {} are near a region of high extinction'.format(nearHighExtinctionRegion))
print('The potential reference point(s) {} are far from a region of high extinction'.format(farHighExtinctionRegion))
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION. --------
