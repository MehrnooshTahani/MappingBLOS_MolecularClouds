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
if abs(regionOfInterest.cloudLatitude) < config.offDiskLatitude:
    Av_threshold = config.onDiskAvThresh
else:
    Av_threshold = config.offDiskAvThresh
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
cloudJeansLength = regionOfInterest.jeanslength  # [pc]
minDiff = np.degrees(np.arctan(cloudJeansLength / cloudDistance))  # [deg]

minDiff_pix = minDiff / abs(hdu.header['CDELT1'])
NDeltNear = config.nearExtinctionMultiplier * math.ceil(minDiff_pix)  # Round up
NDeltFar = config.farExtinctionMultiplier * math.ceil(minDiff_pix)  # Round up
print("\t-A close region around the point has been defined to the suggested {} pixels".format(NDeltNear))
print("\t-A far region around the point has been defined to the suggested {} pixels".format(NDeltFar))

# Choose the minimum extinction value which you want to correspond to an "on" position
highExtinctionThreshold = config.highExtinctionThreshMultiplier * Av_threshold
print("\t-A region of high extinction has been defined to the suggested suggested Av={}".format(highExtinctionThreshold))

# -------- Define the range.

# -------- For each potential reference point
nearHighExtinctionRegion = []
farHighExtinctionRegion = []
for i in list(AllPotentialRefPoints.index):
    idNum = AllPotentialRefPoints['ID#'][i]
    px = AllPotentialRefPoints['Extinction_Index_x'][i]
    py = AllPotentialRefPoints['Extinction_Index_y'][i]

    # ---- Find the extinction range for the given point
    if rjl.nearHighExtinction(px, py, hdu.data, NDeltNear, highExtinctionThreshold):
        nearHighExtinctionRegion.append(i + 1)
    if not rjl.nearHighExtinction(px, py, hdu.data, NDeltFar, highExtinctionThreshold):
        farHighExtinctionRegion.append(i + 1)
    # ---- Find the extinction range for the given point.
# -------- For each potential reference point.
print('The potential reference point(s) {} are near a region of high extinction'.format(nearHighExtinctionRegion))
print('The potential reference point(s) {} are far from a region of high extinction'.format(farHighExtinctionRegion))
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION. --------
