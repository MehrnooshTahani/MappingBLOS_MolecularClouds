"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os

from astropy.wcs import WCS
from astropy.io import fits

import LocalLibraries.config as config
from LocalLibraries.RegionOfInterest import Region

import pandas as pd
import numpy as np

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
saveFilePath_ALlPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_allPotRefPoints + config.cloud + '.txt')
MatchedRMExtincPath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_RMExtinctionMatch + config.cloud + '.txt')
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
# ---- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA
matchedRMExtinctionData = pd.read_csv(MatchedRMExtincPath, sep='\t')
# ---- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA

# ---- FIND ALL POTENTIAL REFERENCE POINTS
# -------- Criterion: Av < threshold
'''
We will only consider points with visual extinction less than the specified threshold value as potential 
reference points
- Here we extract these points and sort the resulting dataframe from smallest to greatest extinction 
'''
# Indices where the extinction is less than the threshold:
ind_extinction = np.where(matchedRMExtinctionData['Extinction_Value'] <= Av_threshold)[0]

# All potential reference points are all reference points with extinction less than the threshold:
AllPotentialRefPoints = matchedRMExtinctionData.loc[ind_extinction].sort_values('Extinction_Value', ignore_index=True)
numAllRefPoints = len(AllPotentialRefPoints)
# -------- Criterion: Av < threshold.
# ---- FIND ALL POTENTIAL REFERENCE POINTS

# ---- SAVE REFERENCE POINT DATA AS A TABLE
if saveFilePath_ALlPotentialRefPoints is not None:
    AllPotentialRefPoints.to_csv(saveFilePath_ALlPotentialRefPoints, index=False)
# ---- SAVE REFERENCE POINT DATA AS A TABLE.
print('Based on this threshold extinction, a total of {} potential reference points were found.'.format(numAllRefPoints))
print(AllPotentialRefPoints)
print('---------------------\n')
# -------- FIND ALL POTENTIAL REFERENCE POINTS. --------