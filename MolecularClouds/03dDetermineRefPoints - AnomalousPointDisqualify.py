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
# -------- Load matched rm and extinction data
MatchedRMExtincPath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_RMExtinctionMatch + config.cloud + '.txt')
# -------- Load matched rm and extinction data.

# -------- DEFINE FILES AND PATHS. --------

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

# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES --------
print('---------------------')
print('We will now check if any of the potential reference points have anomalous rotation measure values.')

# -------- Define "anomalous"
# Load and unpack all the rotation measure data for the region of interest
##
matchedRMExtinctionData = pd.read_csv(MatchedRMExtincPath, sep='\t')

# Choose a rotation measure corresponding to anomalous
rm_avg = np.mean(matchedRMExtinctionData['Rotation_Measure(rad/m2)'])
rm_std = np.std(matchedRMExtinctionData['Rotation_Measure(rad/m2)'])

coeffSTD = config.anomalousSTDNum
rm_upperLimit = rm_avg + coeffSTD * rm_std
rm_lowerLimit = rm_avg - coeffSTD * rm_std
print("\t-Anomalous rotation measure values have been defined to be greater or less"
                                   " than the the suggested {} standard deviations from the mean (rm < {:.2f}rad/m^2 or"
                                   " rm > {:.2f}rad/m^2)".format(coeffSTD, rm_lowerLimit, rm_upperLimit))
# -------- Define "anomalous".

# -------- For each potential reference point
anomalousRMIndex = []
for i in range(len(AllPotenitalRefPoints.AllRefPoints)):
    idNum = AllPotenitalRefPoints.AllRefPoints['ID#'][i]
    if AllPotenitalRefPoints.AllRefPoints['Rotation_Measure(rad/m2)'][i] < rm_lowerLimit or \
            AllPotenitalRefPoints.AllRefPoints['Rotation_Measure(rad/m2)'][i] > rm_upperLimit:
        anomalousRMIndex.append(i + 1)  # To identify points numbered in order of increasing extinction
# -------- For each potential reference point.
print('The potential reference point(s) {} have anomalous rotation measure values'.format(anomalousRMIndex))
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES. --------