"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os
import pandas as pd
import numpy as np

from LocalLibraries.RegionOfInterest import Region
import LocalLibraries.config as config

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

# -------- FIND ALL POTENTIAL REFERENCE POINTS --------
AllPotentialRefPoints = pd.read_csv(saveFilePath_ALlPotentialRefPoints)
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
for i in list(AllPotentialRefPoints.index):
    idNum = AllPotentialRefPoints['ID#'][i]
    if AllPotentialRefPoints['Rotation_Measure(rad/m2)'][i] < rm_lowerLimit or \
            AllPotentialRefPoints['Rotation_Measure(rad/m2)'][i] > rm_upperLimit:
        anomalousRMIndex.append(i + 1)  # To identify points numbered in order of increasing extinction
# -------- For each potential reference point.
print('The potential reference point(s) {} have anomalous rotation measure values'.format(anomalousRMIndex))
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES. --------