"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os

import LocalLibraries.config as config
from LocalLibraries.RegionOfInterest import Region
import LocalLibraries.RefJudgeLib as rjl

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

# -------- CHOOSE THE THRESHOLD EXTINCTION --------
if abs(regionOfInterest.cloudLatitude) < config.offDiskLatitude:
    Av_threshold = config.onDiskAvThresh
else:
    Av_threshold = config.offDiskAvThresh
# -------- CHOOSE THE THRESHOLD EXTINCTION. --------

# -------- FIND ALL POTENTIAL REFERENCE POINTS --------
# ---- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA
matchedRMExtinctionData = pd.read_csv(MatchedRMExtincPath, sep='\t')
# ---- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA

# -------- Criterion: Av < threshold
'''
We will only consider points with visual extinction less than the specified threshold value as potential 
reference points
- Here we extract these points and sort the resulting dataframe from smallest to greatest extinction 
'''

# All potential reference points are all reference points with extinction less than the threshold:
AllPotentialRefPoints, numAllRefPoints = rjl.maskRows(matchedRMExtinctionData, Av_threshold, 'Extinction_Value')
# -------- Criterion: Av < threshold.

# ---- SAVE REFERENCE POINT DATA AS A TABLE
if saveFilePath_ALlPotentialRefPoints is not None:
    AllPotentialRefPoints.to_csv(saveFilePath_ALlPotentialRefPoints, index=False)
# ---- SAVE REFERENCE POINT DATA AS A TABLE.
print('Based on the threshold extinction of {}, a total of {} potential reference points were found.'.format(Av_threshold, numAllRefPoints))
print(AllPotentialRefPoints)
print('---------------------\n')
# -------- FIND ALL POTENTIAL REFERENCE POINTS. --------