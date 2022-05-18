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
saveFigurePath_BLOSvsNRef_AllPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots, 'BLOS_vs_NRef_AllPotentialRefPoints.png')
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

# -------- FIND OPTIMAL NUMBER OF REFERENCE POINTS USING "ALL POTENTIAL REFERENCE POINTS" --------
print('---------------------')
print('By analyzing the stability of calculated BLOS values as a function of number of reference points from 1 to the '
      'total number of reference points ({}):'.format(AllPotenitalRefPoints.numAllRefPoints))
OptimalRefPoints_from_AllPotentialRefPoints = FindOptimalRefPoints(cloudName, AllPotenitalRefPoints.AllRefPoints,
                                                                   saveFigurePath_BLOSvsNRef_AllPotentialRefPoints)

OptimalNumRefPoints_from_AllPotentialRefPoints = OptimalRefPoints_from_AllPotentialRefPoints. \
    Optimal_NumRefPoints_firstMode

chooseOptimalNumRefPoints = input("Given this information, would you like to select the suggested {} reference points? "
                                  "(y/n)".
                                  format(OptimalNumRefPoints_from_AllPotentialRefPoints))
if chooseOptimalNumRefPoints == 'y':
    print('The recommended reference points, numbered in order of increasing extinction, are: {}'.format(
        list([i + 1 for i in range(0, OptimalNumRefPoints_from_AllPotentialRefPoints)])))

if chooseOptimalNumRefPoints == 'n':
    OptimalNumRefPoints_from_AllPotentialRefPoints = int(input('Please enter the number of reference points you '
                                                               'would like to use instead: '))
    print('The recommended reference points, numbered in order of increasing extinction, are: {}'.format(
        list([i + 1 for i in range(0, OptimalNumRefPoints_from_AllPotentialRefPoints)])))

    if OptimalNumRefPoints_from_AllPotentialRefPoints > AllPotenitalRefPoints.numAllRefPoints:
        print('The number of reference points chosen exceeds the total number of potential reference points.  '
              'Using the total number of potential reference points ({})'.format(AllPotenitalRefPoints.numAllRefPoints))
        print('The recommended reference points, numbered in order of increasing extinction, are: {}'.format(
            list([i + 1 for i in range(0, OptimalNumRefPoints_from_AllPotentialRefPoints)])))

print('Please review the BLOS trend stability plot at {} before confirming the number of reference points you would '
      'like to use.'.format(saveFigurePath_BLOSvsNRef_AllPotentialRefPoints))
print('---------------------\n')
# -------- FIND OPTIMAL NUMBER OF REFERENCE POINTS USING "ALL POTENTIAL REFERENCE POINTS". --------