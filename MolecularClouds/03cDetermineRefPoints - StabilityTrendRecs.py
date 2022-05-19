"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os
import pandas as pd

from LocalLibraries.RegionOfInterest import Region
from LocalLibraries.FindOptimalRefPoints import FindOptimalRefPoints
import LocalLibraries.config as config

#Todo: Steps duplicated in f. Unclear if needed. Could be used as an early filter?

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
saveFilePath_ALlPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_allPotRefPoints + config.cloud + '.txt')
saveFigurePath_BLOSvsNRef_AllPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots, 'BLOS_vs_NRef_AllPotentialRefPoints.png')
# -------- DEFINE FILES AND PATHS. --------

# -------- LOAD ALL POTENTIAL REFERENCE POINTS --------
AllPotentialRefPoints = pd.read_csv(saveFilePath_ALlPotentialRefPoints)
print('---------------------\n')
# -------- LOAD ALL POTENTIAL REFERENCE POINTS. --------

# -------- FIND OPTIMAL NUMBER OF REFERENCE POINTS USING "ALL POTENTIAL REFERENCE POINTS" --------
print('---------------------')
print('By analyzing the stability of calculated BLOS values as a function of number of reference points from 1 to the '
      'total number of reference points ({}):'.format(len(AllPotentialRefPoints)))
OptimalRefPoints_from_AllPotentialRefPoints = FindOptimalRefPoints(cloudName, AllPotentialRefPoints,
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

    if OptimalNumRefPoints_from_AllPotentialRefPoints > len(AllPotentialRefPoints):
        print('The number of reference points chosen exceeds the total number of potential reference points.  '
              'Using the total number of potential reference points ({})'.format(len(AllPotentialRefPoints)))
        print('The recommended reference points, numbered in order of increasing extinction, are: {}'.format(
            list([i + 1 for i in range(0, OptimalNumRefPoints_from_AllPotentialRefPoints)])))

print('Please review the BLOS trend stability plot at {} before confirming the number of reference points you would '
      'like to use.'.format(saveFigurePath_BLOSvsNRef_AllPotentialRefPoints))
print('---------------------\n')
# -------- FIND OPTIMAL NUMBER OF REFERENCE POINTS USING "ALL POTENTIAL REFERENCE POINTS". --------