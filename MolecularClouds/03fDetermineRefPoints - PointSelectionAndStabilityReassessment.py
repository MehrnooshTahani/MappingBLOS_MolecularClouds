"""
This is the third stage of the BLOSMapping method where the reference points are determined
"""
import os

import pandas as pd
import numpy as np

from astropy.wcs import WCS
from astropy.io import fits

import matplotlib.pyplot as plt

from LocalLibraries.RegionOfInterest import Region
from LocalLibraries.FindOptimalRefPoints import FindOptimalRefPoints
from LocalLibraries.CalculateB import CalculateB
import LocalLibraries.config as config
import LocalLibraries.RefJudgeLib as rjl

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
saveFilePath_ALlPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_allPotRefPoints + config.cloud + '.txt')
saveFilePath_ReferencePoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_selRefPoints + config.cloud + '.txt')
saveFilePath_ReferenceData = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_refData + config.cloud + '.txt')

saveFigurePath_BLOSvsNRef_AllPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots, 'BLOS_vs_NRef_AllPotentialRefPoints.png')
saveFigurePath_BLOSvsNRef_ChosenPotentialRefPoints = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_plots, 'BLOS_vs_NRef_ChosenRefPoints.png')

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

# -------- LOAD ALL POTENTIAL REFERENCE POINTS --------
AllPotentialRefPoints = pd.read_csv(saveFilePath_ALlPotentialRefPoints)
print('---------------------\n')
# -------- LOAD ALL POTENTIAL REFERENCE POINTS. --------

# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

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

# -------- ASK THE USER WHICH POINTS THEY WANT TO USE AS REFERENCE POINTS --------
chosenRefPoints_Num = [int(item) - 1 for item in input('Please enter the numbers of the reference points you would '
                                                       'like to use as comma separated values').split(',')]
chosenRefPoints = AllPotentialRefPoints.loc[chosenRefPoints_Num].sort_values('Extinction_Value')
#Todo: Here's where it goes.
print(chosenRefPoints)
# -------- ASK THE USER WHICH POINTS THEY WANT TO USE AS REFERENCE POINTS. --------

# -------- FIND REGIONS TO SPLIT THE CLOUD INTO. --------
cloudCenterX, cloudCenterY = rjl.findWeightedCenter(hdu.data, regionOfInterest.xmin, regionOfInterest.xmax, regionOfInterest.ymin, regionOfInterest.ymax)
m, b = rjl.getDividingLine(hdu.data, regionOfInterest.xmin, regionOfInterest.xmax, regionOfInterest.ymin, regionOfInterest.ymax)
mPerp, bPerp = rjl.getPerpendicularLine(cloudCenterX, cloudCenterY, m)
# -------- FIND REGIONS TO SPLIT THE CLOUD INTO. --------

# -------- SORT REF POINTS INTO THESE REGIONS. --------
Q1 = []
Q2 = []
Q3 = []
Q4 = []
for i in range(len(chosenRefPoints)):
    idNum = chosenRefPoints['ID#'][i]
    px = chosenRefPoints['Extinction_Index_x'][i]
    py = chosenRefPoints['Extinction_Index_y'][i]

    # ---- Sort into quadrant
    aboveCloudLine = rjl.isPointAboveLine(px, py, m, b)
    aboveCloudPerpLine = rjl.isPointAboveLine(px, py, mPerp, bPerp)

    if aboveCloudLine and aboveCloudPerpLine:
        Q1.append(i+1)
    elif aboveCloudLine and not aboveCloudPerpLine:
        Q2.append(i+1)
    elif not aboveCloudLine and aboveCloudPerpLine:
        Q3.append(i+1)
    elif not aboveCloudLine and not aboveCloudPerpLine:
        Q4.append(i+1)
    # ---- Sort into quadrant
# -------- SORT REF POINTS INTO THESE REGIONS. --------

# -------- OUTPUT RESULTS. --------
print("The potential reference points, sorted by quadrant, are:")
print("Q1: {}".format(Q1))
print("Q2: {}".format(Q2))
print("Q3: {}".format(Q3))
print("Q4: {}".format(Q4))
# -------- OUTPUT RESULTS. --------

# -------- REASSESS STABILITY --------
'''
We are going to start the plot off with all of the chosen reference points
 - which are in order of increasing extinction.  They may have jumps but this is okay
However, we also want to include points which come after the chosen reference points
The following selects all of the chosen reference points and then adds any of the potential reference points with
extinction greater than the extinction of the last chosen reference point/
'''
RefPoints = chosenRefPoints[:-1].append(AllPotentialRefPoints.set_index('ID#').
                                        loc[list(chosenRefPoints['ID#'])[-1]:].reset_index())\
    .reset_index(drop=True)
# -------- Read the reference point data
numRefPoints = len(RefPoints)
# -------- Read the reference point data.

# -------- Create a table for all blos data
# The rows of this table will represent the number of reference points and the columns of this table will
# represent the individual BLOS points.  Each entry in the table is a calculated BLOS value.
AllData = pd.DataFrame()
# -------- Create a table for all blos data.

# -------- Calculate blos as a function of # ref points
for num in range(numRefPoints):
    # -------- Extract {num} points from the table of potential reference points
    # The extracted potential reference points will be referred to as "candidates"
    candidateRefPoints = RefPoints.loc[:num]
    # -------- Extract {num} points from the table of potential reference points.

    # -------- Use the candidate reference points to calculate BLOS
    B = CalculateB(regionOfInterest.AvFilePath, MatchedRMExtincPath, candidateRefPoints, saveFilePath=None)
    BLOSData = B.BLOSData.set_index('ID#', drop=True)
    # -------- Use the candidate reference points to calculate BLOS

    # -------- Add calculated BLOS to the table of BLOS vs number of reference points
    AllData[str(num + 1)] = BLOSData['Magnetic_Field(uG)']
    # -------- Add calculated BLOS to the table of BLOS vs number of reference points

# -------- CALCULATE BLOS AS A FUNCTION OF # REF POINTS. --------

# -------- FIND OPTIMAL NUM REF POINTS --------
# If a point has been used as a candidate reference point at any time it will not be used to determine the
# optimal number of reference points
DataNoRef = AllData.copy().drop(list(RefPoints['ID#']), errors='ignore')
Identifiers = list(DataNoRef.index)
DataNoRef = DataNoRef.reset_index(drop=True)

# -------- CREATE A FIGURE --------
plt.figure(figsize=(6, 4), dpi=120, facecolor='w', edgecolor='k')

plt.title('Calculated BLOS value as a function of the number of reference points \n ' + cloudName, fontsize=12,
          y=1.08)
plt.xlabel('Number of reference points')
plt.ylabel('Calculated BLOS value ' + r'($\mu G$)')

x = [int(col) for col in DataNoRef.columns]
plt.xticks(x, list(DataNoRef.columns))

cmap = plt.get_cmap('terrain')
colors = [cmap(i) for i in np.linspace(0, 1, len(DataNoRef.index))]

# For each BLOS Point
for i, number in enumerate(DataNoRef.index):
    plt.plot(x, list(DataNoRef.loc[number]), '-o', color=colors[i], markersize=3)

yLower, yUpper = plt.ylim()
plt.vlines(OptimalNumRefPoints_from_AllPotentialRefPoints, yLower, yUpper, color='black', label='Suggested optimal '
                                                                                                'number of reference '
                                                                                                'points')
plt.legend(loc='center right', bbox_to_anchor=(1.1, 0.5), ncol=2, framealpha=1)

plt.savefig(saveFigurePath_BLOSvsNRef_ChosenPotentialRefPoints)
plt.show()
plt.close()
# -------- CREATE A FIGURE. --------
# -------- REASSESS STABILITY. --------

# -------- CALCULATE AND SAVE REFERENCE VALUES --------
cols = ['Number of Reference Points', 'Reference Extinction', 'Reference RM', 'Reference RM AvgErr',
        'Reference RM Std']
referenceData = pd.DataFrame(columns=cols)

referenceData['Number of Reference Points'] = [len(chosenRefPoints)]
referenceData['Reference RM'] = [np.mean(chosenRefPoints['Rotation_Measure(rad/m2)'])]
referenceData['Reference RM AvgErr'] = [np.mean(chosenRefPoints['RM_Err(rad/m2)'])]

# Standard error of the sampled mean:
referenceData['Reference RM Std'] = [np.std(chosenRefPoints['Rotation_Measure(rad/m2)'], ddof=1) /
                                     np.sqrt(len(chosenRefPoints['Rotation_Measure(rad/m2)']))]
referenceData['Reference Extinction'] = [np.mean(chosenRefPoints['Extinction_Value'])]
referenceData.to_csv(saveFilePath_ReferenceData, index=False)
print('Reference values were saved to {}'.format(saveFilePath_ReferenceData))
# -------- CALCULATE AND SAVE REFERENCE VALUES. --------

# -------- SAVE REFERENCE POINTS  --------
chosenRefPoints.to_csv(saveFilePath_ReferencePoints, index=False)
print('Chosen reference points were saved to {}'.format(saveFilePath_ReferencePoints))
# -------- SAVE REFERENCE POINTS. --------
