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
from CalculateB import CalculateB

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
saveFilePath_ALlPotentialRefPoints = os.path.join(currentDir, 'FileOutput/' + cloudName + '/AllPotentialRefPoints'
                                                  + cloudName + '.txt')
saveFilePath_ReferencePoints = os.path.join(currentDir, 'FileOutput/' + cloudName + '/RefPoints' + cloudName + '.txt')
saveFilePath_ReferenceData = os.path.join(currentDir, 'FileOutput/' + cloudName + '/ReferenceData' + cloudName + '.txt')
saveFigurePath_BLOSvsNRef_AllPotentialRefPoints = os.path.join(currentDir, 'FileOutput/' + cloudName +
                                                               '/Plots/BLOS_vs_NRef_AllPotentialRefPoints.png')
saveFigurePath_BLOSvsNRef_ChosenPotentialRefPoints = os.path.join(currentDir, 'FileOutput/' + cloudName +
                                                                  '/Plots/BLOS_vs_NRef_ChosenRefPoints.png')
saveFigureDir_RefPointMap = os.path.join(currentDir, 'FileOutput/' + cloudName + '/Plots/')
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
if abs(regionOfInterest.cloudLatitude) < 15:
    Av_threshold = 2
    print('\t-For clouds that appear near the disk, such as {}, an appropriate threshold value is {}.'
          .format(cloudName, Av_threshold))
else:
    Av_threshold = 1
    print('\t-For clouds that appear off the disk, such as {}, an appropriate threshold value is {}.'
          .format(cloudName, Av_threshold))

chooseAvThreshold = input("Given this information, would you like to set the threshold extinction to the suggested {}? "
                          "(y/n)".format(Av_threshold))
if chooseAvThreshold == 'n':
    Av_threshold = float(input('Please enter the threshold extinction you would like to use instead: '))
# -------- CHOOSE THE THRESHOLD EXTINCTION. --------

# -------- FIND ALL POTENTIAL REFERENCE POINTS --------
AllPotenitalRefPoints = FindAllPotentialReferencePoints(cloudName, Av_threshold,
                                                        saveFilePath=saveFilePath_ALlPotentialRefPoints)
print('Based on this threshold extinction, a total of {} potential reference points were found.'
      .format(AllPotenitalRefPoints.numAllRefPoints))

print(AllPotenitalRefPoints.AllRefPoints)
print('---------------------\n')
# -------- FIND ALL POTENTIAL REFERENCE POINTS. --------

# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS --------
n_AllRef = list(AllPotenitalRefPoints.AllRefPoints['ID#'])
Ra_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Ra(deg)'])
Dec_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Dec(deg)'])
RM_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Rotation_Measure(rad/m2)'])
Av_AllRef = list(AllPotenitalRefPoints.AllRefPoints['Extinction_Value'])
# ---- Convert Ra and Dec of reference points into pixel values of the fits file
x_AllRef = []  # x pixel coordinate of reference
y_AllRef = []  # y pixel coordinate of reference
for i in range(len(Ra_AllRef)):
    pixelRow, pixelColumn = wcs.wcs_world2pix(Ra_AllRef[i], Dec_AllRef[i], 0)
    x_AllRef.append(pixelRow)
    y_AllRef.append(pixelColumn)
# ---- Convert Ra and Dec of reference points into pixel values of the fits file.
# -------- PREPARE TO PLOT ALL POTENTIAL REFERENCE POINTS. --------

# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP --------
fig = plt.figure(figsize=(8, 8), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection=wcs)

plt.title('All Potential Reference Points' + ' in the ' + cloudName + ' region\n', fontsize=12, y=1.08)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
plt.scatter(x_AllRef, y_AllRef, marker='o', facecolor='green', linewidth=.5, edgecolors='black', s=50)

# ---- Annotate the chosen reference points
text = []
for i, number in enumerate(n_AllRef):
    # Each point is labelled in order of increasing extinction value
    # To label with ID number use: txt = ax.text(x_AllRef[i], y_AllRef[i], str(number), size=9, color='w')
    txt = ax.text(x_AllRef[i], y_AllRef[i], str(i + 1), size=9, color='w')
    text.append(txt)
adjustText.adjust_text(text)
# ---- Annotate the chosen reference points

# ---- Style the main axes and their grid
if regionOfInterest.xmax and regionOfInterest.xmin != 'none':
    ax.set_xlim(regionOfInterest.xmin, regionOfInterest.xmax)
if regionOfInterest.ymax and regionOfInterest.ymin != 'none':
    ax.set_ylim(regionOfInterest.ymin, regionOfInterest.ymax)

ra = ax.coords[0]
dec = ax.coords[1]
ra.set_major_formatter('d')
dec.set_major_formatter('d')
ra.set_axislabel('RA (degree)')
dec.set_axislabel('Dec (degree)')

dec.set_ticks(number=10)
ra.set_ticks(number=20)
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_minor_frequency(10)

ra.grid(color='black', alpha=0.5, linestyle='solid')
dec.grid(color='black', alpha=0.5, linestyle='solid')
# ---- Style the main axes and their grid.

# ---- Style the overlay and its grid
overlay = ax.get_coords_overlay('galactic')

overlay[0].set_axislabel('Longitude')
overlay[1].set_axislabel('Latitude')

overlay[0].set_ticks(color='grey', number=20)
overlay[1].set_ticks(color='grey', number=20)

overlay.grid(color='grey', linestyle='solid', alpha=0.7)
# ---- Style the overlay and its grid.

# ---- Style the colour bar
if regionOfInterest.fitsDataType == 'HydrogenColumnDensity':
    cb = plt.colorbar(im, ticklocation='right', fraction=0.02, pad=0.145, format='%.0e')
    cb.ax.set_title('Hydrogen Column Density', linespacing=0.5, fontsize=12)
elif regionOfInterest.fitsDataType == 'VisualExtinction':
    cb = plt.colorbar(im, ticklocation='right', fraction=0.02, pad=0.145)
    cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
# ---- Style the colour bar.

# ---- Display or save the figure
saveFigurePath_RefPointMap = saveFigureDir_RefPointMap + 'RefPointMap_AllPotentialRefPoints.png'
plt.savefig(saveFigurePath_RefPointMap)
plt.show()
plt.close()
# ---- Display or save the figure.
print('Saving the map of all potential reference points to '+saveFigurePath_RefPointMap)
# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP. --------

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

# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION --------
print('---------------------')
print('We will now check if any of the potential reference points are near a region of high extinction.')
# -------- Define the range
# The distance the point can be from a region of high extinction and still be thought to sample the background
cloudDistance = regionOfInterest.distance  # [pc]
cloudJeansLength = 1.  # [pc] CHECK THIS
minDiff = cloudJeansLength / cloudDistance  # [deg]

minDiff_pix = minDiff / abs(hdu.header['CDELT1'])
NDelt = 5 * math.ceil(minDiff_pix)  # Round up
chooseNDelt = input("\t-Would you like the define a region around the given point to the suggested {} pixels? (y/n)".
                    format(NDelt))
if chooseNDelt == 'n':
    NDelt = int(float(input('\tPlease enter the value you would like to use instead: ')))

# Choose the minimum extinction value which you want to correspond to an "on" position
highExtinctionThreshold = 5 * Av_threshold

chooseHighAvThreshold = input("\t-Would you like to define a region of high extinction to the suggested Av={}? (y/n)".
                              format(highExtinctionThreshold))
if chooseHighAvThreshold == 'n':
    highExtinctionThreshold = float(input('\tPlease enter the value you would like to use instead: '))
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
    highExtinction = 0
    for pxx in range(ind_xmin, ind_xmax):
        for pyy in range(ind_ymin, ind_ymax):
            if 0 <= pxx < hdu.data.shape[1] and 0 <= pyy < hdu.data.shape[0]:
                extinction_val = hdu.data[pyy, pxx]
                if extinction_val > highExtinctionThreshold:
                    highExtinction = 1
    if highExtinction == 1:
        nearHighExtinctionRegion.append(i + 1)  # To identify points numbered in order of increasing extinction
    # ---- Cycle through extinction values within the range.
# -------- For each potential reference point.
if len(nearHighExtinctionRegion) != 0:
    print('The potential reference point(s) {} are near a region of high extinction'.format(nearHighExtinctionRegion))

else:
    print('Based on the specified range and definition of high extinction, the none of the potential reference points'
          ' are near a region of high extinction')
print('---------------------\n')
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION. --------


# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES --------
print('---------------------')
print('We will now check if any of the potential reference points have anomalous rotation measure values.')

# -------- Define "anomalous"
# Load and unpack all the rotation measure data for the region of interest
currentDir = os.path.abspath(os.getcwd())
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/' + cloudName + '/MatchedRMExtinction'
                                   + cloudName + '.txt')
matchedRMExtinctionData = pd.read_csv(MatchedRMExtincPath, sep='\t')

# Choose a rotation measure corresponding to anomalous
rm_avg = np.mean(matchedRMExtinctionData['Rotation_Measure(rad/m2)'])
rm_std = np.std(matchedRMExtinctionData['Rotation_Measure(rad/m2)'])

coeff = 3
rm_upperLimit = rm_avg + coeff * rm_std
rm_lowerLimit = rm_avg - coeff * rm_std
chooseAnomalousRMThreshold = input("\t-Would you like to define anomalous rotation measure values to be greater or less"
                                   " than the the suggested {} standard deviations from the mean (rm < {:.2f}rad/m^2 or"
                                   " rm > {:.2f}rad/m^2)? (y/n)".format(coeff, rm_lowerLimit, rm_upperLimit))
if chooseAnomalousRMThreshold == 'n':
    anomalousRMThreshold = float(input('\tPlease enter the value you would like to use instead: '))
# -------- Define "anomalous".

# -------- For each potential reference point
anomalousRMIndex = []
for i in range(len(AllPotenitalRefPoints.AllRefPoints)):
    idNum = AllPotenitalRefPoints.AllRefPoints['ID#'][i]
    if AllPotenitalRefPoints.AllRefPoints['Rotation_Measure(rad/m2)'][i] < rm_lowerLimit or \
            AllPotenitalRefPoints.AllRefPoints['Rotation_Measure(rad/m2)'][i] > rm_upperLimit:
        anomalousRMIndex.append(i + 1)  # To identify points numbered in order of increasing extinction
# -------- For each potential reference point.

if len(anomalousRMIndex) != 0:
    print('The potential reference point(s) {} have anomalous rotation measure values'.format(anomalousRMIndex))

else:
    print('\tBased on the definition of anomalous rotation measure values, none of the potential reference points'
          ' have anomalous rotation measure values')
print('---------------------\n')
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES. --------


# -------- ASK THE USER WHICH POINTS THEY WANT TO USE AS REFERENCE POINTS --------
chosenRefPoints_Num = [int(item) - 1 for item in input('Please enter the numbers of the reference points you would '
                                                       'like to use as comma separated values').split(',')]
chosenRefPoints = AllPotenitalRefPoints.AllRefPoints.loc[chosenRefPoints_Num].sort_values('Extinction_Value')

print(chosenRefPoints)
# -------- ASK THE USER WHICH POINTS THEY WANT TO USE AS REFERENCE POINTS. --------

# -------- REASSESS STABILITY --------
# -------- Load matched rm and extinction data
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/' + cloudName + '/MatchedRMExtinction'
                                   + cloudName + '.txt')
# -------- Load matched rm and extinction data.
'''
We are going to start the plot off with all of the chosen reference points
 - which are in order of increasing extinction.  They may have jumps but this is okay
However, we also want to include points which come after the chosen reference points
The following selects all of the chosen reference points and then adds any of the potential reference points with
extinction greater than the extinction of the last chosen reference point/
'''
RefPoints = chosenRefPoints[:-1].append(AllPotenitalRefPoints.AllRefPoints.set_index('ID#').
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
    B = CalculateB(regionOfInterest.AvFilePath, MatchedRMExtincPath, candidateRefPoints, saveFilePath='none')
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
