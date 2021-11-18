"""
This is the stage of the BLOSMapping method where the optimal number of reference points is determined and magnetic
field values are calculated

THIS IS THE SECOND VERSION OF THIS
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
      'specified threshold. For nearby clouds that appear off the disk, an appropriate threshold value is 1.')
Av_threshold = 1
chooseAvThreshold = input("Given this information, would you like to set the threshold extinction to the suggested {}? "
                          "(y/n)".format(Av_threshold))
if chooseAvThreshold == 'n':
    Av_threshold = float(input('Please enter the threshold extinction you would like to use instead: '))

# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- FIND ALL POTENTIAL REFERENCE POINTS --------
AllPotenitalRefPoints = FindAllPotentialReferencePoints(cloudName, Av_threshold,
                                                        saveFilePath=saveFilePath_ALlPotentialRefPoints)
print('Based on this threshold extinction a total of {} potential reference points were found'
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
fig = plt.figure(figsize=(12, 12), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection=wcs)

plt.title('All Potential Reference Points' + ' in the ' + cloudName + ' region\n', fontsize=12, y=1.08)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
plt.scatter(x_AllRef, y_AllRef, marker='o', facecolor='green', linewidth=.5, edgecolors='black', s=50)

# ---- Annotate the chosen reference points
for i, number in enumerate(n_AllRef):
    # Each point is labelled with its Identification Number
    ax.annotate(number, (x_AllRef[i], y_AllRef[i]), size=9, color='w')
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
cb = plt.colorbar(im, orientation='vertical', ticklocation='right', fraction=0.02, pad=0.145)
cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
# ---- Style the colour bar.

# ---- Display or save the figure
# plt.show()
saveFigurePath_RefPointMap = saveFigureDir_RefPointMap + 'RefPointMap_AllPotentialRefPoints.png'
plt.savefig(saveFigurePath_RefPointMap)
plt.close()
# ---- Display or save the figure.
# -------- CREATE A FIGURE - ALL POTENTIAL REF POINTS MAP. --------

# -------- FIND OPTIMAL NUMBER OF REFERENCE POINTS USING "ALL POTENTIAL REFERENCE POINTS" --------
print('---------------------')
print('By analyzing the stability of calculated BLOS values as a function of number of reference points from 1 to the '
      'total number of reference points ({}):'.format(AllPotenitalRefPoints.numAllRefPoints))
OptimalRefPoints_from_AllPotentialRefPoints = FindOptimalRefPoints(cloudName, AllPotenitalRefPoints.AllRefPoints,
                                                                   'True',
                                                                   saveFigurePath_BLOSvsNRef_AllPotentialRefPoints)

OptimalNumRefPoints_from_AllPotentialRefPoints = OptimalRefPoints_from_AllPotentialRefPoints. \
    Optimal_NumRefPoints_firstMode

chooseOptimalNumRefPoints = input("Given this information, would you like to select the suggested {} reference points? "
                                  "(y/n)".
                                  format(OptimalNumRefPoints_from_AllPotentialRefPoints))
if chooseOptimalNumRefPoints == 'n':
    OptimalNumRefPoints_from_AllPotentialRefPoints = float(input('Please enter the number of reference points you '
                                                                 'would like to use instead: '))
    if OptimalNumRefPoints_from_AllPotentialRefPoints > AllPotenitalRefPoints.numAllRefPoints:
        print('The number of reference points chosen exceeds the total number of potential reference points.  '
              'Using the total number of potential reference points ({})'.format(AllPotenitalRefPoints.numAllRefPoints))

print('Please review the BLOS trend stability plot at {} before confirming the number of reference points you would '
      'like to use.'.format(saveFigurePath_BLOSvsNRef_AllPotentialRefPoints))
print('---------------------\n')
# -------- FIND OPTIMAL NUMBER OF REFERENCE POINTS USING "ALL POTENTIAL REFERENCE POINTS". --------

# -------- EXTRACT POTENTIAL REFERENCE POINTS --------
# Extract {OptimalNumRefPoints_from_AllPotentialRefPoints} points from the table of all potential reference points
chosenRefPoints = AllPotenitalRefPoints.AllRefPoints.loc[:OptimalNumRefPoints_from_AllPotentialRefPoints - 1]
# -------- EXTRACT POTENTIAL REFERENCE POINTS. --------

print(chosenRefPoints)

# -------- PREPARE TO PLOT POTENTIAL REFERENCE POINTS --------
n = list(chosenRefPoints['ID#'])
Ra_ref = list(chosenRefPoints['Ra(deg)'])
Dec_ref = list(chosenRefPoints['Dec(deg)'])
RM_ref = list(chosenRefPoints['Rotation_Measure(rad/m2)'])
Av_ref = list(chosenRefPoints['Extinction_Value'])
# ---- Convert Ra and Dec of reference points into pixel values of the fits file
x_ref = []  # x pixel coordinate of reference
y_ref = []  # y pixel coordinate of reference
for i in range(len(Ra_ref)):
    pixelRow, pixelColumn = wcs.wcs_world2pix(Ra_ref[i], Dec_ref[i], 0)
    x_ref.append(pixelRow)
    y_ref.append(pixelColumn)
# ---- Convert Ra and Dec of reference points into pixel values of the fits file.
# -------- PREPARE TO PLOT reference POINTS. --------

# -------- CREATE A FIGURE - REF POINT MAP --------
fig = plt.figure(figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection=wcs)

plt.title('Selected Reference Points' + ' in the ' + cloudName + ' region\n', fontsize=12, y=1.08)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
plt.scatter(x_ref, y_ref, marker='o', facecolor='green', linewidth=.5, edgecolors='black', s=50)

# ---- Annotate the chosen reference points
for i, number in enumerate(n):
    # Each point is labelled with its Identification Number
    ax.annotate(number, (x_ref[i], y_ref[i]), size=9, color='w')
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
cb = plt.colorbar(im, orientation='vertical', ticklocation='right', fraction=0.02, pad=0.145)
cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
# ---- Style the colour bar.

# ---- Display or save the figure
# plt.show()
saveFigurePath_RefPointMap = saveFigureDir_RefPointMap + 'RefPointMap_NRef' + str(len(chosenRefPoints)) + '.png'
plt.savefig(saveFigurePath_RefPointMap)
plt.close()
# ---- Display or save the figure.
# -------- CREATE A FIGURE - REF POINT MAP. --------

print('Please review the file {} before confirming your final selection of reference points.\n'.
      format(saveFigurePath_RefPointMap))

# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION --------
print('---------------------')
print('We will now check if any of the chosen reference points are near a region of high extinction.')
# -------- Define the range
# The distance the point can be from a region of high extinction and still be thought to sample the background
cloudDistance = regionOfInterest.distance  # [pc]
cloudJeansLength = 1.  # [pc] MEHRNOOSH FIND THIS
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
for i in chosenRefPoints.index:
    px = chosenRefPoints['Extinction_Index_x'][i]
    py = chosenRefPoints['Extinction_Index_y'][i]

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
        nearHighExtinctionRegion.append(i)
    # ---- Cycle through extinction values within the range.
# -------- For each potential reference point.
if len(nearHighExtinctionRegion) != 0:
    print('The potential reference point(s) {} are near a region of high extinction'.format(nearHighExtinctionRegion))

    chooseRemoveNearHighExtinction = input('Based on this information would you like to eliminate any points as'
                                           ' reference points? (y/n)')
    if chooseRemoveNearHighExtinction == 'y':
        wantToRemoveAll = input('\t-Would you like to eliminate all of the suggested points ({})? (y/n)'.
                                format(nearHighExtinctionRegion))
        if wantToRemoveAll == 'y':
            updatedRefPoints = AllPotenitalRefPoints.AllRefPoints.drop(nearHighExtinctionRegion).reset_index(drop=True)
            if len(updatedRefPoints) == 0:
                print('All potential reference point were eliminated. Exiting the program')
                exit()
            # Now redo the stability trend
            saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_NearHighExtn = os.path.join(currentDir,
                                                                                     'FileOutput/' + cloudName + '/Plots/BLOS_vs_NRef_ChosenRefPoints_wo_NearHighExtn.png')
            print('Please review the updated BLOS trend stability plot at  {}'.
                  format(saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_NearHighExtn))
            OptimalRefPoints_from_chosenRefPoints = FindOptimalRefPoints(cloudName, updatedRefPoints,
                                                                         'False',
                                                                         saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_NearHighExtn)

        if wantToRemoveAll == 'n':
            whichToRemove = [int(item) for item in input('\tWhich points would you like to eliminate? (Enter as comma '
                                                         'separated values)').split(',')]
            updatedRefPoints = AllPotenitalRefPoints.AllRefPoints.drop(whichToRemove).reset_index(drop=True)
            if len(updatedRefPoints) == 0:
                print('All potential reference point were eliminated. Exiting the program')
                exit()
            # Now redo the stability trend
            saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_NearHighExtn = os.path.join(currentDir,
                                                                                     'FileOutput/' + cloudName + '/Plots/BLOS_vs_NRef_ChosenRefPoints_wo_NearHighExtn.png')
            print('Please review the updated BLOS trend stability plot at  {}'.
                  format(saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_NearHighExtn))
            OptimalRefPoints_from_chosenRefPoints = FindOptimalRefPoints(cloudName, updatedRefPoints,
                                                                         'False',
                                                                         saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_NearHighExtn)
else:
    print('Based on the specified range and definition of high extinction, the none of the potential reference points'
          ' are near a region of high extinction')
print('---------------------\n')
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS ARE NEAR A REGION OF HIGH EXTINCTION. --------

# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES --------
print('---------------------')
print('We will now check if any of the chosen reference points have anomalous rotation measure values.')

# -------- Define "anomalous"
# Choose a rotation measure corresponding to anomalous
rm_avg = np.mean(chosenRefPoints['Rotation_Measure(rad/m2)'])
rm_std = np.std(chosenRefPoints['Rotation_Measure(rad/m2)'])

coeff = 2
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
for i in chosenRefPoints.index:
    if chosenRefPoints['Rotation_Measure(rad/m2)'][i] < rm_lowerLimit or \
            chosenRefPoints['Rotation_Measure(rad/m2)'][i] > rm_upperLimit:
        anomalousRMIndex.append(i)
# -------- For each potential reference point.

if len(anomalousRMIndex) != 0:
    print('The potential reference point(s) {} have anomalous rotation measure values'.format(anomalousRMIndex))

    chooseRemoveAnomalousRM = input('Based on this information would you like to eliminate any points as'
                                    ' reference points? (y/n)')
    if chooseRemoveAnomalousRM == 'y':
        wantToRemoveAllAnomalousRM = input('\t-Would you like to eliminate all of the suggested points ({})? (y/n)'.
                                           format(anomalousRMIndex))

        if wantToRemoveAllAnomalousRM == 'y':
            updatedRefPoints = AllPotenitalRefPoints.AllRefPoints.drop(anomalousRMIndex).reset_index(drop=True)
            if len(updatedRefPoints) == 0:
                print('All potential reference point were eliminated. Exiting the program')
                exit()
            # Now redo the stability trend
            saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_anomRM = os.path.join(currentDir,
                                                                                     'FileOutput/' + cloudName + '/Plots/BLOS_vs_NRef_ChosenRefPoints_wo_anomRM.png')
            print('Please review the updated BLOS trend stability plot at  {}'.
                  format(saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_anomRM ))
            OptimalRefPoints_from_chosenRefPoints = FindOptimalRefPoints(cloudName, updatedRefPoints,
                                                                         'False',
                                                                         saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_anomRM)

        if wantToRemoveAllAnomalousRM == 'n':
            whichToRemoveAnomalousRM = [int(item) for item in input('\tWhich points would you like to eliminate? (Enter'
                                                                    ' as comma separated values)').split(',')]
            updatedRefPoints = AllPotenitalRefPoints.AllRefPoints.drop(whichToRemoveAnomalousRM).reset_index(drop=True)
            if len(updatedRefPoints) == 0:
                print('All potential reference point were eliminated. Exiting the program')
                exit()
            # Now redo the stability trend
            saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_anomRM = os.path.join(currentDir,
                                                                                     'FileOutput/' + cloudName + '/Plots/BLOS_vs_NRef_ChosenRefPoints_wo_anomRM.png')
            print('Please review the updated BLOS trend stability plot at  {}'.
                  format(saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_anomRM ))
            OptimalRefPoints_from_chosenRefPoints = FindOptimalRefPoints(cloudName, updatedRefPoints,
                                                                         'False',
                                                                         saveFigurePath_BLOSvsNRef_ChosenRefPoints_wo_anomRM)
else:
    print('\tBased on the definition of anomalous rotation measure values, none of the potential reference points'
          ' have anomalous rotation measure values')
print('---------------------\n')
# -------- CHECK TO SEE IF ANY POTENTIAL POINTS HAVE ANOMALOUS RM VALUES. --------

# -------- CALCULATE AND SAVE REFERENCE VALUES --------
cols = ['Number of Reference Points', 'Reference Extinction', 'Reference RM', 'Reference RM AvgErr',
        'Reference RM Std']
referenceData = pd.DataFrame(columns=cols)

referenceData['Number of Reference Points'] = [len(chosenRefPoints)]

referenceData['Reference RM'] = [np.mean(chosenRefPoints['Rotation_Measure(rad/m2)'])]

referenceData['Reference RM AvgErr'] = [np.mean(chosenRefPoints['RM_Err(rad/m2)'])]

# Standard error of the sampled mean:
referenceData['Reference RM Std'] = [np.std(chosenRefPoints['Rotation_Measure(rad/m2)'], ddof=1) / \
                                     np.sqrt(len(chosenRefPoints['Rotation_Measure(rad/m2)']))]

referenceData['Reference Extinction'] = [np.mean(chosenRefPoints['Extinction_Value'])]

referenceData.to_csv(saveFilePath_ReferenceData, index=False)
# -------- CALCULATE AND SAVE REFERENCE VALUES. --------

# -------- SAVE REFERENCE POINTS  --------
chosenRefPoints.to_csv(saveFilePath_ReferencePoints, index=False)
# -------- SAVE REFERENCE POINTS. --------
