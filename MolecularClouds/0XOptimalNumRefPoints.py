"""
This is the stage of the BLOSMapping method where the optimal number of reference points is determined and magnetic
field values are calculated

"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from statistics import mode
from Classes.RegionOfInterest import Region
from Classes.refPointFinder import RefPointFinder
from Classes.CalculateB import CalculateB


# -------- FUNCTION DEFINITION --------
def B2RGB(b):
    """
    Takes BLOS values and assigns them a marker colour and size for use in plotting BLOS data

    :param b: The BLOS value, or list of BLOS values
    :return:  A tuple of (colour, size) corresponding to the rotation measure. Note "colour" is a tuple of (RBG,alpha)
    """
    c = []  # Marker colour
    s = []  # Marker size

    for item in b:
        if abs(item) < 1000:
            s.append(abs(item) / 5)
        if abs(item) >= 1000:
            s.append(1000 / 5)

        alpha = 1  # Optional: set the transparency
        if int(np.sign(item)) == -1:
            c.append((1, 0, 0, alpha))  # Negative rotation measures assigned red
        if int(np.sign(item)) == 1:
            c.append((0, 0, 1, alpha))  # Positive rotation measures assigned blue
        if np.sign(item) == 0:
            c.append((0, 1, 0, alpha))  # Zero-value rotation measures assigned green

    # return the list of RGBA tuples and sizes
    return c, s
# -------- FUNCTION DEFINITION. --------


# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/MatchedRMExtinction'+cloudName+'.txt')
saveFilePath_Fiducial = os.path.join(currentDir, 'FileOutput/'+cloudName+'/FiducialData.txt')
saveFileDir_RefPoints = os.path.join(currentDir, 'FileOutput/'+cloudName+'/RefPoints/')
saveFileDir_BLOS = os.path.join(currentDir, 'FileOutput/'+cloudName+'/BLOS/')
saveFigureDir_RefPointMap = os.path.join(currentDir, 'FileOutput/'+cloudName+'/Plots/')
saveFigurePath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/Plots/BLOS_vs_NRef.png')
# -------- DEFINE FILES AND PATHS. --------

# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- FIND TOTAL NUMBER OF REFERENCE POINTS --------
TotalNRef = RefPointFinder(MatchedRMExtincPath, 0, cloudName, 'none').numAllRefPoints
print('There are {} potential reference points in total'.format(TotalNRef))
# -------- FIND TOTAL NUMBER OF REFERENCE POINTS, --------

# -------- CREATE A TABLE FOR ALL BLOS DATA --------
# The rows of this table will represent the number of reference points and the columns of this table will represent the
# individual BLOS points.  Each entry in the table is a calculated BLOS value.
AllData = pd.DataFrame()
# -------- CREATE A TABLE FOR ALL BLOS DATA. --------

# -------- CREATE A TABLE FOR FIDUCIAL VALUES AS A FUNCTION OF REFERENCE POINTS --------
# The rows of this table will represent the number of reference points and the columns of this table will represent the
# fiducial values
cols = ['Optimal Number of Reference Points', 'Fiducial Extinction', 'Fiducial RM', 'Fiducial RM AvgErr',
        'Fiducial RM Std']
FiducialData = pd.DataFrame(columns=cols)
# -------- CREATE A TABLE FOR FIDUCIAL VALUES AS A FUNCTION OF REFERENCE POINTS. --------

# -------- FIND REF POINTS AND CALCULATE BLOS AS A FUNCTION OF # REF POINTS --------
indRefPoints = []
for NRef in range(1, TotalNRef+1):  # Want to be inclusive of the upper bound

    # -------- Find reference points
    saveFilePath_RefPoints = saveFileDir_RefPoints + 'NRef' + str(NRef) + '.txt'
    RefPointFinder(MatchedRMExtincPath, NRef, cloudName, saveFilePath_RefPoints)
    # -------- Find reference points.

    # -------- READ AND UNPACK REFERENCE POINT DATA --------
    referenceData = pd.read_csv(saveFilePath_RefPoints)

    for identifier in list(referenceData['ID#']):
        if identifier not in indRefPoints:
            indRefPoints.append(identifier)

    n = list(referenceData['ID#'])
    Ra_ref = list(referenceData['Ra(deg)'])
    Dec_ref = list(referenceData['Dec(deg)'])
    RM_ref = list(referenceData['Rotation_Measure(rad/m2)'])
    # -------- READ AND UNPACK REFERENCE POINT DATA. --------

    # -------- PREPARE TO PLOT reference POINTS --------
    # ---- Convert Ra and Dec of reference points into pixel values of the fits file
    x_ref = []  # x pixel coordinate of reference
    y_ref = []  # y pixel coordinate of reference
    for i in range(len(Ra_ref)):
        pixelRow, pixelColumn = wcs.wcs_world2pix(Ra_ref[i], Dec_ref[i], 0)
        x_ref.append(pixelRow)
        y_ref.append(pixelColumn)
    # ---- Convert Ra and Dec of reference points into pixel values of the fits file.
    # -------- PREPARE TO PLOT reference POINTS. --------

    # -------- Calculate BLOS
    saveFilePath_BLOS = saveFileDir_BLOS + 'BLOSNRef' + str(NRef) + '.txt'
    B = CalculateB(regionOfInterest.AvFilePath, MatchedRMExtincPath, saveFilePath_RefPoints, saveFilePath_BLOS)

    FiducialData.loc[NRef] = {'Optimal Number of Reference Points': NRef, 'Fiducial RM': B.fiducialRM,
                              'Fiducial Extinction': B.fiducialExtinction, 'Fiducial RM AvgErr': B.fiducialRMAvgErr,
                              'Fiducial RM Std': B.fiducialRMStd}
    # -------- Calculate BLOS.

    # -------- Read and unpack BLOS
    BLOSData = pd.read_csv(saveFilePath_BLOS)
    # When calculating the optimal num reference points - we want to use the original indexing so we can keep track of
    # how the values of individual points change
    BLOSData = BLOSData.set_index('ID#', drop=True)
    AllData[str(NRef)] = BLOSData['Magnetic_Field(uG)']

    Ra_blos = list(BLOSData['Ra'])
    Dec_blos = list(BLOSData['Dec(deg)'])
    Bval_blos = list(BLOSData['Magnetic_Field(uG)'])
    # -------- Read and unpack BLOS.

    # -------- PREPARE TO PLOT BLOS POINTS --------
    # ---- Convert Ra and Dec of reference points into pixel values of the fits file
    x_blos = []  # x pixel coordinate of reference
    y_blos = []  # y pixel coordinate of reference
    for i in range(len(Ra_blos)):
        pixelRow, pixelColumn = wcs.wcs_world2pix(Ra_blos[i], Dec_blos[i], 0)
        x_blos.append(pixelRow)
        y_blos.append(pixelColumn)
    # ---- Convert Ra and Dec of reference points into pixel values of the fits file.
    color, size = B2RGB(Bval_blos)
    # -------- PREPARE TO PLOT BLOS POINTS. --------

    # -------- CREATE A FIGURE - REF POINT MAP --------
    fig = plt.figure(figsize=(12, 12), dpi=120, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111, projection=wcs)

    plt.title('Reference Rotation Measure Points \n NumPoints = {}'.format(NRef), fontsize=12, y=1.08)
    im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
    plt.scatter(x_ref, y_ref, marker='o', facecolor='green', linewidth=.5, edgecolors='black')
    plt.scatter(x_blos, y_blos, marker='o', s=size, facecolor=color, linewidth=.5, edgecolors='black')

    # ---- Annotate the BLOS Points
    for i, txt in enumerate(n):
        ax.annotate(txt, (x_ref[i], y_ref[i]), size=9, color='w')
    # ---- Annotate the BLOS Points.

    # ---- Style the main axes and their grid
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
    saveFigurePath_RefPointMap = saveFigureDir_RefPointMap + 'RefPointMap_NRef' + str(NRef) + '.png'
    plt.savefig(saveFigurePath_RefPointMap)
    plt.close()
    # ---- Display or save the figure.
    # -------- CREATE A FIGURE - REF POINT MAP. --------
# -------- FIND REF POINTS AND CALCULATE BLOS AS A FUNCTION OF # REF POINTS. --------

# -------- FIND OPTIMAL NUM REF POINTS --------
# If a point has been used as a reference point at any time it will not be used to determine the optimal number of
# reference points
DataNoRef = AllData.copy().drop(indRefPoints, errors='ignore').reset_index(drop=True)

# -------- CREATE A FIGURE --------
fig = plt.figure(figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')

plt.title('Calculated BLOS value as a function of the number of reference points \n '+cloudName, fontsize=12, y=1.08)
plt.xlabel('Number of reference points')
plt.ylabel('Calculated BLOS value (micro Gauss)')

x = [int(col) for col in DataNoRef.columns]
plt.xticks(x, list(DataNoRef.columns))

cmap = plt.get_cmap('terrain')
colors = [cmap(i) for i in np.linspace(0, 1, len(DataNoRef.index))]

# For each BLOS Point
for i, identifier in enumerate(DataNoRef.index):
    plt.plot(x, list(DataNoRef.loc[identifier]), label=str(identifier), color=colors[i],
             marker='${}$'.format(identifier), markevery=int(TotalNRef/2))

plt.legend(loc='center right', bbox_to_anchor=(1.1, 0.5), ncol=2)

plt.savefig(saveFigurePath)
# -------- CREATE A FIGURE. --------

'''
We can now determine the optimal number of reference points using the calculated BLOS values as a function of number of
reference points.

The following algorithm searches for the number of reference points where the trend of calculated BLOS values stabilizes
    - The algorithm checks the values of adjacent pairs of BLOS values
    - If the BLOS values are similar, ie their difference is within a given threshold, then this is called a "run"
    - The algorithm searches for the number of reference points where the longest run occurs 


This algorithm is used on every BLOS point that was not used as a reference point.  The optimal number of reference 
points is taken to be the number of reference points where the longest run occurs most often
'''

threshold = 50
LongestRun_NumRefPoints = []

# For every point within the region of interest that was not used as a reference point:
for index in range(len(DataNoRef)):
    y = list(DataNoRef.loc[index])  # List of BLOS values for the given point
    x = list(np.arange(1, len(y) + 1))  # List of number of reference points

    # -------- STABILITY TREND ALGORITHM --------
    run_length = 1  # To keep track of how many similar values we get in a row
    on_a_run = 0  # To indicate if we are on a run (1) or not (0)
    run_started = []  # To keep track of the number of reference points where the run started
    run = []  # To keep track of the total length of each run

    # For each pair of adjacent BLOS values:
    for i in range(1, len(y)):
        # Find the absolute difference in adjacent BLOS values
        diff = abs(y[i] - y[i - 1])

        if diff <= threshold:
            # If the difference is below the threshold, then we have started, or are in, a "run" of similar values
            run_length += 1
            if on_a_run == 0:
                # We have started a new run
                on_a_run = 1
                run_started.append(x[i - 1])

        if diff > threshold or i == len(y) - 1:
            # If the difference was not below the threshold, or, if it is the end of the list of BLOS values,
            if on_a_run == 1:
                # then if we were on a run, it ends
                run.append(run_length)
                on_a_run = 0
                run_length = 1
    # -------- STABILITY TREND ALGORITHM --------

    if len(run) > 0:
        # We want to find the number of reference points where the longest run occurred
        ind_longestRun = np.where(np.array(run) == max(run))[0][0]
        LongestRun_NumRefPoints.append(run_started[ind_longestRun])

# We will take the optimal number of reference points to be the number of reference points where the longest run
# occurred most often
Optimal_NumRefPoints = mode(LongestRun_NumRefPoints)
print('Most often {} reference points appears as the optimal number of reference points'.format(Optimal_NumRefPoints))

print(FiducialData)

print(FiducialData.loc[Optimal_NumRefPoints])


# Save the optimal number of reference points and the corresponding fiducial data
(pd.DataFrame.from_dict(data=dict(FiducialData.loc[Optimal_NumRefPoints]), orient='index').T.
 to_csv(saveFilePath_Fiducial, index=False))
# -------- FIND OPTIMAL NUM REF POINTS --------
