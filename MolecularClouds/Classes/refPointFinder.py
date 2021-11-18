"""
This file contains the class to find reference points for a given region of interest

THIS IS NEW
"""
import pandas as pd
import numpy as np
import math
from astropy.io import fits
from Classes.RegionOfInterest import Region


# -------- CLASS DEFINITION --------
class RefPointFinder:
    def __init__(self, ExtincRMPath, NRefPoints, cloudName, saveFilePath):
        """
        Takes matched extinction and rotation measure data for the region of interest
        and finds NRefPoints number of reference points

        :param ExtincRMPath: Path to matched extinction and rotation measure data (produced in stage 02)
        :param NRefPoints: Number of reference points to find
        :param saveFilePath: Path to which output is saved. If 'none', no file will be saved
        """

        regionOfInterest = Region(cloudName)

        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA --------
        matchedRMExtinctionData = pd.read_csv(ExtincRMPath, sep='\t')
        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA. --------

        # -------- READ FITS FILE --------
        hdulist = fits.open(regionOfInterest.fitsFilePath)
        hdu = hdulist[0]
        # -------- READ FITS FILE. --------

        # -------- FIND ALL POTENTIAL REFERENCE POINTS --------

        # ---- Criterion 1: Av < 1
        '''
        We will only consider points with visual extinction less than 1 (Av < 1) as reference points
            - Here we extract these points and sort the resulting dataframe from smallest to greatest extinction 
        '''
        # Indices where the extinction is less than one:
        ind_extinction = np.where(matchedRMExtinctionData['Extinction_Value'] < 1.1)[0]

        # All potential reference points are all reference points with extinction less than one:
        AllRefPoints = matchedRMExtinctionData.loc[ind_extinction].sort_values('Extinction_Value', ignore_index=True)
        # ---- Criterion 1: Av < 1.

        # ---- Criterion 2: Remove anomalous rm
        # '''
        # We do not want to consider rm points as reference points if they have extremely anomalous rm values
        #     - Here we remove points with anomalous rm values
        # '''
        # AvgRM = np.mean(AllRefPoints['Rotation_Measure(rad/m2)'])
        # for i in AllRefPoints.index:
        #     if abs(AllRefPoints.loc[i]['Rotation_Measure(rad/m2)']) > 50 * AvgRM:
        #         AllRefPoints = AllRefPoints.drop(i)
        # ---- Criterion 2: Remove anomalous rm.

        # ---- Criterion 3: Remove points near a region of high extinction
        '''
        We do not want to  consider points near a region of high extinction (eg near the cloud itself) as they may not 
        sample the background
        '''
        # -------- Define the range
        # The distance the point can be from a region of high extinction and still be thought to sample the background
        # cloudDistance = regionOfInterest.distance  # [pc]
        # cloudJeansLength = 1.  # [pc] MEHRNOOSH FIND THIS
        # minDiff = cloudJeansLength / cloudDistance  # [deg]
        #
        # minDiff_pix = minDiff / abs(hdu.header['CDELT1'])
        # NDelt = 5 * math.ceil(minDiff_pix)  # Round up
        #
        # # Choose the minimum extinction value which you want to correspond to an "on" position
        # highExtinctionThreshold = 5
        # -------- Define the range.

        # # -------- For each potential reference point
        # for i in AllRefPoints.index:
        #     px = AllRefPoints['Extinction_Index_x'][i]
        #     py = AllRefPoints['Extinction_Index_y'][i]
        #
        #     # ---- Find the extinction range for the given point
        #     ind_xmax = px + NDelt + 1  # add 1 to be inclusive of the upper bound
        #     ind_ymax = py + NDelt + 1  # add 1 to be inclusive of the upper bound
        #     ind_xmin = px - NDelt
        #     ind_ymin = py - NDelt
        #     # ---- Find the extinction range for the given point.
        #
        #     # ---- Cycle through extinction values within the range
        #     # If an extinction value within this range is too high, then it cannot be considered as a reference point
        #     highExtinction = 0
        #     for pxx in range(ind_xmin, ind_xmax):
        #         for pyy in range(ind_ymin, ind_ymax):
        #             if 0 <= pxx < hdu.data.shape[1] and 0 <= pyy < hdu.data.shape[0]:
        #                 extinction_val = hdu.data[pyy, pxx]
        #                 if extinction_val > highExtinctionThreshold:
        #                     highExtinction = 1
        #     if highExtinction == 1:
        #         AllRefPoints = AllRefPoints.drop(i)
        #     # ---- Cycle through extinction values within the range.
        # # -------- For each potential reference point.
        # AllRefPoints = AllRefPoints.reset_index(drop=True)
        # ---- Criterion 3: Remove points near a region of high extinction.

        self.numAllRefPoints = len(AllRefPoints)
        # -------- FIND ALL POTENTIAL REFERENCE POINTS.  --------

        # -------- CREATE A TABLE FOR THE REFERENCE POINTS --------
        ref = pd.DataFrame()
        # -------- CREATE A TABLE FOR THE REFERENCE POINTS --------

        # -------- SELECT REFERENCE POINTS --------
        if NRefPoints != 0:
            '''
            At this point in time we are choosing reference points based on their extinction values only   
            - also this may select multiple points with the same extinction 
            '''
            # Indices of the {NRefPoints} points with the lowest extinction values:
            ind = np.arange(0, NRefPoints)
            '''
            If you want to force the code to use certain points as reference values, set / hardcode
            ind = [ {list of indices} ]
            '''
            # Extract rows corresponding to these indices - these are the reference points
            ref = AllRefPoints.loc[ind]
        # -------- SELECT REFERENCE POINTS. --------

        # -------- SAVE REFERENCE POINT DATA AS A TABLE --------
        if saveFilePath != 'none':
            ref.to_csv(saveFilePath, index=False)
        # -------- SAVE REFERENCE POINT DATA AS A TABLE. --------

# -------- CLASS DEFINITION. --------
