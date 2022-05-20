"""
This file contains the script to find the optimal number of reference points out of a given set.
It is used in the third stage of the BLOS Mapping Method.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from .CalculateB import CalculateB
from .RegionOfInterest import Region
from statistics import mode
import MolecularClouds.LocalLibraries.config as config

#Todo: Only one attribute is used: Optimal_NumRefPoints_firstMode . Consider turning into a function?

def stabilityCheckAlg(DataNoRef):
    '''
    The following algorithm searches for the number of candidate reference points where the trend of calculated BLOS
     values stabilizes
        - The algorithm checks the values of adjacent pairs of BLOS values
        - If the BLOS values are similar, ie their difference is within a given threshold, then this is called a
        "run"
        - The algorithm searches for the number of reference points where the longest run occurs

    This algorithm is used on every BLOS point that was not used as a reference point.  The optimal number of
    reference points is taken to be the number of reference points where the longest run occurs most often

    We repeat this algorithm over all potential threshold values
    :param DataNoRef:
    :return:
    '''
    UpperLimit = max([max(abs(np.diff(list(DataNoRef.loc[number])))) for number in DataNoRef.index])
    LowerLimit = min([min(abs(np.diff(list(DataNoRef.loc[number])))) for number in DataNoRef.index])
    thresholds = np.linspace(LowerLimit, UpperLimit, 500)

    Optimal_NumRefPoints = []

    for threshold in thresholds:

        LongestRun_NumRefPoints = []
        LongestRun_Length = []

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
                    # If the difference is below the threshold, then we have started, or are in, a "run" of similar
                    # values
                    run_length += 1
                    if on_a_run == 0:
                        # We have started a new run
                        on_a_run = 1
                        run_started.append(x[i - 1])

                if diff > threshold or i == len(y) - 1:
                    # If the difference was not below the threshold, or, if it is the end of the list of BLOS values
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
                LongestRun_Length.append(max(run))
        Optimal_NumRefPoints.append(mode(LongestRun_NumRefPoints))
    return Optimal_NumRefPoints

def stabilityTrendGraph(DataNoRef, saveFigurePath):
    Identifiers = list(DataNoRef.index)
    # -------- CREATE A FIGURE --------
    plt.figure(figsize=(6, 4), dpi=120, facecolor='w', edgecolor='k')

    plt.title('Calculated BLOS value as a function of the number of reference points \n ' + config.cloud, fontsize=12,
              y=1.08)
    plt.xlabel('Number of reference points')
    plt.ylabel('Calculated BLOS value ' + r'($\mu G$)')

    x = [int(col) for col in DataNoRef.columns]
    plt.xticks(x, list(DataNoRef.columns))

    cmap = plt.get_cmap('terrain')
    colors = [cmap(i) for i in np.linspace(0, 1, len(DataNoRef.index))]

    # For each BLOS Point
    for i, number in enumerate(DataNoRef.index):
        plt.plot(x, list(DataNoRef.loc[number]), '-o', label=str(Identifiers[i]), color=colors[i], markersize=3)

    # plt.legend(loc='center right', bbox_to_anchor=(1.1, 0.5), ncol=2, framealpha=1, title='Identification Number')

    if saveFigurePath is not None:
        plt.savefig(saveFigurePath)
    plt.show()
    plt.close()
    print('Saving stability trend to ' + saveFigurePath)
    # -------- CREATE A FIGURE. --------

def findTrendData(potentialRefPoints, MatchedRMExtincPath, regionOfInterest):
    # -------- CREATE A TABLE FOR ALL BLOS DATA --------
    # The rows of this table will represent the number of reference points and the columns of this table will
    # represent the individual BLOS points.  Each entry in the table is a calculated BLOS value.
    AllData = pd.DataFrame()
    # -------- CREATE A TABLE FOR ALL BLOS DATA. --------

    # -------- CALCULATE BLOS AS A FUNCTION OF # REF POINTS --------
    for num in range(len(potentialRefPoints)):
        # -------- Extract {num} points from the table of potential reference points
        # The extracted potential reference points will be referred to as "candidates"
        candidateRefPoints = potentialRefPoints.loc[:num]
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
    DataNoRef = AllData.copy().drop(list(potentialRefPoints['ID#']), errors='ignore')
    DataNoRef = DataNoRef.reset_index(drop=True)
    return DataNoRef

# -------- CLASS DEFINITION --------
class FindOptimalRefPoints:
    def __init__(self, cloudName, potentialRefPoints, saveFigurePath):
        """
        This class will determine the optimal number of reference points from a table of potential reference points

        :param cloudName: Name of the region of interest
        :param potentialRefPoints: Table of potential reference points
        :param saveFigurePath: Path to where the Blos vs Nref figure is saved
        """

        regionOfInterest = Region(config.cloud)

        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA --------
        MatchedRMExtincPath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_RMExtinctionMatch + config.cloud + '.txt')
        matchedRMExtinctionData = pd.read_csv(MatchedRMExtincPath, sep='\t')
        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA. --------

        DataNoRef = findTrendData(potentialRefPoints, MatchedRMExtincPath, regionOfInterest)
        DataNoRef.to_csv(
            os.path.join(config.dir_root, config.dir_fileOutput, config.prefix_OptRefPoints + config.cloud + '.txt'))

        stabilityTrendGraph(DataNoRef, saveFigurePath)

        TotalNumPoints = len(matchedRMExtinctionData)
        '''
        We can now determine the optimal number of reference points using the calculated BLOS values as a function of 
        number of candidate reference points.
        '''
        Optimal_NumRefPoints = stabilityCheckAlg(DataNoRef)

        # -------- FIND OPTIMAL NUM REF POINTS --------
        # The number of reference points should be greater than 3 and less than half the total number of points

        Optimal_NumRefPoints_Selection = [value for value in Optimal_NumRefPoints if 4 < value < 0.5 * TotalNumPoints]

        self.Optimal_NumRefPoints_firstMode = mode(Optimal_NumRefPoints_Selection)
        Optimal_NumRefPoints_firstMode_percent = (len(
            np.where(Optimal_NumRefPoints_Selection == self.Optimal_NumRefPoints_firstMode)[0]) / len(
            Optimal_NumRefPoints_Selection)) * 100

        Optimal_NumRefPoints_removeFirstMode = [i for i in Optimal_NumRefPoints_Selection if
                                                i != self.Optimal_NumRefPoints_firstMode]

        if len(Optimal_NumRefPoints_removeFirstMode) != 0:
            self.Optimal_NumRefPoints_secondMode = mode(Optimal_NumRefPoints_removeFirstMode)
            Optimal_NumRefPoints_secondMode_percent = (len(
                np.where(Optimal_NumRefPoints_Selection == self.Optimal_NumRefPoints_secondMode)[0]) / len(
                Optimal_NumRefPoints_Selection)) * 100
            print('\t - {} appears second most often ({:.0f}%) as the optimal number of reference points'.format(
                self.Optimal_NumRefPoints_secondMode, Optimal_NumRefPoints_secondMode_percent))

        self.Optimal_NumRefPoints_max = max(Optimal_NumRefPoints_Selection)
        self.Optimal_NumRefPoints_min = min(Optimal_NumRefPoints_Selection)

        print('\t - {} appears most often ({:.0f}%) as the optimal number of reference points'.format(
            self.Optimal_NumRefPoints_firstMode, Optimal_NumRefPoints_firstMode_percent))

        print('\t - {} is the maximum optimal number of reference points.'.format(
            self.Optimal_NumRefPoints_max))
        print('\t - {} is the minimum optimal number of reference points.'.format(
            self.Optimal_NumRefPoints_min))
        # -------- FIND OPTIMAL NUM REF POINTS. --------

# -------- CLASS DEFINITION. --------
